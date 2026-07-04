"""Tests for building ad-hoc custom databases via ``plannotate makedb``."""

import shutil
from pathlib import Path

import pandas as pd
import pytest
import yaml
from Bio import SeqIO
from typer.testing import CliRunner

from plannotate import _database_builder, _package_data
from plannotate import main as main_module
from plannotate._sqlite import (
    load_descriptions_from_sqlite,
    write_descriptions_to_sqlite,
)
from plannotate.main import app

# shared example inputs for the end-to-end makedb tests
MAKEDB_DATA = Path(__file__).parent / "test_data" / "makedb"
# features present in demo_plasmid.fa that a blastn search reliably recovers (the
# 18 bp His6 repeat is intentionally dropped by overlap resolution)
DEMO_PLASMID_FEATURES = {
    "T7_promoter",
    "lac_operator",
    "RBS",
    "FLAG_tag",
    "SV40_NLS",
    "T7_terminator",
}


def _expected_config(fixture_name: str, out_dir: Path) -> dict:
    """Load an expected-config fixture, substituting the run-specific output dir."""
    text = (MAKEDB_DATA / fixture_name).read_text()
    text = text.replace("__OUTDIR__", str(out_dir.resolve()))
    return yaml.safe_load(text)


def test_normalize_method_accepts_aliases():
    assert _database_builder.normalize_method("BLAST") == "blastn"
    assert _database_builder.normalize_method("nucleotide") == "blastn"
    assert _database_builder.normalize_method("protein") == "diamond"


def test_normalize_method_rejects_unknown():
    with pytest.raises(ValueError, match="Unsupported method"):
        _database_builder.normalize_method("bowtie")


def test_write_descriptions_round_trips_through_loader(tmp_path):
    frame = pd.DataFrame(
        [{"sseqid": "featA", "name": "Feature A", "type": "CDS", "blurb": "a gene"}]
    )
    db_path = tmp_path / "mydb.db"

    count = write_descriptions_to_sqlite("mydb", frame, db_path)

    assert count == 1
    loaded = load_descriptions_from_sqlite(
        "mydb", {"featA"}, {"db_loc": str(tmp_path / "mydb")}
    )
    pd.testing.assert_frame_equal(loaded, frame)


def test_write_descriptions_accepts_column_aliases_and_defaults(tmp_path):
    # accession -> sseqid, Description -> blurb; name/type absent so they default
    frame = pd.DataFrame([{"accession": "featB", "Description": "b summary"}])
    db_path = tmp_path / "aliased.db"

    write_descriptions_to_sqlite("aliased", frame, db_path)

    loaded = load_descriptions_from_sqlite(
        "aliased", None, {"db_loc": str(tmp_path / "aliased")}
    )
    row = loaded.iloc[0]
    assert row["sseqid"] == "featB"
    assert row["name"] == "featB"  # defaulted from the id
    assert row["type"] == "misc_feature"
    assert row["blurb"] == "b summary"


def test_write_descriptions_requires_id_column(tmp_path):
    frame = pd.DataFrame([{"name": "no id", "blurb": "x"}])
    with pytest.raises(ValueError, match="id column"):
        write_descriptions_to_sqlite("bad", frame, tmp_path / "bad.db")


def test_write_descriptions_rejects_colliding_aliases(tmp_path):
    # two columns aliasing to the same canonical field (id and accession -> sseqid)
    # must be rejected up front, not crash opaquely inside to_sql
    frame = pd.DataFrame([{"id": "featA", "accession": "featB", "name": "A"}])
    with pytest.raises(ValueError, match="same field"):
        write_descriptions_to_sqlite("dup", frame, tmp_path / "dup.db")


def test_write_descriptions_dedupes_on_sseqid(tmp_path):
    # duplicate ids would fan a single hit into several annotations at merge time;
    # the first row wins and the rest are dropped
    frame = pd.DataFrame(
        [
            {"sseqid": "featA", "name": "first", "type": "CDS", "blurb": "keep"},
            {"sseqid": "featA", "name": "second", "type": "CDS", "blurb": "drop"},
            {"sseqid": "featB", "name": "other", "type": "CDS", "blurb": "b"},
        ]
    )

    count = write_descriptions_to_sqlite("dedup", frame, tmp_path / "dedup.db")

    assert count == 2
    loaded = load_descriptions_from_sqlite(
        "dedup", None, {"db_loc": str(tmp_path / "dedup")}
    )
    by_id = loaded.set_index("sseqid")
    assert by_id.loc["featA", "name"] == "first"  # first row wins
    assert set(loaded["sseqid"]) == {"featA", "featB"}


def test_build_source_config_uses_default_details(tmp_path):
    config = _database_builder.build_source_config(
        "mydb", "diamond", tmp_path, priority=2
    )
    assert config["method"] == "diamond"
    assert config["location"] == str(tmp_path.resolve())
    assert config["priority"] == 2
    # a descriptions db is always written beside the index -> Default resolves to it
    assert config["details"] == {"default_type": None, "location": "Default"}


def test_descriptions_from_fasta_uses_headers(tmp_path):
    fasta = tmp_path / "feats.fasta"
    fasta.write_text(">featA a helpful summary\nACGTACGT\n>featB\nTTTT\n")

    frame = _database_builder.descriptions_from_fasta(fasta, "blastn")

    by_id = frame.set_index("sseqid")
    # id becomes the name; trailing header text becomes the blurb
    assert by_id.loc["featA", "name"] == "featA"
    assert by_id.loc["featA", "blurb"] == "a helpful summary"
    assert by_id.loc["featB", "blurb"] == ""
    # nucleotide features default to misc_feature; protein would be CDS
    assert set(frame["type"]) == {"misc_feature"}
    assert (
        _database_builder.descriptions_from_fasta(fasta, "diamond")["type"].iloc[0]
        == "CDS"
    )


def test_descriptions_from_fasta_normalizes_diamond_pipe_ids(tmp_path):
    # a protein FASTA with UniProt-style headers: diamond.search reports the hit
    # sseqid as the accession (P12345), so the descriptions must be keyed by it,
    # not the full pipe id, or the per-hit description merge misses.
    fasta = tmp_path / "prot.fasta"
    fasta.write_text(">sp|P12345|GFP_AEQVI Green fluorescent protein\nMSKGEELFTG\n")

    frame = _database_builder.descriptions_from_fasta(fasta, "diamond")

    row = frame.iloc[0]
    assert row["sseqid"] == "P12345"  # matches the normalized diamond hit id
    assert row["name"] == "sp|P12345|GFP_AEQVI"  # original id kept as the label
    assert row["blurb"] == "Green fluorescent protein"


def test_descriptions_from_fasta_keeps_pipe_ids_for_blastn(tmp_path):
    # blastn does not unwrap accessions, so a nucleotide id with a pipe is kept
    fasta = tmp_path / "nucl.fasta"
    fasta.write_text(">gene|001 some gene\nACGTACGT\n")

    frame = _database_builder.descriptions_from_fasta(fasta, "blastn")

    assert frame.iloc[0]["sseqid"] == "gene|001"


def test_descriptions_frame_normalizes_csv_pipe_ids_for_diamond(tmp_path):
    # a CSV whose ids match UniProt-style FASTA headers must be re-keyed to the
    # accession for diamond, since the search adapter reports hits that way.
    csv = tmp_path / "desc.csv"
    csv.write_text("sseqid,name,type,blurb\nsp|P12345|GFP,GFP,CDS,green\n")

    frame = _database_builder._descriptions_frame(
        tmp_path / "unused.fasta", csv, "diamond"
    )

    assert frame.iloc[0]["sseqid"] == "P12345"
    assert frame.iloc[0]["name"] == "GFP"
    assert frame.iloc[0]["blurb"] == "green"


def test_descriptions_frame_keeps_csv_pipe_ids_for_blastn(tmp_path):
    csv = tmp_path / "desc.csv"
    csv.write_text("sseqid,name\ngene|001,my gene\n")

    frame = _database_builder._descriptions_frame(
        tmp_path / "unused.fasta", csv, "blastn"
    )

    assert frame.iloc[0]["sseqid"] == "gene|001"


def test_descriptions_from_fasta_handles_space_after_gt(tmp_path):
    # "> feature" (space after >) parses to id "feature"; the blurb must be empty,
    # not a garbage fragment from slicing the id length off a space-padded header.
    fasta = tmp_path / "spaced.fasta"
    fasta.write_text("> feature\nACGTACGTACGT\n")

    frame = _database_builder.descriptions_from_fasta(fasta, "blastn")

    assert list(frame["sseqid"]) == ["feature"]
    assert frame.iloc[0]["name"] == "feature"
    assert frame.iloc[0]["blurb"] == ""


def test_build_full_config_layers_source_onto_builtins(tmp_path):
    config = _database_builder.build_full_config(
        "mydb", "blastn", tmp_path, 1, include_builtins=True
    )
    # builtins are carried through and the custom source is appended
    assert "snapgene" in config
    assert "swissprot" in config
    assert "mydb" in config
    assert config["snapgene"]["location"] == "Default"


def test_build_full_config_custom_only(tmp_path):
    config = _database_builder.build_full_config(
        "mydb", "blastn", tmp_path, 1, include_builtins=False
    )
    assert set(config) == {"mydb"}


def test_build_full_config_rejects_builtin_name_collision(tmp_path):
    with pytest.raises(ValueError, match="collides with a builtin"):
        _database_builder.build_full_config(
            "snapgene", "blastn", tmp_path, 1, include_builtins=True
        )


def test_full_config_is_loadable_by_get_yaml(tmp_path, monkeypatch):
    # the emitted config for the custom source must survive get_yaml's validation
    frame = pd.DataFrame([{"sseqid": "featA", "name": "A", "type": "CDS"}])
    write_descriptions_to_sqlite("mydb", frame, tmp_path / "mydb.db")
    # give the source a resolvable-looking index path get_yaml can build db_loc from
    (tmp_path / "mydb.dmnd").touch()

    config = {
        "mydb": _database_builder.build_source_config("mydb", "diamond", tmp_path, 1)
    }
    yaml_path = tmp_path / "databases.yml"
    yaml_path.write_text(yaml.safe_dump(config, sort_keys=False))

    loaded = _package_data.get_yaml(yaml_path)
    assert loaded["mydb"]["db_loc"] == str(tmp_path / "mydb")


def test_config_references_builtin_databases(tmp_path):
    # the packaged config uses Default (builtin) locations
    assert _package_data.config_references_builtin_databases(
        _package_data.get_yaml_path()
    )

    custom = tmp_path / "custom.yml"
    custom.write_text(
        yaml.safe_dump(
            {
                "mydb": {
                    "method": "diamond",
                    "location": str(tmp_path),
                    "priority": 1,
                    "details": {"location": None},
                }
            }
        )
    )
    assert not _package_data.config_references_builtin_databases(custom)


class _StubConstruct:
    """Stand-in Construct that skips the real search for CLI plumbing tests."""

    def __init__(self, *args, **kwargs):
        pass

    def to_genbank(self):
        return "LOCUS stub\n//\n"


def test_batch_with_custom_only_config_skips_bundle_gate(monkeypatch, tmp_path):
    # a fully custom config must annotate even without the downloaded bundle
    monkeypatch.setattr(_package_data, "databases_exist", lambda: False)
    monkeypatch.setattr(main_module, "Construct", _StubConstruct)

    custom = tmp_path / "custom.yml"
    custom.write_text(
        yaml.safe_dump(
            {
                "mydb": {
                    "method": "diamond",
                    "location": str(tmp_path),
                    "priority": 1,
                    "details": {"location": None},
                }
            }
        )
    )
    seq = tmp_path / "seq.fa"
    seq.write_text(">s\nATGCATGCATGC\n")

    result = CliRunner().invoke(
        app, ["batch", "-i", str(seq), "-o", str(tmp_path), "-y", str(custom), "-l"]
    )

    assert result.exit_code == 0
    assert "not downloaded" not in result.stdout


def test_batch_with_builtin_config_still_requires_bundle(monkeypatch, tmp_path, caplog):
    monkeypatch.setattr(_package_data, "databases_exist", lambda: False)
    seq = tmp_path / "seq.fa"
    seq.write_text(">s\nATGCATGCATGC\n")

    result = CliRunner().invoke(app, ["batch", "-i", str(seq), "-o", str(tmp_path)])

    assert result.exit_code == 1
    assert "setupdb" in caplog.text


def test_cli_makedb_writes_config(monkeypatch, tmp_path):
    fasta = tmp_path / "seqs.fasta"
    fasta.write_text(">a\nATGC\n")
    captured = {}

    def _fake_make_database(input_file, name, method, out_dir, **kwargs):
        captured["args"] = (input_file, name, method, out_dir, kwargs)
        return {name: {"method": method, "location": str(out_dir)}}

    monkeypatch.setattr(
        main_module._database_builder, "make_database", _fake_make_database
    )

    result = CliRunner().invoke(
        app,
        [
            "makedb",
            "-i",
            str(fasta),
            "-n",
            "mydb",
            "-m",
            "diamond",
            "-o",
            str(tmp_path),
        ],
    )

    assert result.exit_code == 0
    yaml_path = tmp_path / "databases.yml"
    assert yaml_path.is_file()
    written = yaml.safe_load(yaml_path.read_text())
    assert "mydb" in written
    assert captured["args"][1] == "mydb"


def test_cli_makedb_reports_build_errors(monkeypatch, tmp_path, caplog):
    fasta = tmp_path / "seqs.fasta"
    fasta.write_text(">a\nATGC\n")

    def _boom(*args, **kwargs):
        raise RuntimeError("'diamond' was not found on PATH")

    monkeypatch.setattr(main_module._database_builder, "make_database", _boom)

    result = CliRunner().invoke(
        app, ["makedb", "-i", str(fasta), "-n", "mydb", "-m", "diamond"]
    )

    assert result.exit_code == 1
    assert "not found on PATH" in caplog.text


@pytest.mark.integration
def test_makedb_builds_diamond_index_and_annotates(tmp_path):
    if shutil.which("diamond") is None:
        pytest.skip("diamond not installed")

    # a short protein and a nucleotide sequence that translates to it
    protein = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVP"
    fasta = tmp_path / "prot.fasta"
    fasta.write_text(f">gfp_frag\n{protein}\n")
    csv = tmp_path / "desc.csv"
    csv.write_text("sseqid,name,type,blurb\ngfp_frag,GFP fragment,CDS,a fluorophore\n")

    config = _database_builder.make_database(
        fasta,
        "customdb",
        "diamond",
        tmp_path,
        csv=csv,
        include_builtins=False,
    )

    assert (tmp_path / "customdb.dmnd").is_file()
    assert (tmp_path / "customdb.db").is_file()

    yaml_path = tmp_path / "databases.yml"
    yaml_path.write_text(yaml.safe_dump(config, sort_keys=False))

    from plannotate.annotate import annotate

    # back-translate the protein to a DNA query the diamond index can recover
    codon = {
        "A": "GCT",
        "R": "CGT",
        "N": "AAT",
        "D": "GAT",
        "C": "TGT",
        "Q": "CAA",
        "E": "GAA",
        "G": "GGT",
        "H": "CAT",
        "I": "ATT",
        "L": "CTT",
        "K": "AAA",
        "M": "ATG",
        "F": "TTT",
        "P": "CCT",
        "S": "TCT",
        "T": "ACT",
        "W": "TGG",
        "Y": "TAT",
        "V": "GTT",
    }
    dna = "".join(codon[aa] for aa in protein)

    annotations = annotate(dna, yaml_file=yaml_path, linear=True)
    assert (annotations["sseqid"] == "gfp_frag").any()


@pytest.mark.integration
def test_makedb_without_csv_still_annotates(tmp_path):
    # regression: omitting --csv must still produce an annotatable database, with
    # feature names taken from the FASTA headers.
    if shutil.which("makeblastdb") is None:
        pytest.skip("makeblastdb not installed")

    feature = "TAATACGACTCACTATAGGGGAATTGTGAGCGGATAACAATTCCCCTCTAGA"
    fasta = tmp_path / "feats.fasta"
    fasta.write_text(f">T7_lead T7 promoter and operator\n{feature}\n")

    config = _database_builder.make_database(
        fasta, "hdrsdb", "blastn", tmp_path, include_builtins=False
    )
    assert (tmp_path / "hdrsdb.db").is_file()

    yaml_path = tmp_path / "databases.yml"
    yaml_path.write_text(yaml.safe_dump(config, sort_keys=False))

    from plannotate.annotate import annotate

    plasmid = "GCGCGCGCGC" + feature + "GCGCGCGCGC"
    annotations = annotate(plasmid, yaml_file=yaml_path, linear=True)
    hit = annotations[annotations["sseqid"] == "T7_lead"]
    assert not hit.empty
    # name defaults to the FASTA id; blurb comes from the header text
    assert hit.iloc[0]["name"] == "T7_lead"
    assert "promoter" in hit.iloc[0]["blurb"]


@pytest.mark.integration
def test_makedb_diamond_no_csv_pipe_header_annotates(tmp_path):
    # regression (Codex review): a no-CSV protein DB whose header is a UniProt-style
    # pipe id must still attach its synthesized description, because diamond reports
    # the hit sseqid as the accession.
    if shutil.which("diamond") is None:
        pytest.skip("diamond not installed")

    protein = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVP"
    fasta = tmp_path / "prot.fasta"
    fasta.write_text(f">sp|P12345|GFP_AEQVI Green fluorescent protein\n{protein}\n")

    config = _database_builder.make_database(
        fasta, "pipedb", "diamond", tmp_path, include_builtins=False
    )
    yaml_path = tmp_path / "databases.yml"
    yaml_path.write_text(yaml.safe_dump(config, sort_keys=False))

    from plannotate.annotate import annotate

    codon = {
        "A": "GCT",
        "R": "CGT",
        "N": "AAT",
        "D": "GAT",
        "C": "TGT",
        "Q": "CAA",
        "E": "GAA",
        "G": "GGT",
        "H": "CAT",
        "I": "ATT",
        "L": "CTT",
        "K": "AAA",
        "M": "ATG",
        "F": "TTT",
        "P": "CCT",
        "S": "TCT",
        "T": "ACT",
        "W": "TGG",
        "Y": "TAT",
        "V": "GTT",
    }
    dna = "".join(codon[aa] for aa in protein)
    annotations = annotate(dna, yaml_file=yaml_path, linear=True)

    hit = annotations[annotations["sseqid"] == "P12345"]
    assert not hit.empty
    # the synthesized description is attached (would be empty if the merge missed)
    assert hit.iloc[0]["name"] == "sp|P12345|GFP_AEQVI"
    assert "fluorescent" in hit.iloc[0]["blurb"]


@pytest.mark.integration
def test_makedb_diamond_csv_pipe_ids_attach_descriptions(tmp_path):
    # regression (Codex review): a CSV whose ids match UniProt-style FASTA headers
    # must still attach for diamond, whose hits report the unwrapped accession.
    if shutil.which("diamond") is None:
        pytest.skip("diamond not installed")

    protein = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVP"
    fasta = tmp_path / "prot.fasta"
    fasta.write_text(f">sp|P12345|GFP_AEQVI\n{protein}\n")
    csv = tmp_path / "desc.csv"
    csv.write_text(
        "sseqid,name,type,blurb\nsp|P12345|GFP_AEQVI,GFP,CDS,a bright fluorophore\n"
    )

    config = _database_builder.make_database(
        fasta, "csvpipe", "diamond", tmp_path, csv=csv, include_builtins=False
    )
    yaml_path = tmp_path / "databases.yml"
    yaml_path.write_text(yaml.safe_dump(config, sort_keys=False))

    from plannotate.annotate import annotate

    codon = {
        "A": "GCT",
        "R": "CGT",
        "N": "AAT",
        "D": "GAT",
        "C": "TGT",
        "Q": "CAA",
        "E": "GAA",
        "G": "GGT",
        "H": "CAT",
        "I": "ATT",
        "L": "CTT",
        "K": "AAA",
        "M": "ATG",
        "F": "TTT",
        "P": "CCT",
        "S": "TCT",
        "T": "ACT",
        "W": "TGG",
        "Y": "TAT",
        "V": "GTT",
    }
    dna = "".join(codon[aa] for aa in protein)
    annotations = annotate(dna, yaml_file=yaml_path, linear=True)

    hit = annotations[annotations["sseqid"] == "P12345"]
    assert not hit.empty
    # the CSV-supplied name/blurb are attached, not generic fallbacks
    assert hit.iloc[0]["name"] == "GFP"
    assert "fluorophore" in hit.iloc[0]["blurb"]


def test_output_yaml_custom_only_matches_expected(tmp_path):
    config = _database_builder.build_full_config(
        "mylab", "blastn", tmp_path, priority=1, include_builtins=False
    )
    assert config == _expected_config("expected_databases.yml", tmp_path)


def test_output_yaml_diamond_uses_protein_parameters(tmp_path):
    config = _database_builder.build_source_config("mylab", "diamond", tmp_path, 2)
    assert config["method"] == "diamond"
    assert config["priority"] == 2
    assert config["parameters"] == [
        "-k 0",
        "--min-orf 1",
        "--matrix BLOSUM90",
        "--gapopen 10",
        "--gapextend 1",
        "--algo ctg",
        "--id 50",
        "--max-hsps 10",
        "--culling-overlap 200",
        "--seed-cut .001",
        "--comp-based-stats 0",
    ]


def test_output_yaml_layers_onto_builtins_verbatim(tmp_path):
    config = _database_builder.build_full_config(
        "mylab", "blastn", tmp_path, 1, include_builtins=True
    )
    builtins = yaml.safe_load(_package_data.get_yaml_path().read_text())
    # every packaged source is carried through byte-for-byte...
    for source_name, entry in builtins.items():
        assert config[source_name] == entry
    # ...and the one custom source is appended after them
    assert set(config) == set(builtins) | {"mylab"}
    assert list(config)[-1] == "mylab"


def test_output_yaml_round_trips_through_serialization(tmp_path):
    config = _database_builder.build_full_config(
        "mylab", "diamond", tmp_path, 1, include_builtins=True
    )
    written = tmp_path / "databases.yml"
    written.write_text(yaml.safe_dump(config, sort_keys=False))
    assert yaml.safe_load(written.read_text()) == config


@pytest.mark.integration
def test_makedb_end_to_end_with_fixture_files(tmp_path):
    # exercise the shipped example FASTA/CSV/plasmid as real test artifacts
    if shutil.which("makeblastdb") is None:
        pytest.skip("makeblastdb not installed")

    config = _database_builder.make_database(
        MAKEDB_DATA / "features.fasta",
        "mylab",
        "blastn",
        tmp_path,
        csv=MAKEDB_DATA / "descriptions.csv",
        include_builtins=False,
    )

    yaml_path = tmp_path / "databases.yml"
    yaml_path.write_text(yaml.safe_dump(config, sort_keys=False))
    # the written config is exactly the expected custom-only source
    assert yaml.safe_load(yaml_path.read_text()) == _expected_config(
        "expected_databases.yml", tmp_path
    )

    from plannotate.annotate import annotate

    with open(MAKEDB_DATA / "demo_plasmid.fa") as handle:
        plasmid = str(SeqIO.read(handle, "fasta").seq)
    annotations = annotate(plasmid, yaml_file=yaml_path, linear=True)

    assert DEMO_PLASMID_FEATURES <= set(annotations["sseqid"])
    # descriptions from the CSV are attached to the hits
    t7 = annotations[annotations["sseqid"] == "T7_promoter"].iloc[0]
    assert t7["name"] == "T7 promoter"
    assert t7["type"] == "promoter"
    assert "T7 RNA polymerase" in t7["blurb"]
