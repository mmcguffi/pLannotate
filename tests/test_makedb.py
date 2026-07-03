"""Tests for building ad-hoc custom databases via ``plannotate makedb``."""

import shutil

import pandas as pd
import pytest
import yaml
from typer.testing import CliRunner

from plannotate import _database_builder, _package_data
from plannotate import main as main_module
from plannotate._sqlite import (
    load_descriptions_from_sqlite,
    write_descriptions_to_sqlite,
)
from plannotate.main import app


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


def test_build_source_config_with_descriptions_uses_default_details(tmp_path):
    config = _database_builder.build_source_config(
        "mydb", "diamond", tmp_path, priority=2, has_descriptions=True
    )
    assert config["method"] == "diamond"
    assert config["location"] == str(tmp_path.resolve())
    assert config["priority"] == 2
    # a descriptions file is present beside the index -> Default resolves to it
    assert config["details"] == {"default_type": None, "location": "Default"}


def test_build_source_config_without_descriptions_synthesizes(tmp_path):
    config = _database_builder.build_source_config(
        "mydb", "diamond", tmp_path, priority=1, has_descriptions=False
    )
    # protein hits without a descriptions file default to CDS with names from ids
    assert config["details"] == {"default_type": "CDS", "location": None}


def test_build_full_config_layers_source_onto_builtins(tmp_path):
    config = _database_builder.build_full_config(
        "mydb", "blastn", tmp_path, 1, has_descriptions=False, include_builtins=True
    )
    # builtins are carried through and the custom source is appended
    assert "snapgene" in config
    assert "swissprot" in config
    assert "mydb" in config
    assert config["snapgene"]["location"] == "Default"


def test_build_full_config_custom_only(tmp_path):
    config = _database_builder.build_full_config(
        "mydb", "blastn", tmp_path, 1, has_descriptions=False, include_builtins=False
    )
    assert set(config) == {"mydb"}


def test_build_full_config_rejects_builtin_name_collision(tmp_path):
    with pytest.raises(ValueError, match="collides with a builtin"):
        _database_builder.build_full_config(
            "snapgene", "blastn", tmp_path, 1, False, include_builtins=True
        )


def test_full_config_is_loadable_by_get_yaml(tmp_path, monkeypatch):
    # the emitted config for the custom source must survive get_yaml's validation
    frame = pd.DataFrame([{"sseqid": "featA", "name": "A", "type": "CDS"}])
    write_descriptions_to_sqlite("mydb", frame, tmp_path / "mydb.db")
    # give the source a resolvable-looking index path get_yaml can build db_loc from
    (tmp_path / "mydb.dmnd").touch()

    config = {
        "mydb": _database_builder.build_source_config(
            "mydb", "diamond", tmp_path, 1, has_descriptions=True
        )
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
