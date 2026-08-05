"""Unit tests for command-line behavior."""

from Bio import SeqIO
from typer.testing import CliRunner

from plannotate import __version__, _package_data
from plannotate import main as main_module
from plannotate.main import app
from plannotate.models import Construct


def _fake_batch_constructs(records, **kwargs):
    """Build feature-free constructs without running a search (CLI plumbing test)."""
    return [
        Construct(seq=seq, _skip_annotation=True, name=name)
        for name, seq, _prior in records
    ]


def test_cli_help_without_databases():
    result = CliRunner().invoke(app, ["--help"])

    assert result.exit_code == 0
    assert "setupdb" in result.stdout
    assert "batch" in result.stdout
    assert "streamlit" in result.stdout


def test_cli_version(monkeypatch):
    monkeypatch.setattr(
        main_module._package_data,
        "get_database_manifest",
        lambda: {
            "bundle": "plannotate-databases-v2",
            "build_date": "2026-06-20",
            "databases": {"Rfam": {"version": "release 15.1"}},
        },
    )

    result = CliRunner().invoke(app, ["--version"])

    assert result.exit_code == 0
    assert __version__ in result.stdout
    assert "database path" in result.stdout
    assert "plannotate-databases-v2" in result.stdout
    assert "2026-06-20" in result.stdout
    # per-source versions are listed
    assert "Rfam" in result.stdout
    assert "release 15.1" in result.stdout


def test_cli_version_without_databases(monkeypatch):
    def _missing():
        raise FileNotFoundError

    monkeypatch.setattr(main_module._package_data, "get_database_manifest", _missing)

    result = CliRunner().invoke(app, ["--version"])

    assert result.exit_code == 0
    assert __version__ in result.stdout
    assert "not installed" in result.stdout


def test_streamlit_requires_the_server_extra(monkeypatch, caplog):
    monkeypatch.setattr(main_module, "_streamlit_available", lambda: False)

    result = CliRunner().invoke(app, ["streamlit"])

    assert result.exit_code == 1
    assert "plannotate[server]" in caplog.text


def test_streamlit_requires_databases(monkeypatch, caplog):
    monkeypatch.setattr(main_module, "_streamlit_available", lambda: True)
    monkeypatch.setattr(_package_data, "databases_exist", lambda: False)

    result = CliRunner().invoke(app, ["streamlit"])

    assert result.exit_code == 1
    assert "setupdb" in caplog.text


def test_batch_requires_input():
    result = CliRunner().invoke(app, ["batch"])

    assert result.exit_code == 2


def test_setupdb_force_reinstalls_existing_databases(monkeypatch):
    downloads = []
    monkeypatch.setattr(_package_data, "databases_exist", lambda: True)
    monkeypatch.setattr(
        _package_data, "download_databases", lambda: downloads.append(True)
    )

    result = CliRunner().invoke(app, ["setupdb", "--force"])

    assert result.exit_code == 0
    assert downloads == [True]


def test_databases_prints_installed_manifest(monkeypatch):
    monkeypatch.setattr(
        _package_data,
        "get_database_manifest",
        lambda: {"schema_version": 1, "databases": {"Rfam": {"version": "15.1"}}},
    )

    result = CliRunner().invoke(app, ["databases"])

    assert result.exit_code == 0
    assert '"Rfam"' in result.stdout
    assert '"15.1"' in result.stdout


def test_batch_multi_record_writes_one_output_per_record(monkeypatch, tmp_path):
    monkeypatch.setattr(_package_data, "databases_exist", lambda: True)
    monkeypatch.setattr(
        main_module.Construct,
        "annotate_batch",
        staticmethod(_fake_batch_constructs),
    )

    fasta = tmp_path / "multi.fa"
    fasta.write_text(">plasmidA\nACGTACGTACGT\n>plasmidB\nTTTTGGGGCCCC\n")
    output = tmp_path / "out"

    result = CliRunner().invoke(
        app, ["batch", "-i", str(fasta), "-o", str(output), "--csv"]
    )

    assert result.exit_code == 0, result.stdout
    assert (output / "plasmidA_pLann.gbk").exists()
    assert (output / "plasmidB_pLann.gbk").exists()
    assert (output / "plasmidA_pLann.csv").exists()
    assert (output / "plasmidB_pLann.csv").exists()


def _skip_search(monkeypatch):
    """Build constructs on the single-record path without running a search."""

    class _NoSearchConstruct(Construct):
        def __init__(self, **kwargs):
            super().__init__(**{**kwargs, "_skip_annotation": True})

    monkeypatch.setattr(main_module, "Construct", _NoSearchConstruct)


def test_batch_names_locus_after_the_record_not_the_file(monkeypatch, tmp_path):
    monkeypatch.setattr(_package_data, "databases_exist", lambda: True)
    _skip_search(monkeypatch)

    fasta = tmp_path / "some file name.fa"
    fasta.write_text(">plasmidA\nACGTACGTACGT\n")
    output = tmp_path / "out"

    result = CliRunner().invoke(app, ["batch", "-i", str(fasta), "-o", str(output)])

    assert result.exit_code == 0, result.stdout
    # the file name still names the output file, the record still names the locus
    written = output / "some file name_pLann.gbk"
    assert written.exists()
    assert SeqIO.read(written, "genbank").name == "plasmidA"


def test_batch_locus_name_agrees_across_single_and_multi_record(monkeypatch, tmp_path):
    monkeypatch.setattr(_package_data, "databases_exist", lambda: True)
    _skip_search(monkeypatch)
    monkeypatch.setattr(
        main_module.Construct,
        "annotate_batch",
        staticmethod(_fake_batch_constructs),
    )

    single = tmp_path / "single.fa"
    single.write_text(">plasmidA\nACGTACGTACGT\n")
    multi = tmp_path / "multi.fa"
    multi.write_text(">plasmidA\nACGTACGTACGT\n>plasmidB\nTTTTGGGGCCCC\n")

    for source, out in ((single, "one"), (multi, "many")):
        result = CliRunner().invoke(
            app, ["batch", "-i", str(source), "-o", str(tmp_path / out)]
        )
        assert result.exit_code == 0, result.stdout

    alone = SeqIO.read(tmp_path / "one" / "single_pLann.gbk", "genbank")
    shared = SeqIO.read(tmp_path / "many" / "plasmidA_pLann.gbk", "genbank")
    assert alone.name == shared.name == "plasmidA"


def test_batch_fast_rejects_custom_yaml(monkeypatch, tmp_path, caplog):
    monkeypatch.setattr(_package_data, "databases_exist", lambda: True)

    fasta = tmp_path / "in.fa"
    fasta.write_text(">plasmid\nACGTACGTACGT\n")
    custom_yaml = tmp_path / "custom.yml"
    custom_yaml.write_text("{}\n")

    result = CliRunner().invoke(
        app,
        [
            "batch",
            "-i",
            str(fasta),
            "-o",
            str(tmp_path / "out"),
            "--fast",
            "--yaml-file",
            str(custom_yaml),
        ],
    )

    assert result.exit_code == 1
    assert "--fast cannot be combined with --yaml-file" in caplog.text


def test_batch_fast_without_yaml_is_allowed(monkeypatch, tmp_path):
    monkeypatch.setattr(_package_data, "databases_exist", lambda: True)
    monkeypatch.setattr(
        main_module.Construct,
        "annotate_batch",
        staticmethod(_fake_batch_constructs),
    )

    fasta = tmp_path / "multi.fa"
    fasta.write_text(">plasmidA\nACGTACGTACGT\n>plasmidB\nTTTTGGGGCCCC\n")
    output = tmp_path / "out"

    result = CliRunner().invoke(
        app, ["batch", "-i", str(fasta), "-o", str(output), "--fast"]
    )

    assert result.exit_code == 0, result.stdout
    assert (output / "plasmidA_pLann.gbk").exists()


def test_batch_multi_record_deduplicates_colliding_ids(monkeypatch, tmp_path):
    monkeypatch.setattr(_package_data, "databases_exist", lambda: True)
    monkeypatch.setattr(
        main_module.Construct,
        "annotate_batch",
        staticmethod(_fake_batch_constructs),
    )

    fasta = tmp_path / "dup.fa"
    fasta.write_text(">dup\nACGTACGT\n>dup\nTTTTGGGG\n")
    output = tmp_path / "out"

    result = CliRunner().invoke(app, ["batch", "-i", str(fasta), "-o", str(output)])

    assert result.exit_code == 0, result.stdout
    # second record with the same id gets a numeric suffix so it does not overwrite
    assert (output / "dup_pLann.gbk").exists()
    assert (output / "dup_2_pLann.gbk").exists()
