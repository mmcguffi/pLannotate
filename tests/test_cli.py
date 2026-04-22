import tempfile
from pathlib import Path

import pytest
from Bio import SeqIO
from click.testing import CliRunner

pytestmark = pytest.mark.integration


def test_batch():
    from plannotate.pLannotate import main_batch

    runner = CliRunner()
    result = runner.invoke(main_batch)
    assert result.exit_code == 2  # proper exit code for no args


def test_batch_help():
    from plannotate.pLannotate import main_batch

    runner = CliRunner()
    result = runner.invoke(main_batch, ["--help"])
    assert result.exit_code == 0


def test_cli_annotate():
    from plannotate.pLannotate import main_batch

    plasmid = Path("pXampl3.fa")
    with tempfile.TemporaryDirectory() as tmpdir:
        runner = CliRunner()
        result = runner.invoke(
            main_batch,
            [
                "-i",
                f"tests/test_data/{plasmid}",
                "-o",
                tmpdir,
                "-s",
                "",
            ],
        )
        assert result.exit_code == 0
        gbk = SeqIO.read(tmpdir / plasmid.with_suffix(".gbk"), "genbank")
    assert len(gbk.features) > 15


def test_cli_annotate_empty_gbk():
    from plannotate.pLannotate import main_batch

    plasmid = Path("random_dna.fa")
    with tempfile.TemporaryDirectory() as tmpdir:
        runner = CliRunner()
        result = runner.invoke(
            main_batch,
            [
                "-i",
                f"tests/test_data/{plasmid}",
                "-o",
                tmpdir,
                "-s",
                "",
            ],
        )
        assert result.exit_code == 0
        gbk = SeqIO.read(tmpdir / plasmid.with_suffix(".gbk"), "genbank")
    assert len(gbk.features) == 0


def test_cli_annotate_empty_html():
    from plannotate.pLannotate import main_batch

    plasmid = Path("random_dna.fa")
    with tempfile.TemporaryDirectory() as tmpdir:
        runner = CliRunner()
        result = runner.invoke(
            main_batch,
            [
                "-i",
                f"tests/test_data/{plasmid}",
                "-o",
                tmpdir,
                "-s",
                "",
                "-h",
                "-x",
            ],
        )
        assert result.exit_code == 0
        html = tmpdir / plasmid.with_suffix(".html")
        assert html.exists()


def test_cli_save_nan_feature():
    from plannotate.pLannotate import main_batch

    plasmid = Path("nan_feature.fa")
    with tempfile.TemporaryDirectory() as tmpdir:
        runner = CliRunner()
        result = runner.invoke(
            main_batch,
            [
                "-i",
                f"tests/test_data/{plasmid}",
                "-o",
                tmpdir,
                "-s",
                "",
            ],
        )
        assert result.exit_code == 0
        gbk = SeqIO.read(tmpdir / plasmid.with_suffix(".gbk"), "genbank")
    assert len(gbk.features) == 2


def test_bokeh_bakein():
    from plannotate.pLannotate import main_batch

    plasmid = Path("pXampl3.fa")
    with tempfile.TemporaryDirectory() as tmpdir:
        runner = CliRunner()
        cdn_result = runner.invoke(
            main_batch,
            [
                "-i",
                f"tests/test_data/{plasmid}",
                "-o",
                tmpdir,
                "-s",
                ".cdn",
                "-h",
                "-x",
            ],
        )
        inline_result = runner.invoke(
            main_batch,
            [
                "-i",
                f"tests/test_data/{plasmid}",
                "-o",
                tmpdir,
                "-s",
                ".inline",
                "-hf",
                "-x",
            ],
        )
        assert cdn_result.exit_code == 0
        assert inline_result.exit_code == 0
        inline = tmpdir / plasmid.with_suffix(".inline.html")
        cdn = tmpdir / plasmid.with_suffix(".cdn.html")

        assert inline.exists()
        assert cdn.exists()

        assert inline.stat().st_size > cdn.stat().st_size


def test_zero_feature():
    from plannotate.pLannotate import main_batch

    plasmid = Path("nan_feature.fa")
    with tempfile.TemporaryDirectory() as tmpdir:
        runner = CliRunner()
        result = runner.invoke(
            main_batch,
            [
                "-i",
                f"tests/test_data/{plasmid}",
                "-o",
                tmpdir,
                "-s",
                "",
            ],
        )
        assert result.exit_code == 0
        gbk = SeqIO.read(tmpdir / plasmid.with_suffix(".gbk"), "genbank")
    assert len(gbk.features) == 2


def test_annotate_fna(tmp_path):
    from plannotate.pLannotate import main_batch

    input_file = "tests/test_data/pAdDeltaF6.fna"
    arglist = [
        "-i",
        input_file,
        "--output",
        tmp_path,
        "--html",
        "--csv",
        "-f",
        "pAdDeltaF6",
    ]
    result = CliRunner().invoke(main_batch, arglist)
    assert result.exit_code == 0
    gbk = SeqIO.read(tmp_path / "pAdDeltaF6_pLann.gbk", "genbank")
    assert len(gbk.features) == 29
