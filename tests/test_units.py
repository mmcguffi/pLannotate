import os
import os.path as op
import tempfile
from io import StringIO
from pathlib import Path

import pandas as pd
import pytest
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from click.testing import CliRunner

from plannotate import annotate, bokeh_plot, resources, streamlit_app
from plannotate.pLannotate import main_batch

with open("./tests/test_data/RRNB_fragment.txt") as f:
    RRNB = f.read()


def test_yaml():
    db_meta = resources.get_yaml(resources.get_yaml_path())
    snapgene_db = db_meta["snapgene"]
    assert isinstance(snapgene_db, dict)


def test_BLAST():
    db_meta = resources.get_yaml(resources.get_yaml_path())
    snapgene_db = db_meta["snapgene"]
    hits = annotate.BLAST(RRNB, snapgene_db)
    assert not hits.empty


####
def test_get_image():
    name = "icon.png"
    path = ("plannotate", "data", "images", name)
    assert resources.get_image(name) == op.join(
        op.dirname(op.abspath(__package__)), *path
    )


def test_get_template():
    name = "blurb.html"
    path = ("plannotate", "data", "templates", name)
    assert resources.get_template(name) == op.join(
        op.dirname(op.abspath(__package__)), *path
    )


def test_get_fasta():
    # TODO: this is a hack -- fix it
    current_dir = Path((os.path.abspath(__file__))).parent.parent

    fasta_loc = str(current_dir / "plannotate" / "data" / "fastas")
    assert resources.get_example_fastas() == fasta_loc


def test_get_yaml_path():
    path = ("plannotate", "data", "data", "databases.yml")
    assert resources.get_yaml_path() == op.join(
        op.dirname(op.abspath(__package__)), *path
    )


def test_get_yaml():
    yaml = resources.get_yaml(resources.get_yaml_path())

    assert isinstance(yaml, dict)
    assert len(yaml) > 0

    first_key = list(yaml.keys())[0]
    expected_fields = set(
        ("version", "method", "location", "priority", "parameters", "details", "db_loc")
    )
    assert set(yaml[first_key].keys()) == expected_fields


def test_databases_exist():
    assert resources.databases_exist() is True


def test_valid_sequence_correct():
    resources.validate_sequence("ACTG")


def test_valid_sequence_incorrect_base():
    with pytest.raises(ValueError):
        resources.validate_sequence("ACTX")


def test_get_name_ext():
    name, ext = resources.get_name_ext("./a/long/path/test.fasta")
    assert name == "test"
    assert ext == ".fasta"


def test_seq_record():

    NUM_FEATURES = 21
    PLAS_LEN = 6015

    df_path = op.join(__package__, "test_data", "pXampl3.csv")
    df = pd.read_csv(df_path)
    assert len(df) == NUM_FEATURES

    seq_path = op.join(__package__, "test_data", "pXampl3.fa")
    seq = str(SeqIO.read(seq_path, "fasta").seq)
    gbk = resources.get_seq_record(df, seq)

    assert len(gbk.features) == NUM_FEATURES
    assert len(gbk.seq) == PLAS_LEN
    assert "pLannotate" in gbk.annotations["comment"]
    assert gbk.annotations["topology"] == "circular"


def test_get_gbk():
    df_path = op.join(__package__, "test_data", "pXampl3.csv")
    df = pd.read_csv(df_path)

    seq_path = op.join(__package__, "test_data", "pXampl3.fa")
    seq = str(SeqIO.read(seq_path, "fasta").seq)

    gbk_text = resources.get_gbk(df, seq)
    gbk = SeqIO.read(StringIO(gbk_text), "genbank")

    assert type(gbk) is SeqRecord
    assert len(gbk.seq) > 0


def test_get_clean_csv_df():
    df_path = op.join(__package__, "test_data", "pXampl3.csv")
    df = pd.read_csv(df_path)

    df_clean = resources.get_clean_csv_df(df)

    assert len(df.columns) == 28
    assert len(df_clean.columns) == 14


def test_validate_file_fa():
    df_path = op.join(__package__, "test_data", "pXampl3.fa")
    resources.validate_file(df_path, ".fasta")


def test_validate_file_gbk():
    df_path = op.join(__package__, "test_data", "pXampl3_detailed.gbk")
    resources.validate_file(df_path, ".gbk")


####
def test_streamlit_app():
    """this component is hard to test"""
    streamlit_app.run_streamlit(["--yaml-file", resources.get_yaml_path()])


# # runs indefinitely
# def test_streamlit():
#     runner = CliRunner()
#     result = runner.invoke(main_streamlit)
#     assert result.exit_code == 0


def test_batch():
    runner = CliRunner()
    result = runner.invoke(main_batch)
    assert result.exit_code == 2  # proper exit code for no args


def test_batch_help():
    runner = CliRunner()
    result = runner.invoke(main_batch, ["--help"])
    assert result.exit_code == 0


def test_annotate():
    hits = annotate.annotate(RRNB)
    assert hits.iloc[0]["sseqid"] == "rrnB_T1_terminator"
    hits = annotate.annotate(RRNB, linear=True)
    assert hits.iloc[0]["sseqid"] == "rrnB_T1_terminator"


####
def test_get_bokeh():
    df_path = op.join(__package__, "test_data", "pXampl3.csv")
    df = pd.read_csv(df_path)
    bokeh_plot.get_bokeh(df)


def test_cli_annotate():

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


def test_cli_annotate_empty():

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
