import os
import os.path as op
import tempfile
import warnings
from io import StringIO
from pathlib import Path

import pandas as pd
import pytest
from Bio import BiopythonParserWarning, SeqIO
from Bio.SeqRecord import SeqRecord
from click.testing import CliRunner

from plannotate import resources

with open("./tests/test_data/RRNB_fragment.txt") as f:
    RRNB = f.read()

SAM_RIBOSWITCH = (
    "TTTCTATCCAGAGAGGTGGAGGGACTGGCCCTATGAAACCTCGGCAACATTATTGTGCCAATTCCAGCAA"
    "GCGCTAGCTTGAAAGATAGGAA"
)
SAM_RIBOSWITCH_5P_PAD_1 = "A" + SAM_RIBOSWITCH
SAM_RIBOSWITCH_3P_PAD_1 = SAM_RIBOSWITCH + "A"
ISSUE_60_ORIGINAL_SEQUENCE = (
    "cccccacctggcgacaggtgcctctgcggccaaaagccacgtgtataagatacacctgcaaaggcggcaca"
    "accccagtgccacgttgtgagttggatagttgtggaaagagtcaaatggctctcctcaagcgtattcaaca"
    "aggggctgaaggatgcccagaaggtacc"
)


def get_hits_by_db(sequence, db, linear=True):
    from plannotate import annotate

    hits = annotate.annotate(sequence, linear=linear)
    return hits.loc[hits["db"] == db].reset_index(drop=True)


def assert_gbk_round_trips_without_parser_warning(
    gbk_text, expected_start, expected_end
):
    assert "0.." not in gbk_text

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        gbk = SeqIO.read(StringIO(gbk_text), "genbank")

    assert not any(issubclass(w.category, BiopythonParserWarning) for w in caught)
    assert gbk.features[0].location is not None
    assert int(gbk.features[0].location.start) == expected_start
    assert int(gbk.features[0].location.end) == expected_end


def test_yaml():
    db_meta = resources.get_yaml(resources.get_yaml_path())
    snapgene_db = db_meta["snapgene"]
    assert isinstance(snapgene_db, dict)


@pytest.mark.integration
def test_BLAST():
    from plannotate import annotate

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


@pytest.mark.integration
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
@pytest.mark.integration
def test_streamlit_app():
    """this component is hard to test"""
    from plannotate import streamlit_app

    streamlit_app.run_streamlit(["--yaml-file", resources.get_yaml_path()])


# # runs indefinitely
# def test_streamlit():
#     runner = CliRunner()
#     result = runner.invoke(main_streamlit)
#     assert result.exit_code == 0


@pytest.mark.integration
def test_batch():
    from plannotate.pLannotate import main_batch

    runner = CliRunner()
    result = runner.invoke(main_batch)
    assert result.exit_code == 2  # proper exit code for no args


@pytest.mark.integration
def test_batch_help():
    from plannotate.pLannotate import main_batch

    runner = CliRunner()
    result = runner.invoke(main_batch, ["--help"])
    assert result.exit_code == 0


@pytest.mark.integration
def test_annotate():
    from plannotate import annotate

    hits = annotate.annotate(RRNB)
    assert hits.iloc[0]["sseqid"] == "rrnB_T1_terminator"
    hits = annotate.annotate(RRNB, linear=True)
    assert hits.iloc[0]["sseqid"] == "rrnB_T1_terminator"


@pytest.mark.parametrize("linear", [True, False])
@pytest.mark.integration
def test_issue_60_rfam_left_edge_hit_round_trips_after_fix(linear):
    rfam_hits = get_hits_by_db(SAM_RIBOSWITCH, "Rfam", linear=linear)

    assert len(rfam_hits) == 1
    assert rfam_hits.loc[0, "qstart"] == 0
    assert rfam_hits.loc[0, "qend"] == 92
    assert rfam_hits.loc[0, "qseq"] == SAM_RIBOSWITCH

    gbk_text = resources.get_gbk(rfam_hits, SAM_RIBOSWITCH, is_linear=linear)
    assert_gbk_round_trips_without_parser_warning(gbk_text, 0, 92)


@pytest.mark.parametrize(
    ("sequence", "expected_start", "expected_end"),
    [
        (SAM_RIBOSWITCH_3P_PAD_1, 0, 92),
        (SAM_RIBOSWITCH_5P_PAD_1, 1, 93),
    ],
)
@pytest.mark.integration
def test_issue_60_rfam_padding_behaves_by_padding_direction(
    sequence, expected_start, expected_end
):
    rfam_hits = get_hits_by_db(sequence, "Rfam")

    assert len(rfam_hits) == 1
    assert rfam_hits.loc[0, "qstart"] == expected_start
    assert rfam_hits.loc[0, "qend"] == expected_end

    gbk_text = resources.get_gbk(rfam_hits, sequence, is_linear=True)
    assert_gbk_round_trips_without_parser_warning(
        gbk_text, expected_start, expected_end
    )


@pytest.mark.integration
def test_issue_60_original_sequence_snapgene_hit_is_zero_based_in_df_but_valid_in_gbk():
    snapgene_hits = get_hits_by_db(ISSUE_60_ORIGINAL_SEQUENCE, "snapgene")

    assert len(snapgene_hits) == 1
    assert snapgene_hits.loc[0, "qstart"] == 0
    assert snapgene_hits.loc[0, "qend"] == len(ISSUE_60_ORIGINAL_SEQUENCE)

    gbk_text = resources.get_gbk(
        snapgene_hits, ISSUE_60_ORIGINAL_SEQUENCE, is_linear=True
    )
    assert "1..170" in gbk_text

    assert_gbk_round_trips_without_parser_warning(
        gbk_text, 0, len(ISSUE_60_ORIGINAL_SEQUENCE)
    )


####
def test_get_bokeh():
    from plannotate import bokeh_plot

    df_path = op.join(__package__, "test_data", "pXampl3.csv")
    df = pd.read_csv(df_path)
    bokeh_plot.get_bokeh(df)


@pytest.mark.integration
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


@pytest.mark.integration
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


@pytest.mark.integration
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


@pytest.mark.integration
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


@pytest.mark.integration
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


@pytest.mark.integration
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


@pytest.mark.parametrize("ext", ["fasta", "fa", "fas", "fna"])
def test_validate_file_all_fasta_extensions(ext):
    input_file = f"tests/test_data/pAdDeltaF6.{ext}"
    name, ext = resources.get_name_ext(input_file)
    sequence = resources.validate_file(input_file, ext)
    assert len(sequence) == 15420


def test_validate_file_bad_extension():
    input_file = "tests/test_data/pAdDeltaF6.txt"
    name, ext = resources.get_name_ext(input_file)
    with pytest.raises(ValueError, match="must be a FASTA or GenBank file"):
        _ = resources.validate_file(input_file, ext)


@pytest.mark.integration
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
