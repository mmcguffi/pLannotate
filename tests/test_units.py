import os
import tempfile
from io import StringIO
from pathlib import Path

import pandas as pd
import pytest
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from typer.testing import CliRunner

from plannotate import annotate, bokeh_plot, resources
from plannotate.search import blast
from plannotate.main import app
from plannotate.validation import get_name_ext, validate_file, validate_sequence

with open("./tests/test_data/RRNB_fragment.txt") as f:
    RRNB = f.read()


def test_yaml():
    db_meta = resources.get_yaml(resources.get_yaml_path())
    snapgene_db = db_meta["snapgene"]
    assert isinstance(snapgene_db, dict)


def test_BLAST():
    db_meta = resources.get_yaml(resources.get_yaml_path())
    snapgene_db = db_meta["snapgene"]
    hits = blast(RRNB, snapgene_db)
    assert not hits.empty


def test_BLAST_sseqid_split_no_pipe():
    """Test sseqid split when no '|' character is present."""
    df = pd.DataFrame(
        {
            "qstart": [1],
            "qend": [2],
            "sseqid": ["no_pipe"],
            "pident": [99.0],
            "slen": [100],
            "qseq": ["AT"],
            "length": [2],
            "sstart": [1],
            "send": [2],
            "qlen": [100],
            "evalue": [1e-10],
        }
    )

    # Should not raise an exception and should return nan for missing split
    df2 = df.copy()
    df2["sseqid"] = df2["sseqid"].astype(str).str.split("|", n=2).str.get(1)
    assert pd.isna(df2["sseqid"].iloc[0])

    # Note: This tests the edge case where sseqid doesn't contain '|' separators.
    # In diamond BLAST output, sseqid is typically formatted as "db|id|description",
    # but some databases may not follow this format. When splitting by '|' and
    # requesting index 1, pandas returns nan if the split produces fewer than 2 elements.


def test_BLAST_sseqid_split_empty_dataframe():
    """Test sseqid split on empty DataFrame."""
    df_empty = pd.DataFrame(
        columns=[
            "qstart",
            "qend",
            "sseqid",
            "pident",
            "slen",
            "qseq",
            "length",
            "sstart",
            "send",
            "qlen",
            "evalue",
        ]
    )

    # Should not raise an exception and should return empty Series
    df_empty["sseqid"] = df_empty["sseqid"].astype(str).str.split("|", n=2).str.get(1)
    assert df_empty["sseqid"].empty


def test_get_image():
    name = "icon.png"
    project_root = Path(__file__).parent.parent
    expected_path = project_root / "plannotate" / "data" / "images" / name
    assert resources.get_image(name) == expected_path


def test_get_template():
    name = "blurb.html"
    project_root = Path(__file__).parent.parent
    expected_path = project_root / "plannotate" / "data" / "templates" / name
    assert resources.get_template(name) == expected_path


def test_get_fasta():
    # TODO: this is a hack -- fix it
    current_dir = Path((os.path.abspath(__file__))).parent.parent

    fasta_loc = current_dir / "plannotate" / "data" / "fastas"
    assert resources.get_example_fastas() == fasta_loc


def test_get_yaml_path():
    project_root = Path(__file__).parent.parent
    expected_path = project_root / "plannotate" / "data" / "data" / "databases.yml"
    assert resources.get_yaml_path() == expected_path


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
    validate_sequence("ACTG")


def test_valid_sequence_incorrect_base():
    with pytest.raises(ValueError):
        validate_sequence("ACTX")


def test_get_name_ext():
    name, ext = get_name_ext("./a/long/path/test.fasta")
    assert name == "test"
    assert ext == ".fasta"


def test_seq_record():
    from plannotate.models import Construct, df_to_features

    NUM_FEATURES = 21
    PLAS_LEN = 6015

    df_path = Path(__file__).parent / "test_data" / "pXampl3.csv"
    df = pd.read_csv(df_path)
    assert len(df) == NUM_FEATURES

    seq_path = Path(__file__).parent / "test_data" / "pXampl3.fa"
    seq = str(SeqIO.read(seq_path, "fasta").seq)

    # Create a Construct object with the test data
    construct = Construct(seq=seq, linear=False)
    # Replace the features with our test data
    construct.features = df_to_features(df)

    gbk = construct.to_seqrecord()

    assert len(gbk.features) == NUM_FEATURES
    assert len(gbk.seq) == PLAS_LEN
    assert "pLannotate" in str(gbk.annotations["comment"])
    assert gbk.annotations["topology"] == "circular"


def test_get_gbk():
    from plannotate.models import Construct, df_to_features

    df_path = Path(__file__).parent / "test_data" / "pXampl3.csv"
    df = pd.read_csv(df_path)

    seq_path = Path(__file__).parent / "test_data" / "pXampl3.fa"
    seq = str(SeqIO.read(seq_path, "fasta").seq)

    # Create a Construct object with the test data
    construct = Construct(seq=seq, linear=False)
    # Replace the features with our test data
    construct.features = df_to_features(df)

    gbk_text = construct.to_genbank()
    gbk = SeqIO.read(StringIO(gbk_text), "genbank")

    assert isinstance(gbk, SeqRecord)
    assert len(gbk.seq) > 0


def test_get_clean_csv_df():
    from plannotate.models import Construct, df_to_features

    df_path = Path(__file__).parent / "test_data" / "pXampl3.csv"
    df = pd.read_csv(df_path)

    # Create a Construct object with the test data
    seq = "ATCG" * 1500  # 6000 bp sequence
    construct = Construct(seq=seq, linear=False)
    construct.features = df_to_features(df)

    df_clean = construct.to_csv()

    assert len(df.columns) == 28
    assert len(df_clean.columns) == 14


def test_validate_file_fa():
    df_path = Path(__file__).parent / "test_data" / "pXampl3.fa"
    validate_file(str(df_path), ".fasta")


def test_validate_file_gbk():
    df_path = Path(__file__).parent / "test_data" / "pXampl3_detailed.gbk"
    validate_file(str(df_path), ".gbk")


def test_batch():
    runner = CliRunner()
    result = runner.invoke(app, ["batch"])
    assert result.exit_code == 2  # proper exit code for no args


def test_batch_help():
    runner = CliRunner()
    result = runner.invoke(app, ["batch", "--help"])
    assert result.exit_code == 0


def test_annotate():
    hits = annotate.annotate(RRNB)
    assert hits.iloc[0]["sseqid"] == "rrnB_T1_terminator"
    hits = annotate.annotate(RRNB, linear=True)
    assert hits.iloc[0]["sseqid"] == "rrnB_T1_terminator"


def test_get_bokeh():
    df_path = Path(__file__).parent / "test_data" / "pXampl3.csv"
    df = pd.read_csv(df_path)
    bokeh_plot.get_bokeh(df)


def test_cli_annotate():
    plasmid = Path("pXampl3.fa")
    with tempfile.TemporaryDirectory() as tmpdir:
        runner = CliRunner()
        result = runner.invoke(
            app,
            [
                "batch",
                "-i",
                f"tests/test_data/{plasmid}",
                "-o",
                str(tmpdir),
                "-s",
                "",
            ],
        )
        assert result.exit_code == 0
        gbk = SeqIO.read(tmpdir / plasmid.with_suffix(".gbk"), "genbank")
    assert len(gbk.features) > 15


def test_cli_annotate_empty_gbk():
    plasmid = Path("random_dna.fa")
    with tempfile.TemporaryDirectory() as tmpdir:
        runner = CliRunner()
        result = runner.invoke(
            app,
            [
                "batch",
                "-i",
                f"tests/test_data/{plasmid}",
                "-o",
                str(tmpdir),
                "-s",
                "",
            ],
        )
        assert result.exit_code == 0
        gbk = SeqIO.read(tmpdir / plasmid.with_suffix(".gbk"), "genbank")
    assert len(gbk.features) == 0


def test_cli_annotate_empty_html():
    plasmid = Path("random_dna.fa")
    with tempfile.TemporaryDirectory() as tmpdir:
        runner = CliRunner()
        result = runner.invoke(
            app,
            [
                "batch",
                "-i",
                f"tests/test_data/{plasmid}",
                "-o",
                str(tmpdir),
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
    plasmid = Path("nan_feature.fa")
    with tempfile.TemporaryDirectory() as tmpdir:
        runner = CliRunner()
        result = runner.invoke(
            app,
            [
                "batch",
                "-i",
                f"tests/test_data/{plasmid}",
                "-o",
                str(tmpdir),
                "-s",
                "",
            ],
        )
        assert result.exit_code == 0
        gbk = SeqIO.read(tmpdir / plasmid.with_suffix(".gbk"), "genbank")
    assert len(gbk.features) == 2


@pytest.mark.skip(reason="slow, unimportant")
def test_bokeh_bakein():
    plasmid = Path("pXampl3.fa")
    with tempfile.TemporaryDirectory() as tmpdir:
        runner = CliRunner()
        cdn_result = runner.invoke(
            app,
            [
                "batch",
                "-i",
                f"tests/test_data/{plasmid}",
                "-o",
                str(tmpdir),
                "-s",
                ".cdn",
                "-h",
                "-x",
            ],
        )
        inline_result = runner.invoke(
            app,
            [
                "batch",
                "-i",
                f"tests/test_data/{plasmid}",
                "-o",
                str(tmpdir),
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
    plasmid = Path("nan_feature.fa")
    with tempfile.TemporaryDirectory() as tmpdir:
        runner = CliRunner()
        result = runner.invoke(
            app,
            [
                "batch",
                "-i",
                f"tests/test_data/{plasmid}",
                "-o",
                str(tmpdir),
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
    name, ext = get_name_ext(input_file)
    sequence = validate_file(input_file, ext)
    assert len(sequence) == 15420


def test_validate_file_bad_extension():
    input_file = "tests/test_data/pAdDeltaF6.txt"
    _, ext = get_name_ext(input_file)
    with pytest.raises(ValueError, match="must be a FASTA or GenBank file"):
        _ = validate_file(input_file, ext)


@pytest.mark.skip(reason="slow, redundant")
def test_annotate_fna(tmp_path):
    input_file = "tests/test_data/pAdDeltaF6.fna"
    arglist = [
        "batch",
        "-i",
        input_file,
        "--output",
        str(tmp_path),
        "--html",
        "--csv",
        "-f",
        "pAdDeltaF6",
    ]
    result = CliRunner().invoke(app, arglist)
    assert result.exit_code == 0
    gbk = SeqIO.read(tmp_path / "pAdDeltaF6_pLann.gbk", "genbank")
    assert len(gbk.features) == 29


def test_dataframe_to_features():
    """Test conversion of DataFrame to Feature objects."""
    import pandas as pd

    from plannotate.models import df_to_features

    # Create a sample DataFrame with the expected columns
    test_data = {
        "sseqid": ["test1", "test2"],
        "qstart": [0, 100],
        "qend": [50, 150],
        "sstart": [0, 0],
        "send": [50, 50],
        "sframe": [1, -1],
        "score": [95.0, 90.0],
        "evalue": [1e-10, 1e-8],
        "qseq": ["ATCG", "GCTA"],
        "length": [50, 50],
        "slen": [50, 50],
        "pident": [95.0, 90.0],
        "qlen": [1000, 1000],
        "db": ["test_db", "test_db"],
        "Feature": ["Test Feature 1", "Test Feature 2"],
        "Description": ["Test description 1", "Test description 2"],
        "Type": ["CDS", "promoter"],
        "priority": [1, 2],
        "percmatch": [95.0, 90.0],
        "abs percmatch": [95.0, 90.0],
        "pi_permatch": [90.25, 81.0],
        "wiggle": [7, 7],
        "wstart": [7, 107],
        "wend": [43, 143],
        "kind": [1, 1],
        "qstart_dup": [0, 100],
        "qend_dup": [50, 150],
        "fragment": [False, True],
    }

    df = pd.DataFrame(test_data)
    features = df_to_features(df)

    assert len(features) == 2
    assert features[0].feature_name == "Test Feature 1"
    assert features[0].feature_type == "CDS"
    assert features[0].qstart == 0
    assert features[0].qend == 50
    assert features[0].is_forward_strand is True

    assert features[1].feature_name == "Test Feature 2"
    assert features[1].feature_type == "promoter"
    assert features[1].qstart == 100
    assert features[1].qend == 150
    assert features[1].is_reverse_strand is True
    assert features[1].fragment is True


def test_dataframe_to_features_empty():
    """Test conversion of empty DataFrame to Feature objects."""
    import pandas as pd

    from plannotate.models import df_to_features

    empty_df = pd.DataFrame()
    features = df_to_features(empty_df)

    assert len(features) == 0


def test_construct_class_methods():
    """Test Construct class methods for to_genbank, to_csv, and plot."""
    from plannotate.models import Construct, df_to_features

    # Create test data
    test_data = {
        "sseqid": ["test1"],
        "qstart": [0],
        "qend": [50],
        "sstart": [0],
        "send": [50],
        "sframe": [1],
        "score": [95.0],
        "evalue": [1e-10],
        "qseq": ["ATCG"],
        "length": [50],
        "slen": [50],
        "pident": [95.0],
        "qlen": [1000],
        "db": ["test_db"],
        "Feature": ["Test Feature"],
        "Description": ["Test description"],
        "Type": ["CDS"],
        "priority": [1],
        "percmatch": [95.0],
        "abs percmatch": [95.0],
        "pi_permatch": [90.25],
        "wiggle": [7],
        "wstart": [7],
        "wend": [43],
        "kind": [1],
        "qstart_dup": [0],
        "qend_dup": [50],
        "fragment": [False],
    }

    df = pd.DataFrame(test_data)
    seq = "ATCG" * 250  # 1000 bp sequence

    # Create Construct object
    construct = Construct(seq=seq, linear=False)
    construct.features = df_to_features(df)

    # Test to_genbank
    gbk_text = construct.to_genbank()
    assert "LOCUS" in gbk_text
    assert "Test Feature" in gbk_text

    # Test to_csv
    csv_df = construct.to_csv()
    assert not csv_df.empty
    assert "start location" in csv_df.columns

    # Test plot
    plot = construct.plot()
    assert plot is not None

    # Test to_html
    html_content = construct.to_html()
    assert "pLannotate" in html_content
    assert "<!DOCTYPE html>" in html_content

    # Test to_seqrecord
    seq_record = construct.to_seqrecord()
    assert len(seq_record.features) == 1
    assert seq_record.features[0].qualifiers["label"] == "Test Feature"


def test_entry_point_installation():
    """Test that the entry point creates a working command line tool."""
    import subprocess
    import sys

    # Test that the plannotate command exists and can be called
    try:
        result = subprocess.run(
            [sys.executable, "-m", "plannotate.main", "--help"],
            capture_output=True,
            text=True,
            timeout=10,
        )
        assert result.returncode == 0
        assert "Usage:" in result.stdout
    except subprocess.TimeoutExpired:
        pytest.fail("Command timed out")
    except FileNotFoundError:
        pytest.fail("plannotate.main module not found")


def test_diamond_annotation_parsing():
    "diamond output is pretty hacky with pipe splitting"
    INPUT_FILE = "tests/test_data/RNAs.fasta"
    fasta = SeqIO.read(INPUT_FILE, "fasta")

    from plannotate.models import Construct

    # Test in detailed mode to capture ncRNA annotations
    plasmid = Construct(fasta.seq, detailed=True).to_seqrecord()

    feats = [_ for _ in plasmid.features]
    serialized = []
    for feat in feats:
        parts = {
            "start": feat.location.start,  # type: ignore
            "end": feat.location.end,  # type: ignore
            "strand": feat.location.strand,  # type: ignore
            "type": feat.type,
            "is_frag": feat.qualifiers["fragment"],
            "database": feat.qualifiers["database"],
            "name": feat.qualifiers["label"],
        }
        serialized.append(parts)

    new_annos = pd.DataFrame(serialized)
    new_annos = new_annos.sort_values(by=["start", "end"]).reset_index(drop=True)
    old_annos = pd.read_csv("tests/test_data/RNAs_ground-truth-detailed.csv")

    # compare the two DataFrames
    pd.testing.assert_frame_equal(new_annos, old_annos, check_dtype=False)
