import os.path as op
from io import StringIO

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from plannotate import resources


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
