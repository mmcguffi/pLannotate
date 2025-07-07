import shlex
import subprocess
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Dict, List, Union

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from . import resources as rsc
from .filter_annotations import DF_COLS, clean, get_raw_hits
from .infernal import parse_infernal
from .logging_config import get_logger

logger = get_logger(__name__)


# Type definitions for database configuration
DatabaseConfig = Dict[
    str,
    Union[str, int, List[str], Dict[str, Union[str, bool, List[str]]]],
]


def blast(seq: str, db: DatabaseConfig, mode: str = "blastn") -> pd.DataFrame:
    """Run BLASTn search against a nucleotide database."""
    parameters = db["parameters"]
    db_loc = db["db_loc"]
    query = NamedTemporaryFile()
    tmp = NamedTemporaryFile()
    SeqIO.write(SeqRecord(Seq(seq), id="temp"), query.name, "fasta")

    FLAGS = "qstart qend sseqid sframe pident slen qseq length sstart send qlen evalue"
    cmd = (
        f"{mode} -task blastn-short -query {query.name} -out {tmp.name} "
        f'-db {db_loc} {parameters} -outfmt "6 {FLAGS}"'
    )

    _ = subprocess.run(
        shlex.split(cmd),
        shell=False,
        capture_output=True,
        text=True,
    )

    with open(tmp.name, "r") as file_handle:
        align = file_handle.readlines()

    tmp.close()
    query.close()

    inDf = pd.DataFrame([ele.split() for ele in align], columns=FLAGS.split())
    # Convert numeric columns, handling non-numeric values gracefully
    for col in inDf.columns:
        try:
            inDf[col] = pd.to_numeric(inDf[col])
        except (ValueError, TypeError):
            # Keep non-numeric columns as is
            pass

    return inDf


def diamond(seq: str, db: DatabaseConfig) -> pd.DataFrame:
    """Run DIAMOND blastx search against a protein database."""
    parameters = db["parameters"]
    db_loc = db["db_loc"]
    query = NamedTemporaryFile()
    tmp = NamedTemporaryFile()
    SeqIO.write(SeqRecord(Seq(seq), id="temp"), query.name, "fasta")

    FLAGS = "qstart qend sseqid pident slen qseq length sstart send qlen evalue"
    cmd = (
        f"diamond blastx -d {db_loc} -q {query.name} -o {tmp.name} "
        f"{parameters} --outfmt 6 {FLAGS}"
    )
    _ = subprocess.run(
        shlex.split(cmd),
        shell=False,
        capture_output=True,
        text=True,
    )

    with open(tmp.name, "r") as file_handle:
        align = file_handle.readlines()

    tmp.close()
    query.close()

    inDf = pd.DataFrame([ele.split() for ele in align], columns=FLAGS.split())
    # Convert numeric columns, handling non-numeric values gracefully
    for col in inDf.columns:
        try:
            inDf[col] = pd.to_numeric(inDf[col])
        except (ValueError, TypeError):
            # Keep non-numeric columns as is
            pass

    # DIAMOND-specific post-processing
    # Only apply pipe-splitting for databases that use pipe-separated identifiers (like SwissProt)
    if "sseqid" in inDf.columns and not inDf["sseqid"].empty:
        # Check if any sequence IDs contain pipes before attempting to split
        sseqid_str = inDf["sseqid"].astype(str)
        if sseqid_str.str.contains("\\|", regex=True).any():
            inDf["sseqid"] = sseqid_str.str.split("|", n=2).str.get(1)
    inDf["sframe"] = (inDf["qstart"] < inDf["qend"]).astype(int).replace(0, -1)
    inDf["slen"] = inDf["slen"] * 3
    inDf["length"] = abs(inDf["qend"] - inDf["qstart"]) + 1

    return inDf


def infernal(seq: str, db: DatabaseConfig) -> pd.DataFrame:
    """Run Infernal cmscan search against RNA covariance models."""
    parameters = db["parameters"]
    db_loc = db["db_loc"]
    query = NamedTemporaryFile()
    tmp = NamedTemporaryFile()
    SeqIO.write(SeqRecord(Seq(seq), id="temp"), query.name, "fasta")

    FLAGS = "--cut_ga --rfam --noali --nohmmonly --fmt 2"
    cmd = f"cmscan {FLAGS} {parameters} --tblout {tmp.name} --clanin {db_loc} {query.name}"
    _ = subprocess.run(
        shlex.split(cmd),
        shell=False,
        capture_output=True,
        text=True,
    )
    inDf = parse_infernal(tmp.name)

    inDf["qlen"] = len(seq)

    # manually gets DNA sequence from seq(x2)
    if not inDf.empty:
        inDf["qseq"] = inDf.apply(
            lambda x: (seq)[x["qstart"] : x["qend"] + 1].upper(), axis=1
        )

    tmp.close()
    query.close()

    return inDf


def BLAST(seq: str, db: DatabaseConfig) -> pd.DataFrame:
    """Run search against a database using the appropriate method."""
    task = db["method"]

    if task == "blastn":
        return blast(seq, db)
    elif task == "diamond":
        return diamond(seq, db)
    elif task == "infernal":
        return infernal(seq, db)
    else:
        raise ValueError(f"Unknown search method: {task}")


def _is_fragment(feature: pd.Series) -> bool:
    """Determine if a feature is a fragment based on type and match quality."""
    if feature["Type"] == "CDS":
        if feature["pi_permatch"] == 100:
            return False
        elif ((feature["length"] % 3) == 0) & (feature["percmatch"] > 95):
            return False
        else:
            return True
    elif feature["Type"] != "CDS":
        if feature["percmatch"] < 95:
            return True
        else:
            return False
    else:
        logger.error("Fragment error.")
        return False


def annotate(
    seq: str | Seq,
    yaml_file: Path = rsc.get_yaml_path(),
    linear: bool = False,
    is_detailed: bool = False,
) -> pd.DataFrame:
    """Annotate a DNA sequence and return results as DataFrame."""
    # This catches errors in sequence via Biopython
    fileloc = NamedTemporaryFile()
    SeqIO.write(
        SeqRecord(Seq(seq), name="pLannotate", annotations={"molecule_type": "DNA"}),
        fileloc.name,
        "fasta",
    )
    record = list(SeqIO.parse(fileloc.name, "fasta"))
    fileloc.close()

    record = record[0]

    # doubles sequence for origin crossing hits
    if linear is False:
        query = str(record.seq) + str(record.seq)
    elif linear is True:
        query = str(record.seq)

    blastDf = get_raw_hits(query, linear, yaml_file)

    if blastDf.empty:  # if no hits are found
        blastDf = pd.DataFrame(columns=DF_COLS)
        return blastDf

    # this has to re-parse the yaml, so not an elegant solution
    if is_detailed is True:
        blastDf["kind"] = blastDf["Type"]
    else:
        blastDf["kind"] = 1

    blastDf = clean(blastDf)

    if blastDf.empty:  # if no hits are found
        blastDf = pd.DataFrame(columns=DF_COLS)
        return blastDf

    blastDf["fragment"] = blastDf.apply(_is_fragment, axis=1)

    if blastDf.empty:  # if no hits are found
        blastDf = pd.DataFrame(columns=DF_COLS)
        return blastDf

    blastDf["qend"] = blastDf["qend"] + 1  # corrects position for gbk

    # manually gets DNA sequence from inSeq
    # blastDf['qseq'] = inSeq #adds the sequence to the df
    # blastDf['qseq'] = blastDf.apply(lambda x: x['qseq'][x['qstart']:x['qend']+1], axis=1)
    blastDf["qseq"] = blastDf.apply(
        lambda x: (
            str(Seq(x["qseq"]).reverse_complement()) if x["sframe"] == -1 else x["qseq"]
        ),
        axis=1,
    )

    # fill in edge cases (kludge)
    blastDf["Feature"] = blastDf["Feature"].fillna(blastDf["sseqid"])
    blastDf["Description"] = blastDf["Description"].fillna("")
    blastDf["Type"] = blastDf["Type"].fillna("misc_feature")

    return blastDf
