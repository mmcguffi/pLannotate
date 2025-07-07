from pathlib import Path
from tempfile import NamedTemporaryFile

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from . import resources as rsc
from .filter_annotations import DF_COLS, filter_and_clean_hits
from .search import search_all_databases
from .logging_config import get_logger

logger = get_logger(__name__)


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

    blastDf = search_all_databases(query, linear, yaml_file)

    if blastDf.empty:  # if no hits are found
        blastDf = pd.DataFrame(columns=DF_COLS)
        return blastDf

    # this has to re-parse the yaml, so not an elegant solution
    if is_detailed is True:
        blastDf["kind"] = blastDf["Type"]
    else:
        blastDf["kind"] = 1

    blastDf = filter_and_clean_hits(blastDf)

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
