from pathlib import Path
from tempfile import NamedTemporaryFile

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from . import resources as rsc
from .filter_annotations import filter_and_clean_hits
from .logging_config import get_logger
from .search import DF_COLS, search_all_databases

logger = get_logger(__name__)


def _is_fragment(feature: pd.Series) -> bool:
    """Determine if a feature is a fragment based on type and match quality."""
    # If type column doesn't exist, assume it's not a fragment
    if "type" not in feature.index:
        return False

    if feature["type"] == "CDS":
        if feature["pi_permatch"] == 100:
            return False
        elif ((feature["length"] % 3) == 0) & (feature["percmatch"] > 95):
            return False
        else:
            return True
    elif feature["type"] != "CDS":
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
    yaml_file = Path(yaml_file)
    seq = Seq(seq)

    # This catches errors in sequence via Biopython
    fileloc = NamedTemporaryFile()
    SeqIO.write(
        SeqRecord(
            seq,
            name="pLannotate",
            annotations={"molecule_type": "DNA"},
        ),
        fileloc.name,
        "fasta",
    )
    record = list(SeqIO.parse(fileloc.name, "fasta"))
    fileloc.close()

    record = record[0]

    blastDf = search_all_databases(str(seq), linear, yaml_file)
    logger.info(f"Processing {len(blastDf)} raw database hits")

    if blastDf.empty:  # if no hits are found
        blastDf = pd.DataFrame(columns=DF_COLS)
        logger.info("No hits found, returning empty DataFrame")
        return blastDf

    # this has to re-parse the yaml, so not an elegant solution
    logger.debug(
        f"Setting feature kinds for {'detailed' if is_detailed else 'standard'} mode"
    )
    if is_detailed is True:
        blastDf["kind"] = blastDf["type"]
    else:
        blastDf["kind"] = 1

    logger.debug("Starting hit filtering and cleaning...")
    blastDf = filter_and_clean_hits(blastDf, linear)
    logger.info(f"Filtered to {len(blastDf)} high-quality features")

    if blastDf.empty:  # if no hits are found
        blastDf = pd.DataFrame(columns=DF_COLS)
        logger.info("No hits remaining after filtering")
        return blastDf

    logger.debug("Calculating fragment status for features...")
    blastDf["fragment"] = blastDf.apply(_is_fragment, axis=1)
    fragment_count = blastDf["fragment"].sum()
    logger.debug(
        f"Identified {fragment_count} fragments out of {len(blastDf)} features"
    )

    if blastDf.empty:  # if no hits are found
        blastDf = pd.DataFrame(columns=DF_COLS)
        return blastDf

    logger.debug("Adjusting coordinates for GenBank format...")
    blastDf["qend"] = blastDf["qend"] + 1  # corrects position for gbk

    logger.debug("Processing reverse complement sequences...")
    # manually gets DNA sequence from inSeq
    # blastDf['qseq'] = inSeq #adds the sequence to the df
    # blastDf['qseq'] = blastDf.apply(lambda x: x['qseq'][x['qstart']:x['qend']+1], axis=1)
    blastDf["qseq"] = blastDf.apply(
        lambda x: (
            str(Seq(x["qseq"]).reverse_complement()) if x["sframe"] == -1 else x["qseq"]
        ),
        axis=1,
    )

    logger.debug("Filling missing feature information...")
    # fill in edge cases (kludge)
    blastDf["name"] = blastDf["name"].fillna(blastDf["sseqid"])
    blastDf["blurb"] = blastDf["blurb"].fillna("")
    if "type" in blastDf.columns:
        blastDf["type"] = blastDf["type"].fillna("misc_feature")
    else:
        blastDf["type"] = "misc_feature"

    logger.info(f"Annotation complete: {len(blastDf)} features identified")
    return blastDf
