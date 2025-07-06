import os
import subprocess
import sys
from importlib.resources import files
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Any, Dict, Tuple

import yaml
from Bio import SeqIO

from . import __version__ as plannotate_version
from .logging_config import get_logger

logger = get_logger(__name__)

PACKAGE = __package__ or "plannotate"
ROOT_DIR = os.path.dirname(os.path.abspath(__file__))

valid_genbank_exts = [".gbk", ".gb", ".gbf", ".gbff"]
valid_fasta_exts = [".fa", ".fasta", ".fas", ".fna"]
MAX_PLAS_SIZE = 50000

DF_COLS = [
    "sseqid",
    "qstart",
    "qend",
    "sstart",
    "send",
    "sframe",
    "score",
    "evalue",
    "qseq",
    "length",
    "slen",
    "pident",
    "qlen",
    "db",
    "Feature",
    "Description",
    "Type",
    "priority",
    "percmatch",
    "abs percmatch",
    "pi_permatch",
    "wiggle",
    "wstart",
    "wend",
    "kind",
    "qstart_dup",
    "qend_dup",
    "fragment",
]


def get_resource(group: str, name: str) -> Path:
    return Path(str(files(PACKAGE) / f"data/{group}/{name}"))


def get_image(name: str) -> Path:
    return get_resource("images", name)


def get_template(name: str) -> Path:
    return get_resource("templates", name)


def get_example_fastas() -> Path:
    return get_resource("fastas", "")


def get_yaml_path() -> Path:
    return get_resource("data", "databases.yml")


def get_details(name: str) -> Path:
    return get_resource("data", name)


def get_name_ext(file_loc: str) -> Tuple[str, str]:
    base = os.path.basename(file_loc)
    name = os.path.splitext(base)[0]
    ext = os.path.splitext(base)[1]
    return name, ext


def validate_file(file: str, ext: str, max_length: int = MAX_PLAS_SIZE) -> str:
    if ext in valid_fasta_exts:
        # This catches errors on file uploads via Biopython
        temp_fileloc = NamedTemporaryFile()
        record = list(SeqIO.parse(file, "fasta"))
        try:
            record[0].annotations["molecule_type"] = "DNA"
        except IndexError:
            error = (
                "Malformed fasta file --> please submit a fasta file in standard format"
            )
            raise ValueError(error)
        SeqIO.write(record, temp_fileloc.name, "fasta")
        record = list(SeqIO.parse(temp_fileloc.name, "fasta"))
        temp_fileloc.close()

        if len(record) != 1:
            error = "FASTA file contains many entries --> please submit a single FASTA file."
            raise ValueError(error)

    elif ext in valid_genbank_exts:
        temp_fileloc = NamedTemporaryFile()
        try:
            record = list(SeqIO.parse(file, "gb"))[0]
        except IndexError:
            error = "Malformed Genbank file --> please submit a Genbank file in standard format"
            raise ValueError(error)
        # submitted_gbk = record # for combining -- not current imlementated
        SeqIO.write(record, temp_fileloc.name, "fasta")
        record = list(SeqIO.parse(temp_fileloc.name, "fasta"))
        temp_fileloc.close()

    else:
        error = "must be a FASTA or GenBank file"
        raise ValueError(error)

    if len(record) != 1:
        error = (
            "FASTA file contains many entries --> please submit a single FASTA file."
        )
        raise ValueError(error)

    inSeq = str(record[0].seq)

    validate_sequence(inSeq, max_length)

    return inSeq


def validate_sequence(inSeq: str, max_length: int = MAX_PLAS_SIZE) -> None:
    IUPAC = "GATCRYWSMKHBVDNgatcrywsmkhbvdn"
    if not set(inSeq).issubset(IUPAC):
        error = "Sequence contains invalid characters -- must be ATCG and/or valid IUPAC nucleotide ambiguity code"
        raise ValueError(error)

    if len(inSeq) > max_length:
        error = f"Are you sure this is an engineered plasmid? Entry size is too large -- must be {max_length} bases or less."
        raise ValueError(error)


def get_yaml(yaml_file_loc: Path) -> Dict[str, Any]:
    # file_name = get_resource("data", "databases.yml")
    with open(yaml_file_loc, "r") as f:
        dbs = yaml.load(f, Loader=yaml.SafeLoader)

    # collapes list
    for db in dbs.keys():
        blast_database_loc = dbs[db]["location"]
        if blast_database_loc == "Default":
            blast_database_loc = str(files(PACKAGE) / "data/BLAST_dbs")
        try:
            parameters = " ".join(dbs[db]["parameters"])
        except KeyError:
            parameters = ""
        dbs[db]["parameters"] = parameters

        if dbs[db]["method"] == "infernal":
            db_loc = " ".join(
                os.path.join(blast_database_loc, x)
                for x in (f"{db}.clanin", f"{db}.cm")
            )
        else:
            db_loc = os.path.join(blast_database_loc, db)
        # dbs[db]['name'] = db
        dbs[db]["db_loc"] = db_loc
        # db_list.append(dbs[db])
    return dbs


def databases_exist() -> bool:
    return os.path.exists(f"{ROOT_DIR}/data/BLAST_dbs/")


def download_databases() -> None:
    # dynamic version number for the databases
    # this is locked at minor version bumps
    # need to upload a new database into github every minor update
    # patch number bumps just refer to the X.X.0 version
    db_loc = f"https://github.com/mmcguffi/pLannotate/releases/download/v{plannotate_version.rsplit('.', 1)[0]}.0/BLAST_dbs.tar.gz"
    # db_loc = "https://github.com/barricklab/pLannotate/releases/download/v1.1.0/BLAST_dbs.tar.gz"

    # subprocess.call(["wget", "-P", f"{ROOT_DIR}/data/", db_loc])
    subprocess.call(["curl", "-L", "-o", f"{ROOT_DIR}/data/BLAST_dbs.tar.gz", db_loc])

    # check if download was successful
    if not os.path.exists(f"{ROOT_DIR}/data/BLAST_dbs.tar.gz"):
        logger.error(
            "Error downloading databases. Please try again or contact the developer."
        )
        sys.exit()

    logger.info("Download complete.")

    logger.info("Extracting...")
    subprocess.call(
        ["tar", "-xzf", f"{ROOT_DIR}/data/BLAST_dbs.tar.gz", "-C", f"{ROOT_DIR}/data/"]
    )
    logger.info("Extraction complete.")

    logger.info("Removing archive...")
    subprocess.call(["rm", f"{ROOT_DIR}/data/BLAST_dbs.tar.gz"])
    logger.info("Removal complete.")

    logger.info("Done.")
