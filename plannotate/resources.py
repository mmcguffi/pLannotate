import os
import subprocess
import sys
from importlib.resources import files
from pathlib import Path
from typing import Any, Dict

import yaml

from . import __version__ as plannotate_version
from .logging_config import get_logger

logger = get_logger(__name__)

PACKAGE = __package__ or "plannotate"
ROOT_DIR = os.path.dirname(os.path.abspath(__file__))


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


def get_data_directory() -> Path:
    """Get the path to the data directory."""
    return get_resource("data", "")


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
