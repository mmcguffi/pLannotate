import hashlib
import os
import shutil
import sys
import tarfile
import tempfile
from importlib.resources import files
from pathlib import Path
from typing import Any, Dict
from urllib.request import urlopen

import yaml

from . import __version__ as plannotate_version
from .logging_config import get_logger

logger = get_logger(__name__)

PACKAGE = __package__ or "plannotate"
DATABASE_ASSET_NAME = "plannotate-databases-v2.tar.gz"
DATABASE_DIRECTORIES = ("BLAST_dbs", "diamond_dbs", "infernal_dbs")
REQUIRED_DATABASE_FILES = (
    "BLAST_dbs/snapgene.nhr",
    "BLAST_dbs/snapgene.nin",
    "BLAST_dbs/snapgene.nsq",
    "BLAST_dbs/snapgene.ndb",
    "BLAST_dbs/snapgene.nog",
    "BLAST_dbs/snapgene.nos",
    "BLAST_dbs/snapgene.not",
    "BLAST_dbs/snapgene.ntf",
    "BLAST_dbs/snapgene.nto",
    "BLAST_dbs/descriptions.db",
    "diamond_dbs/fpbase.dmnd",
    "diamond_dbs/swissprot.dmnd",
    "diamond_dbs/descriptions.db",
    "infernal_dbs/Rfam.cm",
    "infernal_dbs/Rfam.clanin",
    "infernal_dbs/Rfam.cm.i1f",
    "infernal_dbs/Rfam.cm.i1i",
    "infernal_dbs/Rfam.cm.i1m",
    "infernal_dbs/Rfam.cm.i1p",
)


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
    return Path(str(files(PACKAGE) / "data"))


def get_yaml(yaml_file_loc: Path) -> Dict[str, Any]:
    # file_name = get_resource("data", "databases.yml")
    with open(yaml_file_loc, "r") as f:
        dbs = yaml.load(f, Loader=yaml.SafeLoader)

    # collapes list
    for db in dbs.keys():
        database_loc = dbs[db]["location"]
        method = dbs[db]["method"].lower()
        if database_loc == "Default":
            directory_by_method = {
                "blast": "BLAST_dbs",
                "blastn": "BLAST_dbs",
                "diamond": "diamond_dbs",
                "infernal": "infernal_dbs",
            }
            try:
                database_loc = str(
                    files(PACKAGE) / f"data/{directory_by_method[method]}"
                )
            except KeyError as exc:
                raise ValueError(f"Unsupported database method: {method}") from exc
        try:
            parameters = " ".join(dbs[db]["parameters"])
        except KeyError:
            parameters = ""
        dbs[db]["parameters"] = parameters

        if method == "infernal":
            db_loc = " ".join(
                os.path.join(database_loc, x) for x in (f"{db}.clanin", f"{db}.cm")
            )
        else:
            db_loc = os.path.join(database_loc, db)
        # dbs[db]['name'] = db
        dbs[db]["db_loc"] = db_loc
        # db_list.append(dbs[db])
    return dbs


def databases_exist() -> bool:
    data_directory = get_data_directory()
    return all((data_directory / path).is_file() for path in REQUIRED_DATABASE_FILES)


def _database_asset_url() -> str:
    return os.environ.get(
        "PLANNOTATE_DATABASE_URL",
        (
            "https://github.com/mmcguffi/pLannotate/releases/download/"
            f"v{plannotate_version}/{DATABASE_ASSET_NAME}"
        ),
    )


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _safe_extract(archive: Path, destination: Path) -> None:
    destination = destination.resolve()
    with tarfile.open(archive, "r:gz") as tar:
        for member in tar.getmembers():
            if member.issym() or member.islnk():
                raise ValueError(
                    f"Links are not allowed in the database archive: {member.name}"
                )
            member_path = (destination / member.name).resolve()
            if destination not in member_path.parents and member_path != destination:
                raise ValueError(f"Unsafe path in database archive: {member.name}")
        if sys.version_info >= (3, 12):
            tar.extractall(destination, filter="data")
        else:
            tar.extractall(destination)


def _validate_database_tree(data_directory: Path) -> None:
    missing = [
        path
        for path in REQUIRED_DATABASE_FILES
        if not (data_directory / path).is_file()
    ]
    if missing:
        formatted = ", ".join(missing)
        raise ValueError(f"Database archive is missing required files: {formatted}")


def download_databases() -> None:
    """Download, validate, and install the database bundle for this release."""
    url = _database_asset_url()
    expected_sha256 = os.environ.get("PLANNOTATE_DATABASE_SHA256", "").lower()
    data_directory = get_data_directory()

    logger.info("Downloading databases from %s", url)
    with tempfile.TemporaryDirectory(prefix="plannotate-databases-") as temp_dir:
        temp_path = Path(temp_dir)
        archive_path = temp_path / DATABASE_ASSET_NAME
        with urlopen(url, timeout=60) as response, archive_path.open("wb") as output:
            shutil.copyfileobj(response, output)

        if expected_sha256:
            actual_sha256 = _sha256(archive_path)
            if actual_sha256 != expected_sha256:
                raise ValueError(
                    "Database archive checksum mismatch: "
                    f"expected {expected_sha256}, got {actual_sha256}"
                )

        extracted_data = temp_path / "extracted"
        extracted_data.mkdir()
        _safe_extract(archive_path, extracted_data)
        _validate_database_tree(extracted_data)

        data_directory.mkdir(parents=True, exist_ok=True)
        for directory_name in DATABASE_DIRECTORIES:
            destination = data_directory / directory_name
            if destination.exists():
                shutil.rmtree(destination)
            shutil.copytree(extracted_data / directory_name, destination)

    _validate_database_tree(data_directory)
    logger.info("Database installation complete.")
