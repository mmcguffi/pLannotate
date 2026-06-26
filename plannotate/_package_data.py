"""Manage packaged assets, annotation configuration, and database installation."""

import hashlib
import json
import logging
import os
import shutil
import sys
import tarfile
import tempfile
from importlib.resources import files
from pathlib import Path
from typing import Any
from urllib.request import urlopen

import yaml

from . import __version__ as plannotate_version
from ._tools import methods

logger = logging.getLogger(__name__)

PACKAGE = __package__ or "plannotate"
DATABASE_ASSET_NAME = "plannotate-databases-v2.tar.gz"
DATABASE_MANIFEST_NAME = "database-manifest.json"
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
    DATABASE_MANIFEST_NAME,
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


def get_data_directory() -> Path:
    """Get the path to the data directory."""
    return Path(str(files(PACKAGE) / "data"))


def get_database_manifest() -> dict[str, Any]:
    """Load provenance for the installed database bundle."""
    manifest_path = get_data_directory() / DATABASE_MANIFEST_NAME
    if not manifest_path.is_file():
        raise FileNotFoundError(
            "Database manifest not found; reinstall with 'plannotate setupdb --force'."
        )
    with manifest_path.open() as handle:
        return json.load(handle)


def _normalize_parameters(database_name: str, parameters: object) -> str:
    if isinstance(parameters, str):
        return parameters
    if isinstance(parameters, list) and all(
        isinstance(parameter, str) for parameter in parameters
    ):
        return " ".join(parameters)
    raise ValueError(
        f"Database {database_name!r} parameters must be a string or list of strings"
    )


def _database_paths(
    database_name: str, method: str, configured_location: object
) -> dict[str, str]:
    database_directory = methods.database_directory(method)
    if database_directory is None:
        return {}
    if not isinstance(configured_location, str) or not configured_location:
        raise ValueError(f"Source {database_name!r} location must be a string")
    location = (
        Path(str(files(PACKAGE) / f"data/{database_directory}"))
        if configured_location == "Default"
        else Path(configured_location)
    )
    return methods.database_paths(method, database_name, location)


def get_yaml(yaml_file_loc: str | Path) -> dict[str, dict[str, Any]]:
    """Load and normalize an annotation-source configuration."""
    with Path(yaml_file_loc).open() as handle:
        raw_config = yaml.safe_load(handle)
    if not isinstance(raw_config, dict) or not raw_config:
        raise ValueError("Database configuration must be a non-empty mapping")

    databases: dict[str, dict[str, Any]] = {}
    for database_name, raw_database in raw_config.items():
        if not isinstance(database_name, str) or not isinstance(raw_database, dict):
            raise ValueError("Each database configuration must be a named mapping")
        missing = {"method", "location", "priority", "details"} - raw_database.keys()
        if missing:
            raise ValueError(
                f"Database {database_name!r} is missing: {', '.join(sorted(missing))}"
            )

        database = dict(raw_database)
        method = str(database["method"]).lower()
        methods.database_directory(method)
        database["method"] = method

        location = database["location"]
        priority = database["priority"]
        if isinstance(priority, bool) or not isinstance(priority, int) or priority < 1:
            raise ValueError(
                f"Database {database_name!r} priority must be a positive integer"
            )

        database["parameters"] = _normalize_parameters(
            database_name, database.get("parameters", [])
        )

        details = database["details"]
        if not isinstance(details, dict) or "location" not in details:
            raise ValueError(
                f"Database {database_name!r} details must include a location"
            )
        if details["location"] is not None and not isinstance(details["location"], str):
            raise ValueError(
                f"Database {database_name!r} details location must be a string or null"
            )
        database["details"] = dict(details)
        database.update(_database_paths(database_name, method, location))
        databases[database_name] = database
    return databases


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
            if not (member.isfile() or member.isdir()):
                raise ValueError(
                    f"Special files are not allowed in the database archive: {member.name}"
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


def _validate_database_manifest(data_directory: Path) -> None:
    manifest_path = data_directory / DATABASE_MANIFEST_NAME
    try:
        manifest = json.loads(manifest_path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError("Database archive has an invalid manifest") from exc
    if manifest.get("schema_version") != 1:
        raise ValueError("Database archive has an unsupported manifest schema")
    checksums = manifest.get("files")
    if not isinstance(checksums, dict):
        raise ValueError("Database manifest does not contain file checksums")
    payload_paths = set(REQUIRED_DATABASE_FILES) - {DATABASE_MANIFEST_NAME}
    if set(checksums) != payload_paths:
        raise ValueError("Database manifest file list does not match the archive")
    for relative_path in sorted(payload_paths):
        checksum_entry = checksums[relative_path]
        if not isinstance(checksum_entry, dict):
            raise ValueError(f"Invalid database manifest entry for {relative_path}")
        expected = checksum_entry.get("sha256")
        actual = _sha256(data_directory / relative_path)
        if actual != expected:
            raise ValueError(f"Database checksum mismatch for {relative_path}")


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
        _validate_database_manifest(extracted_data)

        data_directory.mkdir(parents=True, exist_ok=True)
        for directory_name in DATABASE_DIRECTORIES:
            destination = data_directory / directory_name
            if destination.exists():
                shutil.rmtree(destination)
            shutil.copytree(extracted_data / directory_name, destination)
        shutil.copy2(
            extracted_data / DATABASE_MANIFEST_NAME,
            data_directory / DATABASE_MANIFEST_NAME,
        )

    _validate_database_tree(data_directory)
    logger.info("Database installation complete.")
