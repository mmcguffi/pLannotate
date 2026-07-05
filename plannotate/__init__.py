"""Public pLannotate API."""

import logging
from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("plannotate")
except PackageNotFoundError:  # not installed, e.g. running from a source tree
    __version__ = "2.0.0"

from ._rotate import RotationResult, rotate_to_origin
from .gather_databases import build_databases
from .models import Construct, Feature

logging.getLogger(__name__).addHandler(logging.NullHandler())

__all__ = [
    "Construct",
    "Feature",
    "RotationResult",
    "build_databases",
    "rotate_to_origin",
]
