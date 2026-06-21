__version__ = "2.0.0"
import logging

from .logging_config import setup_logging
from .models import Construct

# Set up logging once at module level
setup_logging(level=logging.INFO)

__all__ = ["Construct"]
