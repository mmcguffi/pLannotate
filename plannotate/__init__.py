__version__ = "1.2.2"
import logging

from .logging_config import setup_logging
from .models import Construct

# Set up logging once at module level
setup_logging(level=logging.INFO)

__all__ = ["Construct"]
