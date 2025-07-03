import logging
import sys
from typing import Optional


def setup_logging(
    level: int = logging.INFO,
    format_string: Optional[str] = None,
    stream: Optional[logging.StreamHandler] = None,
) -> logging.Logger:
    """
    Set up logging configuration for pLannotate.
    
    Args:
        level: Logging level (default: INFO)
        format_string: Custom format string for log messages
        stream: Custom stream handler (default: sys.stdout)
    
    Returns:
        Configured logger instance
    """
    if format_string is None:
        format_string = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    
    if stream is None:
        stream = logging.StreamHandler(sys.stdout)
    
    # Configure the root logger
    logging.basicConfig(
        level=level,
        format=format_string,
        handlers=[stream],
        force=True
    )
    
    # Get the pLannotate logger
    logger = logging.getLogger("plannotate")
    logger.setLevel(level)
    
    return logger


def get_logger(name: str = "plannotate") -> logging.Logger:
    """
    Get a logger instance for the specified name.
    
    Args:
        name: Logger name (default: "plannotate")
    
    Returns:
        Logger instance
    """
    return logging.getLogger(name) 