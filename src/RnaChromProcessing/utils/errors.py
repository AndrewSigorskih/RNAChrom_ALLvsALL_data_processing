import logging
import sys

logger = logging.getLogger()

class StageFailedError(RuntimeError):
    pass


def exit_with_error(message: str = '') -> None:
    """Print error and exit"""
    if not message:
        message = 'Empty error message!'
    logger.critical(message)
    sys.exit(1)
