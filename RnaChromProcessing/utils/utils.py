import gzip
import shutil
import logging
import os
import subprocess
import sys

from typing import Any, List, Union

LOGGING_FORMAT = '%(name)s | line %(lineno)-3d | %(levelname)-8s | %(message)s'
logger = logging.getLogger(name=__name__)

def gzip_file(src: str, dst: str, remove_src=True) -> None:
    with open(src, 'rb') as infile, gzip.open(dst, 'wb') as outfile:
        shutil.copyfileobj(infile, outfile)
    if remove_src:
        os.remove(src)

def find_in_list(id: str, lst: List[str]):
    return next(x for x in lst if x.startswith(id))

def exit_with_error(message: str = '') -> None:
    if not message:
        message = 'Empty error message!'
    logger.critical(message)
    sys.exit(1)

def make_directory(path: str, exist_ok: bool = True) -> None:
    try:
        os.makedirs(path, exist_ok=exist_ok)
    except OSError as e:
        exit_with_error(f'Could not create directory: {e}')

def check_file_exists(path: str) -> None:
    if not (os.path.exists(path) and (os.path.getsize(path) > 0)):
        message = f'File {path} does not exist or is empty!'
        exit_with_error(message)

def run_command(cmd: Union[List[str], str],
                **subprocess_args: Any) -> int:
    if isinstance(cmd, list):
        cmd_str = " ".join(cmd)
    else:
        cmd_str = cmd
    logger.debug(f'Running command: {cmd_str}')
    return_code: int = subprocess.run(cmd, **subprocess_args).returncode
    return return_code

def run_get_stdout(cmd: Union[List[str], str],
                **subprocess_args: Any) -> str:
#    if isinstance(cmd, list):
#        cmd_str = " ".join(cmd)
#    else:
#        cmd_str = cmd
#    logger.debug(f'Running command: {cmd_str}')
    result = subprocess.run(cmd, capture_output=True,
                            text=True, **subprocess_args)
    return result.stdout

def configure_logger(debug: bool = False) -> None:
    level: int = logging.DEBUG if debug else logging.INFO
    handler = logging.StreamHandler()
    handler.setLevel(level)
    handler.setFormatter(logging.Formatter(LOGGING_FORMAT))
    logger.setLevel(level)
    logger.addHandler(handler)
