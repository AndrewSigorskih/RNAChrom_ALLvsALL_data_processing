import gzip
import shutil
import logging
import os
import subprocess
import sys

from typing import Any, List, Union

DATE_FORMAT = '%Y-%m-%d %H:%M:%S'
LOGGING_FORMAT = '%(asctime)s | %(levelname)-8s | %(message)s'
logger = logging.getLogger()
VERBOSE = 5

def check_file_exists(path: str) -> None:
    if not (os.path.exists(path) and (os.path.getsize(path) > 0)):
        message = f'File {path} does not exist or is empty!'
        exit_with_error(message)


def configure_logger(logger: logging.Logger, verbose: int = 0) -> None:
    level: int = logging.INFO if not verbose else logging.DEBUG if verbose == 1 else VERBOSE
    logging.addLevelName(VERBOSE, "VERBOSE")
    handler = logging.StreamHandler()
    handler.setLevel(level)
    handler.setFormatter(logging.Formatter(fmt=LOGGING_FORMAT,
                                           datefmt=DATE_FORMAT))
    logger.setLevel(level)
    logger.addHandler(handler)


def exit_with_error(message: str = '') -> None:
    """Print error and exit"""
    if not message:
        message = 'Empty error message!'
    logger.critical(message)
    sys.exit(1)


def find_in_list(id: str, lst: List[str]):
    """A simple helper to find file from
    list of files regardless of extension"""
    try:
        return next(x for x in lst if x.startswith(id))
    except StopIteration:
        return None


def gzip_file(src: str, dst: str, remove_src=True) -> None:
    """Gzip file to specified name"""
    with open(src, 'rb') as infile, gzip.open(dst, 'wb') as outfile:
        shutil.copyfileobj(infile, outfile)
    if remove_src:
        os.remove(src)


def make_directory(path: str, exist_ok: bool = True) -> None:
    try:
        os.makedirs(path, exist_ok=exist_ok)
    except OSError as e:
        exit_with_error(f'Could not create directory: {e}')


def move_exist_ok(src_dir: str, dest_dir: str):
    """Move directory src to dest even if dest_dir
    and some of its contents exist"""
    fileList = os.listdir(src_dir)
    for i in fileList:
        src = os.path.join(src_dir, i)
        dest = os.path.join(dest_dir, i)
        if os.path.exists(dest):
            if os.path.isdir(dest):
                move_exist_ok(src, dest)
                continue
            else:
                os.remove(dest)
        shutil.move(src, dest_dir)


def run_command(cmd: Union[List[str], str],
                **subprocess_args: Any) -> int:
    log_cmd(cmd)
    return_code: int = subprocess.run(cmd, **subprocess_args).returncode
    return return_code


def run_get_stdout(cmd: Union[List[str], str],
                   **subprocess_args: Any) -> str:
    log_cmd(cmd)
    result = subprocess.run(cmd, capture_output=True,
                            text=True, **subprocess_args)
    return result.stdout


def log_cmd(cmd: Union[List[str], str]) -> None:
    if isinstance(cmd, list):
        cmd_str = " ".join(cmd)
    else:
        cmd_str = cmd
    logger.log(VERBOSE, f'Running command: {cmd_str}')
