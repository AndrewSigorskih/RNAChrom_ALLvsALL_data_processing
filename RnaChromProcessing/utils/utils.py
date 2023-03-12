import os
import subprocess
import sys

from typing import Any, Dict, List, Union

def exit_with_error(message: str = '') -> None:
    if not message:
        message = 'Empty error message!'
    print(message, file=sys.stderr)
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
                **subprocess_args: Dict[str, Any]) -> None:
    if isinstance(cmd, list):
        cmd_str = " ".join(cmd)
    else:
        cmd_str = cmd
    #TODO change to logging
    print(f'Running command: {cmd_str}')
    return_code = subprocess.run(cmd, **subprocess_args).returncode
    if return_code != 0:
        msg = f'Failed with code: {return_code}'
        exit_with_error(msg)