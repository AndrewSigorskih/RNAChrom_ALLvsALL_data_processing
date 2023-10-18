import shutil
from functools import partial
from pathlib import Path
from typing import Optional
from typing_extensions import Annotated

from pydantic import BaseModel, Field, field_validator

from .PoolExecutor import PoolExecutor
from ..utils import run_command, validate_tool_path


class StringtieTool(BaseModel):
    stringtie_threads: Optional[int] = 1
    tool_path: Annotated[Optional[Path], Field(validate_default=True)] = None
    
    @field_validator('tool_path')
    @classmethod
    def validate_stringtie_path(cls, val: Optional[Path]) -> Path:
        val_fn = partial(validate_tool_path, tool_name='stringtie')
        return val_fn(val)
    

    def run_stringtie(self,
                      in_file: Path,
                      out_file: Path) -> None:
        pass

    def run_stringtie_merge(self,
                            assembly_list: Path,
                            output_file: Path) -> None:
        cmd = (
            f'{self.tool_path} merge -p {self.stringtie_threads} '
            f'-o {output_file} {assembly_list}'
        )
        run_command(cmd, shell=True)

class StringtiePipeline:
    def __init__(self,
                 work_pth: Path,
                 executor: PoolExecutor,
                 stringtie: StringtieTool):
        self.executor = executor
        self.stringtie_tool = stringtie
