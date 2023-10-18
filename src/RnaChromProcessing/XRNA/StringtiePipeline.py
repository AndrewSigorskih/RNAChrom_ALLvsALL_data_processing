import shutil
from pathlib import Path
from typing import Optional
from typing_extensions import Annotated

from .PoolExecutor import PoolExecutor
from ..utils import exit_with_error, run_command

from pydantic import BaseModel, Field, FieldValidationInfo, field_validator


class StringtieTool(BaseModel):
    stringtie_threads: Optional[int] = 1
    tool_path: Annotated[Optional[Path], Field(validate_default=True)] = None
    
    @field_validator('tool_path')
    @classmethod
    def handle_tool_path(cls, val: Optional[str], info: FieldValidationInfo) -> Path:
        val = val or shutil.which('stringtie')
        if not val:
            exit_with_error('Cannot deduce path to stringtie executable!')
        return Path(val)
    

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
    pass
