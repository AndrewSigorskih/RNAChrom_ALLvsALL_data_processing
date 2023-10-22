from os import remove
from pathlib import Path
from subprocess import DEVNULL
from typing import Literal, List, Optional

from pydantic import BaseModel, PositiveInt

from .basicstage import BasicStage, SamplePair
from ...utils import run_command, validate_tool_path


class _TrimToolParams(BaseModel):
    window: PositiveInt = 5
    qual_th: PositiveInt = 26
    minlen: PositiveInt = 15

class Trim(BasicStage):
    tool: Literal['trimmomatic', 'skip'] = 'trimmomatic'
    tool_path: Optional[Path] = None
    tool_params: Optional[_TrimToolParams] = _TrimToolParams()

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        if self.tool != 'skip':
            self.tool_path = validate_tool_path(self.tool_path, self.tool)
    
    def run(self,
            samples: List[SamplePair]) -> None:
        """Run chosen trimming tool"""
        # choose tool to run
        if self.tool == 'skip':
            func = self._copy_files
        elif (self.tool == 'trimmomatic'):
            func = self._run_trimmomatic
        # prepare filepaths
        dna_outputs = [
            self.stage_dir / sample.dna_file.name for sample in samples
        ]
        rna_outputs = [
            self.stage_dir / sample.rna_file.name for sample in samples
        ]
        # run function
        self.run_function(
            func,
            [sample.dna_file for sample in samples],
            [sample.rna_file for sample in samples],
            dna_outputs, rna_outputs
        )
        # save paths
        for sample, dna_out, rna_out in zip(samples, dna_outputs, rna_outputs):
            sample.set_files(dna_out, rna_out)
    
    def _run_trimmomatic(self,
                         dna_in_file: Path,
                         rna_in_file: Path,
                         dna_out_file: Path,
                         rna_out_file: Path) -> int:
        """run trimmomatic"""
        # make sure conda wrapper works as well
        if self.tool_path.endswith('.jar'):
            tool_alias: str = f'java -jar {self.tool_path}'
        else:  # wrapper
            tool_alias = self.tool_path
        command = (
            f'{tool_alias} PE -phred33 {dna_in_file} '
            f'{rna_in_file} {dna_out_file} {dna_out_file}.unpaired '
            f'{rna_out_file} {rna_out_file}.unpaired '
            f'SLIDINGWINDOW:{self.tool_params.window}:{self.tool_params.qual_th} MINLEN:{self.tool_params.minlen}'
        )
        return_code = run_command(command, shell=True, stdout=DEVNULL, stderr=DEVNULL)
        remove(f'{dna_out_file}.unpaired')
        remove(f'{rna_out_file}.unpaired')
        return return_code
