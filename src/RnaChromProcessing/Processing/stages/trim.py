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
    tool_params: _TrimToolParams = _TrimToolParams()

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        if self.tool != 'skip':
            self.tool_path = validate_tool_path(self.tool_path, self.tool)
    
    def run(self, samples: List[SamplePair]) -> List[SamplePair]:
        """Run chosen trimming tool"""
        # choose tool to run
        if self.tool == 'skip':
            func = self._copy_files
        elif (self.tool == 'trimmomatic'):
            func = self._run_trimmomatic
        # prepare filepaths
        output_samples = self._make_output_samples(samples)
        # run function
        self.run_function(func, samples, output_samples)
        # return results
        return output_samples
    
    def _run_trimmomatic(self,
                         inp_sample: SamplePair,
                         out_sample: SamplePair) -> int:
        """run trimmomatic"""
        # make sure conda wrapper works as well
        if self.tool_path.endswith('.jar'):
            tool_alias: str = f'java -jar {self.tool_path}'
        else:  # wrapper
            tool_alias = self.tool_path
        command = (
            f'{tool_alias} PE -phred33 {inp_sample.dna_file} {inp_sample.rna_file} '
            f'{out_sample.dna_file} {out_sample.dna_file}.unpaired '
            f'{out_sample.rna_file} {out_sample.rna_file}.unpaired '
            f'SLIDINGWINDOW:{self.tool_params.window}:{self.tool_params.qual_th} MINLEN:{self.tool_params.minlen}'
        )
        return_code = run_command(command, shell=True, stdout=DEVNULL, stderr=DEVNULL)
        remove(f'{out_sample.dna_file}.unpaired')
        remove(f'{out_sample.rna_file}.unpaired')
        return return_code
