from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Literal, List, Optional

from pydantic import BaseModel, PositiveInt

from .basicstage import BasicStage, SamplePair
from ...utils import exit_with_error, run_command, validate_tool_path


class _DedupToolParams(BaseModel):
    memlimit: PositiveInt = 4096
    comparison: Literal['tight', 'loose'] = 'loose'

class Dedup(BasicStage):
    tool: Literal['fastuniq', 'fastq-dupaway', 'skip'] = 'fastuniq'
    tool_path: Optional[Path] = None
    tool_params: Optional[_DedupToolParams] = _DedupToolParams()

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        if self.tool != 'skip':
            self.tool_path = validate_tool_path(self.tool_path, self.tool)

    def run(self,
            samples: List[SamplePair]) -> None:
        """Run chosen deduplication tool"""
        # choose tool to run
        if self.tool == 'skip':
            func = self._copy_files
        elif self.tool == 'fastuniq':
            func = self._run_fastuniq
        elif self.tool == 'fastq-dupaway':
            func = self._run_fqdupaway
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

    def _run_fastuniq(self,
                      dna_in_file: Path,
                      rna_in_file: Path,
                      dna_out_file: Path,
                      rna_out_file: Path) -> int:
        """run fastuniq"""
        if (dna_in_file.suffix == '.gz') or (rna_in_file.suffix == '.gz'):
            exit_with_error('Fastuniq tool does not accept gzipped files!')
        with NamedTemporaryFile(mode='w', dir='.') as f:
            tmpfilename = f.name
            # need to flush the file
            with f.file as temp_file:
                print(rna_in_file, dna_in_file, sep='\n', file=temp_file)
            command = [
                self.tool_path, '-i', tmpfilename,
                '-o', rna_out_file, '-p', dna_out_file
            ]
            return_code = run_command(command)
        return return_code
    
    def _run_fqdupaway(self,
                       dna_in_file: Path,
                       rna_in_file: Path,
                       dna_out_file: Path,
                       rna_out_file: Path) -> int:
        """run fastq-dupaway"""
        command = [
            self.tool_path, '-i', rna_in_file, '-u', dna_in_file,
            '-o', rna_out_file, '-p', dna_out_file, '-m', self.tool_params.memlimit,
            '--compare-seq', self.tool_params.comparison
        ]
        return_code = run_command(command)
        return return_code
