import shutil
from os import remove
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional

from pydantic import PositiveInt

from .basicstage import BasicStage, SamplePair
from ...utils import (
    check_file_exists, exit_with_error, run_command, validate_tool_path
)

def _remove_suff(pth: Path) -> str:
    '''remove all suffixes from Path object'''
    return str(pth).removesuffix(''.join(pth.suffixes))


class Align(BasicStage):
    dna_tool: Literal['hisat', 'bwa', 'custom'] = 'hisat'
    rna_tool: Optional[Literal['hisat', 'star', 'bwa', 'custom']] = None

    dna_tool_path: Optional[Path] = None
    rna_tool_path: Optional[Path] = None

    dna_genome_path: Path
    rna_genome_path: Optional[Path] = None

    known_splice: Optional[Path] = None

    tool_threads: PositiveInt = 1

    def __inin__(self, **kwargs):
        super().__init__(**kwargs)
        if not self.rna_tool:
            self.rna_tool = self.dna_tool
        if not self.rna_genome_path:
            self.rna_genome_path = self.dna_genome_path
        if (self.dna_tool == 'custom' and not self.dna_tool_path) or\
              (self.rna_tool == 'custom' and not self.rna_tool_path):
            exit_with_error('Path to custom script is required when using the "custom" option!')
        self.dna_tool_path = validate_tool_path(self.dna_tool_path, self.dna_tool)
        self.rna_tool_path = validate_tool_path(self.rna_tool_path, self.rna_tool)
        check_file_exists(self.dna_genome_path)
        check_file_exists(self.rna_genome_path)
        if self.rna_tool == 'hisat':
            if not self.known_splice:
                exit_with_error('Known splicesite infile is required for hisat RNA alignment!')
            check_file_exists(self.known_splice)
    
    def run(self,
            samples: List[SamplePair]) -> None:
        '''Run chosen alignment tools'''
        # choose dna tools to run
        if self.dna_tool == 'hisat':
            self._align_dna = self._hisat_dna
        elif self.dna_tool == 'bwa':
            self._align_dna = self._bwa
        elif self.dna_tool == 'custom':
            self._align_dna = self._dna_custom
        # choose rna tools to run
        if self.rna_tool == 'hisat':
            self._align_rna = self._hisat_rna
        elif self.rna_tool == 'star':
            self._align_rna = self._star
        elif self.rna_tool == 'bwa':
            self._align_rna = self._bwa
        elif self.rna_tool == 'custom':
            self._align_rna = self._rna_custom
        # prepare filepaths
        dna_outputs = [
            self._stage_dir / f'{_remove_suff(sample.dna_file.name)}.bam' for sample in samples
        ]
        rna_outputs = [
            self._stage_dir / f'{_remove_suff(sample.rna_file.name)}.bam' for sample in samples
        ]
        # run
        self.run_function(
            self._run_one,
            [sample.dna_file for sample in samples],
            [sample.rna_file for sample in samples],
            dna_outputs, rna_outputs
        )
        # outputs
        for sample, dna_out, rna_out in zip(samples, dna_outputs, rna_outputs):
            sample.set_files(dna_out, rna_out)

    def _run_one(self,
                 dna_in_file: Path,
                 rna_in_file: Path,
                 dna_out_file: Path,
                 rna_out_file: Path) -> int:
        ret_1 = self._align_dna(dna_in_file, dna_out_file)
        ret_2 = self._align_rna(rna_in_file, rna_out_file)
        return ret_1 or ret_2
    
    def _hisat_dna(self, dna_in_file: Path, dna_out_file: Path) -> int:
        cmd = (
            f'{self.dna_tool_path} -x {self.dna_genome_path} -p {self.tool_threads} '
            f'--no-unal -k 100 --no-spliced-alignment -U {dna_in_file} | '
            f'samtools view -bSh > {dna_out_file}'
        )
        return run_command(cmd, shell=True)
    
    def _hisat_rna(self, rna_in_file: Path, rna_out_file: Path) -> int:
        cmd = (
            f'{self.rna_tool_path} -x {self.rna_genome_path} -p {self.tool_threads} '
            f'--no-unal -k 100 -known-splicesite-infile {self.known_splice} --dta-cufflinks '
            f'-U {rna_in_file} | samtools view -bSh > {rna_out_file}'
        )
        return run_command(cmd, shell=True)
    
    def _dna_custom(self, dna_in_file: Path, dna_out_file: Path) -> int:
        cmd = [
            self.dna_tool_path, dna_in_file, dna_out_file
        ]
        return run_command(cmd)
    
    def _rna_custom(self, rna_in_file: Path, rna_out_file: Path) -> int:
        cmd = [
            self.rna_tool_path, rna_in_file, rna_out_file
        ]
        return run_command(cmd)

    def _bwa(self, dna_in_file: Path, dna_out_file: Path) -> int:
        pass

    def _star(self, dna_in_file: Path, dna_out_file: Path) -> int:
        pass
