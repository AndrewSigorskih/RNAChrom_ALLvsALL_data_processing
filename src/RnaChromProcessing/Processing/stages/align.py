from logging import getLogger
from pathlib import Path
from typing import List, Literal, Optional

from pydantic import PositiveInt

from .basicstage import BasicStage, SamplePair
from ...utils import (
    check_file_exists, exit_with_error,
    run_command, validate_tool_path
)

logger = getLogger()


class Align(BasicStage):
    tool: Literal['hisat2', 'star', 'bwa', 'custom'] = 'hisat2'
    tool_path: Optional[Path] = None

    dna_genome_path: Optional[Path] = None
    rna_genome_path: Optional[Path] = None

    known_splice: Optional[Path] = None

    tool_threads: PositiveInt = 1

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        # genome path
        if self.tool != 'custom':
            if not self.dna_genome_path:
                exit_with_error('DNA genome path is required if to using custom aligner script!')
            if not self.rna_genome_path:
                logger.debug('RNA genome path not provided, will use DNA genome instead.')
                self.rna_genome_path = self.dna_genome_path
            check_file_exists(self.dna_genome_path)
            check_file_exists(self.rna_genome_path)
        # tool path
        if (self.tool == 'custom' and not self.tool_path) :
            exit_with_error('Path to custom script is required when using the "custom" option!')
        self.tool_path = validate_tool_path(self.tool_path, self.tool)
        # known splice
        if self.tool == 'hisat2':
            if not self.known_splice:
                exit_with_error('Known splicesite infile is required for hisat2 RNA alignment!')
            check_file_exists(self.known_splice)
    
    def run(self, samples: List[SamplePair]) -> List[SamplePair]:
        '''Run chosen alignment tools'''
        # choose tool to run
        if self.tool == 'hisat2':
            func = self._run_hisat
        elif self.tool == 'bwa':
            func = self._run_bwa
        elif self.tool == 'star':
            func = self._run_star
        elif self.tool == 'custom':
            func = self._custom
        # prepare filepaths
        output_samples = self._make_output_samples(samples, new_suff='bam')
        # run
        self.run_function(func, samples, output_samples)
        # return results
        return output_samples

    def _run_hisat(self,
                   inp_sample: SamplePair,
                   out_sample: SamplePair) -> int:
        dna_cmd = (
            f'{self.tool_path} -x {self.dna_genome_path} -p {self.tool_threads} '
            f'--no-spliced-alignment -k 100 --no-unal -U {inp_sample.dna_file} | '  # removed --no-softclip
            f'samtools view -bSh > {out_sample.dna_file}'
        )
        rna_cmd = (
            f'{self.tool_path} -x {self.rna_genome_path} -p {self.tool_threads} -k 100 '
            f'--no-unal --known-splicesite-infile {self.known_splice} --dta-cufflinks '  # removed --no-softclip
            f'-U {inp_sample.rna_file} | samtools view -bSh > {out_sample.rna_file}'
        )
        return_code_1 = run_command(dna_cmd, shell=True)
        return_code_2 = run_command(rna_cmd, shell=True)
        return (return_code_1 or return_code_2)
    
    def _run_bwa(self,
                 inp_sample: SamplePair,
                 out_sample: SamplePair) -> int:
        raise RuntimeError('BWA alignment not implemented yet!')
    
    def _run_star(self,
                  inp_sample: SamplePair,
                  out_sample: SamplePair) -> int:
        raise RuntimeError('STAR alignment not implemented yet!')
