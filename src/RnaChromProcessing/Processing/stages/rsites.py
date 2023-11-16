import gzip
from functools import partial
from mimetypes import guess_type
from typing import Callable, IO, List, Literal, Optional
from pathlib import Path

import pyfastx

from .basicstage import BasicStage, SamplePair
from ...utils import run_command, validate_tool_path


def open_handle(filename: Path) -> Callable[[str], Callable[[str], IO]]:
    encoding = guess_type(filename)[1]
    return partial(gzip.open, mode='rt') if encoding == 'gzip' else open


def output_handle(filename: Path) -> Callable[[str], Callable[[str], IO]]:
    encoding = guess_type(filename)[1]
    return partial(gzip.open, mode='wt') if encoding == 'gzip' else partial(open, mode='w')


def format_fastq(name: str, seq: str, qual: str) -> str:
    return f'@{name}\n{seq}\n+\n{qual}\n'


class Rsites(BasicStage):
    type: Literal['skip', 'grid', 'imargi', 'custom'] = 'skip'
    tool_path: Optional[Path] = None
    rsite_bgn: str = 'AG',
    rsite_end: str = 'CT'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        if self.type == 'custom':
            self.tool_path = validate_tool_path(
                self.tool_path, self.tool_path
            )
    
    def run(self,
            samples: List[SamplePair]) -> None:
        """Manage restriction sites filtration"""
        # choose strategy to run
        if self.type == 'skip':
            func = self._copy_files
        elif self.type == 'imargi':
            func = self._imargi_like
        elif self.type == 'grid':
            func = self._grid_like
        elif self.type == 'custom':
            func = self._custom
        # prepare filepaths
        dna_outputs = [
            self._stage_dir / sample.dna_file.name for sample in samples
        ]
        rna_outputs = [
            self._stage_dir / sample.rna_file.name for sample in samples
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
        
    def _imargi_like(self,
                     dna_in_file: Path,
                     rna_in_file: Path,
                     dna_out_file: Path,
                     rna_out_file: Path) -> int:
        '''Save read pair if DNA read starts with CT or NT\n
           Remove first 2 bases from RNA reads.\n
           Reads in files should be synchronized.'''
        _dna_out = output_handle(dna_out_file)
        _rna_out = output_handle(rna_out_file)
        with _dna_out(dna_out_file) as dna_out_handle,\
             _rna_out(rna_out_file) as rna_out_handle:
            # read = (name, seq, qual)
            for (dna_name, dna_seq, dna_qual), (rna_name, rna_seq, rna_qual) in zip(
                pyfastx.Fastq(dna_in_file, build_index=False, full_name=True),
                pyfastx.Fastq(rna_in_file, build_index=False, full_name=True)
            ):
                if (not dna_seq.startswith('CT')) and (not dna_seq.startswith('NT')):
                    continue
                print(format_fastq(dna_name, dna_seq, dna_qual), file=dna_out_handle, end='')
                print(format_fastq(rna_name, rna_seq[2:], rna_qual[2:]), file=rna_out_handle, end='')
        return 0
    
    def _grid_like(self,
                   dna_in_file: Path,
                   rna_in_file: Path,
                   dna_out_file: Path,
                   rna_out_file: Path) -> int:
        '''Save read pair if DNA reads end with AG\n
           Add CT to the end of selected DNA reads\n
           Assign quality values from terminal AG to novel CT\n
           Reads in files should be synchronized.'''
        n_bases = len(self.rsite_end)
        _dna_out = output_handle(dna_out_file)
        _rna_out = output_handle(rna_out_file)
        with _dna_out(dna_out_file) as dna_out_handle,\
             _rna_out(rna_out_file) as rna_out_handle:
            # read = (name, seq, qual)
            for (dna_name, dna_seq, dna_qual), (rna_name, rna_seq, rna_qual) in zip(
                pyfastx.Fastq(dna_in_file, build_index=False, full_name=True),
                pyfastx.Fastq(rna_in_file, build_index=False, full_name=True)
            ):
                if not dna_seq.endswith(self.rsite_bgn):
                    continue
                print(format_fastq(rna_name, rna_seq, rna_qual), file=rna_out_handle, end='')
                dna_seq += self.rsite_end
                dna_qual += dna_qual[-n_bases:]
                print(format_fastq(dna_name, dna_seq, dna_qual), file=dna_out_handle, end='')
        return 0
    
    def _custom(self,
                dna_in_file: Path,
                rna_in_file: Path,
                dna_out_file: Path,
                rna_out_file: Path) -> int:
        cmd = [self.tool_path, dna_in_file, rna_in_file, dna_out_file, rna_out_file]
        exit_code = run_command(cmd)
        return exit_code
