import gzip
import os
from functools import partial
from mimetypes import guess_type
from typing import Any, Callable, Dict, IO, List

from Bio import SeqIO
import pyfastx

from .basicstage import BasicStage
from ...utils import gzip_file, exit_with_error, run_command


def open_handle(filename: str) -> Callable[[str], Callable[[str], IO]]:
    encoding = guess_type(filename)[1]
    return partial(gzip.open, mode='rt') if encoding == 'gzip' else open


def output_handle(filename: str) -> Callable[[str], Callable[[str], IO]]:
    encoding = guess_type(filename)[1]
    return partial(gzip.open, mode='wt') if encoding == 'gzip' else partial(open, mode='w')


def format_fastq(name: str, seq: str, qual: str) -> str:
    return f'@{name}\n{seq}\n+\n{qual}\n'


class Rsites(BasicStage):
    def __init__(self, 
                 cfg: Dict[str, Any],
                 to_keep: bool,
                 *args, **kwargs):
        super().__init__(*args, **kwargs)
        if (cpus := cfg.get('cpus', None)):
            self.cpus: int = cpus
        self.type: str = cfg.get('type', None)
        self.to_keep = to_keep
        if self.type == 'custom':
            self.tool_path: str = cfg.get('tool_path', None)
            if not self.tool_path:
               exit_with_error('Path to custom tool for rsite procedure is not specified!')
        elif self.type == 'grid':
            self.rsite_bgn: str = cfg.get('rsite_bgn')
            self.rsite_end: str = cfg.get('rsite_end')

    def run(self,
            dna_ids: List[str],
            rna_ids: List[str]):
        """Manage restriction sites filtration"""
        if self.type == 'skip':
            func = self._copy_files if self.to_keep else self._symlink_files
        elif self.type == 'imargi':
            func = self._imargi_like
        elif self.type == 'grid':
            func = self._grid_like
        elif self.type == 'custom':
            func = self._custom
        else:  # unknown rsites managing strategy
            exit_with_error('Unknown restriction-site managing strategy!')
        # run chosen function
        self.run_function(func, dna_ids, rna_ids)

    def x_imargi_like(self,
                     dna_in_file: str,
                     rna_in_file: str,
                     dna_out_file: str,
                     rna_out_file: str) -> int:
        
        '''Save read pair if DNA read starts with CT or NT\n
           Remove first 2 bases from RNA reads.\n
           Reads in files should be synchronized.'''
        tmp_dna_outfile = dna_out_file + '.tmp'
        tmp_rna_outfile = rna_out_file + '.tmp'
        _dna_open = open_handle(dna_in_file)
        _rna_open = open_handle(rna_in_file)
        with open(tmp_dna_outfile, 'w') as dna_out_handle,\
             open(tmp_rna_outfile, 'w') as rna_out_handle,\
             _dna_open(dna_in_file) as dna_in_handle,\
             _rna_open(rna_in_file) as rna_in_handle:
            for dna_read, rna_read in zip(SeqIO.parse(dna_in_handle, format='fastq'),
                                          SeqIO.parse(rna_in_handle, format='fastq')):
                if dna_read.seq.startswith('CT') or dna_read.seq.startswith('NT'):
                    SeqIO.write(dna_read, handle=dna_out_handle, format='fastq')
                    rna_read = rna_read[2:]
                    SeqIO.write(rna_read, handle=rna_out_handle, format='fastq')
        # deal with gz/txt outputs
        if dna_out_file.endswith('.gz'):
            gzip_file(tmp_dna_outfile, dna_out_file)
        else:
            os.rename(tmp_dna_outfile, dna_out_file)
        
        if rna_out_file.endswith('.gz'):
            gzip_file(tmp_rna_outfile, rna_out_file)
        else:
            os.rename(tmp_rna_outfile, rna_out_file)
        return 0

    def x_grid_like(self,
                   dna_in_file: str,
                   rna_in_file: str,
                   dna_out_file: str,
                   rna_out_file: str) -> int:
        '''Save read pair if DNA reads end with AG\n
           Add CT to the end of selected DNA reads\n
           Reads in files should be synchronized.'''
        tmp_dna_outfile = dna_out_file + '.tmp'
        tmp_rna_outfile = rna_out_file + '.tmp'
        _dna_open = open_handle(dna_in_file)
        _rna_open = open_handle(rna_in_file)
        with open(tmp_dna_outfile, 'w') as dna_out_handle,\
             open(tmp_rna_outfile, 'w') as rna_out_handle,\
             _dna_open(dna_in_file) as dna_in_handle,\
             _rna_open(rna_in_file) as rna_in_handle:
            for dna_read, rna_read in zip(SeqIO.parse(dna_in_handle, format='fastq'),
                                          SeqIO.parse(rna_in_handle, format='fastq')):
                if dna_read.seq.endswith('AG'):
                    SeqIO.write(rna_read, handle=rna_out_handle, format='fastq')
                    annotations = dict(dna_read.letter_annotations)
                    annotations['phred_quality'] += annotations['phred_quality'][-2:]
                    dna_read += 'CT'
                    dna_read.letter_annotations = annotations
                    SeqIO.write(dna_read, handle=dna_out_handle, format='fastq')
        # deal with gz/txt outputs
        if dna_out_file.endswith('.gz'):
            gzip_file(tmp_dna_outfile, dna_out_file)
        else:
            os.rename(tmp_dna_outfile, dna_out_file)
        
        if rna_out_file.endswith('.gz'):
            gzip_file(tmp_rna_outfile, rna_out_file)
        else:
            os.rename(tmp_rna_outfile, rna_out_file)
        return 0
    
    def _imargi_like(self,
                     dna_in_file: str,
                     rna_in_file: str,
                     dna_out_file: str,
                     rna_out_file: str) -> int:
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
                   dna_in_file: str,
                   rna_in_file: str,
                   dna_out_file: str,
                   rna_out_file: str) -> int:
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
                dna_in_file: str,
                rna_in_file: str,
                dna_out_file: str,
                rna_out_file: str) -> int:
        cmd = [self.tool_path, dna_in_file, rna_in_file, dna_out_file, rna_out_file]
        exit_code = run_command(cmd)
        return exit_code
