import shutil
import gzip
import os
from functools import partial
from mimetypes import guess_type
from typing import Any, Callable, Dict, IO, List

from Bio import SeqIO

from .basicstage import BasicStage
from ...utils import gzip_file, exit_with_error, run_command


def open_handle(filename: str) -> Callable[[str], IO]:
    encoding = guess_type(filename)[1]
    return partial(gzip.open, mode='rt') if encoding == 'gzip' else open


class Rsites(BasicStage):
    def __init__(self, 
                 cfg: Dict[str, Any],
                 *args, **kwargs):
        super().__init__(*args, **kwargs)
        if (cpus := cfg.get('cpus', None)):
            self.cpus: int = cpus
        self.type: str = cfg.get('type', None)
        if self.type == 'custom':
            self.tool_path: str = cfg.get('tool_path', None)
            if not self.tool_path:
               exit_with_error('Path to custom tool for rsite procedure is not specified!') 

    def run(self,
            dna_ids: List[str],
            rna_ids: List[str]):
        """Manage restriction sites filtration"""
        if self.type == 'skip':
            func = self._copy_files
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

    def _imargi_like(self,
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

    def _grid_like(self,
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

    def _custom(self,
                dna_in_file: str,
                rna_in_file: str,
                dna_out_file: str,
                rna_out_file: str) -> int:
        cmd = [self.tool_path, dna_in_file, rna_in_file, dna_out_file, rna_out_file]
        exit_code = run_command(cmd)
        return exit_code
