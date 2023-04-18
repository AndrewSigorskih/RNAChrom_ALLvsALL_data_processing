import shutil
import os
from typing import Any, Dict, List

import pyfastx

from .basicstage import BasicStage
from ...utils import gzip_file, exit_with_error, run_command


class Rsites(BasicStage):
    def __init__(self, 
                 cfg: Dict[str, Any],
                 *args, **kwargs):
        super().__init__(*args, **kwargs)
        if (cpus := cfg.get('cpus', None)):
            self.cpus: int = cpus
        self.type: str = cfg.get('type', None)

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
        
        tmp_dna_outfile = dna_out_file + '.tmp'
        tmp_rna_outfile = rna_out_file + '.tmp'
        # save DNA reads that start with CT or NT
        dna_reads = pyfastx.Fastq(dna_in_file, build_index=True)
        with open(tmp_dna_outfile, 'w') as f:
            for read in dna_reads:
                if read.seq.startswith('NT') or read.seq.startswith('CT'):
                    print(read.raw, file=f, end='')
        # get RNA reads corresponding to saved DNA reads
        # remove first 2 bases from RNA reads
        dna_reads = pyfastx.Fastq(tmp_dna_outfile, build_index=True)
        rna_reads = pyfastx.Fastq(rna_in_file, build_index=True)
        with open(tmp_rna_outfile, 'w') as f:
            for read in rna_reads:
                if read.id in dna_reads:
                    read.seq = read.seq[2:]
                    print(read.raw, file=f, end='')
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

    def _grid_like(self):
        raise NotImplementedError("GRID procedure is not implemented yet!")

    def _custom(self):
        raise NotImplementedError("Custom rsite-dealing script usage is not implemented yet!")
