import shutil
import os

from tempfile import NamedTemporaryFile
from typing import Any, Dict, List

from .basicstage import BasicStage
from ...utils import exit_with_error, run_command


class Dedup(BasicStage):
    def __init__(self, 
                 cfg: Dict[str, Any],
                 *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.tool: str = cfg.get('tool', None)
        if not self.tool:
            exit_with_error('Deduplication tool not specified!')
        self.params: Dict[str, Any] = cfg.get('params', {})
        self.tool_path: str = cfg.get('tool_path', shutil.which(self.tool))
        if (cpus := cfg.get('cpus', None)):
            self.cpus: int = cpus
        if not self.tool_path:
            exit_with_error(f'Cannot deduce path to {self.tool} executable!')

    def run(self,
            dna_ids: List[str],
            rna_ids: List[str]):
        """Run chosen deduplication tool"""
        # choose tool to run
        if self.tool == 'skip':
            func = self._copy_files
        elif self.tool == 'fastuniq':
            func = self._run_fastuniq
        elif self.tool == 'fastq-dupaway':
            func = self._run_fqdupaway
        else:  # unknown tool
            exit_with_error(f'Unknown deduplication tool: {self.tool}!')
        # run chosen function
        self.run_function(func, dna_ids, rna_ids)

    def _run_fastuniq(self,
                      dna_in_file: str,
                      rna_in_file: str,
                      dna_out_file: str,
                      rna_out_file: str) -> int:
        """run fastuniq"""
        if (dna_in_file.endswith('.gz')) or (rna_in_file.endswith('.gz')):
            exit_with_error('Fastuniq tool does not accept gzipped files!')
        with NamedTemporaryFile(mode='w', delete=False, dir='.') as f:
            tmpfilename = f.name
            print(rna_in_file, dna_in_file, sep='\n', file=f)
        command = (f'{self.tool_path} -i {tmpfilename} '
                    f'-o {rna_out_file} -p {dna_out_file}')
        return_code = run_command(command, shell=True)
        os.remove(tmpfilename)
        return return_code
    
    def _run_fqdupaway(self,
                       dna_in_file: str,
                       rna_in_file: str,
                       dna_out_file: str,
                       rna_out_file: str) -> int:
        memlimit = str(self.params.get('memlimit', '2048'))
        command = (
            f'{self.tool_path} -i {rna_in_file} -u {dna_in_file} '
            f'-o {rna_out_file} -p {dna_out_file} -m {memlimit}'
        )
        return_code = run_command(command, shell=True)
        return return_code
