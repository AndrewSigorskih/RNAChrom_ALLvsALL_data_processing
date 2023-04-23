import shutil

from typing import Any, Dict, List
from subprocess import DEVNULL

from .basicstage import BasicStage
from ...utils import exit_with_error, run_command

class Trim(BasicStage):
    def __init__(self, 
                 cfg: Dict[str, Any],
                 *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.tool: str = cfg.get('tool', None)
        self.params: Dict[str, Any] = cfg.get('params', [])
        self.tool_path: str = cfg.get('tool_path', shutil.which(self.tool))
        if (cpus := cfg.get('cpus', None)):
            self.cpus: int = cpus
        if not self.tool:
            exit_with_error('Deduplication tool not specified!')
        if not self.tool_path:
            exit_with_error(f'Cannot deduce path to {self.tool} executable!')
    
    def run(self,
            dna_ids: List[str],
            rna_ids: List[str]):
        """Run chosen trimming tool"""
        # choose tool to run
        if self.tool == 'skip':
            func = self._copy_files
        elif (self.tool == 'trimmomatic'):
            func = self._run_trimmomatic
        else:  # unknown tool
            exit_with_error(f'Unknown trimming tool: {self.tool}!')
        # run chosen function
        self.run_function(func, dna_ids, rna_ids)
    
    def _run_trimmomatic(self,
                         dna_in_file: str,
                         rna_in_file: str,
                         dna_out_file: str,
                         rna_out_file: str) -> int:
        """run trimmomatic"""
        window = self.params['window']
        qual_th = self.params['qual_th']
        minlen = self.params['minlen']
        # make sure conda wrapper works as well
        if self.tool_path.endswith('.jar'):
            tool_alias: str = f'java -jar {self.tool_path}'
        else:  # wrapper
            tool_alias = self.tool_path
        command = (
            f'{tool_alias} PE -phred33 {dna_in_file} '
            f'{rna_in_file} {dna_out_file} {dna_out_file}.unpaired '
            f'{rna_out_file} {rna_out_file}.unpaired '
            f'SLIDINGWINDOW:{window}:{qual_th} MINLEN:{minlen}'
        )
        return_code = run_command(command, shell=True, stdout=DEVNULL, stderr=DEVNULL)
        return return_code
