import shutil
import os

from typing import Any, Dict, List

from basicstage import BasicStage
from ...utils import exit_with_error, run_command


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
        else:  # unknown rsites mamaging strategy
            exit_with_error('Unknown restriction-site managing strategy!')
        # run chosen function
        self.run_function(func, dna_ids, rna_ids)