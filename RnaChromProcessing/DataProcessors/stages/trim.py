import shutil

from typing import Any, Dict, List

from basicstage import BasicStage
from ...utils import exit_with_error, run_command

class Trim(BasicStage):
    def __init__(self, 
                 cfg: Dict[str, Any],
                 *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.tool: str = cfg.get('tool', None)
        self.params: List[str] = cfg.get('params', [])
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