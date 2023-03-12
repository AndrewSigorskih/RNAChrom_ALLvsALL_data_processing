from typing import Any, Dict, List

from ..utils import exit_with_error, run_command

class Dedup:
    def __init__(self, cfg: Dict[str, Any]):
        self.tool: str = cfg.get('tool', None)
        self.params: List[str] = cfg.get('params', [])
        self.cpus = cfg.get('cpus', None)
        if not self.tool:
            exit_with_error('Deduplication tool not specified!')

    def __call__(self, global_cpus: int):
        if self.tool == 'skip':
            pass

    def _run_fastuniq(self):
        pass