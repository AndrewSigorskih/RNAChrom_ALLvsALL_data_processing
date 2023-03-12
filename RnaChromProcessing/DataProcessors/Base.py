import os
from tempfile import TemporaryDirectory
from typing import Any, Dict, List

from stages import Dedup
from ..utils import exit_with_error, run_command
from ..utils import dedup_default_cfg

# https://stackoverflow.com/questions/74291040/running-multiple-external-commands-in-parallel-in-python

# https://stackoverflow.com/questions/9554544/python-running-command-line-tools-in-parallel

# https://stackoverflow.com/questions/12097406/python-parallel-commands

class BaseProcessor:
    def __init__(self, cfg: Dict[str, Any]):
        # get basic parameters
        self.cpus: int = cfg.get('cpus', 1)
        self.base_dir: str = cfg.get('base_dir', os.getcwd())
        self.input_dir: str = cfg.get('input_dir', None)
        self.output_dir: str = cfg.get('input_dir', None)
        self.rna_ids: List[str] = cfg.get('rna_ids', None)
        self.dna_ids: List[str] = cfg.get('dna_ids', None)
        # spam errors
        self.validate_inputs()
        # create working directory
        self.work_dir = TemporaryDirectory(dir=self.base_dir)
        self.setup_dirs()
        # get stages-specific configs, maybe get default and update???
        self.dupremover: Dedup = Dedup(cfg.get('dedup', dedup_default_cfg))
    def validate_inputs(self):
        if not self.input_dir:
            exit_with_error('Input directory not specified!')
        if not self.output_dir:
            exit_with_error('Output directory not specified!')
        if (not self.rna_ids) or (not self.dna_ids):
            exit_with_error('Input file ids not specified!')
    def setup_dirs(self):
        self.dedup_dir: str = os.path.join(self.work_dir.name, 'dedup')
        self.rsite_dir: str = os.path.join(self.work_dir.name, 'rsites')
        self.trim_dir : str = os.path.join(self.work_dir.name, 'trim')
        self.hisat_dir: str = os.path.join(self.work_dir.name, 'hisat')
        self.bam_dir: str = os.path.join(self.work_dir.name, 'bam')

    def run(self):
        os.chdir(self.workdir)
        self.dupremover()
        ...
        self.chdir(self.base_dir)
