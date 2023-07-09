from abc import ABC, abstractmethod
import os
from tempfile import TemporaryDirectory
from typing import Any, Dict, List

from ..utils import exit_with_error, make_directory

class BaseProcessor(ABC):
    def __init__(self,
                 cfg: Dict[str, Any]):
        # get basic parameters
        self.cpus: int = cfg.get('cpus', 1)
        self.base_dir: str = os.path.abspath(cfg.get('base_dir', os.getcwd()))
        self.input_dir: str = cfg.get('input_dir', None)
        self.output_dir: str = cfg.get('output_dir', None)
        self.dna_ids: List[str] = cfg.get('dna_ids', [])
        self.rna_ids: List[str] = cfg.get('rna_ids', [])
        # spam errors
        self._validate_inputs()
        # create working directory
        os.chdir(self.base_dir)
        self.work_dir = TemporaryDirectory(dir=self.base_dir)

    def _validate_inputs(self):
        """Basic input sanity check"""
        # input dir: check name
        if not self.input_dir:
            exit_with_error('Input directory not specified!')
        self.input_dir = os.path.abspath(self.input_dir)
        # output dir: check name and create if needed
        if not self.output_dir:
            exit_with_error('Output directory not specified!')
        self.output_dir = os.path.abspath(self.output_dir)
        if not os.path.exists(self.output_dir):
            make_directory(self.output_dir)
        # other important inputs
        if (not self.rna_ids) or (not self.dna_ids):
            exit_with_error('Input file ids not specified!')
        if len(self.rna_ids) != len(self.dna_ids):
            exit_with_error('DNA and RNA ids arrays have different length!')

    @abstractmethod
    def run(self):
        pass