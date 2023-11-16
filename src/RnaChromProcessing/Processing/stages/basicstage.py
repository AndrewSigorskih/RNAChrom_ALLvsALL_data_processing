import concurrent.futures
import shutil
import logging
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, List, Optional

from pydantic import BaseModel, PositiveInt

from ...utils.errors import StageFailedError

logger = logging.getLogger()


@dataclass
class SamplePair:
    dna_file: Path
    rna_file: Path

    def set_files(self,
                  new_dna_file: Path,
                  new_rna_file: Path):
        self.dna_file = new_dna_file
        self.rna_file = new_rna_file


class BasicStage(BaseModel):
    cpus: Optional[PositiveInt] = None

    def set_params(self,
                   global_cpus: int,
                   stage_dir: Path) -> None:
        self._stage_dir = stage_dir
        self._stage_dir.mkdir()
        if not self.cpus:
            self.cpus = global_cpus

    def _copy_files(self,
                    dna_in_file: Path,
                    rna_in_file: Path,
                    dna_out_file: Path,
                    rna_out_file: Path) -> int:
        """trivially copy files"""
        shutil.copy(dna_in_file, dna_out_file)
        shutil.copy(rna_in_file, rna_out_file)
        return 0
    
    def _symlink_files(self,
                       dna_in_file: Path,
                       rna_in_file: Path,
                       dna_out_file: Path,
                       rna_out_file: Path) -> int:
        """create symlinks instead of copying"""
        os.symlink(dna_in_file, dna_out_file)
        os.symlink(rna_in_file, rna_out_file)
        return 0
    
    def run_function(self,
                     func: Callable[[str, str, str, str], int],
                     dna_inputs: List[Path],
                     rna_inputs: List[Path],
                     dna_outputs: List[Path],
                     rna_outputs: List[Path],
                     require_zero_code: bool = True) -> None:
        if len(dna_outputs) != len(dna_inputs) != len(rna_outputs) != len(rna_inputs):
            msg = f'Lengths of inputs and outputs lists for {func.__qualname__} do not match!'
            raise StageFailedError(msg)
        
        logger.debug(f'Running function {func.__qualname__} with {self.cpus} threads.')
        with concurrent.futures.ThreadPoolExecutor(max_workers=self.cpus) as executor:
            futures = [
                executor.submit(func, dna_inp, rna_inp, dna_out, rna_out)
                for dna_inp, rna_inp, dna_out, rna_out in 
                zip(dna_inputs, rna_inputs, dna_outputs, rna_outputs)
            ]
            results = [
                future.result() for future in concurrent.futures.as_completed(futures)
            ]
        if require_zero_code and any(x != 0 for x in results):
            msg = f'One or several calls of {func.__qualname__} returned non-zero process exit code!'
            raise StageFailedError(msg)
