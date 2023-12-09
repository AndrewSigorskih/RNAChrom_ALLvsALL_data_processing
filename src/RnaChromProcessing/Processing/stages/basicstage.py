import concurrent.futures
import shutil
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, List, Optional

from pydantic import BaseModel, PositiveInt

from ...utils import remove_suffixes, run_command
from ...utils.errors import StageFailedError

logger = logging.getLogger()


@dataclass(frozen=True)
class SamplePair:
    dna_file: Path
    rna_file: Path


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
                    inp_sample: SamplePair,
                    out_sample: SamplePair) -> int:
        """trivially copy files"""
        shutil.copy(inp_sample.dna_file, out_sample.dna_file)
        shutil.copy(inp_sample.rna_file, out_sample.rna_file)
        return 0
    
    def _custom(self,
                inp_sample: SamplePair,
                out_sample: SamplePair) -> int:
        cmd = [
            self.tool_path, inp_sample.dna_file, inp_sample.rna_file,
            out_sample.dna_file, out_sample.rna_file
        ]
        exit_code = run_command(cmd)
        return exit_code
    
    def _make_output_samples(self,
                             input_samples: List[SamplePair],
                             new_suff: Optional[str] = None) -> List[SamplePair]:
        """Create new SamplePairs pointing to self._stage_dir, 
           changing file extensions if needed."""
        if new_suff:
            return [
                SamplePair(
                    self._stage_dir / f'{remove_suffixes(Path(sample.dna_file.name))}.{new_suff}',
                    self._stage_dir / f'{remove_suffixes(Path(sample.rna_file.name))}.{new_suff}'
                ) for sample in input_samples
            ]
        else:
            return [
                SamplePair(
                    self._stage_dir / sample.dna_file.name,
                    self._stage_dir / sample.rna_file.name
                ) for sample in input_samples
            ]

    def run_function(self,
                     func: Callable[[SamplePair, SamplePair], int],
                     inputs: List[SamplePair],
                     outputs: List[SamplePair],
                     require_zero_code: bool = True) -> None:
        if len(inputs) != len(outputs):
            msg = f'Lengths of inputs and outputs lists for {func.__qualname__} do not match!'
            raise StageFailedError(msg)
        
        logger.debug(f'Running function {func.__qualname__} with {self.cpus} threads.')
        with concurrent.futures.ThreadPoolExecutor(max_workers=self.cpus) as executor:
            futures = [
                executor.submit(func, input_sample, output_sample)
                for input_sample, output_sample in zip(inputs, outputs)
            ]
            results = [
                future.result() for future in concurrent.futures.as_completed(futures)
            ]
        if require_zero_code and any(x != 0 for x in results):
            msg = f'One or several calls of {func.__qualname__} returned non-zero process exit code!'
            raise StageFailedError(msg)
