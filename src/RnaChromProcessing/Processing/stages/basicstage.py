import concurrent.futures
import shutil
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable,  List, Optional

from pydantic import BaseModel, PositiveInt

from ...utils import remove_suffixes
from ...utils.errors import StageFailedError

logger = logging.getLogger()


def _all_equal(iterable: Iterable[int]) -> bool:
    "return True if all elements in iterable are equal"
    iterator = iter(iterable)
    try:
        first = next(iterator)
    except StopIteration:
        return True
    return all(first == x for x in iterator)


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
                    dna_in_file: Path,
                    rna_in_file: Path,
                    dna_out_file: Path,
                    rna_out_file: Path) -> int:
        """trivially copy files"""
        shutil.copy(dna_in_file, dna_out_file)
        shutil.copy(rna_in_file, rna_out_file)
        return 0
    
    def _make_output_samples(self,
                             input_samples: List[SamplePair],
                             new_suff: Optional[str] = None) -> List[SamplePair]:
        """Create new SamplePairs pointing to self._stage_dir, 
           changing file extensions if needed."""
        if new_suff:
            return [
                SamplePair(
                    self._stage_dir / f'{remove_suffixes(sample.dna_file.name)}.{new_suff}',
                    self._stage_dir / f'{remove_suffixes(sample.rna_file.name)}.{new_suff}',
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
                     func: Callable[[Path, Path, Path, Path], int],
                     dna_inputs: List[Path],
                     rna_inputs: List[Path],
                     dna_outputs: List[Path],
                     rna_outputs: List[Path],
                     require_zero_code: bool = True) -> None:
        if not _all_equal((len(dna_inputs), len(dna_outputs), len(rna_outputs), len(rna_inputs))):
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
