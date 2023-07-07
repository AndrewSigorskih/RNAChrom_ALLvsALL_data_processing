import concurrent.futures
import shutil
import logging
import os

from typing import Callable, List

from ...utils import find_in_list

suff_to_filter = ('.unpaired', '.novel_splice')
logger = logging.getLogger()

class StageFailedError(RuntimeError):
    pass

class BasicStage:
    def __init__(self,
                 input_dir: str,
                 output_dir: str,
                 global_cpus: int):
        self.cpus: int = global_cpus # can be changed in subclasses 
        self.input_dir: str = input_dir
        self.output_dir: str = output_dir

    def _copy_files(self,
                    dna_in_file: str,
                    rna_in_file: str,
                    dna_out_file: str,
                    rna_out_file: str) -> int:
        """trivially copy files"""
        shutil.copy(dna_in_file, dna_out_file)
        shutil.copy(rna_in_file, rna_out_file)
        return 0
    
    def _symlink_files(self,
                    dna_in_file: str,
                    rna_in_file: str,
                    dna_out_file: str,
                    rna_out_file: str) -> int:
        """create symlinks instead of copying"""
        os.symlink(dna_in_file, dna_out_file)
        os.symlink(rna_in_file, rna_out_file)
        return 0
    
    def run_function(self,
                     func: Callable[[str, str, str, str], int],
                     dna_ids: List[str],
                     rna_ids: List[str],
                     require_zero_code: bool = True):
        # get filenames
        filenames: List[str] = [x for x in os.listdir(self.input_dir)
                                if not any([x.endswith(y) for y in suff_to_filter])]
        dna_files = [find_in_list(id, filenames) for id in dna_ids]
        rna_files = [find_in_list(id, filenames) for id in rna_ids]
        dna_input_files = [os.path.join(self.input_dir, filename)
                           for filename in dna_files]
        rna_input_files = [os.path.join(self.input_dir, filename)
                           for filename in rna_files]
        dna_output_files = [os.path.join(self.output_dir, filename)
                           for filename in dna_files]
        rna_output_files = [os.path.join(self.output_dir, filename)
                           for filename in rna_files]
        # run func in parallel (threads)
        logger.debug(f'Running function {func.__qualname__} with {self.cpus} threads.')
        with concurrent.futures.ThreadPoolExecutor(max_workers=self.cpus) as executor:
            futures = [executor.submit(func, dna_inp, rna_inp, dna_out, rna_out)
                       for dna_inp, rna_inp, dna_out, rna_out in 
                       zip(dna_input_files, rna_input_files, dna_output_files, rna_output_files)]
            results = [future.result() for future in concurrent.futures.as_completed(futures)]
        if require_zero_code and any([x != 0 for x in results]):
            msg = f'One or several calls of {func.__name__} returned non-zero process exit code!'
            #exit_with_error(msg)
            raise StageFailedError(msg)
