import shutil
import os

import concurrent.futures
from multiprocessing.pool import ThreadPool
from multiprocessing import Pool
from typing import Callable, List

from ...utils import exit_with_error

suff_to_filter = ('.unpaired', '.novel_splice')

def find_in_list(id: str, lst: List[str]):
    return next(x for x in lst if x.startswith(id))

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
    
    def run_function_pool(self,
                          func: Callable[[str, str, str, str], int],
                          dna_ids: List[str],
                          rna_ids: List[str]):
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
        # run func in parallel
        #with ThreadPool(self.cpus) as pool:
        with Pool(self.cpus) as pool:
            results = pool.starmap(func, zip(dna_input_files, rna_input_files,
                                   dna_output_files, rna_output_files))
        if any([x != 0 for x in results]):
            msg = f'One or several calls of {func.__name__} returned non-zero process exit code!'
            exit_with_error(msg)
    
    def run_function(self,
                     func: Callable[[str, str, str, str], int],
                     dna_ids: List[str],
                     rna_ids: List[str]):
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
        # run func in parallel
        # try https://docs.python.org/3/library/concurrent.futures.html#concurrent.futures.ThreadPoolExecutor
        with concurrent.futures.ThreadPoolExecutor(max_workers=self.cpus) as executor:
            futures = [executor.submit(func, dna_inp, rna_inp, dna_out, rna_out)
                       for dna_inp, rna_inp, dna_out, rna_out in 
                       zip(dna_input_files, rna_input_files, dna_output_files, rna_output_files)]
            results = [future.result() for future in concurrent.futures.as_completed(futures)]
        if any([x != 0 for x in results]):
            msg = f'One or several calls of {func.__name__} returned non-zero process exit code!'
            exit_with_error(msg)
            
