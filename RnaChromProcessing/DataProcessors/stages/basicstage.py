import shutil
import os

from multiprocessing.pool import ThreadPool
from multiprocessing import Pool
from typing import Callable, List

from ...utils import exit_with_error


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
    
    def run_function(self,
                     func: Callable[[str, str, str, str], int],
                     dna_ids: List[str],
                     rna_ids: List[str]):
        # get filenames
        filenames: List[str] = os.listdir(self.input_dir)
        dna_input_files = [os.path.join(self.input_dir, filename) 
                           for filename in filenames
                           if any([filename.startswith(id) for id in dna_ids])]
        rna_input_files = [os.path.join(self.input_dir, filename) 
                           for filename in filenames
                           if any([filename.startswith(id) for id in rna_ids])]
        dna_output_files = [filename.replace(self.input_dir, self.output_dir) 
                            for filename in dna_input_files]
        rna_output_files = [filename.replace(self.input_dir, self.output_dir) 
                            for filename in rna_input_files]
        # run func in parallel
        #with ThreadPool(self.cpus) as pool:
        with Pool(self.cpus) as pool:
            results = pool.starmap(func, zip(dna_input_files, rna_input_files,
                                   dna_output_files, rna_output_files))
        if any([x != 0 for x in results]):
            msg = f'One or several calls of {func.__name__} returned non-zero process exit code!'
            exit_with_error(msg)

            
