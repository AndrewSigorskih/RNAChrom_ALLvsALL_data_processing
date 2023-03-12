import shutil
import os

from multiprocessing.pool import ThreadPool
from tempfile import NamedTemporaryFile
from typing import Any, Callable, Dict, List

from ..utils import exit_with_error, run_command

class BasicStage:
    def __init__(self,
                 input_dir: str,
                 output_dir: str):
        self.input_dir: str = input_dir
        self.output_dir: str = output_dir

class Dedup(BasicStage):
    def __init__(self, 
                 cfg: Dict[str, Any],
                 *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.tool: str = cfg.get('tool', None)
        self.params: List[str] = cfg.get('params', [])
        self.tool_path: str = cfg.get('tool_path', shutil.which(self.tool))
        self.cpus: int = cfg.get('cpus', None)
        if not self.tool:
            exit_with_error('Deduplication tool not specified!')
        if not self.tool_path:
            exit_with_error(f'Cannot deduce path to {self.tool} executable!')

    def run(self,
            dna_ids: List[str],
            rna_ids: List[str],
            global_cpus: int):
        """Run chosen deduplication tool"""
        if not self.cpus:
            self.cpus = global_cpus
        # get filenames
        dna_files = (filename for filename in os.listdir(self.input_dir) 
                     if any([filename.startswith(id) for id in dna_ids]))
        rna_files = (filename for filename in os.listdir(self.input_dir) 
                     if any([filename.startswith(id) for id in rna_ids]))
        # choose tool to run
        if self.tool == 'skip':
            func: Callable[[str, str], None] = self._copy_files
        elif self.tool == 'fastuniq':
            func = self._run_fastuniq
        else:  # unknown tool
            exit_with_error('Unknown deduplication tool!')
        # run chosen function
        with ThreadPool(self.cpus) as pool:
            pool.starmap(func, zip(dna_files, rna_files))

    def _copy_files(self, dnafilename: str, rnafilename: str):
        infile1 = os.path.join(self.input_dir, dnafilename)
        infile2 = os.path.join(self.input_dir, rnafilename)
        outfile1 = os.path.join(self.output_dir, dnafilename)
        outfile2 = os.path.join(self.output_dir, rnafilename)
        shutil.copy(infile1, outfile1)
        shutil.copy(infile2, outfile2)

    def _run_fastuniq(self,
                      dnafilename: str,
                      rnafilename: str):
        """run fastuniq"""
        infile1 = os.path.join(self.input_dir, dnafilename)
        infile2 = os.path.join(self.input_dir, rnafilename)
        outfile1 = os.path.join(self.output_dir, dnafilename)
        outfile2 = os.path.join(self.output_dir, rnafilename)
        if (infile1.endswith('.gz')) or (infile2.endswith('.gz')):
            exit_with_error('Fastuniq tool does not accept gzipped files!')
        with NamedTemporaryFile(mode='w', delete=False, dir='.') as f:
            tmpfilename = f.name
            print(infile1, infile2, sep='\n', file=f)
        command = (f'{self.tool_path} -i {tmpfilename} '
                    f'-o {outfile1} -p {outfile2}')
        run_command(command, shell=True)
        os.remove(tmpfilename)