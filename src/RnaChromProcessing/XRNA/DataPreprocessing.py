import shutil
from os import listdir
from pathlib import Path
from typing import List, Optional
from typing_extensions import Annotated

from pydantic import BaseModel, Field, FieldValidationInfo, field_validator

from .PoolExecutor import PoolExecutor
from ..utils import exit_with_error, find_in_list, run_command

PREPROCESS_STAGES = (
    'filter_contacts', 'filter_fastq', 'revc_fastq', 'align', 'sorted_bams'
)

class HisatTool(BaseModel):
    genome_path: Path
    known_splice: Path
    tool_path: Annotated[Optional[Path], Field(validate_default=True)] = None
    hisat_threads: int = 1

    @field_validator('tool_path')
    @classmethod
    def handle_tool_path(cls, val: Optional[str], info: FieldValidationInfo) -> Path:
        val = val or shutil.which('hisat2')
        if not val:
            exit_with_error('Cannot deduce path to hisat2 executable!')
        return Path(val)

    def run(self,
            in_file: Path,
            out_file: Path) -> int:
        """run hisat2"""
        cmd = (
            f'{self.tool_path} -x {self.genome_path} -p {self.hisat_threads} -k 100 '
            f'--no-softclip --known-splicesite-infile {self.known_splice} --dta-cufflinks '
            f'-U {in_file} | samtools view -F 4 -bSh > {out_file}'
        )
        return_code = run_command(cmd, shell=True)
        return return_code


class PreprocessingPipeline:
    def __init__(self,
                 work_pth: Path,
                 executor: PoolExecutor,
                 hisat: HisatTool,
                 file_ids: List[str]):
        self.executor = executor
        self.hisat_tool = hisat
        self.file_ids = file_ids
        for subdir_name in PREPROCESS_STAGES:
            subdir = work_pth / subdir_name
            subdir.mkdir()
            setattr(self, subdir_name, subdir)
    
    def _revc_fastq(self,
                    in_file: Path,
                    out_file: Path) -> int:
        pass

    def run_revc_fastq(self):
        pass

    def run_hisat(self) -> None:
        fnames = listdir(self.revc_fastq)
        inputs = (
            self.revc_fastq / find_in_list(name, fnames) 
            for name in self.file_ids
        )
        outputs = [self.align / f'{name}.bam' for name in self.file_ids]
        self.executor.run_function(
            self.hisat_tool.run,
            inputs, outputs
        )

    def _sort_bams(self,
                   in_file: Path,
                   out_file: Path) -> int:
        cmd = f'samtools sort {in_file} -o {out_file}'
        return_code = run_command(cmd, shell=True)
        return return_code
    
    def run_sort_bams(self) -> None:
        fnames = listdir(self.align)
        inputs = [self.align / name for name in fnames]
        outputs = [self.sorted_bams / name for name in fnames]
        self.executor.run_function(
            self._sort_bams,
            inputs, outputs
        )

    def run(self) -> None:
        pass
