import logging
import shutil
from functools import partial
from os import listdir
from pathlib import Path
from typing import List, Optional
from typing_extensions import Annotated

import pandas as pd
from pydantic import BaseModel, Field, field_validator

from .PoolExecutor import PoolExecutor
from ..utils import (
    exit_with_error, find_in_list, run_command, validate_tool_path
)

PREPROCESS_STAGES = (
    'filter_contacts', 'filter_fastq', 'revc_fastq', 'align', 'sorted_bams'
)

logger = logging.getLogger()


class HisatTool(BaseModel):
    genome_path: Path
    known_splice: Path
    tool_path: Annotated[Optional[Path], Field(validate_default=True)] = None
    hisat_threads: int = 1

    @field_validator('tool_path')
    @classmethod
    def validate_hisat2_path(cls, val: Optional[Path]) -> Path:
        val_fn = partial(validate_tool_path, tool_name='hisat2')
        return val_fn(val)

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
            setattr(self, subdir_name, subdir)
            subdir.mkdir()

    def _rev_compl(self,
                    in_file: Path,
                    out_file: Path) -> int:
        pass

    def run_revc_fastq(self,
                       strand_info_pth: Path) -> None:
        strand_info = pd.read_table(strand_info_pth, sep='\t',
                                    index_col=[0,1])
        if ((unknown_num := (strand_info.strand=='UNKNOWN').sum()) > 0):
            logger.info(f'Discarding {unknown_num} files with unknown true strand..')
        same_strand = strand_info[strand_info.strand == 'SAME'].index.get_level_values(1)
        anti_strand = strand_info[strand_info.strand == 'ANTI'].index.get_level_values(1)

        same_strand = same_strand[same_strand.isin(self.file_ids)]
        anti_strand = anti_strand[anti_strand.isin(self.file_ids)]

        fnames = listdir(self.filter_fastq)
        files_to_move = [find_in_list(name, fnames) for name in same_strand]
        files_to_revc = [find_in_list(name, fnames) for name in anti_strand]

        logger.debug(f'{len(files_to_move)} files with correct strand.')
        self.executor.run_function(
            shutil.move,
            [self.filter_fastq / name for name in files_to_move],
            [self.revc_fastq / name for name in files_to_move],
            require_zero_code=False
        )

        logger.debug(f'{len(files_to_revc)} files with incorrect strand, applying reverse-complement..')
        self.executor.run_function(
            self._rev_compl,
            [self.filter_fastq / name for name in files_to_revc],
            [self.revc_fastq / name for name in files_to_revc]
        )

    def run_hisat(self) -> None:
        fnames = listdir(self.revc_fastq)
        inputs = [
            self.revc_fastq / find_in_list(name, fnames)
            for name in self.file_ids
        ]
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
