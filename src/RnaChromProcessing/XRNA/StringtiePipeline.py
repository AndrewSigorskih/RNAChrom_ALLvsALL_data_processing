import shutil
from functools import partial
from pathlib import Path
from typing import List, Optional, Tuple
from typing_extensions import Annotated

import pandas as pd
from pydantic import BaseModel, Field, field_validator

from .AnnotInfo import AnnotInfo
from .PoolExecutor import PoolExecutor
from ..utils import run_command, validate_tool_path

STRINGTIE_STAGES = (
    'stringtie_raw', 'stringtie_merge', 'bed_transforms', 'stringtie_cov', 'xrna'
)


class StringtieTool(BaseModel):
    stringtie_threads: Optional[int] = 1
    tool_path: Annotated[Optional[Path], Field(validate_default=True)] = None
    
    @field_validator('tool_path')
    @classmethod
    def validate_stringtie_path(cls, val: Optional[Path]) -> Path:
        val_fn = partial(validate_tool_path, tool_name='stringtie')
        return val_fn(val)
    

    def run_stringtie(self,
                      inputs: Tuple[Path, Path],
                      out_file: Path) -> None:
        in_bam, gtf_annot = inputs
        cmd = (
            f'{self.tool_path} -o {out_file} -G {gtf_annot} '
            f'-p {self.stringtie_threads} {in_bam}'
        )
        return_code = run_command(cmd, shell=True)
        return return_code


    def run_stringtie_merge(self,
                            assembly_list: Path,
                            output_file: Path) -> None:
        cmd = (
            f'{self.tool_path} --merge -p {self.stringtie_threads} '
            f'-o {output_file} {assembly_list}'
        )
        return_code = run_command(cmd, shell=True)
        return return_code


class StringtiePipeline:
    def __init__(self,
                 work_pth: Path,
                 executor: PoolExecutor,
                 stringtie: StringtieTool):
        self.executor = executor
        self.stringtie_tool = stringtie
        for subdir_name in STRINGTIE_STAGES:
            subdir = work_pth / subdir_name
            setattr(self, subdir_name, subdir)
            subdir.mkdir()
    
    def run_stringtie(self,
                      input_bams: List[Path],
                      gtf_annot: Path) -> List[Path]:
        inputs = list(
            zip(
                input_bams,
                (gtf_annot for _ in range(len(input_bams)))
            )
        )
        outputs = [
            self.stringtie_raw / f'{file.stem}.gtf' for file in input_bams
        ]
        self.executor.run_function(
            self.stringtie_tool.run_stringtie,
            inputs, outputs
        )
        return outputs
    
    def run_stringtie_merge(self,
                            raw_gtfs: List[Path]) -> None:
        assembly_lst: Path = self.stringtie_merge / 'assembly.lst'
        with open(assembly_lst, 'w') as f:
            print(*raw_gtfs, sep='\n', file=f)
        merged_file: Path = self.stringtie_merge / 'merged.gtf'
        self.stringtie_tool.run_stringtie_merge(assembly_lst, merged_file)
        assembly_lst.unlink()

    def handle_intervals(self,
                         bed_annot: Path) -> None:
        raw_gtf: Path = self.stringtie_merge / 'merged.gtf'
        raw_bed: Path = self.bed_transforms / 'raw.bed'
        nonoverlap_bed: Path = self.bed_transforms / 'non-overlap.bed'
        counts_bed: Path = self.bed_transforms / 'overlaps.bed'
        true_x_bed: Path = self.bed_transforms / 'xrnas.bed'
        # gtf -> sorted bed
        tab = pd.read_csv(
            raw_gtf, sep='\t', header=None, skiprows=2, usecols=[0,2,3,4,5,6,7]
        )
        tab = tab[tab[2] == 'transcript']
        tab = tab[[0,3,4,7,5,6]].sort_values(by=[0,3], inplace=False)
        tab.to_csv(raw_bed, sep='\t', index=False, header=False)
        # sorted bed -> merge overlaps
        cmd = f'bedtools merge -s -c 6 -o distinct -i {raw_bed} > {nonoverlap_bed}'
        run_command(cmd, shell=True)
        # merge overlaps -> intersections w/ annotation
        cmd = (
            f'bedtools coverage -a {nonoverlap_bed} -b {bed_annot} '
            f'-s -counts > {counts_bed}'
        )
        run_command(cmd, shell=True)
        # intersections w/ annotation -> filter out
        

        



    def run(self,
            input_bams: List[Path],
            annot_info: AnnotInfo) -> None:
        raw_gtfs = self.run_stringtie(input_bams, annot_info.gtf_annotation)
        self.run_stringtie_merge(raw_gtfs)
        self.handle_intervals(annot_info.annot_bed)
