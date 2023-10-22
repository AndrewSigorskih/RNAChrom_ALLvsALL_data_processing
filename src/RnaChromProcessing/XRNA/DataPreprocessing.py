import logging
import shutil
from collections import defaultdict
from functools import partial
from os import listdir
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple, Set
from typing_extensions import Annotated

import pandas as pd
from pydantic import BaseModel, Field, field_validator

from .AnnotInfo import AnnotInfo, SampleInfo
from .PoolExecutor import PoolExecutor
from ..utils import (
    find_in_list, run_command, validate_tool_path
)
logger = logging.getLogger()

PREPROCESS_STAGES = (
    'filter_bed', 'filter_fastq', 'revc_fastq', 'align', 'merge_bams', 'sort_bams'
)


def _update_samples_files(file_dir: Path,
                          field_name: str,
                          samples: List[SampleInfo]) -> None:
    files = listdir(file_dir)
    for sample in samples:
        file_path = find_in_list(sample.sample_id, files)
        if not file_path:
            logging.warning(f'Could not find file {sample.sample_id} in {file_dir}!')
            continue
        sample.update_field(field_name, file_path)


def _filter_bed(annot_bed: Path,
                sample: SampleInfo) -> int:
    # require same or opposite strandness
    strand_flag = '-s' if sample.true_strand else '-S'
    counts_file = sample.lst_file.with_suffix('.counts')
    coverage_cmd = (
        f'bedtools coverage -a {sample.bed_file} -b {annot_bed} '
        f'{strand_flag} -counts > {counts_file}'
    )
    return_code = run_command(coverage_cmd, shell=True)
    if return_code != 0:
        return return_code # will be managed by PoolExecutor
    # TODO process by chunks?
    result = pd.read_csv(counts_file, sep='\t', header=None, usecols=[3, 6])
    result = result[result[6] == 0][3]
    #result = '@' + result.apply(str) # for grep!
    result.to_csv(sample.lst_file, header=False, index=False)
    counts_file.unlink()

# https://github.com/lh3/seqtk
def _filter_fq_by_ids(inputs: Tuple[Path, Path],
                      out_file: Path) -> int:
    in_fq, in_lst = inputs
    cmd = f'seqtk subseq {in_fq} {in_lst} > {out_file}'
    return_code = run_command(cmd, shell=True)
    return return_code


def _revc_fastq(sample: SampleInfo,
                out_file: Path) -> int:
    if sample.true_strand:
        shutil.move(sample.fq_file, out_file)
        return_code = 0
    else:
        cmd = f'seqtk seq -r {sample.fq_file} > {out_file}'
        return_code = run_command(cmd, shell=True)
    sample.update_field('fq_file', out_file)
    return return_code

def _merge_bams(in_files: Iterable[Path],
                out_file: Path):
    cmd = (
        f'samtools merge {out_file} '
        f'{"".join(file for file in in_files)}'
    )
    return_code = run_command(cmd)
    return return_code


def _sort_bams(in_file: Path,
               out_file: Path) -> int:
    cmd = ['samtools', 'sort', str(in_file), '-o', str(out_file)]
    return_code = run_command(cmd)
    return return_code


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
                 hisat: HisatTool):
        self.executor = executor
        self.hisat_tool = hisat
        for subdir_name in PREPROCESS_STAGES:
            subdir = work_pth / subdir_name
            setattr(self, subdir_name, subdir)
            subdir.mkdir()

    def run_filter_bed(self,
                       annot_bed: Path,
                       samples_list: List[SampleInfo]) -> None:
        for sample in samples_list:
            sample.update_field(
                'lst_file',
                self.filter_bed / f'{sample.sample_id}.lst'
            )
        annot_inputs = [annot_bed for _ in range(samples_list)]
        self.executor.run_function(
            _filter_bed,
            annot_inputs, samples_list
        )

    def run_filter_fastq(self,
                         samples_list: List[SampleInfo]) -> None:
        inputs = [
            (sample.fq_file, sample.lst_file) for sample in samples_list
        ]
        outputs = [
            self.filter_fastq / f'{sample.sample_id}.fastq' for sample in samples_list
        ]
        self.executor.run_function(
            _filter_fq_by_ids,
            inputs, outputs
        )
        for sample, output in zip(samples_list, outputs):
            sample.update_field('fq_file', output)

    def run_revc_fastq(self,
                       samples_list: List[SampleInfo]) -> None:
        outputs = [
            self.revc_fastq / f'{sample.sample_id}.fastq' for sample in samples_list
        ]
        self.executor.run_function(
            _revc_fastq,
            samples_list, outputs
        )

    def run_hisat(self,
                  samples_list: List[SampleInfo]) -> None:

        outputs = [
            self.align / f'{sample.sample_id}.bam' for sample in samples_list
        ]
        self.executor.run_function(
            self.hisat_tool.run,
            [sample.fq_file for sample in samples_list],
            outputs
        )
        for sample, output in zip(samples_list, outputs):
            sample.update_field('bam_file', output)

    def run_merge_bams(self,
                       replics_dct: Dict[str, List[SampleInfo]]) -> None:
        # maybe too overcomplicated, ensuring same order:
        inputs, outputs = [], []
        for group, samples in replics_dct.items():
            outputs.append(self.merge_bams / f'{group}.bam')
            inputs.append(sample.bam_file for sample in samples)
        self.executor.run_function(
            _merge_bams,
            inputs, outputs,
        )

    def run_sort_bams(self,
                      replics_dct: Dict[str, List[SampleInfo]]) -> None:
        inputs = [self.merge_bams / f'{group}.bam' for group in replics_dct]
        outputs = [self.sort_bams / f'{group}.bam' for group in replics_dct]
        self.executor.run_function(
            _sort_bams,
            inputs, outputs
        )
        self.results = outputs

    def run(self,
            file_ids: Set[str],
            bed_in_dir: Path,
            fq_in_dir: Path,
            annot_info: AnnotInfo) -> List[Path]:
        # here we do not check for inconsistency between fq and bed files
        # since both these dirs came from rnachromprocessing pipeline run
        samples_list = [
            x for x in annot_info.samples if x.sample_id in file_ids          
        ]
        _update_samples_files(bed_in_dir, 'bed_file', samples_list)
        _update_samples_files(fq_in_dir, 'fq_file', samples_list)
        samples_list = [item for item in samples_list if item.inputs_ok()]
        logger.info(f'Started processing {len(samples_list)} files.')

        self.run_filter_bed(annot_info.annot_bed, samples_list)
        self.run_filter_fastq(samples_list)
        self.run_revc_fastq(samples_list)
        self.run_hisat(samples_list)
        # time to merge biological replicas
        replics_dct = defaultdict(list)
        for sample in samples_list:
            replics_dct[sample.sample_group].append(sample)
        self.run_merge_bams(replics_dct)
        self.run_sort_bams(replics_dct)
        return self.results
