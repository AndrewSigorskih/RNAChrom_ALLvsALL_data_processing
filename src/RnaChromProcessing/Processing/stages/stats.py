from collections import defaultdict
from pathlib import Path
from logging import getLogger
from typing import List, Literal, NamedTuple, Optional

import pandas as pd
from pydantic import BaseModel, PositiveInt

from .basicstage import SamplePair
from ...utils import exit_with_error, remove_suffixes, run_get_stdout
from ...utils.PoolExecutor import PoolExecutor


logger = getLogger()

class TablePos(NamedTuple):
    index: int
    column: str

def _get_positions(stage: str, num: int) -> List[TablePos]:
    return [TablePos(i, stage) for i in range(num)]


class StatsCalc(BaseModel):
    cpus: Optional[PositiveInt] = None
    prefix: str = 'stats'
    mode: Literal['skip', 'default'] = 'default'
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._executor = None
        self._result = defaultdict(dict, RNA_ID={}, DNA_ID={})
    
    def set_params(self, global_cpus: int, stage_dir: Path, mism_num: int) -> None:
        self.cpus = self.cpus or global_cpus
        self._mism_num = mism_num
        self._executor = PoolExecutor(self.cpus)
        self._stage_dir = stage_dir
        self._stage_dir.mkdir()

    def _check_ready(self):
        if not self._executor or not self._stage_dir:
            exit_with_error(
                f'{self.__class__.__name__} was not properly instantiated; '
                'set_params method is required to be called before usage!'
            )

    def _count_in_fastq_pair(self, sample: SamplePair, pos: TablePos) -> None:
        """Counts number of reads in fastq files. Performs basic sanity check."""
        res = {}
        rna_name, dna_name = sample.rna_file.name, sample.dna_file.name
        for infile in (sample.dna_file, sample.rna_file):
            cat = 'zcat' if infile.suffix == '.gz' else 'cat'
            cmd = (
                f'{cat} {infile} | wc -l'
            )
            res[infile.name] = int(run_get_stdout(cmd, shell=True))
        if res[rna_name] != res[dna_name]:
            msg = (
                f'Fastq files for the {pos.column} had different line counts: '
                f'{rna_name}: {res[rna_name]}, {dna_name}:{res[dna_name]}. '
                'Any further analysis may be compromised.'
            )
            logger.warning(msg)
        read_num, rem = divmod(res[rna_name], 4)
        if rem != 0:
            msg = (
                f'Fastq files for the {pos.column} had malformed structure! '
                f'Number of lines is not divisible by 4: {res[rna_name]}. '
                'Any further analysis may be compromised.'
            )
            logger.warning(msg)
        self._result[pos.column][pos.index] = read_num

    def _count_in_bam_pair(self, sample: SamplePair, pos: TablePos) -> None:
        """Counts numbers of reads in bam files separately."""
        for infile, data_source in zip(
            (sample.rna_file, sample.dna_file), ('RNA', 'DNA')
        ):
            column = f'{pos.column}_{data_source}'
            cmd = f'samtools view -c {infile}'
            self._result[column][pos.index] = int(run_get_stdout(cmd, shell=True))

    def _count_in_bed_pair(self, sample: SamplePair, pos: TablePos):
        """Counts numbers of reads in bed files separately."""
        for infile, data_source in zip(
            (sample.rna_file, sample.dna_file), ('RNA', 'DNA')
        ):
            column = f'{pos.column}_{data_source}'
            cmd = f'wc -l {infile}'
            self._result[column][pos.index] = int(run_get_stdout(cmd, shell=True))

    def _count_in_contacts_file(self, sample: SamplePair, pos: TablePos) -> None:
        # meta info
        self._result['RNA_ID'][pos.index] = remove_suffixes(sample.rna_file.name)
        self._result['DNA_ID'][pos.index] = remove_suffixes(sample.dna_file.name)
        # read contacts num
        cmd = f'wc -l < {sample.rna_file}'
        self._result[pos.column][pos.index] = int(run_get_stdout(cmd, shell=True)) - 1

    def run (self, stage: str, samples: List[SamplePair]) -> None:
        '''Estimate surviving read pairs number after sertain stages.'''
        if self.mode == 'skip':
            return
        self._check_ready()
        # actual logic
        if stage in ('rsites', 'dedup', 'trim'):
            func = self._count_in_fastq_pair
        elif stage == 'align':
            stage = 'mapped'
            func = self._count_in_bam_pair
        elif stage == 'bam':
            stage = f'mapped_{self._mism_num}mism'
            func = self._count_in_bam_pair
        elif stage == 'bed':
            stage = 'mapped_unique'
            func = self._count_in_bed_pair
        elif stage == 'contacts':
            func = self._count_in_contacts_file
        else:
            return
        # estimate
        positions = _get_positions(stage, len(samples))
        self._executor.run_function(
            func, samples, positions, require_zero_code=False
        )

    def save_result(self, output_dir: Path) -> None:
        output_name = output_dir / f'{self.prefix}.tsv'
        result = pd.DataFrame.from_dict(self._result).sort_index()
        result.to_csv(output_name, sep='\t', index=True, header=True)
