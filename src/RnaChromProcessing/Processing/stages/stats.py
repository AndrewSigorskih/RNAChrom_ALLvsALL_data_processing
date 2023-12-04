from collections import defaultdict
from itertools import repeat
from pathlib import Path
from typing import List, Literal, Optional

import pandas as pd
from pydantic import BaseModel, PositiveInt

from .basicstage import SamplePair
from ...utils import exit_with_error, remove_suffixes, run_command, run_get_stdout
from ...XRNA.PoolExecutor import PoolExecutor


class StatsCalc(BaseModel):
    cpus: Optional[PositiveInt] = None
    prefix: str = 'stats'
    mode: Literal['skip', 'default', 'full'] = 'default'
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._executor = None
        self._result = defaultdict(dict)
    
    def set_params(self, global_cpus: int, stage_dir: Path) -> None:
        self.cpus = self.cpus or global_cpus
        self._executor = PoolExecutor(self.cpus)
        self._stage_dir = stage_dir
        self._stage_dir.mkdir()

    def _count_in_fastq_pair(self, sample: SamplePair, stage: str) -> None:
        """Counts matching IDs in pair of FASTQ files.
        Supported IDs formats:
            * @IDXXX <-> @IDXXX\n
            * @FILE1.IDXXX <-> @FILE2.IDXXX\n
            * @FILE.IDXXX <-> @FILE.IDXXX"""
        rna_id = remove_suffixes(sample.rna_file.name)
        dna_tmp_file = self._stage_dir / sample.dna_file.with_suffix('tmp').name
        rna_tmp_file = self._stage_dir / sample.rna_file.with_suffix('tmp').name
        for infile, outfile in zip(
            (sample.dna_file, sample.rna_file),
            (dna_tmp_file, rna_tmp_file)
        ):
            cat = 'zcat' if infile.suffix == '.gz' else 'cat'
            cmd = (
                f'{cat} {infile} | sed -n  "1~4p" | sed "s/@//g"  | '
                f'cut -d " " -f 1 | cut -d "." -f 2 > {outfile}'
            )
            run_command(cmd, shell=True)
        count_cmd = f"awk 'a[$0]++' {dna_tmp_file} {rna_tmp_file} | wc -l"
        self._result[stage][rna_id] = int(run_get_stdout(count_cmd, shell=True))
        dna_tmp_file.unlink()
        rna_tmp_file.unlink()

    def _count_in_contacts_file(self, sample: SamplePair, stage: str) -> None:
        rna_id = remove_suffixes(sample.rna_file.name)
        count_cmd = f'wc -l < {sample.rna_file}'
        self._result[stage][rna_id] = int(run_get_stdout(count_cmd, shell=True)) - 1
    
    def _count_in_bam_pair(self, sample: SamplePair, mode: str) -> None:
        """Counts matching IDs in pair od BAM files.
        Supported IDs formats:
            * @IDXXX <-> @IDXXX\n
            * @FILE1.IDXXX <-> @FILE2.IDXXX\n
            * @FILE.IDXXX <-> @FILE.IDXXX
        :param: mode either 'mapped' or 'mapped2mism'"""
        pass

    def _count_in_fastqs(self, stage: str, samples: List[SamplePair]) -> None:
        self._executor.run_function(
            self._count_in_fastq_pair, samples, list(repeat(stage, len(samples))),
            require_zero_code=False
        )

    def _count_in_bams(self):
        pass

    def _count_in_contacts(self, stage: str, samples: List[SamplePair]):
        self._executor.run_function(
            self._count_in_contacts_file, samples, list(repeat(stage, len(samples))),
            require_zero_code=False
        )

    def run (self, stage: str, samples: List[SamplePair]) -> None:
        '''Estimate surviving read pairs number after sertain stages.'''
        if self.mode == 'skip':
            return
        if not self._executor or not self._stage_dir:
            exit_with_error(
                f'{self.__class__.__name__} was not properly instantiated; '
                'set_params method is required to be called before usage!'
            )
        # actual logic
        if stage in ('rsites', 'dedup', 'trim'):
            # fastq
            self._count_in_fastqs(stage, samples)
        elif stage == 'align':
            # only if 'full' mode
            if self.mode != 'full':
                return
            pass  #TODO
        elif stage == 'contacts':
            self._count_contacts(stage, samples)
        else:  # skip bam, bed
            return
    
    def save_result(self, output_dir: Path) -> None:
        output_name = output_dir / f'{self.prefix}.tsv'
        result = pd.DataFrame.from_dict(self.result).sort_index()
        result.to_csv(output_name, sep='\t', index=True, header=True)
