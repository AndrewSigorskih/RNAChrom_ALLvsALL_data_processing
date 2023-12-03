from collections import defaultdict
from pathlib import Path
from typing import List, Literal, Optional

import pandas as pd
from pydantic import BaseModel, PositiveInt

from .basicstage import SamplePair
from ...utils import exit_with_error, run_command, run_get_stdout
from ...XRNA.PoolExecutor import PoolExecutor


class StatsCalc(BaseModel):
    cpus: Optional[PositiveInt] = None
    prefix: str = 'stats'
    mode: Literal['skip', 'default', 'full'] = 'default'
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._executor = None
        self._result = defaultdict(dict)
    
    def set_cpus(self, global_cpus: int) -> None:
        self.cpus = self.cpus or global_cpus
        self._executor = PoolExecutor(self.cpus)

    def _count_in_fastqs(self, stage: str, samples: List[SamplePair]) -> None:
        pass

    def _count_in_bams(self):
        pass

    def _count_contacts(self, samples: List[SamplePair]):
        pass

    def run (self, stage: str, samples: List[SamplePair]) -> None:
        '''Estimate surviving read pairs number after sertain stages.'''
        if self.mode == 'skip':
            return
        if not self._executor:
            exit_with_error(
                f'{self.__class__.__name__} was not properly instantiated; '
                f'set_cpus method is required to be called before usage!'
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
            self._count_contacts(samples)
        else:  # skip bam, bed
            return
