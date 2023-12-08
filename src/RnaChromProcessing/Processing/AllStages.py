from logging import getLogger
from os import chdir
from typing import Set

from pydantic import Field

from .Base import BaseProcessor
from .stages import (
    Align, BamFilter, BamToBed, Contacts, Dedup, Rsites, StatsCalc, Trim
)
from ..utils import exit_with_error, move_exist_ok

logger = getLogger()
SUBDIR_LIST = ('rsites', 'dedup', 'trim', 'align', 'bam', 'bed', 'contacts')


class AllStagesProcessor(BaseProcessor):
    rsites: Rsites = Field(default_factory=Rsites)
    dedup: Dedup = Field(default_factory=Dedup)
    trim: Trim = Field(default_factory=Trim)
    align: Align = Field(default_factory=Align)
    bam: BamFilter = Field(default_factory=BamFilter)
    bed: BamToBed = Field(default_factory=BamToBed)
    contacts: Contacts = Field(default_factory=Contacts)

    keep: Set[str] = {'trim', 'bed', 'contacts'}

    stats: StatsCalc = StatsCalc()

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._stage_order = (
            self.rsites, self.dedup, self.trim,
            self.align, self.bam, self.bed, self.contacts
        )
        # init all stages
        for dirname, stage in zip(SUBDIR_LIST, self._stage_order):
            stage.set_params(self.cpus, self._work_pth / dirname)
        self.stats.set_params(self.cpus, self._work_pth / 'stats', self.bam.max_mismatch)

    def save_outputs(self, save_all: bool = False):
        """copy everything needed to out dir"""
        save_list = SUBDIR_LIST if save_all else self.keep
        for to_copy in save_list:
            if to_copy not in SUBDIR_LIST:
                logger.warning(f'Unknown directory to copy: {to_copy}. Skipping..')
                continue
            source_pth = self._work_pth / to_copy
            dest_pth = self.output_dir / to_copy
            dest_pth.mkdir()
            move_exist_ok(source_pth, dest_pth)

    def run(self) -> None:
        samples = self.gather_inputs()
        logger.info(f'Started processing {len(samples)} pairs of files.')
        # run processing
        chdir(self._work_pth)
        try:
            for stagename, stage in zip(SUBDIR_LIST, self._stage_order):
                samples = stage.run(samples)
                self.stats.run(stagename, samples)
        except Exception as _:
            import traceback
            chdir(self.base_dir)
            logger.critical('An error occured during pipeline execution! Trying to save available results..')
            self.save_outputs(save_all=True)
            logger.critical(f'Saved intermediate results to {self.output_dir}.')
            exit_with_error(traceback.format_exc())
        # save outputs
        chdir(self.base_dir)
        self.save_outputs()
        self.stats.save_result(self.output_dir)
        logger.info('Done.')
