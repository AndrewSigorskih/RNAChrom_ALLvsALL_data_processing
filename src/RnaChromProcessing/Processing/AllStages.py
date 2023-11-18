import logging
from os import chdir
from pathlib import Path
from typing import List, Set

from .Base import BaseProcessor
from .stages import (
    Align, BamFilter, BamToBed, Contacts, Dedup, Rsites, Trim
)
from ..utils import exit_with_error, move_exist_ok

logger = logging.getLogger()
SUBDIR_LIST = ('dedup', 'rsites', 'trim', 'hisat', 'bam', 'bed', 'contacts')


class AllStagesProcessor(BaseProcessor):
    rsites: Rsites = Rsites()
    dedup: Dedup = Dedup()
    trim: Trim = Trim()
    align: Align = Align()
    bam: BamFilter = BamFilter()
    bed: BamToBed = BamToBed()
    contacts: Contacts = Contacts()

    keep: Set[str] = {'trim', 'bed', 'contacts'}

    # TODO add stats

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        # init stages

    def save_outputs(self, save_all: bool = False):
        """copy everything needed to out dir"""
        save_list = SUBDIR_LIST if save_all else self.keep
        for to_copy in save_list:
            if to_copy not in SUBDIR_LIST:
                logger.warning(f'Unknown directory to copy: {to_copy}. Skipping..')
                continue
            source_pth = Path(self.work_dir.name) / to_copy
            dest_pth = self.output_dir / to_copy
            dest_pth.mkdir()
            move_exist_ok(source_pth, dest_pth)

    def run(self) -> None:
        samples = self.gather_inputs()
        logger.info(f'Started processing {len(samples)} pairs of files.')
        # run processing
        try:
            for stage in (
                self.rsites, self.dedup, self.trim,
                self.align, self.bam, self.bed
            ):
                samples = stage.run(samples)
            self.contacts.run(samples)
        except Exception as _:
            import traceback
            chdir(self.base_dir)
            logger.critical('An error occured during pipeline execution! Trying to save available results..')
            self.save_outputs(save_all=True)
            logger.critical(f'Saved intermediate results to {self.output_dir}.')
            exit_with_error(traceback.format_exc())
        # TODO calculate stats
        # save outputs
        chdir(self.base_dir)
        self.save_outputs()
        logger.info('Done.')
