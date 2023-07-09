import logging
import os

from typing import Any, Dict

from .Base import BaseProcessor
from .default_configs import dedup_default_cfg, rsites_default_cfg,\
    trim_default_cfg, hisat_default_cfg, contacts_default_cfg
from .stages import AlignedToBed, BamFilter, Contacts,\
                    Dedup, Hisat, Rsites, Trim
from ..utils import exit_with_error, make_directory, move_exist_ok

logger = logging.getLogger()


class SingleStageProcessor(BaseProcessor):
    def __init__(self,
                 cfg: Dict[str, Any],
                 stage: str):
        super().__init__(cfg)
        # workdir
        self.stage_dir: str = os.path.join(self.work_dir.name, stage)
        make_directory(self.stage_dir)
        # setup stage
        self.stage: str = stage
        self.setup_stage(cfg)

    def setup_stage(self, cfg: Dict[str, Any]):
        """For selected stage, create corresponding routine object 
        based on default config and user input"""
        if self.stage == 'rsites':
            rsites_default_cfg.update(cfg.get('rsites', {}))
            self.stage_runner = Rsites(rsites_default_cfg,
                                       self.input_dir, self.stage_dir,
                                       self.cpus)
        elif self.stage == 'dedup':
            dedup_default_cfg.update(cfg.get('dedup', {}))
            self.stage_runner = Dedup(dedup_default_cfg,
                                      self.input_dir, self.stage_dir,
                                      self.cpus)
        elif self.stage == 'trim':
            trim_default_cfg.update(cfg.get('trim', {}))
            self.stage_runner = Trim(trim_default_cfg,
                                     self.input_dir, self.stage_dir,
                                     self.cpus)
        elif self.stage == 'hisat':
            hisat_default_cfg.update(cfg.get('hisat', {}))
            self.stage_runner = Hisat(hisat_default_cfg,
                                      self.input_dir, self.stage_dir,
                                      self.cpus)
        elif self.stage == 'bam':
            self.stage_runner = BamFilter(self.input_dir, self.stage_dir,
                                          self.cpus)
        elif self.stage == 'bed':
            self.stage_runner = AlignedToBed(self.input_dir, self.stage_dir,
                                             self.cpus)
        elif self.stage == 'contacts':
            contacts_default_cfg.update(cfg.get('contacts', {}))
            self.stage_runner = Contacts(contacts_default_cfg,
                                         self.input_dir, self.stage_dir,
                                         self.cpus)
        else:
            exit_with_error(f'Unknown stage: {self.stage}!')
    
    def save_outputs(self):
        """copy results to out dir"""
        dest_pth = os.path.join(self.output_dir, self.stage)
        make_directory(dest_pth)
        move_exist_ok(self.stage_dir, dest_pth)
        
    def run(self):
        """Run selected stage and save result"""
        logger.info(f'Started processing of {len(self.rna_ids)} pairs of files.')
        os.chdir(self.work_dir.name)
        self.stage_runner.run(self.dna_ids, self.rna_ids)
        os.chdir(self.base_dir)
        self.save_outputs()
        logger.info('Done.')
