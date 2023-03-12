import os
from tempfile import TemporaryDirectory
from typing import Any, Dict, List

from stages import Dedup, Rsites, Trim, Hisat
from ..utils import exit_with_error, make_directory, run_command
from ..utils import dedup_default_cfg, rsites_default_cfg, \
                    trim_default_cfg, hisat_default_cfg

class BaseProcessor:
    def __init__(self, cfg: Dict[str, Any]):
        # get basic parameters
        self.cpus: int = cfg.get('cpus', 1)
        self.base_dir: str = cfg.get('base_dir', os.getcwd())
        self.input_dir: str = cfg.get('input_dir', None)
        self.output_dir: str = cfg.get('input_dir', None)
        self.dna_ids: List[str] = cfg.get('dna_ids', [])
        self.rna_ids: List[str] = cfg.get('rna_ids', [])
        self.keep: List[str] = cfg.get('keep', [])
        # spam errors
        self.validate_inputs()
        # create working directory
        self.work_dir = TemporaryDirectory(dir=self.base_dir)
        self.setup_dirs()
        # get stages-specific configs, maybe get default and update???
        self.stages_from_cfgs(cfg)

    def validate_inputs(self):
        if not self.input_dir:
            exit_with_error('Input directory not specified!')
        if not self.output_dir:
            exit_with_error('Output directory not specified!')
        if (not self.rna_ids) or (not self.dna_ids):
            exit_with_error('Input file ids not specified!')
        if not self.keep:
            exit_with_error('Empty keep list!')
    
    def setup_dirs(self):
        self.dedup_dir: str = os.path.join(self.work_dir.name, 'dedup')
        self.rsite_dir: str = os.path.join(self.work_dir.name, 'rsites')
        self.trim_dir : str = os.path.join(self.work_dir.name, 'trim')
        self.hisat_dir: str = os.path.join(self.work_dir.name, 'hisat')
        self.bam_dir: str = os.path.join(self.work_dir.name, 'bam')
        self.bed_dir: str = os.path.join(self.work_dir.name, 'bed')
        self.contacts_dir: str = os.path.join(self.work_dir.name, 'contacts')
        for dirname in (
            self.dedup_dir, self.rsite_dir, self.trim_dir,
            self.hisat_dir, self.bam_dir, self.bed_dir, self.contacts_dir
        ):
            make_directory(dirname)
    
    def stages_from_cfgs(self, cfg: Dict[str, Any]):
        dedup_default_cfg.update(cfg.get('dedup', {}))
        self.dupremover: Dedup = Dedup(dedup_default_cfg,
                                       self.input_dir, self.dedup_dir,
                                       self.cpus)
        rsites_default_cfg.update(cfg.get('rsites', {}))
        self.rsitefilter: Rsites = Rsites(rsites_default_cfg,
                                          self.dedup_dir, self.rsite_dir,
                                          self.cpus)
        trim_default_cfg.update(cfg.get('trim', {}))
        self.trimmer: Trim = Trim(trim_default_cfg,
                                  self.rsite_dir, self.trim_dir,
                                  self.cpus)
        hisat_default_cfg.update(cfg.get('hisat', {}))
        self.aligner: Hisat = Hisat(hisat_default_cfg,
                                    self.trim_dir, self.hisat_dir,
                                    self.cpus)
        # TODO add two more stages

    def run(self):
        os.chdir(self.workdir)
        self.dupremover.run(self.dna_ids, self.rna_ids)
        self.rsitefilter.run(self.dna_ids, self.rna_ids)
        self.trimmer.run(self.dna_ids, self.rna_ids)
        self.aligner.run(self.dna_ids, self.rna_ids)
        # filter bam
        # contacts
        ...
        self.chdir(self.base_dir)
        # copy everythong needed to out dir
