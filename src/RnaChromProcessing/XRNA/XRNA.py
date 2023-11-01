import logging
import os
from pathlib import Path
from tempfile import TemporaryDirectory, mkdtemp
from typing import Optional, Set

from pydantic import BaseModel, Extra, PositiveInt, field_validator

from .AnnotInfo import AnnotInfo
from .DataPreprocessing import HisatTool, PreprocessingPipeline
from .PoolExecutor import PoolExecutor
from .StringtiePipeline import StringtieTool, StringtiePipeline
from ..utils import exit_with_error

logger = logging.getLogger()
logging.getLogger("matplotlib").setLevel(logging.WARNING)


class XRNAProcessor(BaseModel, extra=Extra.allow):
    bed_input_dir: Path
    fq_input_dir: Path
    output_dir: Path
    base_dir: Path = Path('.').resolve()
    cpus: PositiveInt = 1

    rna_ids: Set[str]
    annotation: AnnotInfo

    hisat: HisatTool
    stringtie: StringtieTool

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        # spam errors
        self._validate_inputs()
        # working dir
        os.chdir(self.base_dir)
        #self.work_dir = TemporaryDirectory(dir=self.base_dir)
        self.work_dir = mkdtemp(dir=self.base_dir)
        # other members
        #work_pth = Path(self.work_dir.name)
        work_pth = Path(self.work_dir)
        executor = PoolExecutor(self.cpus)
        self.annotation.prepare_annotation(work_pth)
        self.preprocessing = PreprocessingPipeline(
            work_pth, executor, self.hisat
        )
        self.pipeline = StringtiePipeline(
            work_pth, executor, self.stringtie
        )

    @field_validator('base_dir', 'fq_input_dir', 'bed_input_dir', 'output_dir')
    @classmethod
    def resolve_path(cls, pth: str) -> Path:
        return Path(pth).resolve()
    
    def _validate_inputs(self):
        """Basic input sanity check"""
        if not self.base_dir.exists():
            self.base_dir.mkdir(parents=True)
        if not self.bed_input_dir.exists():
            exit_with_error('Bed input directory does not exist!')
        if not self.fq_input_dir.exists():
            exit_with_error('Fastq input directory does not exist!')
        if not self.output_dir.exists():
            self.output_dir.mkdir(parents=True)   

    def run(self):
        prepared_bams = self.preprocessing.run(
            self.rna_ids, self.bed_input_dir, self.fq_input_dir, self.annotation
        )
        print(f'prepared BAM files: {prepared_bams}')
        self.pipeline.run(prepared_bams, self.annotation)
        logger.info('Done.')
        
        # stringtie pipeline:
            # save outputs
