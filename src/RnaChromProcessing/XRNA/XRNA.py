import logging
import os
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Set

from pydantic import BaseModel, PositiveInt, field_validator

from .AnnotInfo import AnnotInfo
from .DataPreprocessing import HisatTool, PreprocessingPipeline
from .PoolExecutor import PoolExecutor
from .StringtiePipeline import StringtieTool, StringtiePipeline
from ..utils import exit_with_error

logger = logging.getLogger()
logging.getLogger("matplotlib").setLevel(logging.WARNING)


class XRNAProcessor(BaseModel):
    bed_input_dir: Path
    fq_input_dir: Path
    output_dir: Path
    base_dir: Path = Path('.').resolve()
    cpus: PositiveInt = 1

    rna_ids: Set[str]
    annotation: AnnotInfo
    ouputs_prefix: str = 'xrna'
    keep_extras: Set[str] = set()

    hisat: HisatTool
    stringtie: StringtieTool

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        # spam errors
        self._validate_inputs()
        # working dir
        os.chdir(self.base_dir)
        self._work_dir = TemporaryDirectory(dir=self.base_dir)
        # other members
        work_pth = Path(self._work_dir.name)
        executor = PoolExecutor(self.cpus)
        self.annotation.prepare_annotation(work_pth)
        self._preprocessing = PreprocessingPipeline(
            work_pth, executor, self.hisat
        )
        self._pipeline = StringtiePipeline(
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
        # run pipeline
        prepared_bams = self._preprocessing.run(
            self.rna_ids, self.bed_input_dir, self.fq_input_dir, self.annotation
        )
        self._pipeline.run(
            prepared_bams, self.annotation, self.ouputs_prefix
        )
        # save outputs
        self._preprocessing.save_outputs(self.output_dir, self.keep_extras)
        self._pipeline.save_outputs(self.ouputs_prefix, self.output_dir, self.keep_extras)
        logger.info('Done.')
