import logging
import os
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any, Dict, List, Optional

from pydantic import BaseModel, field_validator, PositiveInt

from .DataPreprocessing import HisatTool, PreprocessingPipeline
from .PoolExecutor import PoolExecutor
from .StringtiePipeline import StringtieTool
from ..utils import exit_with_error

logger = logging.getLogger()


class XRNAProcessor(BaseModel):
    cpus: Optional[PositiveInt] = 1
    base_dir: Path
    input_dir: Path
    output_dir: Path

    annot_gtf: Path # take all annot files in one sub config?

    rna_ids: List[str]
    hisat: HisatTool
    stringtie: StringtieTool

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        # spam errors
        self._validate_inputs()
        # working dir
        os.chdir(self.base_dir)
        self.work_dir = TemporaryDirectory(dir=self.base_dir)
        # other members
        work_pth = Path(self.work_dir.name)
        executor = PoolExecutor(self.cpus)
        self.preprocessing = PreprocessingPipeline(
            work_pth, executor, self.hisat, self.rna_ids
        )

    @field_validator('base_dir', 'input_dir', 'output_dir')
    @classmethod
    def resolve_path(cls, pth: str) -> Path:
        return Path(pth).resolve()
    
    def _validate_inputs(self):
        """Basic input sanity check"""
        if not self.base_dir.exists():
            self.base_dir.mkdir(parents=True)
        if not self.input_dir.exists():
            exit_with_error('Input directory does not exist!')
        if not self.output_dir.exists():
            self.output_dir.mkdir(parents=True)   

    def run(self):
        pass


        # reference:
            # save in bed format
        # contacts:
            # filter against reference (bedtools?)
            # save lists of good ids
        # fq inputs:
            # revcompl
            # filter by id lists
            # align using hisat
            # filter bams
        # bam inputs
            # filter by lists of good ids
        # all bams:
            # merge technical replics
            # sort bams
            # split strand and xs tag?????
            # assemble stringtie
            # stringtie merge


    # https://pyfastx.readthedocs.io/en/latest/usage.html#reverse-and-complement-sequence