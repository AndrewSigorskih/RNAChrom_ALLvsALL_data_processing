import logging
import os
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any, Dict, Optional

logger = logging.getLogger()

class XRNAProcessor:
    def __init__(self,
                 cfg: Dict[str, Any]):
        # get basic parameters
        self.cpus: int = cfg.get('cpus', 1)
        self.base_dir: Path = Path(cfg.get('base_dir', os.getcwd())).resolve()
        self.fq_input_dir: Optional[Path] =\
            Path(cfg.get('input_dir')) if 'input_dir' in cfg else None
        self.bam_input_dir: Optional[Path] =\
            Path(cfg.get('bam_input_dir')) if 'bam_input_dir' in cfg else None
        #???????
         # spam errors
        self._validate_inputs()
        # create working directory
        os.chdir(self.base_dir)
        self.work_dir = TemporaryDirectory(dir=self.base_dir)

    def _validate_inputs(self):
        """Basic input sanity check"""
        pass
    
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