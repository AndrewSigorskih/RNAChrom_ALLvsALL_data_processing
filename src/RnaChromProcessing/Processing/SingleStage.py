from logging import getLogger
from os import chdir
from pathlib import Path
from typing import Dict, List, Protocol

from .Base import BaseProcessor
from .stages import (
    Align, BamFilter, BamToBed, Contacts, Dedup, Rsites, Trim
)
from .stages.basicstage import SamplePair
from ..utils import exit_with_error, move_exist_ok

logger = getLogger()


class IStage(Protocol):
    def run(self, samples: List[SamplePair]) -> List[SamplePair]:
        ...
    
    def set_params(self, global_cpus: int, stage_dir: Path) -> None:
        ...

STAGES_MAP: Dict[str, IStage] = {
    'rsites': Rsites, 'dedup': Dedup, 'trim': Trim,
    'align': Align, 'bam': BamFilter, 'bed': BamToBed, 'contacts': Contacts
}


class SingleStageProcessor(BaseProcessor):
    def __init__(self, stagename: str, **kwargs):
        super().__init__(**kwargs)
        stage_conf = kwargs.get(stagename, dict())
        self._stagename = stagename
        self._stage: IStage = STAGES_MAP[stagename](**stage_conf)
        self._stage.set_params(self.cpus, self._work_pth / stagename)

    def save_outputs(self):
        """copy results to out dir"""
        source_pth = self._work_pth / self._stagename
        dest_pth = self.output_dir / self._stagename
        dest_pth.mkdir(exist_ok=True)
        move_exist_ok(source_pth, dest_pth)

    def run(self) -> None:
        samples = self.gather_inputs()
        logger.info(f'Started processing {len(samples)} pairs of files.')
        # run processing
        chdir(self._work_pth)
        try:
            self._stage.run(samples)
        except Exception as _:
            import traceback
            logger.critical(f'An error occured during {self._stagename} stage execution:')
            exit_with_error(traceback.format_exc())
        # save results
        chdir(self.base_dir)
        self.save_outputs()
        logger.info('Done.')
