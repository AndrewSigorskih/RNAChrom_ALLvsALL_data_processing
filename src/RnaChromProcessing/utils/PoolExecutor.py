import concurrent.futures
import logging
from typing import Any, Callable, List, Optional

from .errors import StageFailedError

logger = logging.getLogger()


class PoolExecutor:
    def __init__(self,
                 cpus: int):
        self.cpus = cpus

    def run_function(self,
                     func: Callable[[Any, Any], int],
                     inputs: List[Any],
                     outputs: List[Any],
                     require_zero_code: bool = True,
                     override_cpus: Optional[int] = None) -> None:
        if len(outputs) != len(inputs):
            msg = f'Lengths of inputs and outputs lists for {func.__qualname__} do not match!'
            raise StageFailedError(msg)
        cpus = override_cpus or self.cpus
        logger.debug(f'Running function {func.__qualname__} with {cpus} threads.')
        with concurrent.futures.ThreadPoolExecutor(max_workers=cpus) as executor:
            futures = [
                executor.submit(func, inp, out)
                for inp, out in zip(inputs, outputs)
            ]
            results = [
                future.result() for future in concurrent.futures.as_completed(futures)
            ]
        if require_zero_code and any(x != 0 for x in results):
            msg = f'One or several calls of {func.__qualname__} returned non-zero process exit code!'
            raise StageFailedError(msg)
