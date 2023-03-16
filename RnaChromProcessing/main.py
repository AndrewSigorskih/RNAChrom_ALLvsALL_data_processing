import argparse
import json
import logging

from .DataProcessors import BaseProcessor
from .utils import check_file_exists, configure_logger

logger = logging.getLogger(name=__name__)

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
            formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('-c', '--config', required=True,
                        help='Configuration file.')
    parser.add_argument('-mode', '--run_mode', required=False,
                        choices=('contacts', 'XRNA'), default='contacts',
                        help='''Action to be performed during run:
* contacts: build contacts data from provided raw FASTQ files
* XRNA: attempt to infer X-RNAs from provided contacts and BAM files.''')
    parser.add_argument('-v', '--verbose', default=False, action='store_true',
                        help='Add logging verbosity.')
    return parser.parse_args()

class Program:
    def __init__(self) -> None:
        args: argparse.Namespace = parse_args()
        self.mode: str = args.run_mode
        check_file_exists(args.config)
        with open(args.config, 'r') as f:
            self.config: dict = json.load(f)
        if args.verbose:
            level = logging.DEBUG
        else:
            level = logging.INFO
        configure_logger(level)
        logger.debug(f'Started with arguments: {vars(args)}')

    def run(self) -> None:
        if self.mode == 'XRNA':
            raise NotImplementedError("Warning: XRNA collection is not implemented yet!")
            exit() # for when it will be implemented
        # we are collecting contacts
        # actually run base class or any of its pre-defined subclasses here
        # depending on config
        BaseProcessor(self.config).run()


def main() -> None:
    Program().run()
    

if __name__ == '__main__':
    main()
