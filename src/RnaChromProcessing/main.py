import argparse
import json
import logging

from .DataProcessors import BaseProcessor
from .utils import check_file_exists, configure_logger

logger = logging.getLogger()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
            formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('-c', '--config', required=True,
                        help='Configuration file.')
    parser.add_argument('-v', '--verbose', action='count', default=0,
                        help='''Verbosity level. By default little to none information is printed.
                        Use -v once to increase information logs about each step, and -vv to 
                        print every command that is being run.''')
    return parser.parse_args()

# TODO [DONE] do something about stats being so slow
# TODO [DONE] add cpus argument for contacts building
# TODO [DONE] calculate actual strand of contacts rna parts and store in useful format
# TODO [DONE] remove files that are not needed anymore (stats directory)
# TODO [FAILED] verify symlink creation works correctly
# TODO test GRID data and remove biopython-based outdated rsite solution

def main() -> None:
    args = parse_args()
    configure_logger(logger, args.verbose)
    logger.debug(f'Started with arguments: {vars(args)}')

    check_file_exists(args.config)
    with open(args.config, 'r') as f:
        config: dict = json.load(f)

    BaseProcessor(config).run()
    

if __name__ == '__main__':
    main()
