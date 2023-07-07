import argparse
import json
import logging

from .DataProcessors import AllStagesProcessor, SingleStageProcessor, SUBDIR_LIST
from .utils import check_file_exists, configure_logger

logger = logging.getLogger()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
            formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('-c', '--config', required=True,
                        help='Configuration file.')
    parser.add_argument('-s', '--stage', required=False,
                        choices=SUBDIR_LIST,
                        help='Run only specified stage of the pipeline.')
    parser.add_argument('-v', '--verbose', action='count', default=0,
                        help='''Verbosity level. By default little to none information is printed.
                        Use -v once to increase information logs about each step, and -vv to 
                        print every command that is being run.''')
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    configure_logger(logger, args.verbose)
    logger.debug(f'Started with arguments: {vars(args)}')

    check_file_exists(args.config)
    with open(args.config, 'r') as f:
        config: dict = json.load(f)

    if not args.stage:
        AllStagesProcessor(config).run()
    else:
        SingleStageProcessor(config, args.stage).run()
    

if __name__ == '__main__':
    main()
