import argparse
import logging

from pydantic import ValidationError

#from .DataProcessors import AllStagesProcessor, SingleStageProcessor, SUBDIR_LIST
from .Processing import AllStagesProcessor, SingleStageProcessor, SUBDIR_LIST
from .utils import configure_logger, exit_with_error, load_config

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

    config = load_config(args.config)

    try:
        if not args.stage:
            model = AllStagesProcessor(**config)
        else:
            model = SingleStageProcessor(args.stage, **config)
    except ValidationError as e:
        logger.critical('An error occured during input validation:')
        print(e)
        exit_with_error('Aborting: incorrect input!')

    model.run()
    

if __name__ == '__main__':
    main()
