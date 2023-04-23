import argparse
import logging

from .DataProcessors import StrandCalc
from .utils import check_file_exists, configure_logger

logger = logging.getLogger()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
            formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('-c', '--config', required=True,
                        help='Configuration file.')
    parser.add_argument('-v', '--verbose', default=False, action='store_true',
                        help='Log debug information about each step.')
    return parser.parse_args()
    

def main() -> None:
    args = parse_args()
    configure_logger(logger, args.verbose)
    logger.debug(f'Started with arguments: {vars(args)}')
    check_file_exists(args.config)

    StrandCalc(args.config).run()
    

if __name__ == '__main__':
    main()
