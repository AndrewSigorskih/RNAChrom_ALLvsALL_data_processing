import argparse
import json

from .DataProcessors import BaseProcessor
from .utils import check_file_exists, configure_logger

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
        configure_logger(args.verbose)
        if args.verbose:
            print(f'Started with arguments: {vars(args)}')

    def run(self) -> None:

        if self.mode == 'XRNA':
            raise NotImplementedError("Warning: XRNA collection is not implemented yet!")
            exit() # for when it will be implemented
        # we are collecting contacts
        # actually run base class or any of its pre-defined subclasses here
        # depending on config
        BaseProcessor(self.config).run()

# TODO ?? hide trimmomatic text blanket in tmp log file(s)
# TODO add "grid" and "imargi" rsite-dealing strategies
# TODO add hisat internal cpus argument for even-finer tuning
# TODO add cpus argument for contacts building
# TODO [DONE] add condition on whether exit from run_function if non-zero exitcode returned
# TODO [DONE] make actual exitcode returnal from filter bams and bamtobed
# TODO make statistics calculation as separate step
# TODO calculate actual strand of contacts rna parts and store in useful format


def main() -> None:
    Program().run()
    

if __name__ == '__main__':
    main()
