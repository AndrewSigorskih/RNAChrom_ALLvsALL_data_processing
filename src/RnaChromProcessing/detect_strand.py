import argparse

from .DataProcessors import StrandCalc
from .utils import check_file_exists, configure_logger


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
            formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('-i', '--input-dir', required=True,
                        help='Directory with contacts files.')
    parser.add_argument('-gtf', '--gtf-annotation', required=True,
                        help='Gene annotation in GTF format.')
    parser.add_argument('-g', '--genes', required=True,
                        help='File with list of gene names to test, one per line.')
    parser.add_argument('-o', '--output-dir',
                        help='Directory to store outputs in. Default: current directory')
    parser.add_argument('-p', '--prefix', required=False, default='strand',
                        help='Prefix for output files. Default: "strand"')
    parser.add_argument('-v', '--verbose', default=False, action='store_true',
                        help='Add logging verbosity.')
    return parser.parse_args()
    

def main() -> None:
    args = parse_args()
    configure_logger(args.verbose)
    check_file_exists(args.gtf_annotation)
    check_file_exists(args.genes)
    StrandCalc(input_dir=args.input_dir,
               gtf_annotation=args.gtf_annotation,
               genes=args.genes,
               output_dir=args.output_dir,
               prefix=args.prefix).run()
    

if __name__ == '__main__':
    main()
