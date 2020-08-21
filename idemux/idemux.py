
"""Main module"""

import argparse
import logging
import sys

from idemux import __version__
from idemux.ioutils.parser import parse_sample_sheet
from idemux.processing.demuxer import demux_paired_end

# set logger format
logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    level=logging.INFO)
log = logging.getLogger()


def get_cli_parser():
    """Function that returns an argparse object for command line parsing.
    Returns:
        parser: Returns parser.parse_args() object which contains all the set command line
          options.
    Raises:
        argparse.ArgumentTypeError
    """
    # TODO: update descriptions
    parser = argparse.ArgumentParser(prog='idemux',
                                     description='A tool to demultiplex fastq '
                                                 'files based on Lexogen i7,i5,'
                                                 'i1  barcodes.')
    required_args = parser.add_argument_group('required arguments')
    required_args.add_argument('--r1',
                               type=str,
                               dest='read1',
                               required=True,
                               help='gzipped read 1 fastq file'
                               )
    required_args.add_argument('--r2',
                               type=str,
                               dest='read2',
                               required=True,
                               help='path to gzipped read 2 fastq file'
                               )

    required_args.add_argument('--sample-sheet',
                               type=str,
                               dest='sample_sheet',
                               required=True,
                               help='csv file describing sample names, and barcode '
                                    'combinations '
                               )
    required_args.add_argument('--out',
                               type=str,
                               dest='output_dir',
                               required=True,
                               help='where to write the output files'
                               )
    parser.add_argument('--i5-rc',
                        dest='i5_rc',
                        action='store_true',
                        default=False,
                        help='when the i5 barcode has been sequenced as reverse '
                             'complement. make sure to enter non-reverse complement '
                             'sequences in the barcode file ')
    parser.add_argument('--i1-start',
                        type=int,
                        default=11,
                        dest='i1_start',
                        help='start position of the i1 index (1-based) on read 2 '
                             '(default: 11)')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' +
                                                                     __version__)
    return parser


def main():
    # get the command line arguments
    cli_parser = get_cli_parser()
    # convert to dict so we can use it as **kwargs
    args = vars(cli_parser.parse_args())
    # convert 1 to 0 based indexing
    args['i1_start'] = args['i1_start'] - 1

    # the sample sheet defines sample barcode relations and how the reads should be
    # demultiplexed
    barcode_sample_map, barcodes = parse_sample_sheet(**args)
    # demultiplex and error correct for PE files
    demux_paired_end(barcode_sample_map, barcodes, **args)


if __name__ == "__main__":
    sys.exit(main())
