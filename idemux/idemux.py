import argparse
import logging
from idemux.ioutils.parser import parse_sample_sheet
from idemux.processing.demuxer import demux_paired_end

logging.basicConfig(level=logging.INFO)
log = logging.getLogger()
tool_version = "0.1"


def get_cli_parser():
    """Function that returns an argparse object for command line parsing.
    Returns:
        parser: Returns parser.parse_args() object which contains all the set command line
          options.
    Raises:
        argparse.ArgumentTypeError
    """
    # TODO: update descriptions
    parser = argparse.ArgumentParser(prog='error_correction.py',
                                     description='A tool to demultiplex fastq '
                                                 'files based on Lexogen i7,i5,'
                                                 'i1  barcodes.')
    required_args = parser.add_argument_group('required arguments')
    required_args.add_argument('--r1',
                               type=str,
                               dest='r1',
                               required=True,
                               help='Path to the read 1 fastq file that should be '
                                    'demultiplexed'
                               )
    required_args.add_argument('--r2',
                               type=str,
                               dest='r2',
                               required=True,
                               help='Path to the read 2 fastq file that should be '
                                    'demultiplexed'
                               )

    required_args.add_argument('--sample-sheet',
                               type=str,
                               dest='sample_sheet',
                               required=True,
                               help='Csv file containing sample names, i7, i5 and i1 '
                                    'barcodes'
                               )
    required_args.add_argument('--out',
                                type=str,
                                dest='output_dir',
                                required=True,
                                help='Where to write the output files.'
                                )
    # TODO: check how to do this properly
    parser.add_argument('--version', action='version', version='%(prog)s ' + tool_version)
    return parser


def main():

    # TODO: demultiplexing can be done either on i7 or i5 check what happens in case i5
    #  is specified, and only 1 barcode is in the header
    # get the command line arguments
    cli_parser = get_cli_parser()
    args = cli_parser.parse_args()
    sample_sheet = args.sample_sheet

    # the sample sheet defines sample barcode relations and how the
    # reads should be demultiplexed
    barcode_sample_map, wanted_barcodes, used_lengths = parse_sample_sheet(sample_sheet)

    i7_wanted, i5_wanted, i1_wanted = wanted_barcodes

    # demultiplex and error correct for PE files
    demux_paired_end(args, used_lengths, barcode_sample_map, i7_wanted, i5_wanted,
                     i1_wanted)


if __name__ == "__main__":
    main()
