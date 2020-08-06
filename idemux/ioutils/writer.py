import csv
import gzip
import logging
import os
from contextlib import contextmanager

log = logging.getLogger(__name__)


def write_summary(counter, output_dir):
    """Writes a counter dictionary as tsv file with the name error_correction_stats.tsv.
    Args:
        counter(dict): A dictionary barcodes <sample name, #corrected reads>.
        output_dir (string): The path the file should be written to.
    """
    output_file = os.path.join(output_dir, 'demultipexing_stats.tsv')

    with open(output_file, 'w') as csvfile:
        fieldnames = ['sample', 'written_reads']
        writer = csv.DictWriter(csvfile, delimiter='\t', fieldnames=fieldnames)
        writer.writeheader()
        for file_name, counts in counter.items():
            writer.writerow({'filename': file_name, 'rescued_reads': counts})
    log.info("Run complete! Summary statistics saved to %s", output_file)


@contextmanager
def output_file_handler_pe(barcode_file_map, output_folder):
    """Barcodes exist in a length of either 8, 10 or 12 nt. This function selects the
    right length subset depending on supplied barcodes in the barcode.tsv file. This is
    Args:
        barcode_file_map (dict): A dictionary barcodes barcodes to their respective
          demultiplexed file, <(i7,i5,i1),file_name>.
        output_folder (string): The output folder path.
        max_memory (int): The memory available to the buffered writer. The memory for
          each file corresponds to max_memory // number of files.
    Yield:
        file_handler (dict): A dictionary barcodes barcodes a buffered output file stream.
          <(i7,i5),io_stream>
    Raise:
        Exception: An exception is raised when closing of the handled files is not
          possible
    """
    # TODO: think about adding overwrite protection
    file_handler = {}
    # make output folder if it does not exist
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    # add this one for handling barcodes that dont match any known barcodes
    barcode_file_map["undetermined"] = "undetermined"
    # open buffered writers for output files so we can write reads directly as fastq.gz
    try:
        for barcode, sample_name in barcode_file_map.items():
            # construct the output paths
            mate1_path = out_file_path = os.path.join(output_folder, sample_name,
                                                      "_R1.fastq.gz")
            mate2_path = out_file_path = os.path.join(output_folder, sample_name,
                                                      "_R2.fastq.gz")
            mate1_writer = gzip.open(mate1_path, mode='wt', encoding="utf-8",
                                     compresslevel=4)
            log.debug("File handler opened for mate1: %s", mate1_path)
            mate2_writer = gzip.open(mate2_path, mode='wt', encoding="utf-8",
                                     compresslevel=4)
            log.debug("File handler opened for mate2: %s", mate2_path)
            file_handler[barcode] = mate1_writer, mate2_writer

        file_handler["undetermined"] = mate1_writer, mate2_writer
        yield file_handler
    finally:
        for barcode, io_streams in file_handler.items():
            for stream in io_streams:
                try:
                    stream.close()
                except IOError:
                    log.exception('File handler was unable to close file: "%s"' %
                                  stream.name)
