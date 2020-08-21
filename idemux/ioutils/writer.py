
"""Summary statistics writer"""

import csv
import logging
import os

log = logging.getLogger(__name__)


def write_summary(counter, output_dir):
    """Writes a counter dictionary as tsv file with the name error_correction_stats.tsv.
    Args:
        counter(dict): A dictionary barcodes <sample name, #corrected reads>.
        output_dir (string): The path the file should be written to.
    """
    output_file = os.path.join(output_dir, 'demultipexing_stats.tsv')

    with open(output_file, 'w') as csvfile:
        fieldnames = ['sample_name', 'written_reads']
        writer = csv.DictWriter(csvfile, delimiter='\t', fieldnames=fieldnames)
        writer.writeheader()
        for file_name, counts in counter.items():
            writer.writerow({'sample_name': file_name, 'written_reads': counts})
    log.info("Run complete! Summary statistics saved to %s", output_file)

