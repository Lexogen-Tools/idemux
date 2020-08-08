import gzip
import logging
import pathlib
from collections import defaultdict
from contextlib import ExitStack

log = logging.getLogger(__name__)

class FileHandler(ExitStack):

    def __init__(self, barcode_file_map, output_folder):
        self.barcode_file_map = barcode_file_map
        self.output_folder = pathlib.Path(output_folder)

        super().__init__()
        self.fastq_handler = defaultdict(list)
        self.barcode_file_map["undetermined"] = "undetermined"

    def __enter__(self):
        for barcode, sample_name in self.barcode_file_map.items():
            mate1_path = self.output_folder / "{}{}".format(sample_name, "_R1.fastq.gz")
            mate2_path = self.output_folder / "{}{}".format(sample_name, "_R2.fastq.gz")
            self.fastq_handler[barcode] = self._open_gz_file_handles(
                mate1_path, mate2_path
            )
        return self.fastq_handler

    def _open_gz_file_handles(self, mate1_path, mate2_path):
        out = []
        for mate_path in (mate1_path, mate2_path):
            gz_out = self.enter_context(gzip.open(mate_path, mode='wt',
                                                  encoding="utf-8",
                                                  compresslevel=4))
            out.append(gz_out)
            log.debug("File handler opened for mate: %s", mate_path)
        return out
