
"""File handler module - handles writing of fastq.gz files"""

import gzip
import logging
import os
import pathlib
from collections import defaultdict
from contextlib import ExitStack
import io

log = logging.getLogger(__name__)


class FileHandler(ExitStack):

    def __init__(self, barcode_file_map, output_folder, memory=2 ** 32):
        self.barcode_file_map = barcode_file_map
        self.output_folder = pathlib.Path(output_folder)
        if not os.path.exists(output_folder): # pragma: no cover
            os.makedirs(output_folder)
        self.buffer_per_file = int(memory // (len(barcode_file_map) * 2))
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
            gz_out = io.BufferedWriter(gzip.open(mate_path, mode='wb', compresslevel=4),
                                       buffer_size=self.buffer_per_file)
            out.append(gz_out)
            log.debug("File handler opened for mate: %s", mate_path)
        return out
