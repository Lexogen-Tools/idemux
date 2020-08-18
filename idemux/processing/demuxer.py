import logging
from collections import Counter

from tqdm import tqdm

from idemux.ioutils.file_handler import FileHandler
from idemux.ioutils.parser import get_pe_fastq, peek_into_fastq_files
from idemux.ioutils.writer import write_summary

log = logging.getLogger(__name__)


def process_mate_pair(mate_pair, i7, i5, i1, i1_start, i1_end):
    """process mate pair

    Returns
    -------
    barcodes: ()
    """

    fastq_header = mate_pair[0][0]
    _, _, _barcodes = fastq_header.rpartition(":")
    barcodes = _barcodes[:-1].split("+")
    i7_bc = i5_bc = None

    # when there are 2 barcodes in the fastq header the orientation is i7,i5
    if i7.not_empty and i5.not_empty:
        i7_bc, i5_bc = barcodes
    elif i7.not_empty or i5.not_empty:
        i7_bc, i5_bc = (barcodes[0], None) if i7.not_empty else (None, barcodes[0])

    i7_bc = i7.correction_map.get(i7_bc)
    i5_bc = i5.correction_map.get(i5_bc)

    i1_bc = None

    (m1_hdr, m1_seq, m1_opt, m1_qcs), (m2_hdr, m2_seq, m2_opt, m2_qcs) = mate_pair

    if i1.not_empty:
        i1_bc = m2_seq[i1_start:i1_end]
        _i1_corrected = i1.correction_map.get(i1_bc)

        if _i1_corrected in i1.used_codes:
            m1_hdr = f"{m1_hdr[:-1]}+{i1_bc}\n"

            m2_hdr = f"{m2_hdr[:-1]}+{i1_bc}\n"
            m2_seq = f"{m2_seq[:i1_start]}{m2_seq[i1_end:]}"
            m2_qcs = f"{m2_qcs[:i1_start]}{m2_qcs[i1_end:]}"

            i1_bc = _i1_corrected

    mate_1 = (
        f"{m1_hdr}"
        f"{m1_seq}"
        f"{m1_opt}"
        f"{m1_qcs}"
    )

    mate_2 = (
        f"{m2_hdr}"
        f"{m2_seq}"
        f"{m2_opt}"
        f"{m2_qcs}"
    )

    mate_pair = (mate_1, mate_2)
    return (i7_bc, i5_bc, i1_bc), mate_pair


def demux_paired_end(barcode_sample_map, barcodes, read1, read2, i1_start, output_dir,
                     **kwargs):
    # TODO: add documentation
    # TODO: add logging
    # load the maps that will be used for error correction. As the tool does not allow
    # different length we only need to load the used length
    i7, i5, i1 = barcodes
    i1_end = i1_start + i1.length

    # before doing any processing check if the fastq file is okay.
    peek_into_fastq_files(read1, read2,
                          i7.not_empty, i5.not_empty, i1.not_empty,
                          i7.length, i5.length,
                          i1_start, i1_end)
    read_counter = Counter()
    log.info("Starting demultiplexing")
    # first we need to open the output files the reads should get sorted into
    with FileHandler(barcode_sample_map, output_dir) as file_handler:
        # then we iterate over all the paired end reads
        with get_pe_fastq(read1, read2) as pe_reads:
            for mate_pair in tqdm(pe_reads):
                # here we do the error correction and get obtain the i1 barcode if present
                barcodes, processed_mates = process_mate_pair(mate_pair,
                                                              i7, i5, i1,
                                                              i1_start,
                                                              i1_end)
                # When a barcode combination is unknown/not specified in the sample sheet
                # it will get  sorted into file containing only reads with unknown
                # barcodes.
                fq_out1, fq_out2 = file_handler.get(barcodes,
                                                    file_handler["undetermined"])
                fq_out1.write(processed_mates[0].encode())
                fq_out2.write(processed_mates[1].encode())
                read_counter[barcode_sample_map.get(barcodes, "undetermined")] += 1
    write_summary(read_counter, output_dir)
