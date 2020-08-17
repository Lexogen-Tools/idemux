import logging
from collections import Counter

from tqdm import tqdm

from idemux.ioutils.file_handler import FileHandler
from idemux.ioutils.parser import get_pe_fastq, peek_into_fastq_files
from idemux.ioutils.writer import write_summary

log = logging.getLogger(__name__)

def process_mate_pair(mate_pair,
                      i7_wanted, i5_wanted, i1_wanted,
                      has_i7,
                      map_i7, map_i5, map_i1,
                      i1_start, i1_end):
    """process mate pair

    Returns
    -------
    barcodes: ()
    """

    fastq_header = mate_pair[0][0]
    _, _, _barcodes = fastq_header.rpartition(":")
    barcodes = _barcodes[:-1].split("+")
    # if not barcodes:
    i7_bc = i5_bc = None

    # when there are 2 barcodes in the fastq header the orientation is i7,i5
    if len(barcodes) == 2:
        i7_bc, i5_bc = barcodes

    elif len(barcodes) == 1:
        i7_bc, i5_bc = (barcodes[0], None) if has_i7 else (None, barcodes[0])

    i7_bc = map_i7.get(i7_bc)
    i5_bc = map_i5.get(i5_bc)

    i1_bc = None

    if i7_bc in i7_wanted and i5_bc in i5_wanted:

        i1_bc = mate_pair[1][1][i1_start:i1_end]
        _i1_corrected = map_i1.get(i1_bc)

        if _i1_corrected in i1_wanted:
            i1_bc = _i1_corrected
            (m1_hdr, m1_seq, m1_opt, m1_qcs), (m2_hdr, m2_seq, m2_opt, m2_qcs) = mate_pair

            mate_1 = (
                f"{m1_hdr[:-1]}+{_i1_corrected}\n"
                f"{m1_seq}"
                f"{m1_opt}"
                f"{m1_qcs[:i1_start]}{m1_qcs[i1_end:]}"
            )

            mate_2 = (
                f"{m2_hdr[:-1]}+{_i1_corrected}\n"
                f"{m2_seq[:i1_start]}{m2_seq[i1_end:]}"
                f"{m2_opt}"
                f"{m2_qcs[:i1_start]}{m2_qcs[i1_end:]}"
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
    i7_wanted = i7.used_codes
    i5_wanted = i5.used_codes
    i1_wanted = i1.used_codes
    map_i7 = i7.correction_map
    map_i5 = i5.correction_map
    map_i1 = i1.correction_map
    i1_end = i1_start + i1.length

    # if None is in *_wanted no barcode has been specified
    has_i7 = not i7.empty
    has_i5 = not i5.empty
    has_i1 = not i1.empty

    # before doing any processing check if the fastq file is okay.
    peek_into_fastq_files(read1, read2, has_i7, has_i5, has_i1, i7.length,
                          i5.length, i1_start, i1_end)
    read_counter = Counter()
    log.info("Staring demultiplexing")
    # first we need to open the output files the reads should get sorted into
    with FileHandler(barcode_sample_map, output_dir) as file_handler:
        # then we iterate over all the paired end reads
        with get_pe_fastq(read1, read2) as pe_reads:
            for mate_pair in tqdm(pe_reads):
                # here we do the error correction and get obtain the i1 barcode if present
                barcodes, processed_mates = process_mate_pair(mate_pair,
                                                              i7_wanted,
                                                              i5_wanted,
                                                              i1_wanted,
                                                              has_i7,
                                                              map_i7,
                                                              map_i5,
                                                              map_i1,
                                                              i1_start,
                                                              i1_end)
                # When a barcode combination is unknown/not specified in the sample sheet
                # it will get  sorted into file containing only reads with unknown
                # barcodes.
                fq_out1, fq_out2 = file_handler.get(barcodes,
                                                    file_handler["undetermined"])
                fq_out1.write(processed_mates[0].encode())
                fq_out2.write(processed_mates[1].encode())
    write_summary(read_counter, output_dir)
