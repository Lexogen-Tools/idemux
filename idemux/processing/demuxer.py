import logging
from collections import Counter

from tqdm import tqdm

from idemux.ioutils.file_handler import FileHandler
from idemux.ioutils.parser import load_correction_maps, get_pe_fastq, \
    peek_into_fastq_files
from idemux.ioutils.writer import write_summary

log = logging.getLogger(__name__)


def get_i7_i5_barcodes(mate_pair, correction_maps, has_i7):
    # TODO: add documentation
    # get fastq header of the first mate pair
    fastq_header = mate_pair[0][0]
    barcodes = fastq_header.rpartition(":")[-1].strip().split('+')
    i7_bc = None
    i5_bc = None
    if not barcodes:
        return i7_bc, i5_bc
    # when there are 2 barcodes in the fastq header the orientation is i7,i5
    if len(barcodes) == 2:
        i7_bc = barcodes[0]
        i5_bc = barcodes[1]
    # when there is only 1 barcode present we need to know if its i5 or i7
    elif len(barcodes) == 1:
        if has_i7:
            i7_bc = barcodes[0]
        else:
            i5_bc = barcodes[0]
    i7_bc_corrected = correction_maps.get("i7").get(len(i7_bc)).get(i7_bc)
    i5_bc_corrected = correction_maps.get("i5").get(len(i5_bc)).get(i5_bc)

    return i7_bc_corrected, i5_bc_corrected


def process_mate_pair(mate_pair, i7_wanted, i5_wanted, i1_wanted,
                      has_i7,
                      map_i7, map_i5, map_i1):
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

        i1_bc = mate_pair[1][1][10:22]
        _i1_corrected = map_i1.get(i1_bc)

        if _i1_corrected in i1_wanted:
            i1_bc = _i1_corrected
            (m1_hdr, m1_seq, m1_opt, m1_qcs), (m2_hdr, m2_seq, m2_opt, m2_qcs) = mate_pair

            mate_1 = (
                f"{m1_hdr[:-1]}+{_i1_corrected}\n"
                f"{m1_seq[:10]}{m1_seq[22:]}"
                f"{m1_opt}"
                f"{m1_qcs[:10]}{m1_qcs[22:]}"
            )

            mate_2 = (
                f"{m2_hdr[:-1]}+{_i1_corrected}\n"
                f"{m2_seq[:10]}{m2_seq[22:]}"
                f"{m2_opt}"
                f"{m2_qcs[:10]}{m2_qcs[22:]}"
            )

            mate_pair = (mate_1, mate_2)
    return (i7_bc, i5_bc, i1_bc), mate_pair


def demux_paired_end(args, used_lengths, barcode_sample_map, i7_wanted, i5_wanted,
                     i1_wanted):
    # TODO: add documentation
    # TODO: add logging
    # load the maps that will be used for error correction. As the tool does not allow
    # different length we only need to load the used length
    correction_map = load_correction_maps(used_lengths)
    used_lengths
    # if None is in *_wanted no barcode has been specified
    has_i7 = None not in i7_wanted
    has_i5 = None not in i5_wanted
    # before doing any processing check if the fastq file is okay.
    peek_into_fastq_files(args.r1, args.r2, has_i7, has_i5)
    # TODO: ask if i1 can be reverse complement
    read_counter = Counter()
    log.info("Staring demultiplexing")
    # first we need to open the output files the reads should get sorted into
    with FileHandler(barcode_sample_map, args.output_dir) as file_handler:
        # then we iterate over all the paired end reads
        with get_pe_fastq(args.r1, args.r2) as pe_reads:
            map_i7 = correction_map["i7"][used_lengths["i7"][0]]
            map_i5 = correction_map["i5"][used_lengths["i5"][0]]
            map_i1 = correction_map["i1"][used_lengths["i1"][0]]
            for mate_pair in tqdm(pe_reads):
                # here we do the error correction and get obtain the i1 barcode if present
                barcodes, processed_mates = process_mate_pair(mate_pair,
                                                              i7_wanted,
                                                              i5_wanted,
                                                              i1_wanted,
                                                              has_i7,
                                                              map_i7,
                                                              map_i5,
                                                              map_i1)
                # When a barcode combination is unknown/not specified in the sample sheet
                # it will get  sorted into file containing only reads with unknown
                # barcodes.
                fq_out1, fq_out2 = file_handler.get(barcodes,
                                                    file_handler["undetermined"])
                fq_out1.write(processed_mates[0])
                fq_out2.write(processed_mates[1])
    write_summary(read_counter, args.output_dir)
