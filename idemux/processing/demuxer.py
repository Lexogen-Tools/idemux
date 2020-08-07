import logging
from collections import Counter

from idemux.ioutils.parser import load_correction_maps, get_pe_fastq, \
    peek_into_fastq_files
from idemux.ioutils.writer import output_file_handler_pe, write_summary

log = logging.getLogger(__name__)


def get_i7_i5_barcodes(one_read, correction_maps, has_i7):
    # TODO: add documentation
    fastq_header = one_read[0]
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
    #i5_bc_corrected = correction_maps["i5"].get(i5_bc)

    return i7_bc_corrected, i5_bc_corrected


def process_inline_barcodes(mate_pair, correction_maps, i1_wanted, mate_with_i1=1):
    # TODO: add documentation
    # mapping fq lines to indices
    fq_header_idx = 0
    fq_seq_idx = 1
    fq_seq_scores = 3
    # if there is a inline barcode that we cut out, we need to do this in these lines
    fq_lines_to_shorten = [fq_seq_idx, fq_seq_scores]

    # inline barcodes are 12 nt long
    bc_start_pos = 10
    bc_end_pos = 22
    bc_range = [bc_start_pos, bc_end_pos]

    mate_seq_with_i1 = mate_pair[mate_with_i1][fq_seq_idx]
    i1_bc = mate_seq_with_i1[bc_start_pos:bc_end_pos]
    i1_corrected = correction_maps.get("i1").get(len(i1_bc)).get(i1_bc)
    if i1_corrected in i1_wanted:
        return i1_corrected, update_fq_read_pairs(mate_pair,
                                                  mate_with_i1,
                                                  i1_corrected,
                                                  bc_range,
                                                  fq_lines_to_shorten)
    else:
        # log.warning("Invalid inline barcode present.\n"
        #             "Barcode: %s\n"
        #             "Read header: %s", i1_bc, mate_pair[mate_with_i1][fq_header_idx])
        # log.warning("Read will be sorted into undetermined bin.")
        #
        # # TODO: Fix return values
        return None, mate_pair


def update_fq_read_pairs(mate_pair, mate_with_i1, i1_bc, bc_range, fq_lines_to_mod):
    # if there is an i1 present we need to add it first to the fq header
    updated_mates = [update_fq_header(mate, i1_bc) for mate in mate_pair]
    # then we cut out the barcode from the read. However, when we do that we need to
    # make sure other lines as the quality scores are shortened as well.
    for idx in fq_lines_to_mod:
        updated_mates[mate_with_i1][idx] = remove_substring(mate_pair[mate_with_i1][
                                                                idx], *bc_range)
    return updated_mates


def update_fq_header(fq_read, i1_bc):
    # make a list from the tuples, as they are immutable
    fq_read = list(fq_read)
    new_fq_header = fq_read[0][:-1] + "+" + i1_bc + "\n"
    fq_read[0] = new_fq_header
    # TODO: add documentation
    return fq_read


def remove_substring(s, start, end):
    # # TODO: add documentation
    # if len(s) < end or len(s) < start:
    #     # TODO: add raise
    #     log.exception("String %s is shorter than arguments given for slicing.\n"
    #                   "Start: % d\n"
    #                   "End: %d"
    #                   % s, start, end)
    return s[0:start] + s[end:]


def process_mate_pair(mate_pair, correction_maps, i7_wanted, i5_wanted, i1_wanted, has_i7):
    # TODO: add documentation
    mate_1, mate_2 = mate_pair

    i7_bc, i5_bc = get_i7_i5_barcodes(mate_2, correction_maps, has_i7)
    updated_mate_pair = mate_pair
    # TODO: check if value should ne None or empty string
    i1_bc = None
    if i7_bc in i7_wanted and i5_bc in i5_wanted:
        i1_bc, updated_mate_pair = process_inline_barcodes(mate_pair,
                                                           correction_maps,
                                                           i1_wanted)
    barcodes = (i7_bc, i5_bc, i1_bc)
    return barcodes, updated_mate_pair


def demux_paired_end(args, used_lengths, barcode_sample_map, i7_wanted, i5_wanted,
                     i1_wanted):
    # TODO: add documentation
    # TODO: add logging
    # load the maps that will be used for error correction. As the tool does not allow
    # different length we only need to load the used length
    correction_map = load_correction_maps(used_lengths)

    # if None is in *_wanted no barcode has been specified
    has_i7 = None not in i7_wanted
    has_i5 = None not in i5_wanted
    # before doing any processing check if the fastq file is okay.
    peek_into_fastq_files(args.r1, args.r2, has_i7, has_i5)

    # TODO: ask if i1 can be reverse complement
    read_counter = Counter()
    # first we need to open the output files the reads should get sorted into
    with output_file_handler_pe(barcode_sample_map, args.output_dir) as file_handler:
        # then we iterate over all the paired end reads
        with get_pe_fastq(args.r1, args.r2) as pe_reads:
            for mate_pair in pe_reads:
                # here we do the error correction and get obtain the i1 barcode if present
                barcodes, processed_mates = process_mate_pair(mate_pair,
                                                              correction_map,
                                                              i7_wanted,
                                                              i5_wanted,
                                                              i1_wanted,
                                                              has_i7)
                # When a barcode combination is unknown/not specified in the sample sheet
                # it will get  sorted into file containing only reads with unknown
                # barcodes.
                fq_out1, fq_out2 = file_handler.get(barcodes,
                                                    file_handler["undetermined"])

                fq_out1.writelines(processed_mates[0])
                fq_out2.writelines(processed_mates[1])

                # We want some summary statistics in the end to get an idea how many
                # reads have been sorted into what file/bin. So do some counting for that
                sample_name = barcode_sample_map.get(barcodes, "undetermined")
                read_counter[sample_name] += 1
    write_summary(read_counter, args.output_dir)
