import csv
import gzip
import logging
import sys
from collections import defaultdict, Counter
from contextlib import contextmanager
from itertools import zip_longest

# TODO: add importlib-resources to installer/setup.py dependency
try:
    from importlib import resources
except ImportError:
    import importlib_resources as resources

log = logging.getLogger(__name__)


def parse_sample_sheet(cvs_file):
    # TODO: add documentation
    i7_lengths = set()
    i5_lengths = set()
    i1_lengths = set()

    i7_barcodes = defaultdict(list)
    i5_barcodes = defaultdict(list)
    i1_barcodes = defaultdict(list)

    barcode_sample_map = dict()
    valid_barcodes = get_valid_barcodes()

    sample_counter = Counter()
    # TODO: this is not the right place for the i1 length to be specified. If the i1
    #  length changes in the future this can easily be missed. Maybe move to an ini
    #  file or some singleton?
    # i1 inline barcode is always 12 nt in length
    expected_i1_bc_length = 12
    # expected file header
    sample_sheet_header = ["sample_name", "i7", "i5", "i1"]
    try:
        with open(cvs_file) as sample_sheet:
            csv.register_dialect('strip', skipinitialspace=True)
            lines_read = 0
            reader = csv.DictReader(sample_sheet, dialect='strip')
            for row in reader:
                lines_read += 1
                # check if all the needed columns are present. otherwise we want to
                # disallow further processing
                if not set(sample_sheet_header) <= row.keys():
                    s1 = ",".join(sample_sheet_header)
                    s2 = ",".join(row.keys())
                    error_msg = ("Incorrect sample sheet header. Expected header: %s\n"
                                 "Observed header: %s" % (s1, s2))
                    raise ValueError(error_msg)

                # use get so if there are empty values None is returned
                sample_name = row.get('sample_name')
                i7_bc = row.get('i7')
                i5_bc = row.get('i5')
                i1_bc = row.get('i1')

                # undefined sample names are not allowed
                if sample_name is None:
                    error_msg = ("Incorrect sample name in line %s. Empty values are "
                                 "not allowed as sample names." % lines_read)
                    raise ValueError(error_msg)
                sample_counter[sample_name] += 1

                # make sure the barcodes are actual Lexogen barcodes so we can error
                # correct them later on
                is_valid_barcode(i7_bc, valid_barcodes, "i7", sample_name)
                is_valid_barcode(i5_bc, valid_barcodes, "i5", sample_name)
                is_valid_barcode(i1_bc, valid_barcodes, "i1", sample_name)

                # add barcodes and sample_names to dict so we can do some typ checking
                # later as not all combinations should be allowed
                i7_barcodes[i7_bc].append(sample_name)
                i5_barcodes[i5_bc].append(sample_name)
                i1_barcodes[i1_bc].append(sample_name)
                barcodes = (i7_bc, i5_bc, i1_bc)

                # barcode combinations have to be unique. otherwise we don't know to which
                # sample they belong
                if barcodes in barcode_sample_map:
                    same_barcode = barcode_sample_map[barcodes]
                    error_msg = ("Duplicate barcode combination detected. Sample: %s "
                                 "Barcodes: %s\n Already observed for: %s\n Each "
                                 "barcode combination has to be unique." % (sample_name,
                                                                            barcodes,
                                                                            same_barcode))
                    raise ValueError(error_msg)

                barcode_sample_map[barcodes] = sample_name
                # barcodes of different length are a problem as 8, 10, 12 nucleotide
                # barcodes are a subset of each other. Therefore we only allow one
                # length per barcode type
                i7_lengths.add(len(i7_bc))  # i7 is always required and cannot be none

                if i5_bc is not None:
                    i5_lengths.add(len(i5_bc))
                if i1_bc is not None:
                    i1_lengths.add(len(i1_bc))

        # barcodes for each type need to be of the same length. This is as the
        # sequences of 8 nt long barcodes are contained within 10 nt long barcodes,
        # and the ones of 10 in 12. if we dont enforce the same length per barcode there
        # might be a possibility we cant tell with certainty to which sample a
        # barcoded read belongs.
        has_different_lengths("i7", i7_lengths)
        has_different_lengths("i5", i5_lengths)
        has_different_lengths("i1", i1_lengths)

        observed_i1_length = list(i1_lengths)[0]
        if observed_i1_length != expected_i1_bc_length:
            error_msg = ("i1 inline barcodes should be %s nt long. The i1 indices you "
                         "specified are only %s nt long." % (expected_i1_bc_length,
                                                             observed_i1_length)
                         )
            raise ValueError(error_msg)

        # test if the supplied barcode combinations are valid
        has_valid_barcode_combinations(i7_barcodes, i5_barcodes, i1_barcodes)

        # sample names have to be unique as they determine the outfile name. Otherwise
        # we get problems when we try to write reads belonging to different barcode
        # combinations to one file.
        duplicated_sample_names = [key for key, value in sample_counter.items() if
                                   value > 1]
        if duplicated_sample_names:
            error_msg = ("The sample sheet contains duplicate sample names. Sample "
                         "names have to be unique.\n Sample names: %s" %
                         duplicated_sample_names)
            raise ValueError(error_msg)

    except Exception as e:
        log.exception(e)
        sys.exit(1)
    # we have made it until here until raising an exception. That means the sample sheet
    # information should be okay and we can return required data from the sample sheet.
    i7_used = set(i7_barcodes.keys())
    i5_used = set(i5_barcodes.keys())
    i1_used = set(i1_barcodes.keys())

    used_barcodes = (i7_used, i5_used, i1_used)
    return barcode_sample_map, used_barcodes, {"i7": i7_lengths,
                                               "i5": i5_lengths,
                                               "i1": i1_lengths}


def has_valid_barcode_combinations(i7_barcodes, i5_barcodes, i1_barcodes):
    # TODO: add documentation
    # allowed are:
    # all i7, i5
    # all no i7 all i5
    # all i7  no i5
    # no i7, no i5, all i1

    # forbidden are:
    # some have i7
    # some have i5
    # no i7 no i5, no i1

    # these are some booleans to make the following value checks a bit more human readable
    all_i7 = None not in i7_barcodes
    all_i5 = None not in i5_barcodes
    all_i1 = None not in i1_barcodes

    samples_no_i7, some_i7 = some_dont_have_barcode(i7_barcodes)
    samples_no_i5, some_i5 = some_dont_have_barcode(i5_barcodes)
    samples_no_i1, some_i1 = some_dont_have_barcode(i1_barcodes)

    no_i7 = none_have_barcode(i7_barcodes)
    no_i5 = none_have_barcode(i5_barcodes)
    no_i1 = none_have_barcode(i1_barcodes)

    # demultiplexing on i7 and i5 (and maybe i1). I1 is optional when i7 and i5 are
    # already specified
    #if all([all_i7, all_i5]) and no_i1:
    if all([all_i7, all_i5]):
        return True
    # demultiplexing on i7 and/or i1
    if all([all_i7, no_i5]) and any([all_i1, some_i1]):
        return True
    # demultiplexing on i5 and/or i1
    if all([no_i7, all_i5]) and any([all_i1, some_i1]):
        return True
    # demultiplexing on i1 only
    if all([no_i7, no_i5]) and all_i1:
        return True

    error_messages = []
    # the following things are not allowed because they eventually dont allow assigning
    # reads unambiguously to samples
    if some_i7:
        error_msg = ("Not all samples have an i7 barcode defined. An i7 barcode needs "
                     "to be either specified for all or none of the samples. Samples "
                     "without a barcode: %s" % samples_no_i7)
        error_messages.append(error_msg)
    if some_i5:
        error_msg = ("Not all samples have an i5 barcode defined. An i5 barcode needs "
                     "to be either specified for all or none of the samples. Samples "
                     "without a barcode: %s" % samples_no_i5)
        error_messages.append(error_msg)
    if some_i1:
        error_msg = ("Not all samples have an i7 barcode defined. An i7 barcode needs "
                     "to be either specified for all or none of the samples. Samples "
                     "without a barcode: %s" % samples_no_i1)
        error_messages.append(error_msg)
    # an empty sample sheet does nothing and is therefore disallowed
    if all([no_i7, no_i5, no_i1]):
        error_msg = ("No index sequences have been specified in the sample sheet. "
                     "Please specify some barcodes to run this tool.")
        error_messages.append(error_msg)
    # if we did not return error time!
    raise ValueError("\n".join(error_messages))


def none_have_barcode(barcoded_samples):
    # TODO: add documentation
    return None in barcoded_samples and len(barcoded_samples) == 1


def all_have_barcode(barcoded_samples):
    # TODO: add documentation
    return None not in barcoded_samples


def some_dont_have_barcode(barcoded_samples):
    # TODO: add documentation
    return barcoded_samples.get(None), None in barcoded_samples


def has_different_lengths(barcode_name, barcode_length):
    # TODO: add documentation
    if len(barcode_length) > 1:
        error_msg = ("%s barcodes with different length have been specified. Only one "
                     "length is allowed.\n Observed Length: %s" % barcode_name,
                     barcode_length)
        raise ValueError(error_msg)


def is_valid_barcode(barcode, all_barcodes, barcode_type, sample_name):
    # TODO: add documentation
    barcode_set = all_barcodes[barcode_type].get(len(barcode))
    if barcode not in barcode_set and not None:
        error_msg = ("Barcode %s for sample %s is no valid %s barcode.\n Please check "
                     "your sample sheet." % (barcode, sample_name, barcode_type))
        raise ValueError(error_msg)


@contextmanager
def get_pe_fastq(fq_gz_1, fq_gz_2):
    # TODO: add documentation
    try:
        fq_1_gz_stream = gzip.open(fq_gz_1, 'rt', encoding='utf-8')
        fq_2_gz_stream = gzip.open(fq_gz_2, 'rt', encoding='utf-8')
        fastq_1 = fastq_lines_to_reads(fq_1_gz_stream)
        fastq_2 = fastq_lines_to_reads(fq_2_gz_stream)
        yield zip(fastq_1, fastq_2)
    finally:
        for fq_file in [fq_1_gz_stream, fq_2_gz_stream]:
            log.debug('Trying to close file: %s' % fq_file.name)
            try:
                fq_file.close()
            except IOError:
                log.exception("Could not close the following input file: %s"
                              % fq_file.name)


def fastq_lines_to_reads(fastq_lines):
    """Collect fastq lines and chucks them into reads.
    Args:
        fastq_lines (iterable <str>): An iterable of lines from a fastq file.
    Return:
        read_iterator (iterator): All fastq lines grouped into reads (chunks of 4 lines).
    """
    # see grouper example https://docs.python.org/3/library/itertools.html
    # convert fastq lines into chunked iterator with a length of lines_per_read
    lines_per_read = 4
    chunked_lines = [iter(fastq_lines)] * lines_per_read
    # zip longest calls next when zipping the chunked iterator
    # this leads to all fastq_lines grouped by chunks of size lines_per_read
    return zip_longest(*chunked_lines)


def load_correction_map(barcode_type, barcode_length):
    # TODO: update documentation
    """Reads a tsv barcodes file and returns them as a dict. These dicts are used to
    map erroneous barcodes to their corrected version.
    Args:
        mapping_file (string): A path pointing to a error correction map file.
    Return:
        dict (dict): correction_map : <erroneous_barcode, corrected_barcode>
    """

    package = "idemux.resources.barcodes.%s" % barcode_type
    mapping_file = "correction_map_b96_l%s.tsv" % barcode_length

    mapping = dict()
    log.debug("Loading error correction map from %s", mapping_file)
    with resources.open_text(package, mapping_file) as tsv_file:
        # with open(mapping_file, "r") as tsv_file:
        reader = csv.reader(tsv_file, delimiter='\t')
        for row in reader:
            mapping[row[0]] = row[1]
    return mapping


def load_correction_maps(barcodes):
    # def load_i7_i5_correction_maps(barcode_types=["i7", "i5"], lengths=[8, 10, 12]):
    """Reads i7 and i5 correction maps of the specified types
       Args:
           barcode_types (dict): An dictionary mapping the barcode to an
            iterable of containing the length. Eg. {"i7" : [8,10,12], "i5" :[8] }
       Return:
           dict (dict): correction_map : A nested dict containing the
             correction maps. Format: <barcode_type: length : correction_map>
    """
    barcode_maps = defaultdict(dict)
    # TODO: check if specific case for None needs to be added
    for barcode_type, lengths in barcodes.items():
        for length in lengths:
            correction_map = load_correction_map(barcode_type, length)
        barcode_maps[barcode_type][length] = correction_map
    return barcode_maps


def get_i1_barcodes():
    # TODO update documentation
    """Reads i1 barcode file and returns them as a set. This set contains all
    known/allowed i1 inline barcodes. This set is used to check if the correct i1 barcodes
    have been supplied.
    Args:
        file_path (string): A path pointing the i1 barcode file.
    Return:
        set (String): i1_set : All valid Lexogen i1 inline barcodes.
    """
    package = "idemux.resources.barcodes.i1"
    resource = "i1_barcodes.tsv"
    log.debug("Loading i1 inline barcodes from resource")
    i1_barcodes = set()
    with resources.open_text(package, resource) as fin:
        for line in fin:
            i1_barcodes.add(line.strip())
    return i1_barcodes


def get_valid_barcodes():
    # TODO: add documentation
    barcode_sets = defaultdict(dict)
    all_barcodes = {"i7": [8, 10, 12], "i5": [8, 10, 12], "i1": [12]}
    #i1_length = [12]
    correction_maps = load_correction_maps(all_barcodes)
    for bc_type, bc_length_map in correction_maps.items():
        for length, correction_map in bc_length_map.items():
            barcode_sets[bc_type][length] = set(correction_map.values())
    #barcode_sets["i1"][i1_length] = get_i1_barcodes()
    return barcode_sets


def peek_into_fastq_files(fq_gz_1, fq_gz_2, has_i7, has_i5):
    log.info("Peeking into fastq files to check for barcode formatting errors")
    lines_to_check = 1000
    counter = 0
    try:
        with get_pe_fastq(fq_gz_1, fq_gz_2) as pe_reads:
            for mate_pair in pe_reads:
                check_fastq_headers(mate_pair, has_i7, has_i5)
                counter += 1
                if counter == lines_to_check:
                    break
        log.info("Fastq file formatting seems fine.")
    except Exception as e:
        log.exception(e)
        sys.exit(1)


def check_fastq_headers(mate_pair, has_i7, has_i5):
    # TODO: add documumentation
    # TODO: add read name checking
    header_idx = 0
    mate1, mate2 = mate_pair
    bcs_mate1 = len(mate1[header_idx].rpartition(":")[-1].split('+'))
    bcs_mate2 = len(mate2[header_idx].rpartition(":")[-1].split('+'))

    number_bc_present = [bcs_mate1, bcs_mate2]
    expected_number = sum([has_i7, has_i5])

    example_header_1 = ("@NB502007:379:HM7H2BGXF:1:11101:24585:1069 1:N:0:TCAGGTAANNTT")
    example_header_2 = ("@NB502007:379:HM7H2BGXF:1:11101:24585:1069 "
                        "1:N:0:TCAGGTAANNTT+NANGGNNCNNNN")

    right_number_of_barcodes = [n <= expected_number for n in number_bc_present]
    if not all(right_number_of_barcodes):
        example_header = example_header_2 if expected_number == 2 else example_header_1
        error_msg = ("The fastq file does not contain sufficient barcode information in "
                     "the header.\nExpected number of barcodes: %s\nObserved number of "
                     "barcodes: %s\nPlease check your input file. Your fastq header "
                     "should look similar to this example.\nExample: %s" %
                     (expected_number, number_bc_present, example_header))
        raise ValueError(error_msg)
