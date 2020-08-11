import csv
import gzip
import io
import logging
import sys
from collections import defaultdict, Counter
from contextlib import contextmanager
from itertools import zip_longest

# importlib is > 3.7 Try to import it and if this does not work import the backport
try:
    from importlib import resources
except ImportError:
    import importlib_resources as resources

log = logging.getLogger(__name__)


def parse_sample_sheet(cvs_file):
    """Function to parse the sample_sheet.csv.

    This function takes the path to a idemux sample sheet reads the data and does some
    sanity checks to prevent downstream problems. The csv needs to consist out of 4
    columns specifing a sample name and the respective barcodes. If a sample is only
    i7, i5 or i1 barcoded missing barcodes should be indicated by am empty field.
    In general the sample sheet should be formatted like this and have a header included:

    sample_name,i7,i5,i1
    <sample_name>,<i7_barcode or empty>,<i7_barcode or empty>,<i1_barcode or empty>

    Example:
    These are allowed
        triple indexing (all barcodes specified)
        sample_name,i7,i5,i1
        sample_0,AAAACATGCGTT,CCCCACTGAGTT,AAAACATGCGTT
        sample_1,AAAATCCCAGTT,CCCCTAAACGTT,AAAATCCCAGTT

        dual indexing(either i7+i5, or i7+i1 or i5+i1)
        sample_name,i7,i5,i1
        sample_0,AAAACATGCGTT,,AAAACATGCGTT
        sample_1,AAAATCCCAGTT,,AAAATCCCAGTT

        single indexing(either i7, i5 or i5)
        sample_name,i7,i5,i1
        sample_0,,,AAAACATGCGTT
        sample_1,,,AAAATCCCAGTT

        mixed i1 indexing(some samples have dual or triple indexing. Works only with i1!)
        sample_name,i7,i5,i1
        sample_0,AAAACATGCGTT,CCCCACTGAGTT,AAAACATGCGTT
        sample_1,AAAATCCCAGTT,CCCCTAAACGTT,

        different length for either all i7 (10 nt) or i5 (12 nt)
        sample_name,i7,i5,i1
        sample_0,AAAACATGCG,CCCCACTGAGTT,AAAACATGCGTT
        sample_1,AAAATCCCAG,CCCCTAAACGTT,AAAATCCCAGTT


    These are not allowed and will call sys.exit(1):
        mixed i7/i5 indexing
        sample_name,i7,i5,i1
        sample_0,AAAACATGCGTT,,AAAACATGCGTT
        sample_1,,CCCCTAAACGTT,AAAATCCCAGTT

        duplicate/non-unique barcodes
        sample_name,i7,i5,i1
        sample_0,AAAACATGCGTT,CCCCACTGAGTT,AAAACATGCGTT
        sample_1,AAAACATGCGTT,CCCCACTGAGTT,AAAACATGCGTT

        duplicate sample names
        sample_name,i7,i5,i1
        sample_0,AAAACATGCGTT,CCCCACTGAGTT,AAAACATGCGTT
        sample_0,AAAATCCCAGTT,CCCCTAAACGTT,AAAATCCCAGTT

        different length for one barcode
        sample_name,i7,i5,i1
        sample_0,AAAACATGCGTT,CCCCACTGAGTT,AAAACATGCGTT
        sample_1,AAAATCCCAG,CCCCTAAACGTT,AAAATCCCAGTT

        an i1 barcode shorter than 12 nt
        sample_name,i7,i5,i1
        sample_0,AAAACATGCGTT,CCCCACTGAGTT,AAAACATG
        sample_1,AAAATCCCAG,CCCCTAAACGTT,AAAATCCC

    Arguments:
        cvs_file (str): Path to a sample_sheet.cvs. See examples for formatting

    Returns:
        barcode_sample_map (dict): A dict mapping a tuple of barcodes to sample names
            (i7, i5, i1): sample_name
        used_barcodes (tuple): A tuple of 3 sets, each containing the i7, i5 and i1
            barcodes specified in the sample sheet
        barcode_lengths (dict): A dict mapping the index name to a list barcode
            lengths.  {"i7": list(int), "i5": list(int), "i1": list(int)}

    Except:
        ValueError: Will initiate sys.exit(1)
    """
    # i1 inline barcode is always 12 nt in length
    expected_i1_bc_length = 12
    # we use these to keep track which lengths are beeing used for each barcode
    i7_lengths, i5_lengths, i1_lengths = set(), set(), set()

    # we use these to keep track of which single match to which sample(s)
    i7_barcodes = defaultdict(list)
    i5_barcodes = defaultdict(list)
    i1_barcodes = defaultdict(list)

    # this is for keeping track which unique barcode combination is associated with a
    # sample
    barcode_sample_map = dict()
    # we read all valid Lexogen barcodes here, so we can check if the provided barcodes
    # are valid barcodes. We need to do this as we can only error correct and
    # demultiplex on Lexogen barcodes.
    valid_barcodes = get_valid_barcodes()
    # we use this for checking for duplicated sample names
    sample_count = Counter()
    # expected file header
    sample_sheet_header = ["sample_name", "i7", "i5", "i1"]
    try:
        with open(cvs_file) as sample_sheet:
            # we register a csv dialect to do whitespace trimming for us
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
                sample_count[sample_name] += 1

                # make sure the barcodes are actual Lexogen barcodes so we can error
                # correct them later on
                is_valid_barcode(i7_bc, valid_barcodes, "i7", sample_name)
                is_valid_barcode(i5_bc, valid_barcodes, "i5", sample_name)
                is_valid_barcode(i1_bc, valid_barcodes, "i1", sample_name)

                # add barcodes and sample_names to dict so we can do some value checking
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
        barcode_lengths = {"i7": list(i7_lengths),
                           "i5": list(i5_lengths),
                           "i1": list(i1_lengths)}
        has_different_lengths(barcode_lengths)
        # if me made it here we can be sure he list contains only one barcode
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
        duplicated_sample_names = [k for k, v in sample_count.items() if v > 1]
        if duplicated_sample_names:
            error_msg = ("The sample sheet contains duplicate sample names. Sample "
                         "names have to be unique.\n Sample names: %s" %
                         duplicated_sample_names)
            raise ValueError(error_msg)

    except ValueError as e:
        log.exception(str(e))
        sys.exit(1)
    # we have made it until here until raising an exception. That means the sample sheet
    # information should be okay and we can return required data from the sample sheet.
    i7_used = set(i7_barcodes.keys())
    i5_used = set(i5_barcodes.keys())
    i1_used = set(i1_barcodes.keys())

    used_barcodes = (i7_used, i5_used, i1_used)
    return barcode_sample_map, used_barcodes, {"i7": list(i7_lengths),
                                               "i5": list(i5_lengths),
                                               "i1": list(i1_lengths)}


def has_valid_barcode_combinations(i7_barcodes, i5_barcodes, i1_barcodes):
    """Function that checks if the provided barcodes allow unique sample identification
    for the different usecases.

    Allowed barcoding variants are:
    - all contain an i7 and i5
    - all contain no i7 all an i5
    - all contain an i7 no an i5
    - all contain no i7, no i5, all an i1

    Forbidden barcoding variants are:
    - some have an i7
    - some have an i5
    - no i7 no i5, no i1

    Args:
        i7_barcodes (iterable): i7 barcodes of all samples
        i5_barcodes (iterable): i5 barcodes of all samples
        i1_barcodes (iterable): i1 barcodes of all samples

    Returns:
        is_valid (bool): True when valid barcode combinations are supplied (see above)

    Raises:
        ValueError (err): When supplied an invalid barcode combination (see above)
    """

    # missing values in the csv retruned as None so we want to check fo these
    all_i7 = None not in i7_barcodes
    all_i5 = None not in i5_barcodes
    all_i1 = None not in i1_barcodes

    samples_no_i7, some_i7 = some_dont_have_barcode(i7_barcodes)
    samples_no_i5, some_i5 = some_dont_have_barcode(i5_barcodes)
    samples_no_i1, some_i1 = some_dont_have_barcode(i1_barcodes)

    no_i7 = none_have_barcode(i7_barcodes)
    no_i5 = none_have_barcode(i5_barcodes)
    no_i1 = none_have_barcode(i1_barcodes)

    # allowed cases, we dont actually need these return values, but these make the logic
    # much more obvious and human readable.
    # demultiplexing on i7 and i5 (and maybe i1). I1 is optional when i7 and i5 are
    # already specified
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
    # if we did not return true earlier it is error raising time!
    raise ValueError("\n".join(error_messages))


def none_have_barcode(barcoded_samples):
    """Just a short alias function to check if another value than None is a <barcode:
    sample> map.

     Args:
         barcoded_samples (dict): A dict mapping one type of barcodes (e.g. i7) to a
             sample name

    Returns:
        only_none (bool): Is there another value next to None?
    """
    return None in barcoded_samples and len(barcoded_samples) == 1


def all_have_barcode(barcoded_samples):
    """Just a short alias to function to check if None values are present in a <barcode:
    sample> map.

     Args:
         barcoded_samples (dict): A dict mapping one type of barcode (e.g. i7) to a
             sample name

    Returns:
        bool: Is there no None present?
    """
    return None not in barcoded_samples


def some_dont_have_barcode(barcoded_samples):
    """Just a short alias to function to check if None values are present in a <barcode:
    sample> map and to return the sample names with None values.

     Args:
         barcoded_samples (dict): A dict mapping one type of barcode (e.g. i7) to a
             sample name

    Returns:
        bool: Is there a None in barcoded_samples?
    """
    return barcoded_samples.get(None), None in barcoded_samples


def has_different_lengths(barcode_lengths):
    """Checks if a barcode name maps to a list with more than one element.

     Args:
         barcode_lengths (dict): A dict mapping a barcode name to a list of seen
             barcode lengths e.g. {"i5",[8,12]}.

    Raises:
        ValueError: When the list contains more than one item.
    """
    for name, lengths in barcode_lengths.items():
        if len(lengths) > 1:
            error_msg = ("%s barcodes with different length have been specified. Only "
                         "one length is allowed.\n Observed Lengths: %s" % name,
                         lengths)
            raise ValueError(error_msg)


def is_valid_barcode(barcode, valid_barcodes, barcode_type, sample_name):
    """Checks if the barcode associated with sample is a valid Lexogen barcode.

    Args:
        barcode (byte): A barcode represented as byte string
        valid_barcodes (dict): A nested dict mapping barcode names to length, to valid
            barcodes as byte strings. E.g "i7" : 12 : set(bytes)
        barcode_type (str): Name of the barcode type, e.g "i7"
        sample_name (str): Sample name corresponding to "barcode"

    Raises:
        ValueError: When the barcode is not in valid barcodes
    """
    barcode_set = valid_barcodes[barcode_type].get(len(barcode))
    if barcode not in barcode_set and not None:
        error_msg = ("Barcode %s for sample %s is no valid %s barcode.\n Please check "
                     "your sample sheet." % (barcode, sample_name, barcode_type))
        raise ValueError(error_msg)


@contextmanager
def get_pe_fastq(fq_gz_1, fq_gz_2):
    """Generator that opens two paired fastq.gz files and reads them as byte strings.

    Args:
        fq_gz_1: Path to read1 of the mate pair.
        fq_gz_2: Path to read2 of the mate pair.

    Yields:
        iter,iter: Zipped iterator for read1 and read2
    """
    try:
        input_streams = list()
        input_streams.append(gzip.open(fq_gz_1, 'rt', encoding='utf-8'))
        input_streams.append(gzip.open(fq_gz_2, 'rt', encoding='utf-8'))
        fastq_1 = fastq_lines_to_reads(input_streams[0])
        fastq_2 = fastq_lines_to_reads(input_streams[1])
        yield zip(fastq_1, fastq_2)
    finally:
        for fq_file in input_streams:
            log.debug('Trying to close file: %s' % fq_file.name)
            try:
                fq_file.close()
            except IOError:
                log.exception("Could not close the following input file: %s"
                              % fq_file.name)


def fastq_lines_to_reads(fastq_lines):
    """Collect fastq lines and chucks them into reads. Base on itertools
    grouper example (see https://docs.python.org/3/library/itertools.html).

    Args:
        fastq_lines (iterable <str>): An iterable of lines from a fastq file.

    Return:
        iter: All fastq lines grouped into reads (chunks of 4 lines).
    """
    # convert fastq lines into chunked iterator with a length of lines_per_read
    lines_per_read = 4
    chunked_lines = [iter(fastq_lines)] * lines_per_read
    # zip longest calls next when zipping the chunked iterator
    # this leads to all fastq_lines grouped by chunks of size lines_per_read
    return zip_longest(*chunked_lines)


def load_correction_map(barcode_type, barcode_length):
    """Reads a tsv barcodes file, converts them to byes and returns them as a dict.
    These dicts are used to map erroneous barcodes to their corrected version.
    Args:
        barcode_type (str): A path pointing to a error correction map file.
        barcode_length (int): A path pointing to a error correction map file.

    Return:
        dict (dict): correction_map : <b'erroneous_barcode', b'corrected_barcode'>
    """
    mapping = dict()
    package = "idemux.resources.barcodes.%s" % barcode_type
    mapping_file = "correction_map_b96_l%s.tsv" % barcode_length
    log.debug("Loading error correction map from %s", mapping_file)
    with resources.open_text(package, mapping_file) as tsv_file:
        # with open(mapping_file, "r") as tsv_file:
        reader = csv.reader(tsv_file, delimiter='\t')
        for row in reader:
            mapping[row[0]] = row[1]
    return mapping


def load_correction_maps(barcodes):
    """Reads i7 and i5 correction maps of the specified types.

       Args:
           barcodes (dict): An dictionary mapping the barcode to an iterable of
           lengths. Eg. {"i7" : [8,10,12],"i5" : [8]}
       Return:
           dict (dict): correction_map : A nested dict containing the
           correction maps. Format: <barcode_type: length : correction_map>
    """
    barcode_maps = defaultdict(dict)
    for barcode_type, lengths in barcodes.items():
        for length in lengths:
            correction_map = load_correction_map(barcode_type, length)
            barcode_maps[barcode_type][length] = correction_map
    return barcode_maps


def get_valid_barcodes():
    """Loads all valid Lexogen barcodes.

       Return:
           dict (dict): A nested dict mapping barcode name to, length, to a set of
               valid barcodes. Format: <barcode_type: length : lexogen barcodes>
    """
    barcode_sets = defaultdict(dict)
    all_barcodes = {"i7": [8, 10, 12], "i5": [8, 10, 12], "i1": [12]}
    correction_maps = load_correction_maps(all_barcodes)
    for bc_type, bc_length_map in correction_maps.items():
        for length, correction_map in bc_length_map.items():
            barcode_sets[bc_type][length] = set(correction_map.values())
    return barcode_sets


def peek_into_fastq_files(fq_gz_1, fq_gz_2, has_i7, has_i5):
    """Reads the first 100 lines of paired fastq.gz files and checks if everything is
    okay with the fastq header format.

    Args:
        fq_gz_1 (str): File path of read mate1.
        fq_gz_2 (str): File path of read mate2.
        has_i7 (bool): Did the sample_sheet specify that samples have an i7 index?
        has_i5 (bool): Did the sample_sheet specify that samples have an i5 index?

    Raises:
        ValueError: When the fastq header contains less barcodes than indicated by the
            booleans.
    """
    log.info("Peeking into fastq files to check for barcode formatting errors")
    lines_to_check = 100
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
    """Function to check if the barcodes (i7,i5) specified in the sample sheet are
    as well in the fastq header.

    Args:
        mate_pair (tuple): A tuple of mate_pairs as returned by fastq_lines_to_reads.
        has_i7 (bool): Did the sample_sheet specify that samples have an i7 index?
        has_i5 (bool): Did the sample_sheet specify that samples have an i5 index?

    Raises:
        ValueError: When the fastq header contains less barcodes than indicated by the
            booleans.
    """
    header_idx = 0
    mate1, mate2 = mate_pair
    header_mate_1, header_mate_2 = mate1[header_idx], mate2[header_idx]
    # get the barcodes from the fastq header
    bcs_mate1 = len(header_mate_1.rpartition(":")[-1].split("+"))
    bcs_mate2 = len(header_mate_2.rpartition(":")[-1].split("+"))

    number_bc_present = [bcs_mate1, bcs_mate2]
    expected_number = sum([has_i7, has_i5])
    # this is how a fastq header should look like
    example_header_1 = ("@NB502007:379:HM7H2BGXF:1:11101:24585:1069 1:N:0:TCAGGTAANNTT")
    example_header_2 = ("@NB502007:379:HM7H2BGXF:1:11101:24585:1069 "
                        "1:N:0:TCAGGTAANNTT+NANGGNNCNNNN")

    # check if the header conforms to what was specified in the sample sheet
    right_number_of_barcodes = [n <= expected_number for n in number_bc_present]
    if not all(right_number_of_barcodes):
        example_header = example_header_2 if expected_number == 2 else example_header_1
        error_msg = ("The fastq file does not contain sufficient barcode information in "
                     "the header.\nExpected number of barcodes: %s\nObserved number of "
                     "barcodes: %s\nPlease check your input file. Your fastq header "
                     "should look similar to this example.\nExample: %s\n Observed "
                     "headers: %s" %
                     (expected_number, number_bc_present, example_header,
                      [header_mate_1, header_mate_2])
                     )
        raise ValueError(error_msg)
