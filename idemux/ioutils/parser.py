
"""Parser module - handles all the reading"""

import csv
import gzip
import logging
import sys
from collections import defaultdict, Counter
from contextlib import contextmanager
from itertools import zip_longest

# importlib is > 3.7 Try to import it and if this does not work import the backport
from idemux.ioutils.barcode import Barcode

try: # pragma: no cover
    from importlib import resources
except ImportError: # pragma: no cover
    import importlib_resources as resources

log = logging.getLogger(__name__)


def parse_sample_sheet(sample_sheet, i5_rc, **kwargs):
    # TODO: implement regex check to determine if barcodes are DNA bases
    # TODO: convert strings to uppercase when parsed
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
        sample_sheet (str): Path to a sample_sheet.cvs. See examples for formatting
        i5_rc (bool): Should the i5 indices be reverse complemented?

    Returns:
        dict: A dict mapping a tuple of barcodes to sample names (i7, i5, i1): sample_name
        tuple (Barcodes): A tuple of 3 barcode dataclasses, each containing the i7, i5
            and i1 barcode

    Except:
        ValueError: Will initiate sys.exit(1)
    """
    # we use these to keep track which lengths are beeing used for each barcode
    i7_lengths, i5_lengths, i1_lengths = set(), set(), set()

    # we use these to keep track of which single barcode matches to which sample(s)
    i7_barcodes = defaultdict(list)
    i5_barcodes = defaultdict(list)
    i1_barcodes = defaultdict(list)

    # this is for keeping track which unique barcode combination belongs to a sample
    barcode_sample_map = dict()
    # we use this for checking for duplicated sample names
    sample_count = Counter()
    # expected file header
    sample_sheet_header = ["sample_name", "i7", "i5", "i1"]
    try:
        with open(sample_sheet) as sample_sheet:
            # we register a csv dialect to do whitespace trimming for us
            csv.register_dialect('strip', skipinitialspace=True)
            lines_read = 0
            reader = csv.DictReader(sample_sheet, restval=None, dialect='strip')

            for row in reader:
                lines_read += 1
                # check if all the needed columns are present. otherwise we want to
                # disallow further processing
                if not set(sample_sheet_header) == row.keys():
                    s1 = ",".join(sample_sheet_header)
                    s2 = ",".join(row.keys())
                    error_msg = ("Incorrect sample sheet header. Expected header: %s\n"
                                 "Observed header: %s" % (s1, s2))
                    raise ValueError(error_msg)

                sample_name = row.get('sample_name')
                # undefined sample names are not allowed
                if not sample_name:
                    error_msg = ("Incorrect sample name in line %s. Empty values are "
                                 "not allowed as sample names." % lines_read)
                    raise ValueError(error_msg)
                # we use this to keep track of duplicate sample names
                sample_count[sample_name] += 1
                # as ordered dicts always return an empty string instead of None
                _i7 = row.get('i7')
                i7_bc = _i7 if _i7 else None
                # i5 can be sequenced as reverse complement, translate if needed
                _i5 = row.get('i5') if not i5_rc else reverse_complement(row.get('i5'))
                i5_bc = _i5 if _i5 else None
                _i1 = row.get('i1')
                i1_bc = _i1 if _i1 else None
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
                i7_lengths.add(i7_bc if i7_bc is None else len(i7_bc))
                i5_lengths.add(i5_bc if i5_bc is None else len(i5_bc))
                i1_lengths.add(i1_bc if i1_bc is None else len(i1_bc))
        # putting everything in a data class to make it tidier.
        # TODO: add addtional methods to data class and move this into the loop
        i7 = Barcode('i7', i7_barcodes)
        i5 = Barcode('i5', i5_barcodes, i5_rc)
        i1 = Barcode('i1', i1_barcodes)
        barcodes = (i7, i5, i1)

        # test if the supplied barcode combinations are valid
        has_valid_barcode_combinations(*barcodes)
        # sample names have to be unique as they determine the outfile name. Otherwise
        # we get problems when we try to write reads belonging to different barcode
        # combinations to one file. We do this this way so we throw more accurate errors
        duplicated_sample_names = [k for k, v in sample_count.items() if v > 1]
        if duplicated_sample_names:
            error_msg = ("The sample sheet contains duplicate sample names. Sample "
                         "names have to be unique.\n Sample names: %s" %
                         duplicated_sample_names)
            raise ValueError(error_msg)

    except ValueError as e:
        sys.exit(e)
    # we have made it until here until raising an exception. That means the sample sheet
    # information should be okay and we can return required data from the sample sheet.
    barcodes = [load_correction_map(bc) for bc in barcodes]
    return barcode_sample_map, barcodes


def reverse_complement(sequence):
    """Function that returns the reverse complement of DNA sequence. Accepts A,C,T,G,N
    as input bases.

    Args:
        sequence (str): A DNA string that should be translated to its reverse complement

    Returns:
        str: Reverse complement of the input. Returns None when
            sequence is None

    Raises:
        ValueError: Is raised when the input string contains other letters than A,C,T,G,N
    """
    # if sequence is None, no need to do any work
    if sequence is None:
        return None

    # complements of each base
    base_complements = {'A': 'T',
                        'C': 'G',
                        'G': 'C',
                        'T': 'A',
                        'N': 'N'}
    try:
        # get the reverse complement
        rc_sequence = [base_complements[nt] for nt in sequence[::-1]]
        return "".join(rc_sequence)
    except KeyError:
        # we only get a key error if sequence contains letter that are not covered
        # above
        invalid_bases = [nt for nt in sequence[::-1] if nt not in base_complements]

        error_message = ("The following barcode sequence from the sample sheet "
                         "contains bases that cant be mapped to their reverse "
                         "complement.\nBarcodes: %s\nBases %s" %
                         (sequence, invalid_bases))
        raise ValueError(error_message)


def has_valid_barcode_combinations(i7, i5, i1, *args):
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
        i7 (Barcode): Barcode dataclass of i7 barcodes
        i5 (Barcode): Barcode dataclass of i5 barcodes
        i1 (Barcode): Barcode dataclass of i1 barcodes

    Returns:
        bool: True when valid barcode combinations are supplied (see above)

    Raises:
        ValueError: When supplied an invalid barcode combination (see above)
    """
    # allowed cases, we dont actually need these return values, but these make the logic
    # much more obvious and human readable.
    # demultiplexing on i7 and i5 (and maybe i1). I1 is optional when i7 and i5 are
    # already specified so we dont really need to check for sparsity
    if i7.full and i5.full:
        return True
    # demultiplexing on i7 and/or i1
    if i7.full and i5.empty:
        return True
    # demultiplexing on i5 and/or i1
    if i7.empty and i5.full:
        return True
    # demultiplexing on i1 only
    if i7.empty and i5.empty and i1.full:
        return True

    error_messages = []
    # the following things are not allowed because they eventually dont allow assigning
    # reads unambiguously to samples
    if i7.sparse:
        error_msg = ("Not all samples have an i7 barcode defined. An i7 barcode needs "
                     "to be either specified for all or none of the samples. Samples "
                     "without a barcode: %s" % i7.samples_without_barcodes)
        error_messages.append(error_msg)
    if i5.sparse:
        error_msg = ("Not all samples have an i5 barcode defined. An i5 barcode needs "
                     "to be either specified for all or none of the samples. Samples "
                     "without a barcode: %s" % i5.samples_without_barcodes)
        error_messages.append(error_msg)
    if i1.sparse:
        error_msg = ("Not all samples have an i7 barcode defined. An i7 barcode needs "
                     "to be either specified for all or none of the samples. Samples "
                     "without a barcode: %s" % i1.samples_without_barcodes)
        error_messages.append(error_msg)
    # an empty sample sheet does nothing and is therefore disallowed
    if all([i7.empty, i5.empty, i1.empty]):
        error_msg = ("No index sequences have been specified in the sample sheet. "
                     "Please specify some barcodes to run this tool.")
        error_messages.append(error_msg)
    # if we did not return true earlier it is error raising time!
    raise ValueError("\n".join(error_messages))


@contextmanager
def get_pe_fastq(fq_gz_1, fq_gz_2):
    """Generator that opens two paired fastq.gz files and yields them as utf-8 strings.

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


def get_map_from_resource(package, resource):
    """Load a 2 column tsv file from the specified package and return as dict.
    Args:
        package (str): Package to load resourcess from
        resource (str): Name of the tsv to load

    Return:
        dict(str,str): A dict mapping the first and second column.
    """
    log.debug("Loading error correction map from %s", resource)
    mapping = dict()
    with resources.open_text(package, resource) as tsv_file:
        reader = csv.reader(tsv_file, delimiter='\t')
        for row in reader:
            mapping[row[0]] = row[1]
    return mapping


def load_correction_map(barcode):
    """Adds the right correction map to the provided barcode object returns them as a
    dict. Correction maps are used to map erroneous barcodes to their corrected version.
    Args:
        barcode (Barcode): A Barcode data class object

    Return:
        Barcode: Updated data class with correction map set.
    """
    log.info(f"Trying to find the appropriate barcode set for {barcode.name}...")
    # When there are no barcodes specified, there is nothing to correct.
    if barcode.length == 0:
        log.info(f"No barcodes have been specified for {barcode.name}.")
        barcode.correction_map = {None: None}
        return barcode

    # barcodes can be bough as set of 96 and 348. 348 can also be not fully sequenced
    # and are then ident with the 96 set. Therefore we need to check from which set
    # the provided barcodes originate. We then construct a string describing the
    # proper resource to load.
    for set_size in barcode.get_set_sizes():
        package_str = f"idemux.resources.barcodes.{barcode.name}"
        file_str = f"base_mapping_b{set_size}_l{barcode.length}.tsv"

        corr_map = get_map_from_resource(package_str, file_str)
        _barcode_set = set(corr_map.values())
        # we don't want Nones here as this will mess up the subset testing, as the error
        # correction maps don't contain Nones
        _barcodes_given = barcode.get_used_codes(drop_none=True)
        if _barcodes_given <= _barcode_set:
            log.info(f"Correct set found. Used set is {set_size} barcodes with "
                     f"{barcode.length} nt length.")
            barcode.correction_map = corr_map
            return barcode
    # we want the user to know when non Lexogen barcodes are supplied as no no error
    # correction is happening. However, we suppress it for length 6 as these cant
    # be error corrected and we don't want to incent insecurity when using old 6 nt
    # barcodes
    if barcode.length != 6:
        log.warning(f"No fitting Lexogen barcode set found for {barcode.name}. No "
                    f"error correction will take place for this barcode. Are you using "
                    f"valid Lexogen barodes?")
    # if there is map for error correction, return a map that maps the given barcodes
    # to itself
    _bc_list = list(barcode.used_codes)
    barcode.correction_map = dict(zip(_bc_list, _bc_list))
    return barcode


def peek_into_fastq_files(fq_gz_1, fq_gz_2, has_i7, has_i5, has_i1, i7_length,
                          i5_length, i1_start, i1_end, **kwargs):
    """Reads the first 1000 lines of paired fastq.gz files and checks if everything is
    okay with the fastq header format.

    Args:
        fq_gz_1 (str): File path of read mate1.
        fq_gz_2 (str): File path of read mate2.
        has_i7 (bool): Did the sample_sheet specify that samples have an i7 index?
        has_i5 (bool): Did the sample_sheet specify that samples have an i5 index?
        has_i1 (bool): Did the sample_sheet specify that samples have an i1 index?
        i7_length (int): i7 barcode length, indirectly specified via the sample sheet.
        i5_length (int): i5 barcode length, indirectly specified via the sample sheet.
        i1_start (int): Start position of the i1 inline barcode
        i1_end (int): End position of the i1 inline barcode
    """
    log.info("Peeking into fastq files to check for barcode formatting errors")
    lines_to_check = 1000
    counter = 0
    log.info("Checking fastq input files...")

    with get_pe_fastq(fq_gz_1, fq_gz_2) as pe_reads:
        for mate_pair in pe_reads:
            # check if the formatting is okay
            check_mate_pair(mate_pair, has_i7, has_i5, has_i1, i7_length, i5_length,
                            i1_start, i1_end)
            # TODO: could be replaced with enumerate
            counter += 1
            if counter == lines_to_check:
                break
    log.info("Input file formatting seems fine.")


def check_mate_pair(mate_pair, has_i7, has_i5, has_i1, i7_length, i5_length,
                    i1_start, i1_end):
    """Reads the first 1000 lines of paired fastq.gz files and checks if everything is
       okay with the fastq header format.

       Args:
           mate_pair (tuple): A tuple of mate_pairs as returned by fastq_lines_to_reads.
           has_i7 (bool): Did the sample_sheet specify that samples have an i7 index?
           has_i5 (bool): Did the sample_sheet specify that samples have an i5 index?
           has_i1 (bool): Did the sample_sheet specify that samples have an i1 index?
           i7_length (int): i7 barcode length, indirectly specified via the sample sheet.
           i5_length (int): i5 barcode length, indirectly specified via the sample sheet.
           i1_start (int): Start position of the i1 inline barcode
           i1_end (int): End position of the i1 inline barcode

    Except:
        ValueError: Will initiate sys.exit(1)
       """

    _, mate2 = mate_pair
    try:
        check_fastq_headers(mate_pair, has_i7, has_i5, i7_length, i5_length)
        if has_i1:
            check_mate2_length(mate2, i1_start, i1_end)
    except Exception as e:
        log.info(e)
        sys.exit(e)


def check_mate2_length(mate2, i1_start, i1_end):
    """Check if the mate2 sequence is long enough to contain the specified i1 barcodes.

    Args:
        mate2 (tuple): A fastq read (4 lines).
        i1_start (int): Start position of the i1 inline barcode
        i1_end (int): End position of the i1 inline barcode

    Raises:
        ValueError: When the fastq header contains less barcodes than indicated by the
            booleans.
    """

    seq_idx = 1
    seq = mate2[seq_idx]
    if len(seq[:-1]) < i1_end:
        raise ValueError(f"Mate 2 is too short for the provided i1 barcode settings. "
                         f"According to your settings i1 starts at position {i1_start} "
                         f"and has a length of {i1_end - i1_start}. The sequence of "
                         f"mate 2 is however only {len(seq[:-1])} nt long.")


def check_fastq_headers(mate_pair, has_i7, has_i5, i7_length, i5_length):
    """Function to check if the barcodes (i7,i5) specified in the sample sheet are
    as well in the fastq header.

    Args:
        mate_pair (tuple): A tuple of mate_pairs as returned by fastq_lines_to_reads.
        has_i7 (bool): Did the sample_sheet specify that samples have an i7 index?
        has_i5 (bool): Did the sample_sheet specify that samples have an i5 index?
        i7_length (int): i7 barcode length, indirectly specified via the sample sheet.
        i5_length (int): i5 barcode length, indirectly specified via the sample sheet.

    Raises:
        ValueError: 1) When the fastq header contains less barcodes than indicated by the
            booleans. 2) When the barcodes have a different length than specified in the
            sample sheet.
    """
    header_idx = 0
    m_1, m_2 = mate_pair
    header_mate_1, header_mate_2 = m_1[header_idx], m_2[header_idx]
    # get the barcodes from the fastq header if present.
    # first get the last element. This is either a barcode, or if no barcodes are present
    # a digit describing the sample number
    # see https://support.basespace.illumina.com/articles/descriptive/fastq-files/
    bcs_mate1 = header_mate_1.strip().rpartition(":")[-1]
    # when the last element is only numeric it is not a barcode. if it is not numeric
    # it can be either one or two barcodes. these are normally seperated by a  +
    bcs_mate1 = None if bcs_mate1.isnumeric() else bcs_mate1.split("+")
    bcs_mate2 = header_mate_2.strip().rpartition(":")[-1]
    bcs_mate2 = None if bcs_mate2.isnumeric() else bcs_mate2.split("+")

    # when two mates have have different barcodes, the fastq files ist probably not sorted
    # this will cause trouble and should not be allowed
    if bcs_mate1 != bcs_mate2:
        error_msg = (f"Mate1 and mate2 contain different barcode information. Please "
                     f"make sure the reads in your fastq files are paired.\n"
                     f"Mate1 header: {header_mate_1}\n"
                     f"Mate2 header: {header_mate_2}")
        raise ValueError(error_msg)

    # when there are no barcodes present, set the number to 0
    number_bc_m1 = 0 if bcs_mate1 is None else len(bcs_mate1)
    number_bc_m2 = 0 if bcs_mate2 is None else len(bcs_mate2)

    number_bc_present = [number_bc_m1, number_bc_m2]
    expected_number = sum([has_i7, has_i5])
    # this is how a fastq header should look like
    example_header_1 = ("@NB502007:379:HM7H2BGXF:1:11101:24585:1069 1:N:0:TCAGGTAANNTT")
    example_header_2 = ("@NB502007:379:HM7H2BGXF:1:11101:24585:1069 "
                        "1:N:0:TCAGGTAANNTT+NANGGNNCNNNN")

    # check if the header conforms to what was specified in the sample sheet and has
    # at lest the number of barcodes specified in the sample sheet
    right_number_of_barcodes = [expected_number <= n for n in number_bc_present]
    if not all(right_number_of_barcodes):
        example_header = example_header_2 if expected_number == 2 else example_header_1
        error_msg = (f"The fastq file does not contain sufficient barcode information "
                     f"in the header.\nExpected number of barcodes: {expected_number}\n"
                     f"Observed number of barcodes: {number_bc_present}\n"
                     f"Please check your input file. Your fastq header should look "
                     f"similar to this example.\n"
                     f"Example: {example_header}\n"
                     f"Observed headers: {[header_mate_1, header_mate_2]}")
        raise ValueError(error_msg)

    # check if barcodes specified in the sample sheet and observed in the fastq file
    # are equally long.

    # when there are 2 barcodes in the fastq header the orientation is i7,i5
    if has_i7 and has_i5:
        if len(bcs_mate1[0]) != i7_length or len(bcs_mate1[1]) != i5_length:
            raise ValueError(f"i7 and i5 have a different length than specified in the "
                             f"sample_sheet. "
                             f"Observed length(i7,i5): {len(bcs_mate1[0])}"
                             f",{len(bcs_mate1[1])}\n "
                             f"Expected length(i7,i5): {i7_length},{i5_length}")
    # when there is 1 barcodes in the fastq header we need to check which one it is
    if has_i7 and not has_i5:
        if len(bcs_mate1[0]) != i7_length:
            raise ValueError(f"i7 has a different length than specified in the "
                             f"sample_sheet. "
                             f"Observed length(i7): {len(bcs_mate1[0])}\n"
                             f"Expected length(i7): {i7_length}\n")
    if not has_i7 and has_i5:
        if len(bcs_mate1[0]) != i5_length:
            raise ValueError(f"i5 has a different length than specified in the "
                             f"sample_sheet. "
                             f"Observed length(i5): {len(bcs_mate1[0])}\n"
                             f"Expected length(i5): {i5_length}\n")
