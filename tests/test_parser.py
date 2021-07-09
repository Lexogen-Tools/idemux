#!/usr/bin/env python
"""Tests for `idemux` package."""
from collections import defaultdict, namedtuple

import pytest
from os import path
import pathlib

from idemux.ioutils.barcode import Barcode
from idemux.ioutils.parser import parse_sample_sheet, load_correction_map, \
    reverse_complement, peek_into_fastq_files, get_pe_fastq

RESOURCE = pathlib.Path(path.dirname(__file__)) / "resources"

CSVS_TO_FAIL = RESOURCE.glob("sample_sheet/fail/*.csv")
CSVS_TO_PASS = RESOURCE.glob("sample_sheet/pass/*.csv")

BARCODES_6NT = {"A" * 6: ["test0"], "C" * 6: ["test1"]}
BARCODES_NON_LEX = {"A" * 12: ["test0"], "C" * 12: ["test1"]}
BARCODES_NONE = {None: ["test0", "test1"]}

BARCODES_TO_TEST = [Barcode("i7", BARCODES_6NT),
                    Barcode("i7", BARCODES_NON_LEX),
                    Barcode("i7", BARCODES_NONE)]

FOR_RC_PASS = [("AAAACATCGTTN", "NAACGATGTTTT"),
               (None, None),
               ("", "")]
FOR_RC_FAIL = ["AAAACATCGTTNXXXXX", 10]

FQ_RES = RESOURCE / "fastq"
INDEX_LENGTH = 12
I1_START = 10
Fq_data = namedtuple('Fq_data', ['condition',
                                 'fq_gz_1', 'fq_gz_2',
                                 'has_i7', 'has_i5', 'has_i1',
                                 'i7_length', 'i5_length', 'i1_start', 'i1_end'])
FQ_TO_PASS = [
    # Here we test a few default usecases that need idemux need to be able to handle
    # i7,i5,i1 specified and present
    Fq_data("i7(12)nt, i5(12)nt, i1(12)nt",
            FQ_RES / "i7_i5_i1_read_1.fastq.gz", FQ_RES / "i7_i5_i1_read_2.fastq.gz",
            True, True, True,
            INDEX_LENGTH, INDEX_LENGTH, I1_START, I1_START + INDEX_LENGTH),
    # i7,i5 not specified but present
    Fq_data("i7(0)nt, i5(0)nt, i1(12)nt",
            FQ_RES / "i7_i5_i1_read_1.fastq.gz", FQ_RES / "i7_i5_i1_read_2.fastq.gz",
            False, False, True,
            None, None, I1_START, I1_START + INDEX_LENGTH),
    # i7 not specified but present
    Fq_data("i7(0)nt, i5(12)nt, i1(12)nt",
            FQ_RES / "i7_i5_i1_read_1.fastq.gz", FQ_RES / "i7_i5_i1_read_2.fastq.gz",
            False, True, True,
            None, INDEX_LENGTH, I1_START, I1_START + INDEX_LENGTH),
    # i5 not specified but present
    Fq_data("i7(12)nt, i5(0)nt, i1(12)nt",
            FQ_RES / "i7_i5_i1_read_1.fastq.gz", FQ_RES / "i7_i5_i1_read_2.fastq.gz",
            True, False, True,
            INDEX_LENGTH, None, I1_START, I1_START + INDEX_LENGTH),
    # short i7 but full i5, i1 present
    Fq_data("i7(10)nt, i5(12nt), i1(12)nt",
            FQ_RES / "i7-2_i5_i1_read_1.fastq.gz", FQ_RES / "i7-2_i5_i1_read_2.fastq.gz",
            True, True, True,
            INDEX_LENGTH - 2, INDEX_LENGTH, I1_START, I1_START + INDEX_LENGTH),
]

FQ_TO_FAIL = [
    # Here we mostly check of length of the by the user specified barcodes have the same
    # length as in the fastq file. We need to to this and otherwise throw an error
    # as we otherwise cant error correct due to ambiguity. All of the conditions below
    # need to fail.
    # i7 in the fastq file is longer than specified.
    Fq_data("too long i7",
            FQ_RES / "i7+2_i5_i1_read_1.fastq.gz", FQ_RES / "i7_i5_i1_read_2.fastq.gz",
            True, True, True,
            INDEX_LENGTH, INDEX_LENGTH, I1_START, I1_START + INDEX_LENGTH),
    # i5 in the fastq file is longer than specified
    Fq_data("too long i5",
            FQ_RES / "i7_i5+2_i1_read_1.fastq.gz", FQ_RES / "i7_i5+2_i1_read_2.fastq.gz",
            True, True, True,
            INDEX_LENGTH, INDEX_LENGTH, I1_START, I1_START + INDEX_LENGTH),
    # i1 index the fastq file shorter than specified. Needs to raise an error
    Fq_data("too short i1 sequence",
            FQ_RES / "i7_i5_noi1_read_1.fastq.gz", FQ_RES / "i7_i5_noi1_read_2.fastq.gz",
            True, True, True,
            INDEX_LENGTH, INDEX_LENGTH, I1_START, I1_START + INDEX_LENGTH),
    # mate1 and mate 2 have different barcode length. Most likely the user provided
    # files from to different runs. Needs to raise an error
    Fq_data("different barcode headers",
            FQ_RES / "i7+2_i5_i1_read_1.fastq.gz", FQ_RES / "i7_i5_i1_read_2.fastq.gz",
            True, True, True,
            INDEX_LENGTH, INDEX_LENGTH, I1_START, I1_START + INDEX_LENGTH),
    # i7 barcode in the fq file is shorter than specified. Needs to raise error.
    Fq_data("too short i7",
            FQ_RES / "i7-2_i5_i1_read_1.fastq.gz", FQ_RES / "i7-2_i5_i1_read_2.fastq.gz",
            True, True, True,
            INDEX_LENGTH, INDEX_LENGTH, I1_START, I1_START + INDEX_LENGTH)
]


@pytest.mark.parametrize('csv_file',
                         [pytest.param(f, id=f.name) for f in CSVS_TO_FAIL])
def test_parse_sample_sheet_to_fail(csv_file):
    with pytest.raises(SystemExit) as e:
        parse_sample_sheet(csv_file, False)


@pytest.mark.parametrize('csv_file', [pytest.param(f, id=f.name) for f in CSVS_TO_PASS])
def test_parse_sample_sheet_to_pass(csv_file):
    parse_sample_sheet(csv_file, False)


@pytest.mark.parametrize('bc_to_load',
                         [pytest.param(f, id=f.__repr__()) for f in BARCODES_TO_TEST])
def test_load_correction_maps(bc_to_load):
    bc = load_correction_map(bc_to_load)
    corr_map = bc.correction_map
    assert list(corr_map.keys()) == list(corr_map.values())


@pytest.mark.parametrize("test_in, expected", FOR_RC_PASS)
def test_reverse_complement_pass(test_in, expected):
    res = reverse_complement(test_in)
    assert res == expected


@pytest.mark.parametrize('sequences_fail', FOR_RC_FAIL)
def test_reverse_complement_fail(sequences_fail):
    with pytest.raises((ValueError, TypeError)) as e:
        reverse_complement(sequences_fail)


@pytest.mark.parametrize('fq_pass',
                         [pytest.param(f._asdict(), id=f.condition) for f in FQ_TO_PASS])
def test_peek_into_fastq_files_fq_pass(fq_pass):
    peek_into_fastq_files(**fq_pass)


@pytest.mark.parametrize('fq_fail',
                         [pytest.param(f._asdict(), id=f.condition) for f in FQ_TO_FAIL])
def test_peek_into_fastq_files_fq_fail(fq_fail):
    with pytest.raises(SystemExit) as e:
        peek_into_fastq_files(**fq_fail)


def test_fq_gz_parser(one_paired_read):
    test_r_1, test_r_2 = one_paired_read
    fq_1 = FQ_RES / "i7_i5_i1_read_1.fastq.gz"
    fq_2 = FQ_RES / "i7_i5_i1_read_2.fastq.gz"
    with get_pe_fastq(fq_1, fq_2) as pe_reads:
        for reads in pe_reads:
            mate1, mate2 = reads
            assert mate1 == test_r_1
            assert mate2 == test_r_2

@pytest.fixture
def one_paired_read():
    read1 = ("@read_name:AANACATGCGTT+CGCCACTGAGTT\n",
              "AAAACATGCGTTCCCCACTGAGTTAAAACATGCGTT\n",
              "+\n",
              "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE\n")
    read2 = ("@read_name:AANACATGCGTT+CGCCACTGAGTT\n",
              "TTTTAGCATGAAAACATGCGTT\n",
              "+\n",
              "EEEEEEEEEEEEEEEEEEEEEE\n")
    return read1, read2
