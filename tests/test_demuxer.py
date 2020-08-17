#!/usr/bin/env python
"""Tests for `idemux` package."""
from idemux.ioutils.parser import get_pe_fastq, parse_sample_sheet
from idemux.processing.demuxer import process_mate_pair

import pathlib
from os import path

import pytest

from idemux.ioutils.parser import parse_sample_sheet, get_pe_fastq
from idemux.processing.demuxer import process_mate_pair

I1_START = 10

I7_POS = slice(0, 12)
I5_POS = slice(12, 24)
I1_POS = slice(24, 36)

@pytest.fixture
def demux_i7_i5_i1():
    res = pathlib.Path(path.dirname(__file__)) / "resources" / "end_to_end"
    read_1 = res / "i7_i5_i1_read_1.fastq.gz"
    read_2 = res / "i7_i5_i1_read_2.fastq.gz"
    csv = res / "i7_i5_i1_sample_sheet.csv"
    return read_1, read_2, csv


@pytest.fixture
def demux_i7_i1():
    res = pathlib.Path(path.dirname(__file__)) / "resources" / "end_to_end"
    read_1 = res / "i7_i1_read_1.fastq.gz"
    read_2 = res / "i7_i1_read_2.fastq.gz"
    csv = res / "i7_i1_sample_sheet.csv"
    return read_1, read_2, csv


@pytest.fixture
def demux_i5_i1():
    res = pathlib.Path(path.dirname(__file__)) / "resources" / "end_to_end"
    read_1 = res / "i7_i1_read_1.fastq.gz"
    read_2 = res / "i7_i1_read_2.fastq.gz"
    csv = res / "i7_i1_sample_sheet.csv"
    return read_1, read_2, csv


@pytest.fixture
def demux_i1():
    # demulitplexing only on i1 one. doesnt matter of the file
    # should not matter if the header contains i7 and i5 barcodes
    res = pathlib.Path(path.dirname(__file__)) / "resources" / "end_to_end"
    read_1 = res / "i7_i5_i1_read_1.fastq.gz"
    read_2 = res / "i7_i5_i1_read_1.fastq.gz"
    csv = res / "i1_sample_sheet.csv"
    return read_1, read_2, csv


@pytest.fixture
def demux_i7_i5():
    res = pathlib.Path(path.dirname(__file__)) / "resources" / "end_to_end"
    read_1 = res / "i7_i5_read_1.fastq.gz"
    read_2 = res / "i7_i5_read_2.fastq.gz"
    csv = res / "i7_i1_sample_sheet.csv"
    return read_1, read_2, csv


def test_demux_i7_i5_i1(demux_i7_i5_i1):
    """Testing expected vs calcualted barcodes from paired-end fastq files

    * reads in read1 fastq file look as follows:
    <read_name>:<correctable_i7>+<correctable_i5>
    <correct_i7><correct_i5><correct_i1>
    +
    E{36}

    * reads in read2 fastq file look as follows:
    <read_name>:<correctable_i7>+<correctable_i5>
    <rand_10nt><correctable_i1>
    +
    E{22}

    """

    read1, read2, csv = demux_i7_i5_i1
    barcode_sample_map, barcodes = parse_sample_sheet(csv, i5_rc=False)

    i7, i5, i1 = barcodes
    i1_end = I1_START + i1.length

    # then we iterate over all the paired end reads
    with get_pe_fastq(read1, read2) as pe_reads:
        for mate_pair in pe_reads:
            # here we do the error correction and get obtain the i1 barcode if present
            corrected_bc, processed_mates = process_mate_pair(mate_pair,
                                                              i7.used_codes,
                                                              i5.used_codes,
                                                              i1.used_codes,
                                                              i7.not_empty,
                                                              i7.correction_map,
                                                              i5.correction_map,
                                                              i1.correction_map,
                                                              I1_START,
                                                              i1_end)

            i7_corr, i5_corr, i1_corr = corrected_bc
            m1, m2 = processed_mates

            m1_seq = m1.split("\n")[1]
            expected_i7 = m1_seq[I7_POS]
            expected_i5 = m1_seq[I5_POS]
            expected_i1 = m1_seq[I1_POS]

            assert i7_corr == expected_i7
            assert i5_corr == expected_i5
            assert i1_corr == expected_i1


def test_demux_i7_i1(demux_i7_i1):
    read1, read2, csv = demux_i7_i1
    barcode_sample_map, barcodes = parse_sample_sheet(csv, i5_rc=False)

    i7, i5, i1 = barcodes
    i1_start = 10
    i1_end = I1_START + i1.length

    # then we iterate over all the paired end reads
    with get_pe_fastq(read1, read2) as pe_reads:
        for mate_pair in pe_reads:
            # here we do the error correction and get obtain the i1 barcode if present
            corrected_bc, processed_mates = process_mate_pair(mate_pair,
                                                              i7.used_codes,
                                                              i5.used_codes,
                                                              i1.used_codes,
                                                              i7.not_empty,
                                                              i7.correction_map,
                                                              i5.correction_map,
                                                              i1.correction_map,
                                                              I1_START,
                                                              i1_end)

            i7_corr, i5_corr, i1_corr = corrected_bc
            m1, m2 = processed_mates

            m1_seq = m1.split("\n")[1]
            expected_i7 = m1_seq[I7_POS]
            expected_i5 = None
            expected_i1 = m1_seq[I1_POS]

            assert i7_corr == expected_i7
            assert i5_corr == expected_i5
            assert i1_corr == expected_i1


def test_demux_i1(demux_i1):
    read1, read2, csv = demux_i1
    barcode_sample_map, barcodes = parse_sample_sheet(csv, i5_rc=False)

    i7, i5, i1 = barcodes
    i1_end = I1_START + i1.length

    # then we iterate over all the paired end reads
    with get_pe_fastq(read1, read2) as pe_reads:
        for mate_pair in pe_reads:
            # here we do the error correction and get obtain the i1 barcode if present
            corrected_bc, processed_mates = process_mate_pair(mate_pair,
                                                              i7.used_codes,
                                                              i5.used_codes,
                                                              i1.used_codes,
                                                              i7.not_empty,
                                                              i7.correction_map,
                                                              i5.correction_map,
                                                              i1.correction_map,
                                                              I1_START,
                                                              i1_end)

            i7_corr, i5_corr, i1_corr = corrected_bc
            m1, m2 = processed_mates

            m1_seq = m1.split("\n")[1]
            expected_i7 = None
            expected_i5 = None
            expected_i1 = m1_seq[I1_POS]

            assert i7_corr == expected_i7
            assert i5_corr == expected_i5
            assert i1_corr == expected_i1
