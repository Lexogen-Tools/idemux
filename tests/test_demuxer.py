#!/usr/bin/env python
"""Tests for `idemux` package."""
import csv
import pathlib
import pytest
from os import path
from idemux.ioutils.parser import parse_sample_sheet, get_pe_fastq
from idemux.processing.demuxer import process_mate_pair, demux_paired_end

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
    # demulitplexing only on i1 one.
    # should not matter if the header contains i7 and i5 barcodes
    res = pathlib.Path(path.dirname(__file__)) / "resources" / "end_to_end"
    read_1 = res / "i7_i5_i1_read_1.fastq.gz"
    read_2 = res / "i7_i5_i1_read_2.fastq.gz"
    csv = res / "i1_sample_sheet.csv"
    return read_1, read_2, csv


@pytest.fixture
def demux_i7_i5():
    res = pathlib.Path(path.dirname(__file__)) / "resources" / "end_to_end"
    read_1 = res / "i7_i5_read_1.fastq.gz"
    read_2 = res / "i7_i5_read_2.fastq.gz"
    csv = res / "i7_i5_sample_sheet.csv"
    return read_1, read_2, csv


def demux_loop(read1, read2, barcodes, i1_start, i1_end):
    i7, i5, i1 = barcodes
    with get_pe_fastq(read1, read2) as pe_reads:
        for mate_pair in pe_reads:
            # here we do the error correction and get obtain the i1 barcode if present
            corrected_bc, processed_mates = process_mate_pair(mate_pair,
                                                              i7, i5, i1,
                                                              i1_start,
                                                              i1_end)
            yield corrected_bc, processed_mates


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

    _, _, i1 = barcodes

    i1_end = I1_START + i1.length

    for corrected_bc, processed_mates in demux_loop(read1, read2,
                                                    barcodes,
                                                    I1_START, i1_end):
        i7_corr, i5_corr, i1_corr = corrected_bc
        m1, m2 = processed_mates

        m1_split = m1.split("\n")
        m2_split = m2.split("\n")

        m1_seq = m1_split[1]
        m1_qs = m1_split[3]

        m2_seq = m2_split[1]
        m2_qs = m2_split[3]

        expected_i7 = m1_seq[I7_POS]
        expected_i5 = m1_seq[I5_POS]
        expected_i1 = m1_seq[I1_POS]

        assert i7_corr == expected_i7
        assert i5_corr == expected_i5
        assert i1_corr == expected_i1

        assert len(m1_seq) == len(m1_qs)
        assert len(m2_seq) == len(m2_qs)


def test_demux_i7_i1(demux_i7_i1):
    read1, read2, csv = demux_i7_i1
    barcode_sample_map, barcodes = parse_sample_sheet(csv, i5_rc=False)

    _, _, i1 = barcodes

    i1_end = I1_START + i1.length

    for corrected_bc, processed_mates in demux_loop(read1, read2,
                                                    barcodes,
                                                    I1_START, i1_end):
        i7_corr, i5_corr, i1_corr = corrected_bc

        m1, m2 = processed_mates
        m1_split = m1.split("\n")
        m2_split = m2.split("\n")

        m1_seq = m1_split[1]
        m1_qs = m1_split[3]

        m2_seq = m2_split[1]
        m2_qs = m2_split[3]

        expected_i7 = m1_seq[I7_POS]
        expected_i5 = None
        expected_i1 = m1_seq[I1_POS]

        assert i7_corr == expected_i7
        assert i5_corr == expected_i5
        assert i1_corr == expected_i1

        assert len(m1_seq) == len(m1_qs)
        assert len(m2_seq) == len(m2_qs)


def test_demux_i7_i5(demux_i7_i5):
    read1, read2, csv = demux_i7_i5
    barcode_sample_map, barcodes = parse_sample_sheet(csv, i5_rc=False)

    _, _, i1 = barcodes

    i1_end = I1_START + i1.length

    for corrected_bc, processed_mates in demux_loop(read1, read2,
                                                    barcodes,
                                                    I1_START, i1_end):
        i7_corr, i5_corr, i1_corr = corrected_bc

        m1, m2 = processed_mates
        m1_split = m1.split("\n")
        m2_split = m2.split("\n")

        m1_seq = m1_split[1]
        m1_qs = m1_split[3]

        m2_seq = m2_split[1]
        m2_qs = m2_split[3]

        expected_i7 = m1_seq[I7_POS]
        expected_i5 = m1_seq[I5_POS]
        expected_i1 = None

        assert i7_corr == expected_i7
        assert i5_corr == expected_i5
        assert i1_corr == expected_i1

        assert len(m1_seq) == len(m1_qs)
        assert len(m2_seq) == len(m2_qs)


def test_demux_i1(demux_i1):
    read1, read2, csv = demux_i1
    barcode_sample_map, barcodes = parse_sample_sheet(csv, i5_rc=False)

    _, _, i1 = barcodes
    i1_end = I1_START + i1.length

    for corrected_bc, processed_mates in demux_loop(read1, read2,
                                                    barcodes,
                                                    I1_START, i1_end):
        i7_corr, i5_corr, i1_corr = corrected_bc

        m1, m2 = processed_mates
        m1_split = m1.split("\n")
        m2_split = m2.split("\n")

        m1_seq = m1_split[1]
        m1_qs = m1_split[3]

        m2_seq = m2_split[1]
        m2_qs = m2_split[3]

        expected_i7 = None
        expected_i5 = None
        expected_i1 = m1_seq[I1_POS]

        assert i7_corr == expected_i7
        assert i5_corr == expected_i5
        assert i1_corr == expected_i1

        assert len(m1_seq) == len(m1_qs)
        assert len(m2_seq) == len(m2_qs)


def test_demux_paired_end(demux_i7_i5_i1, tmp_path):
    expected_reads = 100

    read1, read2, csv_file = demux_i7_i5_i1
    barcode_sample_map, barcodes = parse_sample_sheet(csv_file, i5_rc=False)
    demux_paired_end(barcode_sample_map, barcodes, read1, read2, I1_START, tmp_path)

    stats_file = pathlib.Path(tmp_path / "demultipexing_stats.tsv")
    with open(stats_file, 'r') as stats:
        csv.register_dialect('strip', skipinitialspace=True)
        reader = csv.DictReader(stats, delimiter='\t', dialect='strip')
        for row in reader:
            n_reads = int(row.get("written_reads"))
            assert n_reads == expected_reads
