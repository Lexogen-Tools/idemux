#!/usr/bin/env python
"""Tests for `idemux` package."""
from idemux.ioutils.parser import get_pe_fastq
from idemux.processing.demuxer import process_mate_pair

#
# def test_process_mate_pair(read1, read2):
#     with get_pe_fastq(read1, read2) as pe_reads:
#         for mate_pair in (pe_reads):
#                 barcodes, processed_mates = process_mate_pair(mate_pair,
#                                                                i7_wanted,
#                                                                i5_wanted,
#                                                                i1_wanted,
#                                                                has_i7,
#                                                                map_i7,
#                                                                map_i5,
#                                                                map_i1,
#                                                                i1_start,
#                                                                i1_end)
#
#
#
#                 i7, i5, i1 = barcodes
#                i7_expected = processed_mates[0][:12]
#                i5_expected = processed_mates[0][12:24]
#                i1_expected = processed_mates[24:36]
#
#                 assert i7==i7_expected
#                 assert i5==i5_expected
#                 assert i1==i1_expected