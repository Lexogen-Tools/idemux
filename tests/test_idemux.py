#!/usr/bin/env python
"""Tests for `idemux` package."""
import logging
import pytest
from idemux.ioutils.barcode import Barcode
from idemux.ioutils.parser import parse_sample_sheet


def test_barcode_class():
    with pytest.raises(ValueError):
        Barcode("i7", {7}, {"AAAAACGGCAGG"})
#
#     with pytest.raises(ValueError):
#         Barcode("i7", {None}, {"AAAAACGGCAGG"})
#
#     with pytest.raises(ValueError):
#         Barcode("i7", {10, 12}, {"AAAAACGGCAGG"})
#
#     assert Barcode("i7", {10}, {"AAAAACGGCAGG":}).get_set_sizes() == [96, 384]
#     assert Barcode("i7", {8}, {"AAAAACGGCAGG":}).get_set_sizes() == [96]
#     logging.getLogger().info("Testing barcode data classes - Done")
