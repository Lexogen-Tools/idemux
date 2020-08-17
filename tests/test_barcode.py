#!/usr/bin/env python
"""Tests for `idemux` package."""
from collections import defaultdict

import pytest

from idemux.ioutils.barcode import Barcode


@pytest.fixture
def invalid_lengths():
    vals = defaultdict(list)
    vals[("A" * 7)].append("test0")
    vals["C" * 7].append("test1")
    return vals


@pytest.fixture
def none_values():
    vals = defaultdict(list)
    vals[None].append("test0")
    return vals


@pytest.fixture
def mixed_none_values():
    vals = defaultdict(list)
    vals[None].append("test0")
    vals["C" * 10].append("test1")
    return vals


@pytest.fixture
def different_lengths():
    vals = defaultdict(list)
    vals[("A" * 10)].append("test0")
    vals["C" * 12].append("test1")
    return vals


@pytest.fixture
def correct_length():
    vals = defaultdict(list)
    vals[("A" * 12)].append("test0")
    return vals


def test_invalid_lengths(invalid_lengths):
    with pytest.raises(ValueError):
        Barcode("invalid length", invalid_lengths),


def test_different_lengths(different_lengths):
    with pytest.raises(ValueError):
        Barcode("different lengths", different_lengths)


def test_barcode_rc(correct_length):
    bc = Barcode("i5", correct_length, True)
    assert bc.reverse_complement
    assert bc.name == "i5_rc"


def test_is_full(correct_length):
    bc = Barcode("i5", correct_length)
    assert bc.full


def test_is_sparse(mixed_none_values):
    bc = Barcode("i5", mixed_none_values)
    assert bc.sparse


def test_is_empty(none_values):
    bc = Barcode("i5", none_values)
    assert bc.empty
    assert bc.length == 0


def test_length_getter(correct_length):
    bc = Barcode("i5", correct_length, True)
    assert bc.length == 12
