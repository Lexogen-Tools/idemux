#!/usr/bin/env python
"""Tests for `idemux` package."""

import pytest
from idemux.ioutils.parser import peek_into_fastq_files


@pytest.fixture
def response():
    """Sample pytest fixture.

    See more at: http://doc.pytest.org/en/latest/fixture.html
    """
    # import requests
    # return requests.get('https://github.com/audreyr/cookiecutter-pypackage')


def test_content(response):
    """Sample pytest test function with the pytest fixture as an argument."""
    # from bs4 import BeautifulSoup
    # assert 'GitHub' in BeautifulSoup(response.content).title.string


def test_fastq_headers(response):

    peek_into_fastq_files("test_1.fastq.gz",
                          "test_2.fastq.gz",
                          True,
                          True)
    assert True


def test_fastq_header_check(response):
    pass
