"""Top-level package for idemux."""

__author__ = """Falko Hofmann"""
__email__ = 'falko.hofmann@lexogen.com'

try:
    from idemux._version import __version__
except ImportError:
    __version__ = "not-installed"
