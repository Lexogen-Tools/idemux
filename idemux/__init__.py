"""Top-level package for idemux."""

__author__ = """Falko Hofmann"""
__email__ = 'falko.hofmann@lexogen.com'

try:
    from idemux._version import version
    __version__ = version
except ImportError:
    __version__ = "not-installed"
