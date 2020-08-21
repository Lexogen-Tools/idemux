
"""Barcode data class - stores data during processing"""

from collections import defaultdict
from dataclasses import dataclass, field
import logging

log = logging.getLogger(__name__)


@dataclass
class Barcode:
    name: str
    _sample_map: defaultdict = field(default_factory=lambda: defaultdict(list))
    reverse_complement: bool = False
    correction_map: dict = None
    length: int = 0
    empty = True
    not_empty = False
    used_codes: set = field(default_factory=set)
    _observed_lengths: set = field(default_factory=set)
    _lengths_96_barcodes = {8, 10, 12}
    _lengths_384_barcodes = {10, 12}
    _allowed_lengths = {6, 8, 10, 12}


    def __post_init__(self):
        # if the barcode is reverse complement add _rc to the name
        # we need this later to load the correction maps
        if self.reverse_complement:
            self.name = f"{self.name}_rc"
        self._check_length()
        self.used_codes = set(self._sample_map.keys())
        self.empty = None in self.used_codes and len(self.used_codes) == 1
        self.not_empty = not self.empty

    def _check_length(self):
        # barcodes for each type need to be of the same length. This is as the
        # sequences of 8 nt long barcodes are contained within 10 nt long barcodes,
        # and the ones of 10 in 12. if we dont enforce the same length per barcode there
        # might be a possibility we cant tell with certainty to which sample a
        # barcoded read belongs.
        for barcode, sample_name in self._sample_map.items():
            if barcode is not None:
                if len(barcode) in self._allowed_lengths:
                    self._observed_lengths.add(len(barcode))
                    self.length = len(barcode)
                    if len(self._observed_lengths) > 1:
                        raise ValueError(f"{self.name} barcodes with a different length "
                                         f"have been observed for {sample_name}. Barcodes"
                                         f" need to have the same length for all "
                                         f"samples.\nObserved barcode length:"
                                         f" {len(barcode)} \nPrevious observed length:"
                                         f" {self._observed_lengths.pop()}")
                else:
                    raise ValueError(f"{self.name} barcodes are {len(barcode)} nt long."
                                     f"{self.name} barcodes are only allowed to be"
                                     f" {self._allowed_lengths} nt long.")

    def get_set_sizes(self):
        set_sizes = []
        if self.length in self._allowed_lengths:
            if self.length in self._lengths_96_barcodes:
                set_sizes.append(96)
            if self.length in self._lengths_384_barcodes:
                set_sizes.append(384)
        return set_sizes

    def get_used_codes(self, drop_none=False):
        _return_val = self.used_codes.copy()
        if drop_none:
            _return_val.discard(None)
        return _return_val

    # attribute getters
    @property
    def samples_without_barcodes(self):
        return self._sample_map.get(None)

    @property
    def sparse(self):
        return None in self.used_codes and len(self.used_codes) > 1

    @property
    def full(self):
        return None not in self.used_codes
