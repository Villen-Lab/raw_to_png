"""Shared test fixtures: a fake Thermo reader so tests run without .NET / fisher_py."""

from dataclasses import dataclass

import pytest


@dataclass
class FakeMsOrder:
    """Stand-in for fisher_py's MsOrderType enum: a plain enum exposing .value.

    Crucially NOT an int — this is what makes ``int(ms_order)`` fail in the real
    library, so the production code must use ``.value``.
    """

    value: int


@dataclass
class FakeFilter:
    ms_order: FakeMsOrder


@dataclass
class FakeStats:
    tic: float
    base_peak_intensity: float


@dataclass
class FakeRunHeader:
    first_spectrum: int
    last_spectrum: int


class FakeReader:
    """Minimal reader implementing the surface extract_chromatograms() uses.

    ``scans`` is a list of (ms_level, rt, tic, bpc) tuples, 1-indexed by scan no.
    """

    def __init__(self, scans):
        self._scans = scans
        self.disposed = False

    @property
    def run_header_ex(self):
        return FakeRunHeader(first_spectrum=1, last_spectrum=len(self._scans))

    def _scan(self, scan_no):
        return self._scans[scan_no - 1]

    def get_filter_for_scan_number(self, scan_no):
        return FakeFilter(ms_order=FakeMsOrder(self._scan(scan_no)[0]))

    def get_scan_stats_for_scan_number(self, scan_no):
        _, _, tic, bpc = self._scan(scan_no)
        return FakeStats(tic=tic, base_peak_intensity=bpc)

    def retention_time_from_scan_number(self, scan_no):
        return self._scan(scan_no)[1]

    def select_instrument(self, device, n):
        pass

    def dispose(self):
        self.disposed = True


@pytest.fixture
def mixed_scans():
    """Interleaved MS1/MS2 scans: 3 MS1, 2 MS2."""
    return [
        # (ms_level, rt, tic, bpc)
        (1, 0.10, 1000.0, 500.0),
        (2, 0.11, 50.0, 25.0),
        (1, 0.20, 2000.0, 900.0),
        (2, 0.21, 60.0, 30.0),
        (1, 0.30, 1500.0, 700.0),
    ]
