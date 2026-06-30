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
    base_peak_mass: float


@dataclass
class FakeCentroidStream:
    """Stand-in for fisher_py's centroid stream: parallel m/z / intensity / charge.

    ``charges`` may hold 0 for ions whose charge the instrument left unassigned,
    mirroring real MS1 data; ``_base_peak_charge`` must read that as "unknown".
    """

    masses: list
    intensities: list
    charges: list


@dataclass
class FakeRunHeader:
    first_spectrum: int
    last_spectrum: int


class FakeReader:
    """Minimal reader implementing the surface extract_chromatograms() uses.

    ``scans`` is a list of (ms_level, rt, tic, bpc, bpm[, bp_charge]) tuples,
    1-indexed by scan number. ``bp_charge`` is optional (defaults to 0 =
    unassigned) and only consulted via the centroid stream when charge is wanted.
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
        s = self._scan(scan_no)
        return FakeStats(tic=s[2], base_peak_intensity=s[3], base_peak_mass=s[4])

    def get_centroid_stream(self, scan_no, include_reference):
        s = self._scan(scan_no)
        charge = s[5] if len(s) > 5 else 0
        # Single-peak spectrum: the base peak is the only (hence tallest) ion.
        return FakeCentroidStream(
            masses=[s[4]], intensities=[s[3]], charges=[charge])

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
        # (ms_level, rt, tic, bpc, base_peak_mass)
        (1, 0.10, 1000.0, 500.0, 400.1234),
        (2, 0.11, 50.0, 25.0, 410.5),
        (1, 0.20, 2000.0, 900.0, 524.2671),
        (2, 0.21, 60.0, 30.0, 300.0),
        (1, 0.30, 1500.0, 700.0, 612.3456),
    ]


@pytest.fixture
def charged_scans():
    """MS1 scans carrying an assigned base-peak charge (and one left unassigned)."""
    return [
        # (ms_level, rt, tic, bpc, base_peak_mass, bp_charge)
        (1, 0.10, 1000.0, 500.0, 400.1234, 2),
        (1, 0.20, 2000.0, 900.0, 524.2671, 3),
        (1, 0.30, 1500.0, 700.0, 612.3456, 0),  # 0 => unassigned -> nan
    ]
