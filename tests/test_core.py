"""Tests for chromatogram extraction and PNG rendering (reader fully mocked)."""

import numpy as np

from raw_to_png import core
from raw_to_png.core import extract_chromatograms, render_png, top_peaks

from conftest import FakeReader


def test_extract_ms1_only(mixed_scans):
    rt, tic, bpc = extract_chromatograms(FakeReader(mixed_scans), ms_level=1)
    # Only the 3 MS1 scans should survive.
    assert np.allclose(rt, [0.10, 0.20, 0.30])
    assert np.allclose(tic, [1000.0, 2000.0, 1500.0])
    assert np.allclose(bpc, [500.0, 900.0, 700.0])


def test_extract_ms2_only(mixed_scans):
    rt, _, _ = extract_chromatograms(FakeReader(mixed_scans), ms_level=2)
    assert np.allclose(rt, [0.11, 0.21])


def test_extract_all_scans(mixed_scans):
    """ms_level=None keeps every scan regardless of MS order."""
    rt, _, _ = extract_chromatograms(FakeReader(mixed_scans), ms_level=None)
    assert rt.size == len(mixed_scans)


def test_extract_no_matching_scans():
    only_ms1 = [(1, 0.1, 100.0, 50.0)]
    rt, tic, bpc = extract_chromatograms(FakeReader(only_ms1), ms_level=2)
    assert rt.size == tic.size == bpc.size == 0


def test_uses_enum_value_not_int(mixed_scans):
    """Guard against regressing to int(ms_order), which raises TypeError on the
    real (non-int) MsOrderType enum that FakeMsOrder mimics."""
    # FakeMsOrder is not an int; int() on it would raise, so a passing run proves
    # the production code compares .value.
    extract_chromatograms(FakeReader(mixed_scans), ms_level=1)


def test_top_peaks_picks_highest_separated():
    rt = np.arange(0, 10, dtype=float)
    y = np.array([1, 9, 2, 3, 8, 1, 7, 2, 6, 1], dtype=float)
    # min_sep_frac=0 means no separation constraint -> pure top-n by height.
    idx = top_peaks(rt, y, n=3, min_sep_frac=0.0)
    assert list(idx) == sorted([1, 4, 6])  # the 9, 8, 7
    # Returned indices are sorted by retention time.
    assert list(idx) == sorted(idx)


def test_top_peaks_enforces_separation():
    rt = np.arange(0, 10, dtype=float)  # span 9
    # Two tall adjacent scans (idx 4,5) belong to one peak; separation should
    # keep only the taller of the pair, then move to the next distinct peak.
    y = np.array([1, 1, 1, 1, 10, 9, 1, 1, 8, 1], dtype=float)
    idx = top_peaks(rt, y, n=2, min_sep_frac=0.3)  # min_sep = 2.7 min
    assert list(idx) == [4, 8]  # not 4 and 5


def test_top_peaks_fewer_than_n():
    rt = np.array([0.0, 1.0])
    y = np.array([5.0, 3.0])
    assert list(top_peaks(rt, y, n=5)) == [0, 1]


def test_top_peaks_empty():
    empty = np.asarray([])
    assert top_peaks(empty, empty, n=5).size == 0


def test_stagger_levels_stacks_close_labels():
    rt = np.arange(0, 100, dtype=float)  # span 99, x_thresh = 0.06*99 ~ 5.9
    # Two clusters: {10,12} are close (stack), then 60 is isolated (reset), and
    # {62} is close to 60 (stack).
    idx = np.array([10, 12, 60, 62])
    assert core._stagger_levels(rt, idx) == [0, 1, 0, 1]


def test_stagger_levels_all_isolated():
    rt = np.arange(0, 100, dtype=float)
    idx = np.array([5, 40, 80])
    assert core._stagger_levels(rt, idx) == [0, 0, 0]


def test_stagger_levels_empty():
    rt = np.arange(0, 10, dtype=float)
    assert core._stagger_levels(rt, np.array([], dtype=int)) == []


def test_stagger_levels_respects_threshold():
    rt = np.arange(0, 100, dtype=float)  # span 99
    idx = np.array([10, 12, 14])  # 2 apart
    # Tiny threshold (0.01*99 < 2): none count as close -> all baseline row.
    assert core._stagger_levels(rt, idx, x_thresh_frac=0.01) == [0, 0, 0]
    # Wide threshold (0.5*99): each is close to the prior -> stack up.
    assert core._stagger_levels(rt, idx, x_thresh_frac=0.5) == [0, 1, 2]


def test_render_png_writes_file(tmp_path, monkeypatch):
    rt = np.linspace(0, 1, 10)
    monkeypatch.setattr(core, "get_chromatograms",
                        lambda path, ms_level=1: (rt, rt * 2, rt * 3))
    out = render_png("sample_file.raw", tmp_path, ms_level=1, dpi=72)
    assert out is not None
    assert out.name == "sample_file.png"
    assert out.exists()


def test_render_png_skips_when_empty(tmp_path, monkeypatch):
    empty = np.asarray([])
    monkeypatch.setattr(core, "get_chromatograms",
                        lambda path, ms_level=1: (empty, empty, empty))
    out = render_png("empty.raw", tmp_path, ms_level=1)
    assert out is None
    assert list(tmp_path.iterdir()) == []
