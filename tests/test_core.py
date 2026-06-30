"""Tests for chromatogram extraction and PNG rendering (reader fully mocked)."""

import numpy as np

from raw_to_png import core
from raw_to_png.core import extract_chromatograms, render_png

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
