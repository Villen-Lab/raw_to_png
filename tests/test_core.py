"""Tests for chromatogram extraction and PNG rendering (reader fully mocked)."""

import numpy as np

from raw_to_png import core
from raw_to_png.core import extract_chromatograms, render_png, top_peaks

from conftest import FakeReader


def test_extract_ms1_only(mixed_scans):
    rt, tic, bpc, bpm, charge = extract_chromatograms(FakeReader(mixed_scans), ms_level=1)
    # Only the 3 MS1 scans should survive.
    assert np.allclose(rt, [0.10, 0.20, 0.30])
    assert np.allclose(tic, [1000.0, 2000.0, 1500.0])
    assert np.allclose(bpc, [500.0, 900.0, 700.0])
    assert np.allclose(bpm, [400.1234, 524.2671, 612.3456])
    # Charge is opt-in: defaults to all-nan without want_charge.
    assert np.all(np.isnan(charge))


def test_extract_ms2_only(mixed_scans):
    rt, *_ = extract_chromatograms(FakeReader(mixed_scans), ms_level=2)
    assert np.allclose(rt, [0.11, 0.21])


def test_extract_all_scans(mixed_scans):
    """ms_level=None keeps every scan regardless of MS order."""
    rt, *_ = extract_chromatograms(FakeReader(mixed_scans), ms_level=None)
    assert rt.size == len(mixed_scans)


def test_extract_no_matching_scans():
    only_ms1 = [(1, 0.1, 100.0, 50.0, 400.0)]
    rt, tic, bpc, bpm, charge = extract_chromatograms(FakeReader(only_ms1), ms_level=2)
    assert rt.size == tic.size == bpc.size == bpm.size == charge.size == 0


def test_extract_charge_when_requested(charged_scans):
    """want_charge resolves the base-peak charge; an unassigned (0) charge -> nan."""
    *_, charge = extract_chromatograms(FakeReader(charged_scans), ms_level=1,
                                       want_charge=True)
    assert np.allclose(charge[:2], [2.0, 3.0])
    assert np.isnan(charge[2])


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
    # {62} is close to 60 (stack). Pass an explicit (positive) threshold since the
    # default label spacing is now negative — stacking off — by design.
    idx = np.array([10, 12, 60, 62])
    assert core._stagger_levels(rt, idx, x_thresh_frac=0.06) == [0, 1, 0, 1]


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


def _fake_chromatograms(rt):
    """A get_chromatograms stand-in returning the full 5-array bundle."""
    return rt, rt * 2, rt * 3, 400.0 + rt * 100, np.full(rt.size, 2.0)


def test_render_png_writes_file(tmp_path, monkeypatch):
    rt = np.linspace(0, 1, 10)
    monkeypatch.setattr(core, "get_chromatograms",
                        lambda path, ms_level=1, want_charge=False: _fake_chromatograms(rt))
    out = render_png("sample_file.raw", tmp_path, ms_level=1, dpi=72)
    assert out is not None
    assert out.name == "sample_file.png"
    assert out.exists()


def test_render_png_with_mz_and_charge(tmp_path, monkeypatch):
    """show_mz/show_charge must render without error and forward want_charge."""
    rt = np.linspace(0, 1, 10)
    captured = {}

    def fake(path, ms_level=1, want_charge=False):
        captured["want_charge"] = want_charge
        return _fake_chromatograms(rt)

    monkeypatch.setattr(core, "get_chromatograms", fake)
    out = render_png("sample_file.raw", tmp_path, ms_level=1, dpi=72,
                     show_mz=True, show_charge=True)
    assert out is not None and out.exists()
    assert captured["want_charge"] is True  # show_charge drives the extra read


def test_render_png_skips_when_empty(tmp_path, monkeypatch):
    empty = np.asarray([])
    monkeypatch.setattr(
        core, "get_chromatograms",
        lambda path, ms_level=1, want_charge=False: (empty, empty, empty, empty, empty))
    out = render_png("empty.raw", tmp_path, ms_level=1)
    assert out is None
    assert list(tmp_path.iterdir()) == []


def test_peak_label_intensity_only():
    """show_rt=False isolates the scientific-notation intensity, single line."""
    label = core._peak_label(3.8e10, show_rt=False)
    assert "\n" not in label
    assert "times10" in label  # mathtext scientific notation


def test_peak_label_default_includes_rt():
    """Retention time labels alongside intensity by default (2 lines, bare value)."""
    label = core._peak_label(3.8e10, rt=12.3)
    lines = label.split("\n")
    assert len(lines) == 2
    assert "times10" in lines[0]
    assert lines[1] == "12.30"  # bare minutes, no "RT"/"min" text


def test_peak_label_rt_only():
    label = core._peak_label(3.8e10, rt=12.3, show_intensity=False)
    assert "\n" not in label
    assert label == "12.30"


def test_peak_label_rt_unassigned_renders_question_mark():
    """Missing rt (None, the default) renders as '?', consistent with mz/charge."""
    label = core._peak_label(3.8e10, show_intensity=False)
    assert label == "?"


def test_peak_label_rt_decimals():
    """rt_decimals controls the retention-time precision; default stays at 2."""
    assert core._peak_label(3.8e10, rt=12.3456, show_intensity=False) == "12.35"
    assert core._peak_label(
        3.8e10, rt=12.3456, show_intensity=False, rt_decimals=0) == "12"


def test_peak_label_lines_rt_font_matches_mz():
    """The retention-time line uses the smaller m/z font size, not the intensity one."""
    (_, intensity_fs), (_, rt_fs) = core._peak_label_lines(3.8e10, rt=12.3)
    assert rt_fs == core.MZ_LABEL_FONTSIZE
    assert intensity_fs == core.PEAK_LABEL_FONTSIZE
    assert rt_fs < intensity_fs


def test_peak_label_with_mz_and_charge():
    label = core._peak_label(3.8e10, mz=524.2671, charge=2.0,
                             show_mz=True, show_charge=True, show_rt=False)
    lines = label.split("\n")
    assert len(lines) == 3
    assert "524.2671" in lines[1]
    assert lines[2] == "z=2"


def test_peak_label_lines_mz_font_is_smaller():
    """The m/z line carries a smaller font than the intensity / charge lines."""
    (_, intensity_fs), (_, mz_fs), (_, charge_fs) = core._peak_label_lines(
        3.8e10, mz=524.2671, charge=2.0, show_mz=True, show_charge=True, show_rt=False)
    assert mz_fs == core.MZ_LABEL_FONTSIZE
    assert intensity_fs == charge_fs == core.PEAK_LABEL_FONTSIZE
    assert mz_fs < intensity_fs


def test_peak_label_mz_decimals():
    """mz_decimals controls the m/z precision; default stays at 4."""
    assert "524.2671" in core._peak_label(3.8e10, mz=524.26713, show_mz=True)
    assert "524.27" in core._peak_label(
        3.8e10, mz=524.26713, show_mz=True, mz_decimals=2)
    # No more decimals than requested (2 dp -> not the 4-dp rendering).
    assert "524.2671" not in core._peak_label(
        3.8e10, mz=524.26713, show_mz=True, mz_decimals=2)


def test_peak_label_mz_only_without_intensity():
    """show_intensity=False and show_rt=False drop everything but m/z."""
    label = core._peak_label(3.8e10, mz=524.2671, show_intensity=False, show_mz=True,
                             show_rt=False)
    assert "\n" not in label
    assert "times10" not in label  # no intensity line
    assert "524.2671" in label


def test_peak_label_nothing_selected_is_empty():
    """All fields off -> empty string, so the caller draws just the dot."""
    assert core._peak_label(3.8e10, show_intensity=False, show_rt=False) == ""


def test_peak_label_unassigned_renders_question_mark():
    """Opted-in but missing values stay as aligned '?' rows, not dropped lines."""
    label = core._peak_label(3.8e10, mz=np.nan, charge=np.nan,
                             show_mz=True, show_charge=True, show_rt=False)
    lines = label.split("\n")
    assert len(lines) == 3
    assert lines[1].endswith("?")
    assert lines[2] == "z=?"
