"""Tests for file discovery and CLI orchestration (rendering mocked out)."""

from raw_to_png import cli
from raw_to_png.cli import build_parser, find_raw_files, main


def test_find_raw_files_sorted_unique(tmp_path):
    for n in ("b.raw", "a.raw", "c.RAW"):
        (tmp_path / n).touch()
    found = find_raw_files(tmp_path)
    names = [p.name for p in found]
    # No duplicates even though *.raw and *.RAW are globbed separately,
    # and the result is sorted.
    assert names == sorted(set(names))
    assert len(names) == len(set(names))
    assert "c.RAW" in names


def test_find_raw_files_empty(tmp_path):
    assert find_raw_files(tmp_path) == []


def test_main_empty_dir_returns_1(tmp_path, capsys):
    out_dir = tmp_path / "out"
    rc = main(["--in", str(tmp_path), "--out", str(out_dir)])
    assert rc == 1
    assert "No .RAW files found" in capsys.readouterr().out


def test_main_ms_level_0_maps_to_none(tmp_path, monkeypatch):
    """--ms-level 0 is the 'all scans' sentinel and must reach render as None."""
    (tmp_path / "x.raw").touch()
    captured = {}

    def fake_render(raw, out_dir, ms_level=1, dpi=200, n_peak_labels=5, label_spacing=0.05,
                    show_mz=False, show_charge=False, show_intensity=True, show_rt=True,
                    label_rotation=-60, mz_decimals=4, rt_decimals=2,
                    max_label_corner="left", tic_ymax=None, bpc_ymax=None):
        captured["ms_level"] = ms_level
        return None

    monkeypatch.setattr(cli, "render_png", fake_render)
    rc = main(["--in", str(tmp_path), "--out", str(tmp_path / "out"), "--ms-level", "0"])
    assert rc == 0
    assert captured["ms_level"] is None


def test_main_one_bad_raw_does_not_abort_batch(tmp_path, monkeypatch):
    """Fail-soft: an exception on one file is caught; the batch still completes."""
    for n in ("good1.raw", "bad.raw", "good2.raw"):
        (tmp_path / n).touch()
    seen = []

    def fake_render(raw, out_dir, ms_level=1, dpi=200, n_peak_labels=5, label_spacing=0.05,
                    show_mz=False, show_charge=False, show_intensity=True, show_rt=True,
                    label_rotation=-60, mz_decimals=4, rt_decimals=2,
                    max_label_corner="left", tic_ymax=None, bpc_ymax=None):
        seen.append(raw.name)
        if raw.name == "bad.raw":
            raise RuntimeError("boom")
        return None

    monkeypatch.setattr(cli, "render_png", fake_render)
    rc = main(["--in", str(tmp_path), "--out", str(tmp_path / "out")])
    assert rc == 0
    assert seen == ["bad.raw", "good1.raw", "good2.raw"]  # all attempted, sorted


def test_parser_defaults():
    args = build_parser().parse_args(["--in", "i", "--out", "o"])
    assert args.ms_level == 1
    assert args.dpi == 200
    assert args.top_peaks == 5
    assert args.label_spacing == -0.01  # negative: stagger off, slant separates
    assert args.label_rotation == -60
    assert args.mz_decimals == 4
    assert args.rt_decimals == 2
    # New annotation toggles default off, leaving the existing plot unchanged.
    assert args.show_mz is False
    assert args.show_charge is False
    assert args.hide_intensity is False
    # Retention time is labelled by default (opt-out, not opt-in).
    assert args.hide_rt is False


def test_main_show_mz_and_charge_passed_through(tmp_path, monkeypatch):
    (tmp_path / "x.raw").touch()
    captured = {}

    def fake_render(raw, out_dir, ms_level=1, dpi=200, n_peak_labels=5, label_spacing=0.05,
                    show_mz=False, show_charge=False, show_intensity=True, show_rt=True,
                    label_rotation=-60, mz_decimals=4, rt_decimals=2,
                    max_label_corner="left", tic_ymax=None, bpc_ymax=None):
        captured["show_mz"] = show_mz
        captured["show_charge"] = show_charge
        return None

    monkeypatch.setattr(cli, "render_png", fake_render)
    main(["--in", str(tmp_path), "--out", str(tmp_path / "out"),
          "--show-mz", "--show-charge"])
    assert captured["show_mz"] is True
    assert captured["show_charge"] is True


def test_main_hide_intensity_inverts_to_show_intensity(tmp_path, monkeypatch):
    """--hide-intensity reaches render as show_intensity=False (default True)."""
    (tmp_path / "x.raw").touch()
    captured = {}

    def fake_render(raw, out_dir, ms_level=1, dpi=200, n_peak_labels=5, label_spacing=0.05,
                    show_mz=False, show_charge=False, show_intensity=True, show_rt=True,
                    label_rotation=-60, mz_decimals=4, rt_decimals=2,
                    max_label_corner="left", tic_ymax=None, bpc_ymax=None):
        captured["show_intensity"] = show_intensity
        return None

    monkeypatch.setattr(cli, "render_png", fake_render)
    main(["--in", str(tmp_path), "--out", str(tmp_path / "out"),
          "--show-mz", "--hide-intensity"])
    assert captured["show_intensity"] is False


def test_main_hide_rt_inverts_to_show_rt(tmp_path, monkeypatch):
    """--hide-rt reaches render as show_rt=False (default True)."""
    (tmp_path / "x.raw").touch()
    captured = {}

    def fake_render(raw, out_dir, ms_level=1, dpi=200, n_peak_labels=5, label_spacing=0.05,
                    show_mz=False, show_charge=False, show_intensity=True, show_rt=True,
                    label_rotation=-60, mz_decimals=4, rt_decimals=2,
                    max_label_corner="left", tic_ymax=None, bpc_ymax=None):
        captured["show_rt"] = show_rt
        return None

    monkeypatch.setattr(cli, "render_png", fake_render)
    main(["--in", str(tmp_path), "--out", str(tmp_path / "out"), "--hide-rt"])
    assert captured["show_rt"] is False


def test_main_rt_decimals_passed_through(tmp_path, monkeypatch):
    (tmp_path / "x.raw").touch()
    captured = {}

    def fake_render(raw, out_dir, ms_level=1, dpi=200, n_peak_labels=5, label_spacing=0.05,
                    show_mz=False, show_charge=False, show_intensity=True, show_rt=True,
                    label_rotation=-60, mz_decimals=4, rt_decimals=2,
                    max_label_corner="left", tic_ymax=None, bpc_ymax=None):
        captured["rt_decimals"] = rt_decimals
        return None

    monkeypatch.setattr(cli, "render_png", fake_render)
    main(["--in", str(tmp_path), "--out", str(tmp_path / "out"), "--rt-decimals", "1"])
    assert captured["rt_decimals"] == 1


def test_main_top_peaks_passed_through(tmp_path, monkeypatch):
    (tmp_path / "x.raw").touch()
    captured = {}

    def fake_render(raw, out_dir, ms_level=1, dpi=200, n_peak_labels=5, label_spacing=0.05,
                    show_mz=False, show_charge=False, show_intensity=True, show_rt=True,
                    label_rotation=-60, mz_decimals=4, rt_decimals=2,
                    max_label_corner="left", tic_ymax=None, bpc_ymax=None):
        captured["n_peak_labels"] = n_peak_labels
        return None

    monkeypatch.setattr(cli, "render_png", fake_render)
    main(["--in", str(tmp_path), "--out", str(tmp_path / "out"), "--top-peaks", "3"])
    assert captured["n_peak_labels"] == 3


def test_main_label_rotation_passed_through(tmp_path, monkeypatch):
    (tmp_path / "x.raw").touch()
    captured = {}

    def fake_render(raw, out_dir, ms_level=1, dpi=200, n_peak_labels=5, label_spacing=0.05,
                    show_mz=False, show_charge=False, show_intensity=True, show_rt=True,
                    label_rotation=-60, mz_decimals=4, rt_decimals=2,
                    max_label_corner="left", tic_ymax=None, bpc_ymax=None):
        captured["label_rotation"] = label_rotation
        return None

    monkeypatch.setattr(cli, "render_png", fake_render)
    main(["--in", str(tmp_path), "--out", str(tmp_path / "out"), "--label-rotation", "0"])
    assert captured["label_rotation"] == 0


def test_main_decimals_passed_through(tmp_path, monkeypatch):
    (tmp_path / "x.raw").touch()
    captured = {}

    def fake_render(raw, out_dir, ms_level=1, dpi=200, n_peak_labels=5, label_spacing=0.05,
                    show_mz=False, show_charge=False, show_intensity=True, show_rt=True,
                    label_rotation=-60, mz_decimals=4, rt_decimals=2,
                    max_label_corner="left", tic_ymax=None, bpc_ymax=None):
        captured["mz_decimals"] = mz_decimals
        return None

    monkeypatch.setattr(cli, "render_png", fake_render)
    main(["--in", str(tmp_path), "--out", str(tmp_path / "out"), "--decimals", "2"])
    assert captured["mz_decimals"] == 2


def test_main_max_label_corner_passed_through(tmp_path, monkeypatch):
    (tmp_path / "x.raw").touch()
    captured = {}

    def fake_render(raw, out_dir, **kwargs):
        captured["max_label_corner"] = kwargs["max_label_corner"]

    monkeypatch.setattr(cli, "render_png", fake_render)
    main(["--in", str(tmp_path), "--out", str(tmp_path / "out"), "--max-label-corner", "right"])
    assert captured["max_label_corner"] == "right"


def test_max_label_corner_defaults_left_and_rejects_bad_value():
    import pytest

    assert build_parser().parse_args(["--in", "i", "--out", "o"]).max_label_corner == "left"
    with pytest.raises(SystemExit):
        build_parser().parse_args(["--in", "i", "--out", "o", "--max-label-corner", "middle"])


def test_main_ymax_passed_through(tmp_path, monkeypatch):
    (tmp_path / "x.raw").touch()
    captured = {}

    def fake_render(raw, out_dir, **kwargs):
        captured.update(tic_ymax=kwargs["tic_ymax"], bpc_ymax=kwargs["bpc_ymax"])

    monkeypatch.setattr(cli, "render_png", fake_render)
    main(["--in", str(tmp_path), "--out", str(tmp_path / "out"), "--tic-ymax", "3e10"])
    assert captured == {"tic_ymax": 3e10, "bpc_ymax": None}


def test_ymax_rejects_non_positive_and_non_numeric():
    import pytest

    for bad in ("0", "-1e9", "abc"):
        with pytest.raises(SystemExit):
            build_parser().parse_args(["--in", "i", "--out", "o", "--tic-ymax", bad])


def test_peak_labels_alias_still_accepted():
    """--peak-labels is kept as a deprecated alias and maps to the same dest."""
    args = build_parser().parse_args(["--in", "i", "--out", "o", "--peak-labels", "7"])
    assert args.top_peaks == 7


def test_main_label_spacing_passed_through(tmp_path, monkeypatch):
    (tmp_path / "x.raw").touch()
    captured = {}

    def fake_render(raw, out_dir, ms_level=1, dpi=200, n_peak_labels=5, label_spacing=0.05,
                    show_mz=False, show_charge=False, show_intensity=True, show_rt=True,
                    label_rotation=-60, mz_decimals=4, rt_decimals=2,
                    max_label_corner="left", tic_ymax=None, bpc_ymax=None):
        captured["label_spacing"] = label_spacing
        return None

    monkeypatch.setattr(cli, "render_png", fake_render)
    main(["--in", str(tmp_path), "--out", str(tmp_path / "out"), "--label-spacing", "0.1"])
    assert captured["label_spacing"] == 0.1
