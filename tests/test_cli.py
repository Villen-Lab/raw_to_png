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

    def fake_render(raw, out_dir, ms_level=1, dpi=200, n_peak_labels=5, label_spacing=0.05):
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

    def fake_render(raw, out_dir, ms_level=1, dpi=200, n_peak_labels=5, label_spacing=0.05):
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
    assert args.peak_labels == 5
    assert args.label_spacing == 0.05


def test_main_peak_labels_passed_through(tmp_path, monkeypatch):
    (tmp_path / "x.raw").touch()
    captured = {}

    def fake_render(raw, out_dir, ms_level=1, dpi=200, n_peak_labels=5, label_spacing=0.05):
        captured["n_peak_labels"] = n_peak_labels
        return None

    monkeypatch.setattr(cli, "render_png", fake_render)
    main(["--in", str(tmp_path), "--out", str(tmp_path / "out"), "--peak-labels", "3"])
    assert captured["n_peak_labels"] == 3


def test_main_label_spacing_passed_through(tmp_path, monkeypatch):
    (tmp_path / "x.raw").touch()
    captured = {}

    def fake_render(raw, out_dir, ms_level=1, dpi=200, n_peak_labels=5, label_spacing=0.05):
        captured["label_spacing"] = label_spacing
        return None

    monkeypatch.setattr(cli, "render_png", fake_render)
    main(["--in", str(tmp_path), "--out", str(tmp_path / "out"), "--label-spacing", "0.1"])
    assert captured["label_spacing"] == 0.1
