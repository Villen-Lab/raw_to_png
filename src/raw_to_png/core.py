"""Core logic: read chromatograms from a .RAW file and render them to a PNG.

The reader-touching work is split in two so the scan-iteration logic can be
unit-tested with a mock reader on any platform (the real Thermo RawFileReader
only runs on Windows / .NET):

* ``extract_chromatograms(reader, ms_level)`` — pure logic over an *already-open*
  reader object. Mockable; no fisher_py import needed.
* ``get_chromatograms(raw_path, ms_level)`` — opens the .RAW via fisher_py,
  delegates to ``extract_chromatograms``, and disposes the reader. fisher_py is
  imported lazily inside this function so the module imports cleanly where
  fisher_py / .NET is unavailable (e.g. Linux CI).
"""

from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")  # headless: no display needed
import matplotlib.pyplot as plt  # noqa: E402  (must follow matplotlib.use)

# Plot colours, matched to the original prototype.
TIC_COLOR = "#1f4e79"  # dark blue
BPC_COLOR = "#7d2e2e"  # dark red


def extract_chromatograms(reader, ms_level=1):
    """Return ``(rt, tic, bpc)`` numpy arrays from an already-open reader.

    Parameters
    ----------
    reader
        An object exposing the Thermo RawFileReader surface used here:
        ``run_header_ex.first_spectrum`` / ``.last_spectrum``,
        ``get_filter_for_scan_number(n).ms_order`` (a ``MsOrderType`` enum),
        ``get_scan_stats_for_scan_number(n).tic`` / ``.base_peak_intensity``,
        and ``retention_time_from_scan_number(n)``.
    ms_level
        MS order to keep (1 = MS1). ``None`` keeps every scan.

    Notes
    -----
    TIC and base-peak values are read from each scan's *stored statistics*
    (what FreeStyle plots), not recomputed from centroids — a deliberate
    fidelity choice.

    ``ms_order`` is a plain ``Enum`` (``MsOrderType``), **not** an ``IntEnum``,
    so it must be compared via ``.value``; ``int(ms_order)`` raises ``TypeError``.
    MS1 is ``MsOrderType.Ms`` (value 1).
    """
    first = reader.run_header_ex.first_spectrum
    last = reader.run_header_ex.last_spectrum

    rt, tic, bpc = [], [], []
    for scan_no in range(first, last + 1):
        if ms_level is not None:
            ms_order = reader.get_filter_for_scan_number(scan_no).ms_order
            if ms_order.value != ms_level:
                continue
        stats = reader.get_scan_stats_for_scan_number(scan_no)
        rt.append(reader.retention_time_from_scan_number(scan_no))
        tic.append(stats.tic)
        bpc.append(stats.base_peak_intensity)

    return np.asarray(rt), np.asarray(tic), np.asarray(bpc)


def get_chromatograms(raw_path, ms_level=1):
    """Open one .RAW file via fisher_py and return ``(rt, tic, bpc)`` arrays.

    rt  : retention time (minutes) per included scan
    tic : total ion current per scan
    bpc : base-peak intensity per scan
    """
    # Imported lazily: fisher_py pulls in .NET assemblies that are unavailable
    # on platforms used only to import this module (e.g. CI).
    from fisher_py.raw_file_reader import RawFileReaderAdapter
    from fisher_py.data import Device

    reader = RawFileReaderAdapter.file_factory(str(raw_path))
    try:
        reader.select_instrument(Device.MS, 1)
        return extract_chromatograms(reader, ms_level=ms_level)
    finally:
        reader.dispose()


def render_png(raw_path, out_dir, ms_level=1, dpi=200):
    """Render a 2-panel (TIC / base peak) PNG for one .RAW file.

    Returns the output path, or ``None`` if no scans matched ``ms_level``.
    """
    name = Path(raw_path).stem
    rt, tic, bpc = get_chromatograms(raw_path, ms_level=ms_level)

    if rt.size == 0:
        print(f"  ! no MS{ms_level} scans found in {name}, skipping")
        return None

    fig, (ax_tic, ax_bpc) = plt.subplots(
        2, 1, figsize=(10, 6), sharex=True, constrained_layout=True
    )
    fig.suptitle(name, fontsize=13, fontweight="bold")

    ax_tic.plot(rt, tic, linewidth=0.7, color=TIC_COLOR)
    ax_tic.set_ylabel("TIC intensity")
    ax_tic.set_title("Total Ion Chromatogram", fontsize=10, loc="left")
    ax_tic.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

    ax_bpc.plot(rt, bpc, linewidth=0.7, color=BPC_COLOR)
    ax_bpc.set_ylabel("Base peak intensity")
    ax_bpc.set_xlabel("Retention time (min)")
    ax_bpc.set_title("Base Peak Chromatogram", fontsize=10, loc="left")
    ax_bpc.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

    for ax in (ax_tic, ax_bpc):
        ax.margins(x=0)
        ax.grid(True, alpha=0.25)

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)  # gitignored; absent on fresh clone
    out_path = out_dir / f"{name}.png"
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)
    print(f"  -> {out_path}")
    return out_path
