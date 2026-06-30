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

# How many of the most abundant peaks to label, and how far apart (as a fraction
# of the retention-time span) two labelled peaks must be so we don't tag several
# adjacent scans of the same peak.
DEFAULT_PEAK_LABELS = 5
PEAK_MIN_SEP_FRAC = 0.01

# How close (as a fraction of the RT span) two peak labels may be before the
# later one is staggered onto a higher row to avoid overlap. ~one label width.
DEFAULT_LABEL_SPACING = 0.05


def top_peaks(rt, y, n=DEFAULT_PEAK_LABELS, min_sep_frac=PEAK_MIN_SEP_FRAC):
    """Return indices of the ``n`` most abundant, well-separated peaks.

    Greedy: take the highest point, then the next-highest that is at least
    ``min_sep_frac`` of the RT span away from every point already chosen, and so
    on. This keeps the picks on *distinct* chromatographic peaks rather than
    several scans straddling one apex. Returned indices are sorted by RT.
    """
    if y.size == 0:
        return np.array([], dtype=int)
    span = float(rt.max() - rt.min())
    min_sep = min_sep_frac * span if span > 0 else 0.0

    chosen = []
    for idx in np.argsort(y)[::-1]:  # highest intensity first
        if len(chosen) >= n:
            break
        if all(abs(rt[idx] - rt[c]) >= min_sep for c in chosen):
            chosen.append(int(idx))
    return np.array(sorted(chosen), dtype=int)


def _sci_label(v):
    """Format an intensity as mathtext scientific notation, e.g. 3.80x10^10."""
    if v <= 0:
        return "0"
    exp = int(np.floor(np.log10(v)))
    mant = v / 10**exp
    return rf"${mant:.2f}\times10^{{{exp}}}$"


def _stagger_levels(rt, idx, x_thresh_frac=DEFAULT_LABEL_SPACING):
    """Assign a stacking level to each labelled peak so close labels don't overlap.

    ``idx`` is the peak indices sorted by RT. A peak within ``x_thresh_frac`` of
    the RT span of the previous labelled peak is bumped one level higher than it,
    so a horizontal cluster of labels stacks 0, 1, 2, ... up the page. Isolated
    peaks reset to level 0. Larger ``x_thresh_frac`` => labels spread apart more
    readily (stack sooner); smaller => labels stay on the baseline row longer.
    """
    if len(idx) == 0:
        return []
    x_span = float(rt.max() - rt.min())
    x_thresh = x_thresh_frac * x_span
    levels = [0]
    for k in range(1, len(idx)):
        close = (rt[idx[k]] - rt[idx[k - 1]]) < x_thresh
        levels.append(levels[-1] + 1 if close else 0)
    return levels


def _annotate_top_peaks(ax, rt, y, color, n=DEFAULT_PEAK_LABELS,
                        label_spacing=DEFAULT_LABEL_SPACING):
    """Mark and label the n most abundant peaks with their intensity.

    Labels stay horizontally centred over their peak (so the dot identifies which
    peak each belongs to) but are staggered vertically when peaks crowd together;
    ``label_spacing`` controls how readily that staggering kicks in.
    """
    if n <= 0 or y.size == 0:
        return
    idx = top_peaks(rt, y, n=n)
    ax.scatter(rt[idx], y[idx], s=14, color=color, zorder=5)

    levels = _stagger_levels(rt, idx, x_thresh_frac=label_spacing)
    # Anchor every label in a horizontal cluster (a run of increasing levels) to
    # that cluster's tallest apex, so the point-offset rows separate the labels
    # cleanly even when the peaks themselves differ in height.
    baselines = [float(v) for v in y[idx]]
    start = 0
    for k in range(1, len(idx) + 1):
        if k == len(idx) or levels[k] == 0:
            base = max(baselines[start:k])
            for j in range(start, k):
                baselines[j] = base
            start = k

    for pos, (i, level) in enumerate(zip(idx, levels)):
        ax.annotate(
            _sci_label(y[i]),
            xy=(rt[i], baselines[pos]),
            xytext=(0, 4 + 13 * level),  # 13 pt ~ one label row per level
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=8,
            color=color,
        )
    # Headroom so the topmost (possibly staggered) label isn't clipped.
    headroom = 1.18 + 0.12 * (max(levels) if levels else 0)
    ax.set_ylim(top=float(y.max()) * headroom)


def _annotate_max(ax, y, color, label):
    """Write the max intensity in a large bold label in the panel's top-left."""
    if y.size == 0:
        return
    ax.text(
        0.012, 0.95, f"{label}: {_sci_label(float(y.max()))}",
        transform=ax.transAxes, ha="left", va="top",
        fontsize=13, fontweight="bold", color=color,
        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=color, alpha=0.85),
    )


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


def render_png(raw_path, out_dir, ms_level=1, dpi=200, n_peak_labels=DEFAULT_PEAK_LABELS,
               label_spacing=DEFAULT_LABEL_SPACING):
    """Render a 2-panel (TIC / base peak) PNG for one .RAW file.

    The ``n_peak_labels`` most abundant peaks in each panel are marked and
    labelled with their intensity in scientific notation (0 disables labelling).
    ``label_spacing`` (fraction of the RT span) tunes how close two labels may be
    before the later one staggers onto a higher row.

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

    # Label the most abundant peaks in each panel...
    _annotate_top_peaks(ax_tic, rt, tic, TIC_COLOR, n=n_peak_labels,
                        label_spacing=label_spacing)
    _annotate_top_peaks(ax_bpc, rt, bpc, BPC_COLOR, n=n_peak_labels,
                        label_spacing=label_spacing)
    # ...and the panel maximum, large, in the top-left corner.
    _annotate_max(ax_tic, tic, TIC_COLOR, "Max TIC")
    _annotate_max(ax_bpc, bpc, BPC_COLOR, "Max base peak")

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)  # gitignored; absent on fresh clone
    out_path = out_dir / f"{name}.png"
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)
    print(f"  -> {out_path}")
    return out_path
