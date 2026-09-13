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

import math
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")  # headless: no display needed
import matplotlib.pyplot as plt

# Plot colours, matched to the original prototype.
TIC_COLOR = "#1f4e79"  # dark blue
BPC_COLOR = "#7d2e2e"  # dark red

# How many of the most abundant peaks to label, and how far apart (as a fraction
# of the retention-time span) two labelled peaks must be so we don't tag several
# adjacent scans of the same peak.
DEFAULT_PEAK_LABELS = 5
PEAK_MIN_SEP_FRAC = 0.01

# How close (as a fraction of the RT span) two peak labels may be before the
# later one is staggered onto a higher row to avoid overlap. The default is
# negative so staggering never triggers (rt gaps are >= 0): rotated labels (see
# DEFAULT_LABEL_ROTATION) clear each other horizontally, so vertical stacking is
# off by default. Raise it (e.g. 0.05) to bring stacking back for upright labels.
DEFAULT_LABEL_SPACING = -0.01

# Angle (degrees, counter-clockwise) the peak labels are drawn at; slanting them
# lets adjacent labels sit side by side without overlapping or needing to stack.
DEFAULT_LABEL_ROTATION = -60

# Label font sizes (points). The m/z line is deliberately a couple of points
# smaller than the intensity / charge lines so it reads as secondary detail.
PEAK_LABEL_FONTSIZE = 8
MZ_LABEL_FONTSIZE = 6

# Perpendicular gap (points) between a label's stacked lines. Kept small so the
# compact-font lines (retention time, m/z) read as one tight block rather than
# drifting apart; the intensity line is only 8 pt so this still clears it.
LABEL_LINE_GAP = 7

# Horizontal nudge (points) applied to every peak label so the slanted block sits
# more centred over its peak — a touch to the right — instead of leaning fully to
# its left. Raise it to shift labels further right, lower (or negative) to shift
# them back over/left of the peak.
LABEL_DX = 8

# Decimal places shown for the m/z label (the only multi-decimal field; intensity
# is scientific notation, charge an integer). High-res Thermo data warrants 4.
DEFAULT_MZ_DECIMALS = 4

# Decimal places shown for the retention-time label (minutes).
DEFAULT_RT_DECIMALS = 2

# Which top corner ("left" or "right") of each panel holds the large panel-max label.
MAX_LABEL_CORNERS = ("left", "right")
DEFAULT_MAX_LABEL_CORNER = "left"


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


def _peak_label_lines(intensity, mz=None, charge=None, rt=None,
                      show_intensity=True, show_mz=False, show_charge=False, show_rt=True,
                      mz_decimals=DEFAULT_MZ_DECIMALS, rt_decimals=DEFAULT_RT_DECIMALS):
    """Build a peak's label as a list of ``(text, fontsize)`` lines.

    Each of intensity / retention time / ``m/z`` / charge is its own opt-in line,
    so a label can carry any combination (e.g. m/z only, with the intensity line
    suppressed via ``show_intensity=False``). Retention time, like intensity, is
    on by default — it is the peak's own position and is always known, unlike
    m/z / charge which are opt-in extras. Each line carries its own font size —
    the m/z and retention-time lines are smaller (``MZ_LABEL_FONTSIZE``) than the
    intensity / charge lines — which is why the lines are drawn as separate
    annotations rather than one ``\\n``-joined string (a single matplotlib Text
    can't mix sizes). ``mz_decimals``/``rt_decimals`` set the m/z and
    retention-time decimal places. The retention-time line is the bare value in
    minutes (no ``RT``/``min`` text). ``mz``/``charge``/``rt`` that are missing or
    unassigned (``None``/non-finite) render as ``?`` rather than being dropped, so
    an opted-in column stays visually aligned across peaks. For
    the TIC panel the supplied ``mz``/``charge`` are the scan's *base peak* (its
    tallest ion), since a TIC value itself has no single m/z. Empty when nothing
    selected.
    """
    lines = []
    if show_intensity:
        lines.append((_sci_label(intensity), PEAK_LABEL_FONTSIZE))
    if show_rt:
        ok = rt is not None and np.isfinite(rt)
        text = f"{rt:.{rt_decimals}f}" if ok else "?"
        lines.append((text, MZ_LABEL_FONTSIZE))
    if show_mz:
        ok = mz is not None and np.isfinite(mz) and mz > 0
        text = rf"$m/z$ {mz:.{mz_decimals}f}" if ok else r"$m/z$ ?"
        lines.append((text, MZ_LABEL_FONTSIZE))
    if show_charge:
        ok = charge is not None and np.isfinite(charge) and charge > 0
        lines.append((f"z={int(charge)}" if ok else "z=?", PEAK_LABEL_FONTSIZE))
    return lines


def _peak_label(intensity, mz=None, charge=None, rt=None,
                show_intensity=True, show_mz=False, show_charge=False, show_rt=True,
                mz_decimals=DEFAULT_MZ_DECIMALS, rt_decimals=DEFAULT_RT_DECIMALS):
    """The ``\\n``-joined text of a peak label (font sizes dropped).

    A convenience over ``_peak_label_lines`` for callers that only need the text.
    """
    return "\n".join(text for text, _ in _peak_label_lines(
        intensity, mz=mz, charge=charge, rt=rt, show_intensity=show_intensity,
        show_mz=show_mz, show_charge=show_charge, show_rt=show_rt,
        mz_decimals=mz_decimals, rt_decimals=rt_decimals))


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
                        label_spacing=DEFAULT_LABEL_SPACING,
                        mz=None, charge=None,
                        show_intensity=True, show_mz=False, show_charge=False,
                        show_rt=True,
                        rotation=DEFAULT_LABEL_ROTATION,
                        mz_decimals=DEFAULT_MZ_DECIMALS,
                        rt_decimals=DEFAULT_RT_DECIMALS):
    """Mark the n most abundant peaks and label each with the selected fields.

    A peak's label is built from any combination of intensity / retention time /
    m/z / charge (see ``_peak_label_lines``); with all four off the peaks are
    still dotted but unlabelled. Labels are anchored at their peak (so the dot
    identifies which peak each belongs to) and slanted by ``rotation`` degrees so
    neighbours sit side by side. ``label_spacing`` still drives the vertical
    stagger, but the default is negative (off) because the slant already
    separates crowded labels; raise it to re-enable stacking.

    Each label's lines are drawn as separate annotations (so the m/z line can use a
    smaller font); they are offset along the rotation-perpendicular so they read as
    one slanted block. ``mz``/``charge`` are arrays parallel to ``rt``; the peak's
    own retention time (in ``rt``) is used directly for the ``show_rt`` line.
    """
    if n <= 0 or y.size == 0:
        return
    idx = top_peaks(rt, y, n=n)
    ax.scatter(rt[idx], y[idx], s=14, color=color, zorder=5)

    levels = _stagger_levels(rt, idx, x_thresh_frac=label_spacing)

    n_lines = int(show_intensity) + int(show_rt) + int(show_mz) + int(show_charge)
    row_pt = 13 * n_lines  # one stacked label's height ~ one row per text line
    # Each label starts at its own peak's tip and rises clear of the trace. The
    # text is anchored at the apex with its bottom there (va="bottom"); the
    # horizontal alignment is chosen from the slant so the body extends into the
    # upper half-plane for any rotation — left when the text reads upward (e.g.
    # +90 deg, vertical), right when it reads downward (e.g. the -60 default),
    # centred when horizontal. Successive lines step along +n, the perpendicular
    # 90 deg counter-clockwise of the reading direction, so they fan outward as a
    # tidy block (the first line sits nearest the dot).
    theta = math.radians(rotation)
    step = (-math.sin(theta), math.cos(theta))
    sin_t = math.sin(theta)
    ha = "left" if sin_t > 1e-9 else "right" if sin_t < -1e-9 else "center"
    line_gap = LABEL_LINE_GAP  # perpendicular gap between a label's stacked lines
    for pos, (i, level) in enumerate(zip(idx, levels)):
        lines = _peak_label_lines(
            y[i],
            mz=mz[i] if mz is not None else None,
            charge=charge[i] if charge is not None else None,
            rt=rt[i],
            show_intensity=show_intensity, show_mz=show_mz, show_charge=show_charge,
            show_rt=show_rt,
            mz_decimals=mz_decimals, rt_decimals=rt_decimals,
        )
        if not lines:  # nothing selected: leave just the dot
            continue
        base_dy = 4 + row_pt * level
        for li, (text, fs) in enumerate(lines):
            ax.annotate(
                text,
                xy=(rt[i], float(y[i])),  # anchor at this peak's own tip
                xytext=(LABEL_DX + step[0] * line_gap * li,
                        base_dy + step[1] * line_gap * li),
                textcoords="offset points",
                ha=ha,
                va="bottom",
                rotation=rotation,
                rotation_mode="anchor",
                fontsize=fs,
                color=color,
            )
    # Headroom so the topmost label isn't clipped. Rotated labels rise off the
    # peak, and the closer to vertical the taller they stand, so scale the reserve
    # with |sin(rotation)|; the per-level stagger step adds more on top.
    max_level = max(levels) if levels else 0
    if rotation % 180 != 0:
        base_head = 0.30 + 0.40 * abs(sin_t)
    else:
        base_head = 0.18 * max(n_lines, 1)
    headroom = 1.0 + base_head + 0.14 * n_lines * max_level
    ax.set_ylim(top=float(y.max()) * headroom)


def _annotate_max(ax, y, color, label, corner=DEFAULT_MAX_LABEL_CORNER):
    """Write the max intensity in a large bold label in the panel's top ``corner``."""
    if corner not in MAX_LABEL_CORNERS:
        raise ValueError(f"corner must be one of {MAX_LABEL_CORNERS}, got {corner!r}")
    if y.size == 0:
        return
    x = 0.012 if corner == "left" else 0.988
    ax.text(
        x, 0.95, f"{label}: {_sci_label(float(y.max()))}",
        transform=ax.transAxes, ha=corner, va="top",
        fontsize=13, fontweight="bold", color=color,
        bbox={"boxstyle": "round,pad=0.3", "fc": "white", "ec": color, "alpha": 0.85},
    )


def _base_peak_charge(reader, scan_no):
    """Best-effort charge state of a scan's base peak, from its centroid stream.

    Returns ``np.nan`` when the centroid stream is unavailable, empty, or carries
    no assigned charge — common for MS1 / profile data, where Thermo leaves the
    charge unset. This is why ``--show-charge`` is a lower-confidence option: the
    value simply does not exist for many scans, and the label renders ``z=?``.

    The reader call is wrapped defensively: a single scan that fails to yield a
    centroid stream must not abort extraction of the whole run.
    """
    try:
        stream = reader.get_centroid_stream(scan_no, False)
    except Exception:  # noqa: BLE001  (best-effort: one bad scan must not abort the run)
        return np.nan
    intensities = getattr(stream, "intensities", None)
    charges = getattr(stream, "charges", None)
    if intensities is None or charges is None:
        return np.nan
    intensities = np.asarray(intensities, dtype=float)
    charges = np.asarray(charges, dtype=float)
    if intensities.size == 0 or charges.size != intensities.size:
        return np.nan
    z = charges[int(np.argmax(intensities))]
    return float(z) if np.isfinite(z) and z > 0 else np.nan


def extract_chromatograms(reader, ms_level=1, want_charge=False):
    """Return ``(rt, tic, bpc, bpm, charge)`` numpy arrays from an open reader.

    Parameters
    ----------
    reader
        An object exposing the Thermo RawFileReader surface used here:
        ``run_header_ex.first_spectrum`` / ``.last_spectrum``,
        ``get_filter_for_scan_number(n).ms_order`` (a ``MsOrderType`` enum),
        ``get_scan_stats_for_scan_number(n).tic`` / ``.base_peak_intensity`` /
        ``.base_peak_mass``, ``retention_time_from_scan_number(n)``, and — only
        when ``want_charge`` — ``get_centroid_stream(n, False)``.
    ms_level
        MS order to keep (1 = MS1). ``None`` keeps every scan.
    want_charge
        When True, also resolve each scan's base-peak charge via its centroid
        stream (see ``_base_peak_charge``). This costs one extra reader call per
        kept scan, so it is opt-in; when False, ``charge`` is all-``nan``.

    Returns
    -------
    rt, tic, bpc, bpm, charge : np.ndarray
        Per kept scan: retention time (min), total ion current, base-peak
        intensity, base-peak m/z, and base-peak charge (``nan`` where unknown).

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

    rt, tic, bpc, bpm, charge = [], [], [], [], []
    for scan_no in range(first, last + 1):
        if ms_level is not None:
            ms_order = reader.get_filter_for_scan_number(scan_no).ms_order
            if ms_order.value != ms_level:
                continue
        stats = reader.get_scan_stats_for_scan_number(scan_no)
        rt.append(reader.retention_time_from_scan_number(scan_no))
        tic.append(stats.tic)
        bpc.append(stats.base_peak_intensity)
        bpm.append(stats.base_peak_mass)
        charge.append(_base_peak_charge(reader, scan_no) if want_charge else np.nan)

    return (np.asarray(rt), np.asarray(tic), np.asarray(bpc),
            np.asarray(bpm), np.asarray(charge))


def get_chromatograms(raw_path, ms_level=1, want_charge=False):
    """Open one .RAW file via fisher_py and return the chromatogram arrays.

    rt     : retention time (minutes) per included scan
    tic    : total ion current per scan
    bpc    : base-peak intensity per scan
    bpm    : base-peak m/z per scan
    charge : base-peak charge per scan (``nan`` unless ``want_charge``)
    """
    # Imported lazily: fisher_py pulls in .NET assemblies that are unavailable
    # on platforms used only to import this module (e.g. CI).
    from fisher_py.data import Device
    from fisher_py.raw_file_reader import RawFileReaderAdapter

    reader = RawFileReaderAdapter.file_factory(str(raw_path))
    try:
        reader.select_instrument(Device.MS, 1)
        return extract_chromatograms(reader, ms_level=ms_level, want_charge=want_charge)
    finally:
        reader.dispose()


def render_png(raw_path, out_dir, ms_level=1, dpi=200, n_peak_labels=DEFAULT_PEAK_LABELS,
               label_spacing=DEFAULT_LABEL_SPACING, show_mz=False, show_charge=False,
               show_intensity=True, show_rt=True, label_rotation=DEFAULT_LABEL_ROTATION,
               mz_decimals=DEFAULT_MZ_DECIMALS, rt_decimals=DEFAULT_RT_DECIMALS,
               max_label_corner=DEFAULT_MAX_LABEL_CORNER, tic_ymax=None, bpc_ymax=None):
    """Render a 2-panel (TIC / base peak) PNG for one .RAW file.

    The ``n_peak_labels`` most abundant peaks in each panel are marked; ``0``
    disables marking and labelling entirely. Labels are slanted by
    ``label_rotation`` degrees; ``label_spacing`` (fraction of the RT span) tunes
    how close two labels may be before the later one staggers onto a higher row
    (negative => stacking off, relying on the slant to separate them).

    Each marked peak's label is composed of independently selectable fields:
    ``show_intensity`` (the intensity in scientific notation, on by default),
    ``show_rt`` (the peak's own retention time in minutes, on by default),
    ``show_mz`` (the scan's base-peak m/z), and ``show_charge`` (the base-peak
    charge, best-effort — ``z=?`` where the .RAW assigns none). Turn ``show_intensity``
    and/or ``show_rt`` off to label the peaks with m/z and/or charge alone. For the
    TIC panel the m/z and charge refer to that scan's tallest ion, since a TIC value
    has no single m/z. ``mz_decimals``/``rt_decimals`` set how many decimal places
    the m/z and retention time are shown to. ``max_label_corner`` ("left" or
    "right") picks which top corner holds each panel's large max-intensity label.
    ``tic_ymax``/``bpc_ymax`` fix the top of the TIC / base-peak y-axis (e.g. 3e10)
    so runs share a scale; ``None`` (default) autoscales with label headroom.
    Signal above a fixed limit is clipped at the top of the panel.

    Returns the output path, or ``None`` if no scans matched ``ms_level``.
    """
    for flag, ymax in (("tic_ymax", tic_ymax), ("bpc_ymax", bpc_ymax)):
        if ymax is not None and not ymax > 0:
            raise ValueError(f"{flag} must be a positive number, got {ymax!r}")

    name = Path(raw_path).stem
    rt, tic, bpc, bpm, charge = get_chromatograms(
        raw_path, ms_level=ms_level, want_charge=show_charge)

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

    # Label the most abundant peaks in each panel. Both panels carry the same
    # per-scan base-peak m/z and charge (for TIC that's the scan's tallest ion).
    _annotate_top_peaks(ax_tic, rt, tic, TIC_COLOR, n=n_peak_labels,
                        label_spacing=label_spacing, mz=bpm, charge=charge,
                        show_intensity=show_intensity, show_mz=show_mz,
                        show_charge=show_charge, show_rt=show_rt,
                        rotation=label_rotation,
                        mz_decimals=mz_decimals, rt_decimals=rt_decimals)
    _annotate_top_peaks(ax_bpc, rt, bpc, BPC_COLOR, n=n_peak_labels,
                        label_spacing=label_spacing, mz=bpm, charge=charge,
                        show_intensity=show_intensity, show_mz=show_mz,
                        show_charge=show_charge, show_rt=show_rt,
                        rotation=label_rotation,
                        mz_decimals=mz_decimals, rt_decimals=rt_decimals)
    # ...and the panel maximum, large, in the chosen top corner.
    _annotate_max(ax_tic, tic, TIC_COLOR, "Max TIC", corner=max_label_corner)
    _annotate_max(ax_bpc, bpc, BPC_COLOR, "Max base peak", corner=max_label_corner)
    # A fixed y-axis top overrides the autoscaled headroom set while labelling.
    if tic_ymax is not None:
        ax_tic.set_ylim(top=tic_ymax)
    if bpc_ymax is not None:
        ax_bpc.set_ylim(top=bpc_ymax)

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)  # gitignored; absent on fresh clone
    out_path = out_dir / f"{name}.png"
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)
    print(f"  -> {out_path}")
    return out_path
