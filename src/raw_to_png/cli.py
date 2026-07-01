"""Command-line interface: batch-render PNGs for every .RAW in a directory."""

import argparse
import sys
from pathlib import Path

from raw_to_png.core import (
    DEFAULT_LABEL_ROTATION,
    DEFAULT_LABEL_SPACING,
    DEFAULT_MZ_DECIMALS,
    DEFAULT_PEAK_LABELS,
    DEFAULT_RT_DECIMALS,
    render_png,
)


def find_raw_files(in_dir):
    """Return the sorted, de-duplicated list of .RAW files in ``in_dir``.

    Globs both ``*.raw`` and ``*.RAW`` then de-dupes via ``set()``: a
    case-insensitive filesystem (Windows) would otherwise list each file twice.
    """
    in_dir = Path(in_dir)
    raw_files = list(in_dir.glob("*.raw")) + list(in_dir.glob("*.RAW"))
    return sorted(set(raw_files))


def build_parser():
    p = argparse.ArgumentParser(
        prog="raw-to-png",
        description="Render TIC + base-peak chromatogram PNGs from Thermo .RAW files.",
    )
    p.add_argument("--in", dest="in_dir", required=True,
                   help="Directory containing .RAW files")
    p.add_argument("--out", dest="out_dir", required=True,
                   help="Directory to write PNGs (created if absent)")
    p.add_argument("--ms-level", type=int, default=1,
                   help="MS order to plot (1=MS1, default). Use 0 for all scans.")
    p.add_argument("--dpi", type=int, default=200, help="Output PNG resolution")
    p.add_argument("--top-peaks", "--peak-labels", dest="top_peaks", type=int,
                   default=DEFAULT_PEAK_LABELS, metavar="N",
                   help="Mark the N most abundant peaks per panel (default 5; use "
                        "0 to mark none). What each marked peak is labelled with is "
                        "set by --show-mz / --show-charge / --hide-intensity. "
                        "(--peak-labels is a deprecated alias.)")
    p.add_argument("--label-spacing", type=float, default=DEFAULT_LABEL_SPACING,
                   help="How close (fraction of the retention-time span) two peak "
                        "labels may be before the later one is staggered onto a "
                        f"higher row (default {DEFAULT_LABEL_SPACING}). Negative "
                        "keeps every label on the baseline row (the slant separates "
                        "them); raise it to stack crowded labels vertically.")
    p.add_argument("--label-rotation", type=float, default=DEFAULT_LABEL_ROTATION,
                   metavar="DEG",
                   help="Angle (degrees, counter-clockwise) to slant peak labels "
                        f"(default {DEFAULT_LABEL_ROTATION}). Use 0 for upright "
                        "labels.")
    p.add_argument("--show-mz", action="store_true",
                   help="Also label each peak with its scan's base-peak m/z "
                        "(for the TIC panel, the scan's tallest ion).")
    p.add_argument("--show-charge", action="store_true",
                   help="Also label each peak with its scan's base-peak charge "
                        "state. Best-effort: shows 'z=?' when the .RAW assigns no "
                        "charge (common for MS1) and costs an extra read per scan.")
    p.add_argument("--decimals", dest="mz_decimals", type=int,
                   default=DEFAULT_MZ_DECIMALS, metavar="N",
                   help="Decimal places shown for m/z values "
                        f"(default {DEFAULT_MZ_DECIMALS}). Only affects --show-mz.")
    p.add_argument("--rt-decimals", dest="rt_decimals", type=int,
                   default=DEFAULT_RT_DECIMALS, metavar="N",
                   help="Decimal places shown for the peak retention-time label "
                        f"(default {DEFAULT_RT_DECIMALS}). Only affects the "
                        "retention-time line (on unless --hide-rt is given).")
    p.add_argument("--hide-intensity", action="store_true",
                   help="Drop the intensity line from peak labels, leaving only "
                        "--show-mz / --show-charge. Use e.g. --top-peaks 5 "
                        "--show-mz --hide-intensity to label N peaks with m/z "
                        "alone. (The top-left panel-max label is unaffected.)")
    p.add_argument("--hide-rt", action="store_true",
                   help="Drop the retention-time (minutes) line from peak labels. "
                        "Retention time is labelled by default alongside intensity.")
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)

    in_dir = Path(args.in_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    # 0 is the sentinel for "all scans" -> None disables MS-level filtering.
    ms_level = None if args.ms_level == 0 else args.ms_level

    raw_files = find_raw_files(in_dir)
    if not raw_files:
        print(f"No .RAW files found in {in_dir}")
        return 1

    print(f"Found {len(raw_files)} .RAW file(s)")
    for raw in raw_files:
        print(f"Processing {raw.name}")
        try:
            render_png(raw, out_dir, ms_level=ms_level, dpi=args.dpi,
                       n_peak_labels=args.top_peaks,
                       label_spacing=args.label_spacing,
                       show_mz=args.show_mz, show_charge=args.show_charge,
                       show_intensity=not args.hide_intensity,
                       show_rt=not args.hide_rt,
                       label_rotation=args.label_rotation,
                       mz_decimals=args.mz_decimals,
                       rt_decimals=args.rt_decimals)
        except Exception as e:  # fail-soft: one bad .RAW must not abort the batch
            print(f"  ! failed on {raw.name}: {e}")

    print("Done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
