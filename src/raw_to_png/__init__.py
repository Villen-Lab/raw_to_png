"""raw_to_png — render TIC and base-peak chromatogram PNGs from Thermo .RAW files."""

from raw_to_png.core import (
    extract_chromatograms,
    get_chromatograms,
    render_png,
    top_peaks,
)

__version__ = "0.1.0"

__all__ = [
    "extract_chromatograms",
    "get_chromatograms",
    "render_png",
    "top_peaks",
    "__version__",
]
