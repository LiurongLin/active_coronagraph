from __future__ import annotations

import os
from pathlib import Path


PACKAGE_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = PACKAGE_DIR.parents[1]


def _path_from_env(name: str, default: Path) -> Path:
    value = os.environ.get(name)
    if value:
        return Path(value).expanduser().resolve()
    return default


def data_dir() -> Path:
    """Return the directory containing bundled or user-provided runtime data."""
    return _path_from_env("ACTIVE_CORONAGRAPH_DATA_DIR", PROJECT_ROOT / "data")


def output_dir() -> Path:
    """Return the base directory used for generated simulation outputs."""
    return _path_from_env("ACTIVE_CORONAGRAPH_OUTPUT_DIR", Path.cwd())


def pixel_noise_output_dir() -> Path:
    """Return the base directory used by the pixel-noise averaging workflow."""
    return _path_from_env("ACTIVE_CORONAGRAPH_PIXEL_NOISE_OUTPUT_DIR", output_dir())


def phase_screen_dir() -> Path:
    """Return the directory containing phase-screen FITS cubes."""
    return data_dir() / "phase_screen"


def binary_mask_dir() -> Path:
    """Return the directory containing binary-mask FITS files."""
    return data_dir() / "binary_mask"


def legacy_phase_screen_file() -> Path:
    """Return the optional legacy single phase-screen cube path."""
    return _path_from_env(
        "ACTIVE_CORONAGRAPH_PHASE_SCREEN_FILE",
        phase_screen_dir() / "TROIA_phase_screens_new.fits",
    )
