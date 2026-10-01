from __future__ import annotations

"""Utility functions for running and post-processing coronagraph simulations.

This module keeps the original public API (function names and parameters)
but improves style, safety, and maintainability:
- Safer filesystem handling with a temporary working-directory context manager.
- More robust path construction via pathlib (no stringly-typed paths).
- Clearer variable names; avoid shadowing built-ins like ``dir``.
- Type hints and richer docstrings.
- Optional logging hooks and consistent directory creation.
- Input validation (e.g., matching broadband lists and weights).

NOTE: We purposely keep wildcard imports from user modules (basic, make_plot)
since their exported symbols are not enumerated here and the original code
relies on them. If you know the exact names, prefer importing explicitly.
"""

from contextlib import contextmanager
from pathlib import Path
from typing import Iterable, List, Optional, Sequence

import numpy as np
from astropy.io import fits

# External project modules (kept as wildcard to preserve behavior)
# External project modules (imported explicitly to avoid wildcard pollution)
from . import coronagraphs_polished
from . import make_plot
from .config import output_dir



# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

@contextmanager
def pushd(path: Path):
    """Temporarily change the working directory.

    Parameters
    ----------
    path : Path
        Directory to switch into for the duration of the context.
    """
    import os

    old = Path.cwd()
    path.mkdir(parents=True, exist_ok=True)
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(old)


def _root_folder(obstruction: bool, noise_level: Optional[float], folder_name: Optional[str]) -> Path:
    """Compute the root output folder path for a run."""
    base = f"ideal_coro_2rd_mirror_{obstruction}"
    if noise_level is not None:
        base = f"{base}_noise_{noise_level:g}"
    if folder_name:
        return output_dir() / f"{base}{folder_name}"
    return output_dir() / base


def _run_dir(
    *,
    lyot_stop: Optional[str],
    name: str,
    charge: Optional[int],
    wavelength: Optional[float],
    nsamp: int,
    fpm_sam: int,
    binary: bool,
    obstruction: bool,
    greyscale: Optional[int | bool],
    offset: Optional[Sequence[float]],
    phase_shift: Optional[Sequence[float]] = (0.0, 0.0),
    phase_shift_unit: str = "array",
) -> Path:
    """Build the per-run directory name.

    Mirrors the original naming scheme while avoiding shadowing ``dir``.
    """
    # Common prefix
    prefix = f"lyot_{lyot_stop}_{name}_{charge}"

    if offset is not None:
        if greyscale is not None:
            tail = (
                f"sam_{nsamp}_fpm_sam={fpm_sam}_binary_{binary}_"
                f"obstruction_{obstruction}_greyscale_{greyscale}_offset_True"
            )
        else:
            tail = (
                f"sam_{nsamp}_fpm_sam={fpm_sam}_binary_{binary}_"
                f"obstruction_{obstruction}_offset_True"
            )
    else:
        if greyscale is not None:
            tail = (
                f"lambda_{wavelength}_fpm_sam={fpm_sam}_binary_{binary}_"
                f"obstruction_{obstruction}_greyscale_{greyscale}"
            )
        else:
            tail = (
                f"lambda_{wavelength}_fpm_sam={fpm_sam}_binary_{binary}_"
                f"obstruction_{obstruction}"
            )

    shift_suffix = ""
    if phase_shift is not None:
        dy, dx = float(phase_shift[0]), float(phase_shift[1])
        if dy != 0.0 or dx != 0.0:
            shift_suffix = f"_phase_shift_{dy:g}_{dx:g}_{phase_shift_unit}"

    return Path(f"{prefix}_{tail}{shift_suffix}")


# -----------------------------------------------------------------------------
# Public API
# -----------------------------------------------------------------------------

def run(
    dim: int,
    name: str,
    wavelength: float,
    nsamp: int,
    fpm_sam: int,
    charge: Optional[int] = None,
    lyot_stop: Optional[str] = None,
    binary: bool = True,
    obstruction: bool = True,
    spider: bool = False,
    get_lyot: bool = False,
    rotate: bool = False,
    greyscale: Optional[int | bool] = True,
    offset: Optional[Sequence[float]] = None,
    roddier_d: float = 1.06,
    do_res: bool = False,
    cal_factor: float = 1.0,
    phase=None,
    folder_name: Optional[str] = None,
    noise_level: Optional[float] = None,
    fill_factor: bool = False,
    ghost: bool = False,
    phase_shift: Sequence[float] = (0.0, 0.0),
    phase_shift_unit: str = "array",
    coro_shift: bool = True,
    mask_opt_params: Optional[dict] = None,
):
    """Set up folders and execute a coronagraph simulation or resolution-energy run.

    Parameters mirror the original function. See original docstring for details.
    """
    root = _root_folder(obstruction, noise_level, folder_name)
    run_dir = _run_dir(
        lyot_stop=lyot_stop,
        name=name,
        charge=charge,
        wavelength=wavelength,
        nsamp=nsamp,
        fpm_sam=fpm_sam,
        binary=binary,
        obstruction=obstruction,
        greyscale=greyscale,
        offset=offset,
        phase_shift=phase_shift,
        phase_shift_unit=phase_shift_unit,
    )

    with pushd(root / run_dir):
        if do_res:
            # Keep behavior identical to original
            res_ene_value = coronagraphs_polished.res_ene(
                dim,
                name,
                wavelength,
                nsamp,
                fpm_sam,
                lyot_stop=lyot_stop,
                charge=charge,
                binary=binary,
                obstruction=obstruction,
                get_lyot=get_lyot,
                rotate=rotate,
                greyscale=greyscale,
                offset=offset,
                roddier_d=roddier_d,
                phase_shift=phase_shift,
                phase_shift_unit=phase_shift_unit,
                coro_shift=coro_shift,
            )
            return res_ene_value

        # Regular simulation
        result = coronagraphs_polished.simple_coro(
            dim,
            name,
            wavelength,
            nsamp,
            fpm_sam,
            lyot_stop=lyot_stop,
            charge=charge,
            binary=binary,
            obstruction=obstruction,
            spider=spider,
            get_lyot=get_lyot,
            rotate=rotate,
            greyscale=greyscale,
            offset=offset,
            roddier_d=roddier_d,
            cal_factor=cal_factor,
            phase=phase,
            noise_level=noise_level,
            fill_factor=fill_factor,
            ghost = ghost,
            phase_shift=phase_shift,
            phase_shift_unit=phase_shift_unit,
            coro_shift=coro_shift,
            mask_opt_params=mask_opt_params,
        )
        return result


def get_through_put(
    name: str,
    nsamp: int,
    fpm_sam: int,
    offset: Sequence[float],
    ang: float,
    wl: Optional[float] = None,
    charge: Optional[int] = None,
    lyot_stop: Optional[str] = None,
    binary: bool = False,
    obstruction: bool = True,
    get_lyot: bool = True,
    greyscale: Optional[int | bool] = True,
    bb: bool = False,
):
    """Compute throughput at an offset angle and save outputs into the standard tree."""
    root = _root_folder(obstruction, None, None)

    # Build container directory for throughput runs
    if bb:
        top = root / f"lyot_{lyot_stop}_tp_bb_wl_{wl}_ang_{ang}"
    else:
        top = root / f"lyot_{lyot_stop}_tp/"

    # Per-run directory
    run_dir = (
        f"lyot_{lyot_stop}_{name}_{charge}_sam_{nsamp}_fpm_sam={fpm_sam}_"
        f"binary_{binary}_obstruction_{obstruction}_greyscale_{greyscale}_{offset[0] * 2}"
    )

    with pushd(top / run_dir):
        coronagraphs_polished.simple_coro(
            100,
            name,
            wl,
            nsamp,
            fpm_sam,
            lyot_stop=lyot_stop,
            charge=charge,
            binary=binary,
            obstruction=obstruction,
            get_lyot=get_lyot,
            rotate=False,
            greyscale=greyscale,
            offset=offset,
        )


def combine_throughput(
    name: str,
    nsamp: int,
    fpm_sam: int,
    offset: Sequence[float],
    ang: float,
    charge: Optional[int] = None,
    lyot_stop: Optional[str] = None,
    binary: bool = True,
    obstruction: bool = True,
    greyscale: Optional[int | bool] = True,
    bb: bool = False,
    wl: Optional[float] = None,
):
    """Aggregate throughput values over offsets and save a combined file.

    Fixes a path selection bug from the original implementation: the broadband
    (``bb=True``) case now correctly uses the "bb_wl_{wl}" container; the
    monochromatic case uses "lyot_{lyot_stop}_tp/".
    """
    root = _root_folder(obstruction, None, None)

    if bb:
        base = root / f"lyot_{lyot_stop}_tp_bb_wl_{wl}_ang_{ang}"
    else:
        base = root / f"lyot_{lyot_stop}_tp"

    throughput_list: List[np.ndarray] = []

    for off in offset:
        run_dir = (
            f"lyot_{lyot_stop}_{name}_{charge}_sam_{nsamp}_fpm_sam={fpm_sam}_"
            f"binary_{binary}_obstruction_{obstruction}_greyscale_{greyscale}_{off * 2}"
        )
        with pushd(base / run_dir if not bb else base):
            if bb:
                tp = np.loadtxt(f"tp_bb_{name}_{charge}_{off * 2}.txt")
            else:
                tp = np.loadtxt("throughput.txt")
            throughput_list.append(tp)

    out_name = (
        f"lyot_{lyot_stop}_{name}_{charge}_sam_{nsamp}_fpm_sam={fpm_sam}_"
        f"binary_{binary}_obstruction_{obstruction}_greyscale_{greyscale}"
    )

    if bb:
        np.savetxt(base / f"{out_name}_bb.txt", throughput_list)
    else:
        np.savetxt(base / f"{out_name}_wl_{wl}.txt", throughput_list)


def broadband_combine_final_psf(
    wavelength: Sequence[float],
    fpm_sam: int,
    nsamp: int,
    weights: Sequence[float],
    name: str,
    charge: Optional[int],
    obstruction: bool,
    lyot_stop: Optional[str],
    noise_level: Optional[float] = None,
    folder_name: str | None = "",
):
    """Weighted-average the final PSF over multiple wavelengths and save FITS.

    Returns
    -------
    np.ndarray
        The weighted-average PSF array.
    """
    if len(wavelength) != len(weights):
        raise ValueError("wavelength and weights must have the same length")

    root = _root_folder(obstruction, noise_level, folder_name)

    base_dir = (
        root
        / _run_dir(
            lyot_stop=lyot_stop,
            name=name,
            charge=charge,
            wavelength=1.0,
            nsamp=nsamp,
            fpm_sam=fpm_sam,
            binary=False,
            obstruction=obstruction,
            greyscale=8,
            offset=None,
        )
    )

    base_path = (
        base_dir
        / f"final_focal_plane_coro_mask={name}_{charge}_sample={nsamp}_fpm_sam={fpm_sam}_obstruction_{obstruction}.fits"
    )
    output_path = base_path.with_name(
        f"final_focal_plane_coro_mask={name}_{charge}_sample={nsamp}_fpm_sam={fpm_sam}_obstruction_{obstruction}_bb.fits"
    )

    arrays = []
    for lam in wavelength:
        lam_path = Path(str(base_path).replace("lambda_1.0", f"lambda_{lam}"))
        with fits.open(lam_path) as hdul:
            arrays.append(hdul[0].data)

    arrays = np.asarray(arrays)
    weights = np.asarray(weights)

    weighted_average = np.tensordot(weights, arrays, axes=1) / weights.sum()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    hdu = fits.PrimaryHDU(weighted_average)
    hdu.writeto(output_path, overwrite=True)

    return weighted_average


def broadband_combine_tp(
    name: str,
    charge: Optional[int],
    obstruction: bool,
    wavelength: Sequence[float],
    lyot_stop: Optional[str],
    greyscale: int | bool,
    ang: float,
):
    """Combine saved throughput curves across wavelengths using fixed weights."""
    root = _root_folder(obstruction, None, None)
    bb_tp: List[np.ndarray] = []

    # NOTE: Weights preserved from original implementation
    wl_weight = np.asarray([0.5, 1.0, 1.0, 1.0, 0.5])

    base_path = (
        root
        / f"lyot_{lyot_stop}_tp_bb_wl_1.0_ang_{ang}/"
        / f"lyot_{lyot_stop}_{name}_{charge}_sam_10_fpm_sam=10_binary_False_"
        f"obstruction_True_greyscale_{greyscale}_wl_1.0.txt"
    )

    for lam in wavelength:
        file_path = Path(str(base_path).replace("wl_1.0", f"wl_{lam}"))
        with open(file_path, "r"):
            bb_tp.append(np.loadtxt(file_path))

    bb_tp = np.asarray(bb_tp)
    if wl_weight.shape[0] != bb_tp.shape[0]:
        raise ValueError("Weight count must match number of wavelength files.")

    tp = (wl_weight @ bb_tp) / 4.0

    out = (
        root
        / f"lyot_{lyot_stop}_{name}_{charge}_sam_10_fpm_sam=10_binary_False_"
        f"obstruction_True_greyscale_8_bb_ang_{ang}.txt"
    )
    np.savetxt(out, tp)


def psf_into_fits(wl: float, name: str, charge: Optional[int], offset: Sequence[float]):
    """Stack PSFs over offsets into a FITS cube and save via ``wfits``."""
    root = Path(f"ideal_coro_2rd_mirror_True/lyot_None_1_1_tp_bb_wl_{wl}")
    psf_cube: List[np.ndarray] = []

    for off in offset:
        txt = (
            root
            / f"lyot_None_1_1_{name}_{charge}_sam_10_fpm_sam=10_binary_False_"
            f"obstruction_True_greyscale_8_{off}/final_psf.txt"
        )
        psf_cube.append(np.loadtxt(txt))

    make_plot.wfits(psf_cube, root / f"psf_cube_{name}_{charge}.fits")


def combine_phase_screen_psf(
    name: str,
    charge: Optional[int],
    obstruction: bool,
    folder_name: str = "_phase_screen_lag2",
    n_realizations: int = 100,
):
    """Average PSFs over multiple phase-screen realizations and save a FITS file."""
    base = _root_folder(obstruction, None, folder_name)

    # Where the averaged file will be written
    output_path = (
        base
        / f"lyot_None_1_1_{name}_{charge}_lambda_1.0_fpm_sam=10_binary_False_"
        f"obstruction_{obstruction}_greyscale_8/"
        f"final_focal_plane_coro_mask={name}_{charge}_sample=10_fpm_sam=10_obstruction_{obstruction}.fits"
    )

    arrays: List[np.ndarray] = []
    for i in range(n_realizations):
        file_path = (
            base
            / f"{i}/"
            / f"lyot_None_1_1_{name}_{charge}_lambda_1.0_fpm_sam=10_binary_False_"
            f"obstruction_{obstruction}_greyscale_8/"
            f"final_focal_plane_coro_mask={name}_{charge}_sample=10_fpm_sam=10_obstruction_{obstruction}.fits"
        )
        if not file_path.exists():
            continue
        with fits.open(file_path) as hdul:
            arrays.append(hdul[0].data)

    if not arrays:
        raise FileNotFoundError(
            f"No phase-screen PSFs found for {name=} {charge=} in {base}"
        )

    mean_psf = np.mean(arrays, axis=0)[0]
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Save with project-provided convenience writer
    make_plot.wfits([mean_psf], output_path)
