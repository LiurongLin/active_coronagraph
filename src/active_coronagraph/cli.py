from __future__ import annotations
"""
Main CLI entrypoint for coronagraph simulations.

This refactor keeps the original behavior but improves safety and readability:
- Fixes indentation/syntax errors in `main_func` and `res_ene_collect`.
- Adds type hints and docstrings.
- Adds "res_ene" to --function choices (it was referenced but missing).
- Removes fragile `os.chdir("..")` calls (use functions' own path handling).
- Avoids out-of-range indexes in the tp_bb branch.
- Simplifies imports (kept wildcard project imports to avoid breaking behavior).
"""

from typing import List, Optional, Sequence, Tuple
import argparse
import math
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from astropy.io import fits

# Scientific / plotting (kept for compatibility with downstream calls)

from .basic import *            # noqa: F401,F403
from .make_plot import *        # noqa: F401,F403
from .coronagraphs_polished import *     # noqa: F401,F403  (simple_coro, res_ene, etc.)
from .main_functions_polished import *   # noqa: F401,F403  (run, broadband_combine_final_psf, etc.)
from .main_functions_polished import _root_folder, _run_dir
from .new_mask import directional_throughput_mono, throughput_map_mono, throughput_1d_mono
from .config import (
    data_dir as configured_data_dir,
    output_dir as configured_output_dir,
    phase_screen_dir,
)
import os
import hcipy





# -----------------------------------------------------------------------------
# Default parameters (preserved from original file)
# -----------------------------------------------------------------------------

dim = 100
default_name = ["vortex", "vortex", "vortex", "vortex", "FQPM", "roddier", "dual_zone", "ACM", "ACM"]
nsamp = 10
fpm_sam = 10
default_charge = [8, 2, 4, 6, None, None, None, 2, 4]
obstruction = False
high_sam = True
get_lyot = False
lyot_stop_default = [
    "large_cross", "small_cross", "zernike_cross",
    "None_1_1", "None", "None_1_1.1", "None_1_1.3", "None_1_1.4", "None_1_1.5"
]
lyot_stop = lyot_stop_default[0]
fpm_sam_list = [5, 10]
cal_factor = 1.0
noise_level = 0
spider = False
ghost = False
gs = None # greyscale
phase_shift = (0.0, 0.0)
phase_shift_unit = "array"
coro_shift = True

# Wavelength ratios are lambda/lambda_0. Monochromatic workflows default to 1.0;
# broadband workflows default to a 20% band centered on lambda_0.
MONO_WAVELENGTH = [1.0]
BROADBAND_20_PERCENT_WAVELENGTH = [0.9, 0.95, 1.0, 1.05, 1.1]
BROADBAND_FUNCTIONS = {"broadband_combine", "broadband_kit", "tp_bb"}
wavelength = MONO_WAVELENGTH.copy()
weight = np.ones(len(wavelength))
weight[[0, -1]] = 0.5
ang = 45
folder_name = ""
phase_screen_folder = phase_screen_dir()

contrast_ind = np.arange(1, 60, 2)
roddier_d_list = np.arange(1, 1.5, 0.02)
name = default_name[6]
charge = default_charge[6]
# cus_name = "psf_ghost.fits"
cus_name = None

# Throughput offsets
offset_30_1 = np.arange(0, 6, 0.1) * (np.sqrt(3) / 2)
offset_30_2 = np.arange(0, 6, 0.1) * (1 / 2)
offset_15_1 = np.arange(0, 6, 0.1) * math.cos(15)
offset_15_2 = np.arange(0, 6, 0.1) * math.sin(15)
offset_0_1 = np.arange(0, 6, 0.1)
offset_0_2 = np.zeros_like(offset_0_1)


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

def _call_preserving_cwd(func, *args, **kwargs):
    """Call a legacy helper that may change cwd, then restore the caller cwd."""
    old = Path.cwd()
    try:
        return func(*args, **kwargs)
    finally:
        os.chdir(old)


def _call_from_output_dir_preserving_cwd(func, *args, **kwargs):
    """Call a legacy output reader from the configured output directory."""
    old = Path.cwd()
    try:
        out_dir = configured_output_dir()
        out_dir.mkdir(parents=True, exist_ok=True)
        os.chdir(out_dir)
        return func(*args, **kwargs)
    finally:
        os.chdir(old)


def parse_charge(charge_str: str) -> Optional[int]:
    """argparse type for vortex charge that accepts 'None'."""
    if charge_str.lower() == 'none':
        return None
    try:
        return int(charge_str)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"Invalid value for charge: {charge_str}") from exc


def str2bool(v) -> bool:
    """argparse boolean parser supporting many common truthy/falsey spellings."""
    if isinstance(v, bool):
        return v
    s = v.lower()
    if s in ('yes', 'true', 't', 'y', '1'):
        return True
    if s in ('no', 'false', 'f', 'n', '0'):
        return False
    raise argparse.ArgumentTypeError('Boolean value expected.')


def trapezoid_weights(wavelengths: Sequence[float]) -> np.ndarray:
    """Return endpoint-half trapezoid weights for a wavelength grid."""
    weights = np.ones(len(wavelengths), dtype=float)
    if len(weights) > 1:
        weights[[0, -1]] = 0.5
    return weights


def resolve_wavelengths(function_name: str, requested: Optional[Sequence[float]]) -> list[float]:
    """Resolve CLI wavelength ratios for monochromatic and broadband workflows."""
    if requested is not None:
        if len(requested) == 0:
            raise argparse.ArgumentTypeError("--wavelength requires at least one value.")
        return [float(wl) for wl in requested]
    if function_name in BROADBAND_FUNCTIONS:
        return BROADBAND_20_PERCENT_WAVELENGTH.copy()
    return MONO_WAVELENGTH.copy()


# -----------------------------------------------------------------------------
# Workflows
# -----------------------------------------------------------------------------

def main_func(mask_name: str, charge_val: Optional[int], phase: Optional[np.ndarray] = None, folder: str = folder_name) -> None:
    """Run `run()` over all wavelengths for a given mask + charge."""
    for wl in wavelength:
        run(dim, mask_name, float(wl), nsamp, fpm_sam,
            lyot_stop=lyot_stop, charge=charge_val, binary=False,
            obstruction=obstruction, spider=spider, get_lyot=get_lyot,
            rotate=False, greyscale=gs, cal_factor=cal_factor, phase=phase,
            folder_name=folder, noise_level = noise_level, ghost = ghost,
            phase_shift=phase_shift, phase_shift_unit=phase_shift_unit,
            coro_shift=coro_shift, )
        contrast = _call_from_output_dir_preserving_cwd(
            annuli_mean,
            mask_name, nsamp, fpm_sam, lyot_stop, charge_val, obstruction,
            wl=float(wl), rotate=False, greyscale=gs, broadband=False,
            noise_level=noise_level, folder_name=folder,
        )
        r_mid = (np.arange(len(contrast)) + 0.5) / 10.0
        plot_dir = _root_folder(obstruction, noise_level, folder) / _run_dir(
            lyot_stop=lyot_stop,
            name=mask_name,
            charge=charge_val,
            wavelength=float(wl),
            nsamp=nsamp,
            fpm_sam=fpm_sam,
            binary=False,
            obstruction=obstruction,
            greyscale=gs,
            offset=None,
            phase_shift=phase_shift,
            phase_shift_unit=phase_shift_unit,
        )
        plot_dir.mkdir(parents=True, exist_ok=True)
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.plot(r_mid, contrast)
        ax.set_yscale("log")
        ax.set_xlabel(r"$\lambda$/D")
        ax.set_ylabel("Contrast")
        ax.set_title(f"{mask_name} pupil size {dim}, LS 0.95 annular contrast")
        fig.tight_layout()
        fig.savefig(
            plot_dir / (
                f"final_focal_plane_coro_mask={mask_name}_{charge_val}_sample={nsamp}_"
                f"fpm_sam={fpm_sam}_obstruction_{obstruction}_wl_{float(wl)}_pupil_{dim}_annuli_mean.png"
            ),
            dpi=200,
        )
        plt.close(fig)

# def ghost(mask_name: str, charge_val: Optional[int], phase: Optional[np.ndarray] = None, folder: str = folder_name) -> None:


def broadband_final_psf(mask_name: str, charge_val: Optional[int], folder: str = folder_name) -> None:
    """Run the broadband PSF workflow and save the combined product."""
    if len(wavelength) < 2:
        print(
            "Warning: broadband workflow is using fewer than two wavelengths; "
            f"resolved wavelength list is {wavelength}."
        )

    print(f"Running broadband PSF workflow over wavelength ratios: {wavelength}")
    print(f"Using trapezoid weights: {weight.tolist()}")

    # broadband_combine_final_psf() reads greyscale_8 directories, so generate
    # the wavelength-specific products with greyscale=8 here.
    for wl in wavelength:
        run(
            dim,
            mask_name,
            float(wl),
            nsamp,
            fpm_sam,
            lyot_stop=lyot_stop,
            charge=charge_val,
            binary=False,
            obstruction=obstruction,
            spider=spider,
            get_lyot=get_lyot,
            rotate=False,
            greyscale=8,
            cal_factor=cal_factor,
            phase=None,
            folder_name=folder,
            noise_level=noise_level,
            ghost=ghost,
            phase_shift=phase_shift,
            phase_shift_unit=phase_shift_unit,
            coro_shift=coro_shift,
        )

    broadband_combine_final_psf(
        wavelength,
        fpm_sam,
        nsamp,
        weight,
        mask_name,
        charge_val,
        obstruction,
        lyot_stop,
        noise_level=noise_level,
        folder_name=folder,
    )

    annuli_mean(
        mask_name,
        nsamp,
        fpm_sam,
        lyot_stop,
        charge_val,
        obstruction,
        wl=1.0,
        rotate=False,
        greyscale=8,
        broadband=True,
        noise_level=noise_level,
        folder_name=folder,
    )


def _phase_screen_case_tag(path: Path) -> str:
    stem = path.stem
    marker = "jitter"
    suffix = "percentLamdaOverD"
    if marker in stem and suffix in stem:
        start = stem.index(marker) + len(marker)
        end = stem.index(suffix, start)
        return f"{stem[start:end]}%"
    return stem


def compare_phase_screen_cases(mask_name: str, charge_val: Optional[int]) -> None:
    """Compare the no-screen baseline with all phase-screen cubes in `phase_screen/`."""
    local_obstruction = obstruction
    base_compare_folder = "_phase_screen_compare"
    baseline_folder = f"{base_compare_folder}/no_screen"
    start_dir = Path.cwd()

    if not phase_screen_folder.exists():
        raise FileNotFoundError(f"Phase-screen folder not found: {phase_screen_folder}")

    os.chdir(start_dir)
    try:
        main_func(mask_name, charge_val, phase=None, folder=baseline_folder)
        baseline_curve = _call_from_output_dir_preserving_cwd(
            annuli_mean,
            mask_name, nsamp, fpm_sam, lyot_stop, charge_val, local_obstruction,
            wl=1.0, rotate=False, greyscale=gs, broadband=False,
            noise_level=noise_level, folder_name=baseline_folder,
        )
        os.chdir(start_dir)

        curves: list[tuple[str, np.ndarray]] = [("no_screen", np.asarray(baseline_curve))]

        for fits_path in sorted(phase_screen_folder.glob("*.fits")):
            case_tag = _phase_screen_case_tag(fits_path)
            case_folder = f"{base_compare_folder}/{case_tag}"

            with fits.open(fits_path) as hdul:
                phase_cube = hdul[0].data

            for i, phase_screen in enumerate(phase_cube):
                os.chdir(start_dir)
                main_func(mask_name, charge_val, phase=phase_screen, folder=f"{case_folder}/{i}")

            os.chdir(start_dir)
            combine_phase_screen_psf(
                mask_name, charge_val, local_obstruction,
                folder_name=case_folder, n_realizations=len(phase_cube),
                lyot_stop=lyot_stop, nsamp=nsamp, fpm_sam=fpm_sam,
                greyscale=gs, noise_level=noise_level,
            )

            case_curve = _call_from_output_dir_preserving_cwd(
                annuli_mean,
                mask_name, nsamp, fpm_sam, lyot_stop, charge_val, local_obstruction,
                wl=1.0, rotate=False, greyscale=gs, broadband=False,
                noise_level=noise_level, folder_name=case_folder,
            )
            os.chdir(start_dir)
            curves.append((case_tag, np.asarray(case_curve)))

        r_mid = (np.arange(len(curves[0][1])) + 0.5) / 10.0
        out_dir = _root_folder(local_obstruction, noise_level, base_compare_folder)
        out_dir.mkdir(parents=True, exist_ok=True)

        fig, ax = plt.subplots(figsize=(8, 6))
        for label, curve in curves:
            ax.plot(r_mid, curve, label=label)
        ax.set_yscale("log")
        ax.set_xlabel(r"$\lambda$/D")
        ax.set_ylabel("Median contrast")
        ax.set_title(f"Phase-screen comparison: {mask_name} {charge_val}")
        ax.legend()
        fig.tight_layout()
        fig.savefig(out_dir / f"phase_screen_comparison_{mask_name}_{charge_val}.png", dpi=200)
        plt.close(fig)
    finally:
        os.chdir(start_dir)


def res_ene_collect(mask_name: str, charge_val: Optional[int], phase: Optional[np.ndarray] = None, folder: str = folder_name) -> None:
    """Collect residual energy across wavelengths using `run(..., do_res=True)` and save to txt."""
    res_ene_list: list[float] = []
    for wl in wavelength:
        metric = run(dim, mask_name, float(wl), nsamp, fpm_sam,
                     lyot_stop=lyot_stop, charge=charge_val, binary=False,
                     obstruction=obstruction, get_lyot=get_lyot, rotate=False,
                     greyscale=None, cal_factor=cal_factor, phase=phase,
                     folder_name=folder, noise_level=noise_level, do_res=True,
                     phase_shift=phase_shift, phase_shift_unit=phase_shift_unit,
                     coro_shift=coro_shift)
        res_ene_list.append(float(metric))
    np.savetxt(f"res_ene_{mask_name}_{charge_val}_100_dense.txt", res_ene_list)


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description='Run coronagraph pipelines with different parameters.')
    p.add_argument('--function',
                   choices=['main_func', 'res_ene', 'annuli_mean', 'summary_plot', 'broadband_combine',
                            'tp_bb', 'plot_tp', 'fpm_plot', 'broadband_kit', 'add_fill_factor',
                            'ghost_im', 'compare_phase_screens', 'tp_directional_mono',
                            'tp_map_mono', 'tp_1d_mono'],
                   required=True, help='Which workflow to run.')
    p.add_argument('--name', nargs='*', default=default_name,
                   help='List of phase mask names to use.')
    p.add_argument('--charge', nargs='*', type=parse_charge, default=default_charge,
                   help='List of vortex charges to use.')
    p.add_argument('--lyot', nargs='?', type=str2bool, default=False,
                   help='Whether to plot Lyot planes (if supported by downstream calls).')
    p.add_argument('--wavelength', nargs='*', type=float, default=None,
                   help=(
                       'Wavelength list (ratio λ/λ0). Defaults to [1.0] for '
                       'monochromatic workflows and [0.9, 0.95, 1.0, 1.05, 1.1] '
                       'for broadband workflows.'
                   ))
    p.add_argument('--phase-shift', nargs=2, type=float, metavar=('DY', 'DX'), default=(0.0, 0.0),
                   help='Shift PSF center by focal-plane pixels (DY DX) via entrance-pupil phase ramp.')
    p.add_argument('--phase-shift-unit', choices=['array', 'fpm'], default='array',
                   help='Unit of --phase-shift: array pixels or fpm pixels.')
    p.add_argument('--coro-shift', type=str2bool, default=coro_shift,
                   help='Whether to apply the half-pixel coronagraph centering shift.')
    p.add_argument('--dv-sigma-ld', type=float, default=1.0,
                   help='Double-vortex sigma in lambda/D (used when --name double_vortex).')
    p.add_argument('--dv-angle-deg', type=float, default=0.0,
                   help='Double-vortex azimuthal angle in degrees (used when --name double_vortex).')
    p.add_argument('--dv-separation-ld', type=float, default=2.0,
                   help='Double-vortex separation in lambda/D (used when --name double_vortex).')
    p.add_argument('--dv-pa-deg', type=float, default=0.0,
                   help='Double-vortex position angle in degrees (used when --name double_vortex).')
    p.add_argument('--phase-map-fits', type=str, default=None,
                   help='Path to one specific binary FITS mask. Used with --name binary_mask/binary_vortex.')
    p.add_argument('--tp-map-max-ld', type=float, default=20.0,
                   help='Upper bound of the 2D throughput map in lambda/D. The map samples [0, max) on each axis.')
    p.add_argument('--tp-map-min-ld', type=float, default=0.0,
                   help='Lower bound of the 2D throughput map in lambda/D.')
    p.add_argument('--tp-map-step-ld', type=float, default=1.0,
                   help='2D throughput map sampling in lambda/D.')
    p.add_argument('--tp-1d-max-ld', type=float, default=20.0,
                   help='Upper bound of the 1D throughput sweep in lambda/D.')
    p.add_argument('--tp-1d-step-ld', type=float, default=0.1,
                   help='1D throughput sampling in lambda/D.')
    p.add_argument('--data-dir', type=Path, default=None,
                   help='Directory containing runtime data. Defaults to ACTIVE_CORONAGRAPH_DATA_DIR or ./data.')
    p.add_argument('--phase-screen-folder', type=Path, default=None,
                   help='Directory containing phase-screen FITS cubes.')
    return p


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    global wavelength
    global weight
    global phase_shift
    global phase_shift_unit
    global coro_shift
    global phase_screen_folder
    try:
        wavelength = resolve_wavelengths(args.function, args.wavelength)
    except argparse.ArgumentTypeError as exc:
        parser.error(str(exc))
    weight = trapezoid_weights(wavelength)
    phase_shift = (float(args.phase_shift[0]), float(args.phase_shift[1]))
    phase_shift_unit = args.phase_shift_unit
    coro_shift = args.coro_shift
    if args.data_dir is not None:
        import os

        os.environ["ACTIVE_CORONAGRAPH_DATA_DIR"] = str(args.data_dir.expanduser().resolve())
    phase_screen_folder = (
        args.phase_screen_folder.expanduser().resolve()
        if args.phase_screen_folder is not None
        else phase_screen_dir()
    )

    if args.function == 'main_func':
        for mask_name, ch in zip(args.name, args.charge):
            main_func(mask_name, ch)

    elif args.function == 'res_ene':
        for mask_name, ch in zip(args.name, args.charge):
            res_ene_collect(mask_name, ch)

    elif args.function == 'annuli_mean':
        for mask_name, ch in zip(args.name, args.charge):
            for wl in wavelength:
                annuli_mean(mask_name, nsamp, fpm_sam, lyot_stop, ch, obstruction,
                            wl=float(wl), rotate=False, greyscale=8, broadband=False,
                            noise_level=noise_level, folder_name=folder_name, cus_name = cus_name)

    elif args.function == 'summary_plot':
        make_summary_plot("contrast", 0, diff_lyot=False, diff_fpm_sam=False,
                          diff_nsamp=True, diff_bb=False, obstruction=False,
                          greyscale=False, broadband=False)

    elif args.function == 'broadband_combine':
        for mask_name, ch in zip(args.name, args.charge):
            broadband_final_psf(mask_name, ch)

    elif args.function == 'broadband_kit':
        # original range was np.arange(7,10,1) which exceeds default_name length; clamp safely
        for i in tqdm(range(7, min(10, len(default_name)))):
            broadband_final_psf(default_name[i], default_charge[i], folder=folder_name)

    elif args.function == 'add_fill_factor':
        for mask_name, ch in zip(args.name, args.charge):
            run(dim, mask_name, 1.0, 100, 10, lyot_stop=lyot_stop, charge=ch, binary=False,
                obstruction=False, get_lyot=get_lyot, rotate=False, greyscale=8,
                cal_factor=cal_factor, phase=None, folder_name="_fill_factor",
                fill_factor=True, noise_level=noise_level, phase_shift=phase_shift,
                phase_shift_unit=phase_shift_unit, coro_shift=coro_shift)
            annuli_mean(mask_name, 100, 10, lyot_stop, ch, False, wl=1.0,
                        rotate=False, greyscale=8, broadband=False,
                        noise_level=noise_level, folder_name="_fill_factor")

    elif args.function == 'compare_phase_screens':
        for mask_name, ch in zip(args.name, args.charge):
            compare_phase_screen_cases(mask_name, ch)

    elif args.function == 'tp_bb':
        # Safe indices: original had [6,7,8,9] but lists have len 8 -> clamp to valid range
        indices = [i for i in [6, 7, 8, 9] if i < len(default_name)]
        for i in indices:
            ang = 45 if default_name[i] == "FQPM" else 0
            # Combine throughput over offsets for a subset of wavelengths (kept behavior but safe length)
            for j in range(min(5, len(wavelength))):
                combine_throughput(default_name[i], nsamp, fpm_sam,
                                   np.arange(0, 6, 0.1) * math.cos(math.radians(ang)),
                                   ang, charge=default_charge[i], lyot_stop=lyot_stop,
                                   binary=False, obstruction=True, greyscale=8,
                                   bb=False, wl=float(wavelength[j]))

        for i in indices:
            ang = 45 if default_name[i] == "FQPM" else 0
            broadband_combine_tp(default_name[i], default_charge[i], True,
                                  wavelength, lyot_stop, 8, ang)

    elif args.function == 'tp_directional_mono':
        # Uses wl=1.0 by default for monochromatic throughput.
        wl = float(wavelength[0])
        for mask_name, ch in zip(args.name, args.charge):
            out_files = directional_throughput_mono(
                mask_name,
                ch,
                wl=wl,
                max_ld=10.0,
                step_ld=0.1,
                dv_sigma_ld=args.dv_sigma_ld,
                dv_angle_deg=args.dv_angle_deg,
                dv_separation_ld=args.dv_separation_ld,
                dv_pa_deg=args.dv_pa_deg,
                phase_map_fits=args.phase_map_fits,
            )
            if isinstance(out_files, (list, tuple)):
                for out_file in out_files:
                    print(f"Saved directional throughput to: {out_file}")
            else:
                print(f"Saved directional throughput to: {out_files}")

    elif args.function == 'tp_map_mono':
        wl = float(wavelength[0])
        for mask_name, ch in zip(args.name, args.charge):
            out_files = throughput_map_mono(
                mask_name,
                ch,
                wl=wl,
                min_ld=args.tp_map_min_ld,
                max_ld=args.tp_map_max_ld,
                step_ld=args.tp_map_step_ld,
                dv_sigma_ld=args.dv_sigma_ld,
                dv_angle_deg=args.dv_angle_deg,
                dv_separation_ld=args.dv_separation_ld,
                dv_pa_deg=args.dv_pa_deg,
                phase_map_fits=args.phase_map_fits,
            )
            if isinstance(out_files, (list, tuple)):
                for out_file in out_files:
                    print(f"Saved throughput map to: {out_file}")
            else:
                print(f"Saved throughput map to: {out_files}")

    elif args.function == 'tp_1d_mono':
        wl = float(wavelength[0])
        for mask_name, ch in zip(args.name, args.charge):
            out_files = throughput_1d_mono(
                mask_name,
                ch,
                wl=wl,
                max_ld=args.tp_1d_max_ld,
                step_ld=args.tp_1d_step_ld,
                dv_sigma_ld=args.dv_sigma_ld,
                dv_angle_deg=args.dv_angle_deg,
                dv_separation_ld=args.dv_separation_ld,
                dv_pa_deg=args.dv_pa_deg,
                phase_map_fits=args.phase_map_fits,
            )
            if isinstance(out_files, (list, tuple)):
                for out_file in out_files:
                    print(f"Saved 1D throughput to: {out_file}")
            else:
                print(f"Saved 1D throughput to: {out_files}")

    elif args.function == 'fpm_plot':
        make_fpm_image([10, 100], 8, mode="paper")

    elif args.function == 'plot_tp':
        # Placeholder hook if you re-enable plot_tp in your project.
        print("plot_tp is currently commented out in the original file.")

    elif args.function == 'ghost_im':
        for wl in wavelength:
            ghost_im(dim = dim, nsamp = nsamp, fpm_sam = fpm_sam, wavelength = float(wl), lyot_sum = None, phase = None,
            shift = True, save = False, obstruction = False, spider = False)

    else:
        raise ValueError(f"Unknown function: {args.function}")


if __name__ == '__main__':
    main()
