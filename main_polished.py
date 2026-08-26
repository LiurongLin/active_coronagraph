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

from basic import *            # noqa: F401,F403
from make_plot import *        # noqa: F401,F403
from coronagraphs_polished import *     # noqa: F401,F403  (simple_coro, res_ene, etc.)
from main_functions_polished import *   # noqa: F401,F403  (run, broadband_combine_final_psf, etc.)
from main_functions_polished import _root_folder, _run_dir
from new_mask import directional_throughput_mono, throughput_map_mono, throughput_1d_mono
import inspect
import os, glob
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

# Broadband setup
# wavelength = np.array(100) / np.arange(90, 114, 2)  # array around 1.0
# wavelength = [0.9, 0.95, 1.0, 1.05, 1.1] # the ratio between lambda and lambda_0
wavelength = [1.0]
weight = np.ones(len(wavelength))
weight[[0, -1]] = 0.5
ang = 45
folder_name = ""
phase_screen_folder = Path("phase_screen")
pixel_noise_mean_output_dir = Path("/media/liurong/My Passport/PLACID")

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
        contrast = annuli_mean(
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


def _combine_pixel_noise_annuli_mean(
    mask_name: str,
    charge_val: Optional[int],
    n_realizations: int = 1,
    source_prefix: str = "_pixel_noise",
    output_folder: str = "_pixel_noise_mean",
):
    """Average annular-contrast curves over repeated pixel-noise realizations."""
    curve_name = (
        f"final_focal_plane_coro_mask={mask_name}_{charge_val}_sample={nsamp}_"
        f"fpm_sam={fpm_sam}_obstruction_{obstruction}_bb.txt"
    )

    arrays: list[np.ndarray] = []
    for i in range(n_realizations):
        file_path = (
            _root_folder(obstruction, noise_level, f"{source_prefix}_{i}")
            / f"lyot_{lyot_stop}_{mask_name}_{charge_val}_lambda_1.0_fpm_sam={fpm_sam}_"
              f"binary_False_obstruction_{obstruction}_greyscale_8"
            / curve_name
        )
        if not file_path.exists():
            continue
        arrays.append(np.loadtxt(file_path, dtype=float))

    if not arrays:
        raise FileNotFoundError(
            f"No pixel-noise annuli-mean files found for {mask_name=} {charge_val=}"
        )

    mean_curve = np.mean(np.asarray(arrays), axis=0)
    noise_folder = f"noise_{noise_level:g}" if noise_level is not None else "noise_none"
    run_dir = pixel_noise_mean_output_dir / noise_folder / output_folder / (
        f"lyot_{lyot_stop}_{mask_name}_{charge_val}_lambda_1.0_fpm_sam={fpm_sam}_"
        f"binary_False_obstruction_{obstruction}_greyscale_8"
    )
    output_path = run_dir / curve_name
    run_dir.mkdir(parents=True, exist_ok=True)
    np.savetxt(output_path, mean_curve)
    return output_path


def add_pixel_noise(mask_name: str, charge_val: Optional[int], phase: Optional[np.ndarray] = None):
    n_realizations = 1
    for i in range(n_realizations):
        folder_n = f"_pixel_noise_{i}"
        main_func(mask_name, charge_val, folder = folder_n)
        broadband_combine_final_psf(wavelength, fpm_sam, nsamp, weight,
                                    mask_name, charge_val,
                                    obstruction, lyot_stop, noise_level=noise_level,
                                    folder_name=folder_n)

        annuli_mean(mask_name, nsamp, fpm_sam, lyot_stop, charge_val, obstruction,
                    wl=1.0, rotate=False, greyscale=8, broadband=True,
                    noise_level=noise_level, folder_name=folder_n)

    _combine_pixel_noise_annuli_mean(mask_name, charge_val, n_realizations=n_realizations)

    for i in range(n_realizations):
        folder_n = f"_pixel_noise_{i}"
        fits_dir = os.path.join(_root_folder(obstruction, noise_level, folder_n), "")  # same as used by run()
        for f in glob.glob(os.path.join(fits_dir, "**", "final_focal_plane_coro_mask=*.fits"), recursive=True):
            try:
                os.remove(f)
                print(f"Deleted temporary file: {f}")
            except OSError as e:
                print(f"Warning: could not delete {f} — {e}")


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
    local_obstruction = True
    base_compare_folder = "_phase_screen_compare"
    baseline_folder = f"{base_compare_folder}/no_screen"
    repo_root = Path(__file__).resolve().parent

    if not phase_screen_folder.exists():
        raise FileNotFoundError(f"Phase-screen folder not found: {phase_screen_folder}")

    os.chdir(repo_root)
    try:
        main_func(mask_name, charge_val, phase=None, folder=baseline_folder)
        baseline_curve = annuli_mean(
            mask_name, nsamp, fpm_sam, lyot_stop, charge_val, local_obstruction,
            wl=1.0, rotate=False, greyscale=8, broadband=False,
            folder_name=baseline_folder,
        )
        os.chdir(repo_root)

        curves: list[tuple[str, np.ndarray]] = [("no_screen", np.asarray(baseline_curve))]

        for fits_path in sorted(phase_screen_folder.glob("*.fits")):
            case_tag = _phase_screen_case_tag(fits_path)
            case_folder = f"{base_compare_folder}/{case_tag}"

            with fits.open(fits_path) as hdul:
                phase_cube = hdul[0].data

            for i, phase_screen in enumerate(phase_cube):
                os.chdir(repo_root)
                main_func(mask_name, charge_val, phase=phase_screen, folder=f"{case_folder}/{i}")

            os.chdir(repo_root)
            combine_phase_screen_psf(
                mask_name, charge_val, local_obstruction,
                folder_name=case_folder, n_realizations=len(phase_cube),
            )

            case_curve = annuli_mean(
                mask_name, nsamp, fpm_sam, lyot_stop, charge_val, local_obstruction,
                wl=1.0, rotate=False, greyscale=8, broadband=False,
                folder_name=case_folder,
            )
            os.chdir(repo_root)
            curves.append((case_tag, np.asarray(case_curve)))

        r_mid = (np.arange(len(curves[0][1])) + 0.5) / 10.0
        out_dir = _root_folder(local_obstruction, None, base_compare_folder)
        out_dir.mkdir(parents=True, exist_ok=True)

        fig, ax = plt.subplots(figsize=(8, 6))
        for label, curve in curves:
            ax.plot(r_mid, curve, label=label)
        ax.set_yscale("log")
        ax.set_xlabel(r"$\\lambda$/D")
        ax.set_ylabel("Median contrast")
        ax.set_title(f"Phase-screen comparison: {mask_name} {charge_val}")
        ax.legend()
        fig.tight_layout()
        fig.savefig(out_dir / f"phase_screen_comparison_{mask_name}_{charge_val}.png", dpi=200)
        plt.close(fig)
    finally:
        os.chdir(repo_root)


def res_ene_collect(mask_name: str, charge_val: Optional[int], phase: Optional[np.ndarray] = None, folder: str = folder_name) -> None:
    """Collect residual energy across wavelengths using `run(..., do_res=True)` and save to txt."""
    res_ene_list: list[float] = []
    for wl in wavelength:
        metric = run(mask_name, float(wl), nsamp, fpm_sam,
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
                            'tp_bb', 'plot_tp', 'fpm_plot', 'broadband_kit', 'add_phase_screen', 'add_fill_factor',
                            'pixel_noise', 'ghost_im', 'compare_phase_screens', 'tp_directional_mono',
                            'tp_map_mono', 'tp_1d_mono'],
                   required=True, help='Which workflow to run.')
    p.add_argument('--name', nargs='*', default=default_name,
                   help='List of phase mask names to use.')
    p.add_argument('--charge', nargs='*', type=parse_charge, default=default_charge,
                   help='List of vortex charges to use.')
    p.add_argument('--lyot', nargs='?', type=str2bool, default=False,
                   help='Whether to plot Lyot planes (if supported by downstream calls).')
    p.add_argument('--wavelength', nargs='*', default=wavelength,
                   help='Wavelength list (ratio λ/λ0).')
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
    return p


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    global phase_shift
    global phase_shift_unit
    global coro_shift
    phase_shift = (float(args.phase_shift[0]), float(args.phase_shift[1]))
    phase_shift_unit = args.phase_shift_unit
    coro_shift = args.coro_shift

    # Set globals based on args when appropriate (kept minimal to preserve behavior)
    # For simplicity we leave module-level defaults unless you want full arg plumbing.

    if args.function == 'main_func':
        for mask_name, ch in zip(args.name, args.charge):
            main_func(mask_name, ch)


    elif args.function == 'pixel_noise':
        print(simple_coro.__module__)
        print(inspect.getsourcefile(simple_coro))
        for mask_name, ch in zip(args.name, args.charge):
            add_pixel_noise(mask_name, ch)

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
            broadband_combine_final_psf(wavelength, fpm_sam, nsamp, weight,
                                        mask_name, ch, obstruction, lyot_stop,
                                        noise_level=noise_level)

    elif args.function == 'broadband_kit':
        # original range was np.arange(7,10,1) which exceeds default_name length; clamp safely
        for i in tqdm(range(7, min(10, len(default_name)))):
            main_func(default_name[i], default_charge[i])
            broadband_combine_final_psf(wavelength, fpm_sam, nsamp, weight,
                                        default_name[i], default_charge[i],
                                        obstruction, lyot_stop, noise_level=noise_level,
                                        folder_name=folder_name)
            annuli_mean(default_name[i], nsamp, fpm_sam, lyot_stop, default_charge[i],
                obstruction, wl=1.0, rotate=False, greyscale=8, broadband=True,
                noise_level=noise_level, folder_name=folder_name)

    elif args.function == 'add_fill_factor':
        for mask_name, ch in zip(args.name, args.charge):
            run(mask_name, 1.0, 100, 10, lyot_stop=lyot_stop, charge=ch, binary=False,
                obstruction=False, get_lyot=get_lyot, rotate=False, greyscale=8,
                cal_factor=cal_factor, phase=None, folder_name="_fill_factor",
                fill_factor=True, noise_level=noise_level, phase_shift=phase_shift,
                phase_shift_unit=phase_shift_unit, coro_shift=coro_shift)
            annuli_mean(mask_name, 100, 10, lyot_stop, ch, False, wl=1.0,
                        rotate=False, greyscale=8, broadband=False,
                        noise_level=noise_level, folder_name="_fill_factor")

    elif args.function == 'add_phase_screen':
        # This branch previously used os.chdir(".."); we avoid that.
        local_obstruction = True
        for mask_name, ch in zip(args.name, args.charge):
            for i in range(100):
                try:
                    with fits.open("2)seeing=1,vmag=8,ZA=30,lag=2/TROIA_phase_screens_new.fits") as hdul:
                        phase_screen = hdul[0].data[i]
                except FileNotFoundError:
                    print("FileNotFoundError: TROIA_phase_screens_new.fits")
                    break  # nothing to process
                main_func(mask_name, ch, phase=phase_screen, folder=f"_phase_screen_lag2/{i}")
            combine_phase_screen_psf(mask_name, ch, local_obstruction)
            annuli_mean(mask_name, nsamp, fpm_sam, lyot_stop, ch, local_obstruction,
                        wl=1.0, rotate=False, greyscale=8, broadband=False,
                        folder_name="_phase_screen_lag2")

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
        wl = float(args.wavelength[0]) if len(args.wavelength) > 0 else 1.0
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
        wl = float(args.wavelength[0]) if len(args.wavelength) > 0 else 1.0
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
        wl = float(args.wavelength[0]) if len(args.wavelength) > 0 else 1.0
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
