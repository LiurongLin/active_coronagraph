from __future__ import annotations
"""
Coronagraph utilities: pupil/stop generation, focal-plane mask creation,
and end-to-end coronagraph simulation pipelines.

This refactor keeps the original public API (function names & parameters)
but improves readability, safety, and maintainability:

- Remove duplicate/unused imports; organize remaining imports.
- Add docstrings and type hints.
- Use consistent naming (e.g., `fpm_sam` throughout).
- Safer math and array handling; avoid shadowing built-ins.
- Fix small bugs (e.g., undefined `mask_params` -> `opt_params` in `focal_mask_new`).
- Keep wildcard imports from project modules to avoid breaking behavior.
"""

# --- External project / scientific stack ---
from typing import Optional, Sequence, Tuple, List
import numpy as np
from astropy.io import fits  # noqa: F401  (used indirectly by wfits)
from mpl_toolkits.axes_grid1 import make_axes_locatable  # noqa: F401 (left for backward compat)

# Project modules (wildcard preserved: functions are used by name in this file)
from .basic import *          # noqa: F401,F403  (Isum, phi_ramp, round_to_even, etc.)
from .phase_masks import *    # noqa: F401,F403  (Phase_masks)
from .new_mask import double_vortex_phase_mask

import matplotlib.pyplot as plt  # only used for optional viz in `focal_mask_new`


# -----------------------------------------------------------------------------
# Pupil and spiders
# -----------------------------------------------------------------------------

def make_entrance_pupil(dim: int, d: Optional[float] = None, nor: bool = True) -> np.ndarray:
    """
    Build a binary entrance pupil of size (dim x dim).

    Parameters
    ----------
    dim : int
        Output array size (pixels).
    d : float, optional
        Diameter of the clear aperture in pixels. If None, defaults to `dim`.
    nor : bool, default True
        If True, normalize total intensity to 100 using project helper `Isum`.

    Returns
    -------
    np.ndarray (complex)
        Complex field of the entrance pupil (unit amplitude, zero phase by default).
    """
    if d is None:
        d = float(dim)

    im = np.zeros((dim, dim))
    ampli = circle_mask(im, (dim - 1) // 2, (dim - 1) // 2, d / 2)
    phase = np.zeros_like(im)
    A = ampli * np.exp(1j * phase)

    if nor:
        # Normalize intensity to 100 intensity units using project helper Isum(A)
        A = A / np.sqrt(Isum(A) / 100.0)

    return A


def casted_spider(p1: Sequence[float], p2: Sequence[float], spider_width: float):
    """Return a float64-casted spider callable for supersampling evaluation."""
    return lambda grid: make_spider(p1, p2, spider_width)(grid).astype(np.float64)


def spider_arms(dims: int, pupil_grid, normalized: bool = True) -> np.ndarray:
    """
    Create a four-arm spider obscuration pattern evaluated on `pupil_grid`.

    Parameters
    ----------
    dims : int
        Linear dimension of the pupil (pixels), used when `normalized=True` for scaling.
    pupil_grid : Grid
        hcipy grid where the aperture is evaluated.
    normalized : bool, default True
        If True, convert physical units to pixel units using `dims` as pupil diameter.

    Returns
    -------
    np.ndarray
        Product of four supersampled spiders, reshaped to the grid.
    """
    pupil_diameter = 4.0  # meters
    spider_width = 1.958e-2  # meters
    spider_offset = np.array([0.0, 0.37251])  # meters

    if normalized:
        scale = dims / pupil_diameter
        spider_width *= scale
        spider_offset = (spider_offset * scale).tolist()

    # Four mirror edges (45-degree arms)
    s = dims / (2 * np.sqrt(2))
    mirror_edge1 = (s,  s)
    mirror_edge2 = (-s, s)
    mirror_edge3 = (s, -s)
    mirror_edge4 = (-s, -s)

    spider1 = casted_spider(spider_offset, mirror_edge1, spider_width)
    spider2 = casted_spider(spider_offset, mirror_edge2, spider_width)
    spider3 = casted_spider([-v for v in spider_offset], mirror_edge3, spider_width)
    spider4 = casted_spider([-v for v in spider_offset], mirror_edge4, spider_width)

    spider1 = evaluate_supersampled(spider1, pupil_grid, 8)
    spider2 = evaluate_supersampled(spider2, pupil_grid, 8)
    spider3 = evaluate_supersampled(spider3, pupil_grid, 8)
    spider4 = evaluate_supersampled(spider4, pupil_grid, 8)

    return spider1 * spider2 * spider3 * spider4


# -----------------------------------------------------------------------------
# Apertures & Lyot stop
# -----------------------------------------------------------------------------

def make_aperture(
    dim: int,
    oversample: int,
    obstruction: bool = True,
    spider: bool = False,
    obstruction_ratio: float = 0.25,
    lyot: bool = False,
    lyot_fraction: Optional[float] = None,
    nsamp: int = 10,
    nor: bool = True,
    phase: Optional[np.ndarray] = None,
) -> np.ndarray:
    """
    Build an (optionally obstructed) circular aperture or a Lyot stop.

    Notes
    -----
    - If `lyot` is True, the pupil is upsampled to `dims = 100 * nsamp` and
      the physical pupil diameter is scaled by `lyot_fraction` (expected in (0,1]).
    - Uses supersampling (`evaluate_supersampled`) to reduce sampling artifacts.
    """
    if lyot:
        dims = dim * nsamp
        # scale physical diameter by lyot_fraction (fallback to 1.0 if not given)
        lyot_fraction = 1.0 if lyot_fraction is None else lyot_fraction
        effective_dim = dim * lyot_fraction
    else:
        dims = dim
        effective_dim = dim

    pupil_grid = make_pupil_grid(dims, dims)

    # Build geometric aperture (binary real function on the grid)
    if obstruction:
        aperture = make_obstructed_circular_aperture(effective_dim, obstruction_ratio)
    else:
        aperture = make_circular_aperture(effective_dim)

    # Supersample to mitigate staircasing
    aperture = evaluate_supersampled(aperture, pupil_grid, oversample)

    # Optionally apply spiders
    if spider:
        aperture *= spider_arms(dims, pupil_grid)

    A = np.reshape(aperture, (dims, dims))

    if lyot:
        # Return binary Lyot stop
        return A

    # Complex field with optional phase
    if phase is None:
        phase = np.zeros_like(A)
    A = A * np.exp(1j * phase)

    return A / np.sqrt(Isum(A) / 1.0) if nor else A


def make_lyot_stop(
    dim: int,
    fraction: float = 0.95,
    sec_fraction: float = 1.1,
    obstruction: bool = True,
    nsamp: int = 10,
) -> np.ndarray:
    """Simulate the Lyot stop component for the coronagraph."""
    return make_aperture(
        dim,
        oversample=8,
        obstruction=obstruction,
        obstruction_ratio=0.25 * sec_fraction / fraction,
        lyot=True,
        lyot_fraction=fraction,
        nsamp=nsamp,
        nor=False,
    )


# -----------------------------------------------------------------------------
# Focal-plane mask generation
# -----------------------------------------------------------------------------

def focal_mask_new(
    name: str,
    dim: int,
    nsamp: int,
    fpm_sam: int,
    wavelength: float,
    charge: Optional[int],
    opt_params: Optional[dict] = None,
    greyscale: Optional[int] = None,
    greyscale_method: str = "legacy",
    cal_factor: float = 1.0,
    fill_factor: bool = False,
    noise_level: Optional[float] = None,
) -> np.ndarray:
    """
    Build and resample a focal-plane phase mask (FQPM/ACM/vortex/roddier/dual_zone/double_vortex).

    Parameters
    ----------
    name : {"FQPM", "ACM", "vortex", "roddier", "dual_zone"}
    dim : int
        Field-of-view in units of λ/D.
    nsamp : int
        Sampling per λ/D on the numerical grid.
    fpm_sam : int
        Sampling per λ/D of the SLM grid from which the phase is generated.
    wavelength : float
        Wavelength scaling factor used in the final mask phase.
    charge : int or None
        Charge for vortex/ACM masks.
    opt_params : dict
        Optional parameters: {"rotate": bool, "roddier_d": float, "plot_diff": bool}.
    greyscale : int or None
        If set, quantize phase to 2**greyscale levels (with `cal_factor` scale).
    cal_factor : float
        Phase calibration (multiplicative) before quantization.
    fill_factor : bool
        If True, zero every 10th row/column to simulate dead lines.
    noise_level : float or None
        If set, add Gaussian noise with std = noise_level * 2π.

    Returns
    -------
    np.ndarray (complex)
        Complex-valued focal-plane mask.
    """
    if opt_params is None:
        opt_params = {"rotate": False, "roddier_d": 1.06, "vortex_center": 0}

    dims = dim * nsamp
    phase_mask = Phase_masks(fpm_sam=fpm_sam, pupil_size=dim)

    # Base phase (radians in [0, 2π] typically)
    if name == "FQPM":
        phase = phase_mask.FQPM(center=opt_params.get("vortex_center", 0)).astype(float, copy=True)
    elif name == "ACM":
        phase = phase_mask.ACM(charge)
    elif name == "vortex":
        phase = phase_mask.vortex(
            charge,
            opt_params.get("rotate", False),
            center=opt_params.get("vortex_center", 0),
        )
    elif name == "roddier":
        phase = phase_mask.roddier(opt_params.get("roddier_d", 1.06))
    elif name == "dual_zone":
        phase = phase_mask.dual_zone(a1=0.515, a2=0.705, z1=0.47, z2=0.92)
    elif name == "dual_zone_old":
        phase = phase_mask.dual_zone_old(a1=0.515, a2=0.705, z1=0.47, z2=0.92)
    elif name == "double_vortex":
        phase = double_vortex_phase_mask(
            dim=dim,
            fpm_sam=fpm_sam,
            charge=charge,
            separation_ld=opt_params.get("dv_separation_ld", 2.0),
            pa_deg=opt_params.get("dv_pa_deg", 0.0),
            field_rotation_deg=opt_params.get("dv_field_rotation_deg", 0.0),
            invert_field_rotation=opt_params.get("dv_invert_field_rotation", False),
            angle_deg=opt_params.get("dv_angle_deg", 0.0),
            flux_ratio=opt_params.get("dv_flux_ratio", 1.0),
            sigma_ld=opt_params.get("dv_sigma_ld", 1.0),
            alpha=opt_params.get("dv_alpha", 1.0),
            r0_ld=opt_params.get("dv_r0_ld", 0.5),
            flip_lr=opt_params.get("dv_flip_lr", False),
        )
    elif name == "binary_fits":
        phase_map_fits = opt_params.get("phase_map_fits")
        if not phase_map_fits:
            raise ValueError("binary_fits requires opt_params['phase_map_fits'].")
        phase = fits.getdata(phase_map_fits).astype(float, copy=False)
        expected_shape = (int(dim * fpm_sam), int(dim * fpm_sam))
        if phase.shape != expected_shape:
            raise ValueError(
                f"binary_fits mask shape {phase.shape} does not match expected {expected_shape} "
                f"for dim={dim}, fpm_sam={fpm_sam}."
            )
    else:
        raise ValueError(f"Unknown focal mask type: {name!r}")

    # Resample phase to the numerical grid
    scale = nsamp / fpm_sam
    phase_rs = np.repeat(np.repeat(phase, scale, axis=0), scale, axis=1)

    # Roddier: ensure binary {0, π}
    if name == "roddier":
        thres = 0.8
        phase_rs = np.where(phase_rs > thres, np.pi, 0.0)

    # Optional grayscale quantization (apply after calibration)
    if greyscale is not None:
        max_phase = 2 * np.pi * cal_factor
        n_levels = 2 ** greyscale

        if greyscale_method == "new":
            phase_scaled = np.mod(phase_rs * cal_factor, max_phase)
            idx = np.floor(phase_scaled / max_phase * n_levels).astype(int)
            idx = np.clip(idx, 0, n_levels - 1)
            levels = np.arange(n_levels) * max_phase / n_levels
            phase_rs = levels[idx]
        elif greyscale_method == "legacy":
            bins = np.linspace(0.0, max_phase, n_levels, endpoint=True)
            idx = np.digitize(phase_rs * cal_factor, bins) - 1
            idx = np.clip(idx, 0, bins.size - 1)
            phase_rs = bins[idx]
        else:
            raise ValueError(
                f"Unknown greyscale_method: {greyscale_method!r}. "
                "Use 'legacy' or 'new'."
            )

    # Normalize/shift phase following the original logic
    m = np.ones((dims, dims), dtype=float)
    phase_rs = phase_rs - np.pi
    phase_rs = np.angle(m * np.exp(1j * phase_rs)) + np.pi
    phase_rs = phase_rs * (1.0 / wavelength) - np.pi

    # Fill factor emulation: zero out stripes on *amplitude* mask
    if fill_factor:
        m[::10, :] = 0.0
        m[:, ::10] = 0.0

    # Optional plot
    if opt_params.get("plot_diff", False):
        plt.figure(figsize=(6, 5))
        im = plt.imshow(phase_rs, cmap="jet")
        plt.colorbar(im)
        plt.title(f"Phase Mask: {name} (Charge={charge if charge is not None else 'N/A'})")
        plt.tight_layout()
        plt.show()

    # Add Gaussian phase noise if requested
    if noise_level is not None:
        # rng = np.random.default_rng()  # Automatically seeded from system entropy
        noise = np.random.normal(0.0, noise_level*2*np.pi, phase_rs.shape)
        phase_rs = phase_rs+noise

    m_vor = m * np.exp(1j * phase_rs)

    # Save the phase (angle) for reference
    wfits(
        [np.angle(m_vor) + np.pi],
        f"m_vor_{name}_{charge}_nsamp={nsamp}_fpm_sam={fpm_sam}_rotate={opt_params.get('rotate', False)}"
        f"_greyscale={greyscale}_gmethod={greyscale_method}_wl_{wavelength}_cal_{cal_factor}.fits",
    )

    return m_vor


# -----------------------------------------------------------------------------
# Core propagation
# -----------------------------------------------------------------------------

def coro_high_sam(
    A: np.ndarray,
    m: np.ndarray,
    N: np.ndarray,
    dim: int,
    nsamp: int,
    fpm_sam: int,
    charge: Optional[int] = None,
    lyot_sum: Optional[float] = None,
    shift: bool = True,
    save: bool = False,
    get_lyot: bool = False,
    ghost: bool = False,
):
    """
    High-sampling coronagraph propagation with optional Lyot normalization.

    If `get_lyot` is True, returns the integrated intensity after Lyot plane.
    Otherwise returns (final_focal_plane_field, res_ene).
    """
    if shift:
        print("shift")
        # Half-pixel phase ramp to center sampling symmetrically in pupil
        ph = np.zeros_like(A, dtype=float)
        X, Y = np.meshgrid(
            np.linspace(-np.pi, np.pi, A.shape[0], endpoint=False),
            np.linspace(-np.pi, np.pi, A.shape[1], endpoint=False),
        )
        A = A * np.exp(1j * (ph - 0.5 / nsamp * X - 0.5 / nsamp * Y))

    # First focal plane
    A = FFT(A, nsamp * dim)

    # Apply focal mask (multiplicative in focal plane) then propagate back
    A = m * A
    A = IFFT(A, 1)

    if save and lyot_sum is not None:
        # Save a small cut-out profile of the Lyot plane normalized by `lyot_sum`
        lyot_plane = (np.abs(A) ** 2 / lyot_sum)
        wfits(lyot_plane, f"lyot_plane_profile_{nsamp}_{charge}_{fpm_sam}.fits")

    # Plane after Lyot stop
    A_pre_lyot_sum = np.sum(np.abs(A) ** 2)
    A = N * A
    res_ene = float(np.sum(np.abs(A) ** 2) / A_pre_lyot_sum)

    if save and lyot_sum is not None:
        after_lyot = (np.abs(A) ** 2 / lyot_sum)
        wfits(after_lyot, f"after_lyot_plane_profile_{nsamp}_{charge}_{fpm_sam}.fits")

    if get_lyot:
        return float(np.sum(np.abs(A) ** 2))

    # Final focal plane
    A = FFT(A, nsamp * dim)
    return A, res_ene



def ghost_im(
    dim: int,
    nsamp: int,
    fpm_sam: int,
    wavelength: float = 1,
    charge: Optional[int] = None,
    lyot_sum: Optional[float] = None,
    phase: Optional[np.ndarray] = None,
    shift: bool = True,
    save: bool = False,
    obstruction: bool = False,
    spider: bool = False,
):
    """
    High-sampling coronagraph propagation with optional Lyot normalization.

    If `get_lyot` is True, returns the integrated intensity after Lyot plane.
    Otherwise returns (final_focal_plane_field, res_ene).
    """

    lyot_stop_para = {
        "small_cross": [0.98, 1.2],
        "large_cross": [0.95, 1.3],
        "zernike_cross": [0.98, 1.05],
        "None": [1.0, 1.2],
        "None_1_1": [1.0, 1.0],
        "None_1_1.1": [1.0, 1.1],
        "None_1_1.3": [1.0, 1.3],
        "None_1_1.4": [1.0, 1.4],
        "None_1_1.5": [1.0, 1.5],
    }
    A = make_aperture(round_to_even(dim / wavelength), 8,
                      obstruction=obstruction, spider=spider,
                      obstruction_ratio=0.25, phase=phase)
    A_1 = A
    N = make_lyot_stop(int(dim / wavelength),
                               fraction=lyot_stop_para["None_1_1"][0],
                               sec_fraction=lyot_stop_para["None_1_1"][1],
                               obstruction=obstruction)
    if shift:
        # Half-pixel phase ramp to center sampling symmetrically in pupil
        ph = np.zeros_like(A, dtype=float)
        X, Y = np.meshgrid(
            np.linspace(-np.pi, np.pi, A.shape[0], endpoint=False),
            np.linspace(-np.pi, np.pi, A.shape[1], endpoint=False),
        )
        A = A * np.exp(1j * (ph - 0.5 / nsamp * X - 0.5 / nsamp * Y))

    # First focal plane
    A = FFT(A, nsamp * 100)

    # refractive ghost
    A = np.sqrt(0.005) * A
    # second pupil plane
    A = IFFT(A, 1)

    if save and lyot_sum is not None:
        # Save a small cut-out profile of the Lyot plane normalized by `lyot_sum`
        lyot_plane = (np.abs(A) ** 2 / lyot_sum)[400:600, 400:600]
        wfits(lyot_plane, f"lyot_plane_profile_{nsamp}_{charge}_{fpm_sam}_ghost.fits")

    A = N * A

    if save and lyot_sum is not None:
        after_lyot = (np.abs(A) ** 2 / lyot_sum)[400:600, 400:600]
        wfits(after_lyot, f"after_lyot_plane_profile_{nsamp}_{charge}_{fpm_sam}_ghost.fits")

    # Final focal plane
    A = FFT(A, nsamp * 100)
    m_vor_ones = np.ones_like(N)
    wo_mask_peak = np.max(
        np.abs(
            coro_high_sam(A_1, m_vor_ones, N, 100, nsamp, fpm_sam,
                          charge=charge, lyot_sum=None, shift=shift, save=False)[0]
        ) ** 2)
    wfits([np.abs(A)**2/wo_mask_peak], f"psf_ghost.fits")
    return A


# -----------------------------------------------------------------------------
# User-facing pipelines
# -----------------------------------------------------------------------------

def simple_coro(
    dim: int,
    name: str,
    wavelength: float,
    nsamp: int,
    fpm_sam: int,
    lyot_stop: Optional[str] = None,
    offset: Optional[Sequence[float]] = None,
    charge: Optional[int] = None,
    binary: bool = False,
    obstruction: bool = True,
    spider: bool = False,
    get_lyot: bool = True,
    rotate: bool = False,
    greyscale: bool | int = True,
    greyscale_method: str = "legacy",
    roddier_d: float = 1.06,
    cal_factor: float = 1.0,
    phase: Optional[np.ndarray] = None,
    noise_level: Optional[float] = None,
    fill_factor: bool = False,
    ghost: bool = False,
    phase_shift: Sequence[float] = (0.0, 0.0),
    phase_shift_unit: str = "array",
    vortex_center: float = 0,
    coro_shift: bool = True,
    mask_opt_params: Optional[dict] = None,
):
    """
    End-to-end coronagraph simulation with optional Lyot normalization.

    Returns either the normalized final focal plane (if not throughput mode),
    or a scalar throughput value when `offset` is provided.
    """
    # Entrance pupil (binary or smooth edge via analytic aperture)
    if binary:
        A = make_entrance_pupil(dim)
    else:
        A = make_aperture(round_to_even(dim / wavelength), 8,
                          obstruction=obstruction, spider=spider,
                          obstruction_ratio=0.25, phase=phase)

    # Optional off-axis planet via phase ramp
    A_1 = A

    # Shift PSF center by adding a phase ramp at the entrance pupil.
    dy, dx = float(phase_shift[0]), float(phase_shift[1])
    if phase_shift_unit == "array":
        dy_arr, dx_arr = dy, dx
    elif phase_shift_unit == "fpm":
        scale = nsamp / float(fpm_sam)
        dy_arr, dx_arr = dy * scale, dx * scale
    else:
        raise ValueError(f"Unknown phase_shift_unit: {phase_shift_unit!r}. Use 'array' or 'fpm'.")

    if dy_arr != 0.0 or dx_arr != 0.0:
        # `phi_ramp` expects (x, y) while `phase_shift` is (dy, dx).
        # Divide by nsamp: nsamp pixels per lambda/D in the focal plane.
        A_1 = A_1 * np.exp(1j * phi_ramp(A_1, dx_arr / nsamp, dy_arr / nsamp))

    if offset is not None:
        scale = round_to_even(dim / wavelength) / dim
        A_1 = A_1 * np.exp(1j * phi_ramp(A_1, offset[0] * scale, offset[1] * scale))

    # Lyot stop sizing presets
    lyot_stop_para = {
        "small_cross":  [0.98, 1.2],
        "large_cross":  [0.95, 1.3],
        "zernike_cross": [0.98, 1.05],
        "None":         [1.0, 1.2],
        "None_1_1":     [1.0, 1.0],
        "None_1_1.1":   [1.0, 1.1],
        "None_1_1.3":   [1.0, 1.3],
        "None_1_1.4":   [1.0, 1.4],
        "None_1_1.5":   [1.0, 1.5],
    }

    N_vor = make_lyot_stop(
        dim,
        fraction=lyot_stop_para[lyot_stop][0] if lyot_stop in lyot_stop_para else 1.0,
        sec_fraction=lyot_stop_para[lyot_stop][1] if lyot_stop in lyot_stop_para else 1.2,
        obstruction=obstruction,
        nsamp=nsamp,
    )

    # Lyot normalization (compute sum without phase mask)
    save = False
    lyot_sum = None
    shift = bool(coro_shift)
    print(shift)
    if get_lyot:
        m_vor_ones = np.ones_like(N_vor)
        lyot_sum = coro_high_sam(A, m_vor_ones, N_vor, dim, nsamp, fpm_sam,
                                  charge=charge, lyot_sum=None, shift=shift, save=False, get_lyot=True, ghost=ghost)
        save = True

    # Build focal-plane phase mask
    local_opt_params = {"rotate": rotate, "roddier_d": roddier_d, "vortex_center": vortex_center}
    if mask_opt_params:
        local_opt_params.update(mask_opt_params)

    m_vor = focal_mask_new(
        name, dim, nsamp, fpm_sam, wavelength, charge=charge,
        opt_params=local_opt_params,
        greyscale=greyscale if isinstance(greyscale, int) else (8 if greyscale else None),
        greyscale_method=greyscale_method,
        cal_factor=cal_factor, fill_factor=fill_factor, noise_level=noise_level,
    )

    # Final propagation (returns final field + resolution energy)
    final_field, _res_ene = coro_high_sam(
        A_1, m_vor, N_vor, dim, nsamp, fpm_sam, charge = charge,
        lyot_sum=lyot_sum, shift=shift, save=save, ghost = ghost
    )

    # Throughput mode: compare off-axis energy to no-mask energy
    if not get_lyot:
        m_vor_ones = np.ones_like(m_vor)
        if fill_factor:
            m_vor_ones[::10, :] = 0.0
            m_vor_ones[:, ::10] = 0.0

        if offset is not None:
            N_vor_nor = make_lyot_stop(dim,
                                       fraction=lyot_stop_para["None_1_1"][0],
                                       sec_fraction=lyot_stop_para["None_1_1"][1],
                                       obstruction=obstruction,
                                       nsamp=nsamp)
            wo_mask = coro_high_sam(A, m_vor_ones, N_vor_nor, dim, nsamp, fpm_sam,
                                    charge=charge,
                                    lyot_sum=None, shift=shift, save=False)[0]
            wo_mask_sum = np.sum(np.abs(wo_mask) ** 2)
            throughput = np.sum(np.abs(final_field) ** 2) / wo_mask_sum
            np.savetxt("throughput.txt", [throughput])
            return float(throughput)
        else:
            wo_mask_peak = np.max(
                np.abs(
                    coro_high_sam(A_1, m_vor_ones, N_vor, dim, nsamp, fpm_sam,
                                  charge=charge, lyot_sum=None, shift=shift, save=False)[0]
                ) ** 2
            )
            norm_final = (np.abs(final_field) ** 2) / wo_mask_peak
            wfits([norm_final], f"final_focal_plane_coro_mask={name}_{charge}_sample={nsamp}_fpm_sam={fpm_sam}_obstruction_{obstruction}.fits")
            return norm_final

    return final_field


def combine_psf_with_ghost(
    dim: int,
    name: str,
    wavelength: float,
    nsamp: int,
    fpm_sam: int,
    charge: Optional[int] = None,
    binary: bool = False,
    obstruction: bool = True,
    spider: bool = False,
    lyot_stop: Optional[str] = None,
    offset: Optional[Sequence[float]] = None,
    rotate: bool = False,
    greyscale: bool | int = True,
    roddier_d: float = 1.06,
    cal_factor: float = 1.0,
    phase: Optional[np.ndarray] = None,
    noise_level: Optional[float] = None,
    fill_factor: bool = False,
    phase_shift: Sequence[float] = (0.0, 0.0),
    phase_shift_unit: str = "array",
    # Ghost-specific options
    ghost_scale: float = 1.0,
    normalize_by_main_peak: bool = True,
    save: bool = True,
    out_name: str = "psf_combined.fits",
) -> np.ndarray:
    """
    Build a combined PSF intensity image that adds the coronagraphic (main) PSF
    and the refractive ghost PSF.

    This function calls `simple_coro` to get the main coronagraphic complex focal
    field and `ghost_im` to get the ghost complex focal field. It then forms
    intensities, optionally scales the ghost intensity, sums them, optionally
    normalizes by the main PSF peak, and saves the result using `wfits`.

    Parameters
    ----------
    (Most arguments mirror `simple_coro` and `ghost_im`.)
    ghost_scale : float
        Multiplicative factor applied to the ghost intensity (not amplitude).
        Use 1.0 to keep the ghost as implemented in `ghost_im` (which already
        applies an internal amplitude scaling). Values <1 reduce the ghost.
    normalize_by_main_peak : bool
        If True (default), divide the combined image by the peak of the main
        PSF (so the on-axis peak of the main PSF becomes 1.0). If False, the
        raw intensities are returned.
    save : bool
        If True, write the combined PSF to `out_name` via `wfits`.
    out_name : str
        Filename for the saved FITS.

    Returns
    -------
    np.ndarray
        2D array with the combined intensity image.
    """
    # Request complex focal fields from existing helpers. We force get_lyot=True
    # to obtain the complex final field from `simple_coro`.
    main_field = simple_coro(
        dim=dim,
        name=name,
        wavelength=wavelength,
        nsamp=nsamp,
        fpm_sam=fpm_sam,
        lyot_stop=lyot_stop,
        offset=offset,
        charge=charge,
        binary=binary,
        obstruction=obstruction,
        spider=spider,
        get_lyot=True,
        rotate=rotate,
        greyscale=greyscale,
        roddier_d=roddier_d,
        cal_factor=cal_factor,
        phase=phase,
        noise_level=noise_level,
        fill_factor=fill_factor,
        ghost=False,
        phase_shift=phase_shift,
        phase_shift_unit=phase_shift_unit,
    )

    ghost_field = ghost_im(
        dim=dim,
        nsamp=nsamp,
        fpm_sam=fpm_sam,
        wavelength=wavelength,
        charge=charge,
        lyot_sum=None,
        phase=phase,
        shift=True,
        save=False,
        obstruction=obstruction,
        spider=spider,
    )

    # Compute intensities
    I_main = np.abs(main_field) ** 2
    I_ghost = np.abs(ghost_field) ** 2

    # Apply ghost scale (interpreted on intensity)
    I_ghost = I_ghost * ghost_scale

    combined = I_main + I_ghost

    if normalize_by_main_peak:
        peak = np.max(I_main)
        if peak > 0:
            combined = combined / peak

    if save:
        # wrap in list to preserve compatibility with wfits calls elsewhere
        wfits([combined], out_name)

    return combined


def res_ene(
    dim: int,
    name: str,
    wavelength: float,
    nsamp: int,
    fpm_sam: int,
    lyot_stop: Optional[str] = None,
    offset: Optional[Sequence[float]] = None,
    charge: Optional[int] = None,
    binary: bool = False,
    obstruction: bool = True,
    get_lyot: bool = True,
    rotate: bool = False,
    greyscale: bool | int = True,
    roddier_d: float = 1.06,
    phase_shift: Sequence[float] = (0.0, 0.0),
    phase_shift_unit: str = "array",
    coro_shift: bool = True,
) -> float:
    """
    Compute the resolution/energy metric (Lyot-normalized) for a given setup.
    """
    # Pupil
    if binary:
        A = make_entrance_pupil(dim)
    else:
        A = make_aperture(round_to_even(dim / wavelength), 8, obstruction=obstruction, obstruction_ratio=0.25)

    # Optional PSF-center shift via pupil ramp
    A_1 = A
    dy, dx = float(phase_shift[0]), float(phase_shift[1])
    if phase_shift_unit == "array":
        dy_arr, dx_arr = dy, dx
    elif phase_shift_unit == "fpm":
        scale = nsamp / float(fpm_sam)
        dy_arr, dx_arr = dy * scale, dx * scale
    else:
        raise ValueError(f"Unknown phase_shift_unit: {phase_shift_unit!r}. Use 'array' or 'fpm'.")

    if dy_arr != 0.0 or dx_arr != 0.0:
        A_1 = A_1 * np.exp(1j * phi_ramp(A_1, dx_arr / nsamp, dy_arr / nsamp))

    # Optional planet offset
    if offset is not None:
        scale = round_to_even(dim / wavelength) / dim
        A_1 = A_1 * np.exp(1j * phi_ramp(A_1, offset[0] * scale, offset[1] * scale))

    lyot_stop_para = {"small_cross": [0.98, 1.2], "large_cross": [0.95, 1.3], "zernike_cross": [0.98, 1.05],
                      "None": [1.0, 1.2], "None_1_1": [1.0, 1.0], "None_1_1.1": [1.0, 1.1],
                      "None_1_1.3": [1.0, 1.3], "None_1_1.4": [1.0, 1.4], "None_1_1.5": [1.0, 1.5]}

    N_vor = make_lyot_stop(round_to_even(dim / wavelength),
                           fraction=lyot_stop_para[lyot_stop][0] if lyot_stop in lyot_stop_para else 1.0,
                           sec_fraction=lyot_stop_para[lyot_stop][1] if lyot_stop in lyot_stop_para else 1.2,
                           obstruction=obstruction, nsamp=nsamp)

    # Lyot normalization (no phase mask)
    m_vor_ones = np.ones_like(N_vor)
    lyot_sum = coro_high_sam(A, m_vor_ones, N_vor, dim, nsamp, fpm_sam,
                             charge=charge,
                             lyot_sum=None, shift=coro_shift, save=False, get_lyot=True)

    # Chosen focal mask (grayscale default: if True, use 8 levels for backward compat)
    m_vor = focal_mask_new(name, dim, nsamp, fpm_sam, wavelength, charge=charge,
                           opt_params={"rotate": rotate, "roddier_d": roddier_d},
                           greyscale=greyscale if isinstance(greyscale, int) else (8 if greyscale else None))

    # Compute resolution energy
    _final, metric = coro_high_sam(A_1, m_vor, N_vor, dim, nsamp, fpm_sam,
                                   charge=charge,
                                   lyot_sum=lyot_sum, shift=coro_shift, save=True)
    return float(metric)
