from __future__ import annotations

from pathlib import Path

from astropy.io import fits

try:
    from tqdm import tqdm
except ImportError:  # pragma: no cover - fallback for environments without tqdm
    tqdm = None


def load_phase_map(path, dim, fpm_sam):
    """Load a user-provided 2D FITS phase map in radians on the FPM grid."""
    if path is None:
        raise ValueError("new_mask requires a phase_map_fits path.")
    phase = fits.getdata(Path(path).expanduser().resolve())
    expected = (int(dim * fpm_sam), int(dim * fpm_sam))
    if phase.shape != expected:
        raise ValueError(f"new_mask phase map shape {phase.shape} must be {expected}.")
    if phase.dtype.kind not in "biuf" or not np.all(np.isfinite(phase)):
        raise ValueError("new_mask phase map must contain finite real numeric values in radians.")
    return np.array(phase, dtype=float, copy=True)


def _mask_job(mask_name, charge_val, phase_map_fits, dim, fpm_sam):
    """Describe one analytic or user-provided mask for a throughput workflow."""
    phase = None
    params = None
    output_name = mask_name
    if mask_name == "new_mask":
        phase = load_phase_map(phase_map_fits, dim, fpm_sam)
        path = Path(phase_map_fits).expanduser().resolve()
        params = {"phase_map_fits": str(path)}
        output_name = f"new_mask_{path.stem}"
        charge_val = None
    return {
        "runtime_mask_name": mask_name,
        "output_mask_name": output_name,
        "charge": charge_val,
        "mask_opt_params": params,
        "phase_preview": phase,
    }


def directional_throughput_mono(
    mask_name,
    charge_val,
    wl=1.0,
    max_ld=10.0,
    step_ld=0.1,
    dim=100,
    nsamp=10,
    fpm_sam=10,
    lyot_stop="large_cross",
    obstruction=False,
    spider=False,
    greyscale=None,
    cal_factor=1.0,
    folder_name="",
    noise_level=0,
    ghost=False,
    phase_shift=(0.0, 0.0),
    phase_shift_unit="array",
    coro_shift=True,
    phase_map_fits=None,
):
    """Measure monochromatic throughput from -max_ld..+max_ld along 3 directions."""
    from main_functions_polished import run, _root_folder

    repo_root = Path(__file__).resolve().parent


    radii = np.arange(-max_ld, max_ld + 0.5 * step_ld, step_ld)
    inv_sqrt2 = 1.0 / np.sqrt(2.0)
    # IMPORTANT: `simple_coro(offset=...)` expects (x, y), not (dy, dx).
    offsets = {
        "horizontal": [(float(r), 0.0) for r in radii],                    # x sweep
        "vertical": [(0.0, float(r)) for r in radii],                      # y sweep
        "diagonal": [(float(r * inv_sqrt2), float(r * inv_sqrt2)) for r in radii],
    }
    colors = {
        "horizontal": "tab:blue",
        "vertical": "tab:orange",
        "diagonal": "tab:green",
    }
    out_dir = _root_folder(obstruction, noise_level, folder_name) / "throughput_directional_mono"
    out_dir.mkdir(parents=True, exist_ok=True)

    out_files = []
    jobs = [_mask_job(mask_name, charge_val, phase_map_fits, dim, fpm_sam)]
    total_steps = len(jobs) * sum(len(pairs) for pairs in offsets.values())

    def iter_sweeps():
        for job in jobs:
            for direction, pairs in offsets.items():
                direction_total = len(pairs)
                for step_idx, (ox, oy) in enumerate(pairs, start=1):
                    yield job, direction, direction_total, step_idx, ox, oy

    progress = None
    if tqdm is not None:
        progress = tqdm(total=total_steps, desc="tp_directional_mono", unit="eval")
    else:
        print(
            f"[tp_directional_mono] Starting {total_steps} throughput evaluations "
            f"across {len(jobs)} mask(s). tqdm not available.",
            flush=True,
        )

    results_by_job = {
        job["output_mask_name"]: {key: [] for key in offsets}
        for job in jobs
    }

    for job, direction, direction_total, step_idx, ox, oy in iter_sweeps():
        tp = run(
            dim,
            job["runtime_mask_name"],
            wl,
            nsamp,
            fpm_sam,
            charge=job["charge"],
            lyot_stop=lyot_stop,
            binary=False,
            obstruction=obstruction,
            spider=spider,
            get_lyot=False,
            rotate=False,
            greyscale=greyscale,
            offset=(ox, oy),
            cal_factor=cal_factor,
            phase=None,
            folder_name=folder_name,
            noise_level=noise_level,
            fill_factor=False,
            ghost=ghost,
            phase_shift=phase_shift,
            phase_shift_unit=phase_shift_unit,
            coro_shift=coro_shift,
            mask_opt_params=job["mask_opt_params"],
        )
        results_by_job[job["output_mask_name"]][direction].append(float(tp))
        if progress is not None:
            progress.set_postfix_str(
                f"mask={job['output_mask_name']} dir={direction} {step_idx}/{direction_total}"
            )
            progress.update(1)
        else:
            print(
                f"[tp_directional_mono] mask={job['output_mask_name']} "
                f"direction={direction} step={step_idx}/{direction_total}",
                flush=True,
            )

    if progress is not None:
        progress.close()
    else:
        print("[tp_directional_mono] Completed.", flush=True)

    for job in jobs:
        results = results_by_job[job["output_mask_name"]]

        charge_tag = "none" if job["charge"] is None else str(job["charge"])
        out_file = out_dir / (
            f"{job['output_mask_name']}_{charge_tag}_wl_{wl}_tp_0_to_{max_ld}_step_{step_ld}.txt"
        )

        table = np.column_stack([
            radii,
            np.asarray(results["horizontal"]),
            np.asarray(results["vertical"]),
            np.asarray(results["diagonal"]),
        ])
        header = (
            "sep_lambda_over_D throughput_horizontal "
            "throughput_vertical throughput_diagonal"
        )
        np.savetxt(out_file, table, header=header)

        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(radii, results["horizontal"], label="Horizontal", color=colors["horizontal"])
        ax.plot(radii, results["vertical"], label="Vertical", color=colors["vertical"])
        ax.plot(radii, results["diagonal"], label="Diagonal", color=colors["diagonal"])
        ax.set_xlabel("Separation (lambda/D), signed")
        ax.set_ylabel("Throughput")
        ax.set_title(f"Monochromatic throughput: {job['output_mask_name']} wl={wl}")
        ax.grid(True, alpha=0.3)
        ax.legend()
        fig.tight_layout()
        fig.savefig(
            out_dir / f"{job['output_mask_name']}_{charge_tag}_wl_{wl}_tp_0_to_{max_ld}.png",
            dpi=200,
        )
        plt.close(fig)

        phase = job["phase_preview"]
        if phase is None:
            phase = np.zeros((int(dim * fpm_sam), int(dim * fpm_sam)))

        n = phase.shape[0]
        c = (n - 1) / 2.0
        fig2, ax2 = plt.subplots(figsize=(8, 7))
        im = ax2.imshow(phase, cmap="gray", origin="upper", vmin=0.0, vmax=2 * np.pi)
        cbar = plt.colorbar(im, ax=ax2)
        cbar.set_label("Phase [rad]")

        for key, offs in offsets.items():
            xs = [c + ox * fpm_sam for ox, _ in offs]
            ys = [c + oy * fpm_sam for _, oy in offs]
            ax2.plot(xs, ys, color=colors[key], lw=1.8, label=f"{key} trace")

        ax2.scatter([c], [c], c="red", s=25, marker="x", label="center")
        ax2.set_title(f"Offset traces on phase mask: {job['output_mask_name']}")
        ax2.set_xlabel("x [pixel]")
        ax2.set_ylabel("y [pixel]")
        ax2.legend(loc="upper right", fontsize=8, framealpha=0.85)
        fig2.tight_layout()
        fig2.savefig(
            out_dir / f"{job['output_mask_name']}_{charge_tag}_wl_{wl}_traces_pm{max_ld}.png",
            dpi=220,
        )
        plt.close(fig2)
        out_files.append(out_file)

    if len(out_files) == 1:
        return out_files[0]
    return out_files


def throughput_map_mono(
    mask_name,
    charge_val,
    wl=1.0,
    min_ld=0.0,
    max_ld=20.0,
    step_ld=1.0,
    dim=100,
    nsamp=10,
    fpm_sam=10,
    lyot_stop="large_cross",
    obstruction=False,
    spider=False,
    greyscale=None,
    cal_factor=1.0,
    folder_name="",
    noise_level=0,
    ghost=False,
    phase_shift=(0.0, 0.0),
    phase_shift_unit="array",
    coro_shift=True,
    phase_map_fits=None,
):
    """Measure monochromatic throughput on an x/y lambda/D grid."""
    from main_functions_polished import run, _root_folder

    repo_root = Path(__file__).resolve().parent


    x_coords = np.arange(min_ld, max_ld + 0.5 * step_ld, step_ld, dtype=float)
    y_coords = np.arange(min_ld, max_ld + 0.5 * step_ld, step_ld, dtype=float)
    out_dir = _root_folder(obstruction, noise_level, folder_name) / "throughput_map_mono"
    out_dir.mkdir(parents=True, exist_ok=True)

    jobs = [_mask_job(mask_name, charge_val, phase_map_fits, dim, fpm_sam)]
    total_steps = len(jobs) * len(x_coords) * len(y_coords)
    progress = None
    if tqdm is not None:
        progress = tqdm(total=total_steps, desc="tp_map_mono", unit="eval")
    else:
        print(
            f"[tp_map_mono] Starting {total_steps} throughput evaluations "
            f"across {len(jobs)} mask(s). tqdm not available.",
            flush=True,
        )

    out_files = []
    for job in jobs:
        throughput_map = np.zeros((len(y_coords), len(x_coords)), dtype=float)
        for iy, oy in enumerate(y_coords, start=1):
            for ix, ox in enumerate(x_coords, start=1):
                tp = run(
                    dim,
                    job["runtime_mask_name"],
                    wl,
                    nsamp,
                    fpm_sam,
                    charge=job["charge"],
                    lyot_stop=lyot_stop,
                    binary=False,
                    obstruction=obstruction,
                    spider=spider,
                    get_lyot=False,
                    rotate=False,
                    greyscale=greyscale,
                    offset=(float(ox), float(oy)),
                    cal_factor=cal_factor,
                    phase=None,
                    folder_name=folder_name,
                    noise_level=noise_level,
                    fill_factor=False,
                    ghost=ghost,
                    phase_shift=phase_shift,
                    phase_shift_unit=phase_shift_unit,
                    coro_shift=coro_shift,
                    mask_opt_params=job["mask_opt_params"],
                )
                throughput_map[iy - 1, ix - 1] = float(tp)
                if progress is not None:
                    progress.set_postfix_str(
                        f"mask={job['output_mask_name']} y={iy}/{len(y_coords)} x={ix}/{len(x_coords)}"
                    )
                    progress.update(1)
                else:
                    print(
                        f"[tp_map_mono] mask={job['output_mask_name']} "
                        f"y={iy}/{len(y_coords)} x={ix}/{len(x_coords)}",
                        flush=True,
                    )

        charge_tag = "none" if job["charge"] is None else str(job["charge"])
        out_file = out_dir / (
            f"{job['output_mask_name']}_{charge_tag}_wl_{wl}_tp_map_{min_ld}_to_{max_ld}_step_{step_ld}.txt"
        )
        header = (
            f"rows=y_lambda_over_D cols=x_lambda_over_D "
            f"x_count={len(x_coords)} y_count={len(y_coords)} "
            f"x_min={x_coords[0]:g} x_max={x_coords[-1]:g} "
            f"y_min={y_coords[0]:g} y_max={y_coords[-1]:g} step={step_ld:g}"
        )
        np.savetxt(out_file, throughput_map, header=header)

        extent = [
            x_coords[0] - 0.5 * step_ld,
            x_coords[-1] + 0.5 * step_ld,
            y_coords[-1] + 0.5 * step_ld,
            y_coords[0] - 0.5 * step_ld,
        ]
        fig, ax = plt.subplots(figsize=(7, 6))
        im = ax.imshow(throughput_map, origin="upper", extent=extent, aspect="equal")
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label("Throughput")
        ax.set_xlabel("x offset [lambda/D]")
        ax.set_ylabel("y offset [lambda/D]")
        ax.set_title(
            f"Monochromatic throughput map: {job['output_mask_name']} "
            f"wl={wl} ({len(x_coords)}x{len(y_coords)})"
        )
        fig.tight_layout()
        fig.savefig(
            out_dir / f"{job['output_mask_name']}_{charge_tag}_wl_{wl}_tp_map_{min_ld}_to_{max_ld}.png",
            dpi=220,
        )
        plt.close(fig)

        phase = job["phase_preview"]
        if phase is None:
            phase = np.zeros((int(dim * fpm_sam), int(dim * fpm_sam)))

        n = phase.shape[0]
        c = (n - 1) / 2.0
        xx = c + x_coords * fpm_sam
        yy = c + y_coords * fpm_sam
        fig2, ax2 = plt.subplots(figsize=(8, 8))
        im2 = ax2.imshow(phase, cmap="gray", origin="upper")
        cbar2 = plt.colorbar(im2, ax=ax2)
        cbar2.set_label("Phase [rad]")
        for xpix in xx:
            ax2.axvline(xpix, color="tab:red", alpha=0.18, lw=0.6)
        for ypix in yy:
            ax2.axhline(ypix, color="tab:blue", alpha=0.18, lw=0.6)
        grid_x, grid_y = np.meshgrid(xx, yy)
        ax2.scatter(grid_x.ravel(), grid_y.ravel(), s=4, c="yellow", alpha=0.35)
        ax2.scatter([c], [c], c="cyan", s=30, marker="x", label="mask center")
        ax2.set_title(
            f"Throughput map sampling on mask: {job['output_mask_name']} "
            f"[{min_ld:g}, {max_ld:g}] step {step_ld:g}"
        )
        ax2.set_xlabel("x [pixel]")
        ax2.set_ylabel("y [pixel]")
        ax2.legend(loc="upper right", fontsize=8, framealpha=0.85)
        fig2.tight_layout()
        fig2.savefig(
            out_dir / f"{job['output_mask_name']}_{charge_tag}_wl_{wl}_tp_map_{min_ld}_to_{max_ld}_on_mask.png",
            dpi=220,
        )
        plt.close(fig2)
        out_files.append(out_file)

    if progress is not None:
        progress.close()
    else:
        print("[tp_map_mono] Completed.", flush=True)

    if len(out_files) == 1:
        return out_files[0]
    return out_files


def throughput_1d_mono(
    mask_name,
    charge_val,
    wl=1.0,
    max_ld=20.0,
    step_ld=0.1,
    dim=100,
    nsamp=10,
    fpm_sam=10,
    lyot_stop="large_cross",
    obstruction=False,
    spider=False,
    greyscale=None,
    cal_factor=1.0,
    folder_name="",
    noise_level=0,
    ghost=False,
    phase_shift=(0.0, 0.0),
    phase_shift_unit="array",
    coro_shift=True,
    phase_map_fits=None,
):
    """Measure monochromatic throughput from 0..max_ld along +x."""
    from main_functions_polished import run, _root_folder

    repo_root = Path(__file__).resolve().parent


    x_coords = np.arange(0.0, max_ld + 0.5 * step_ld, step_ld, dtype=float)
    out_dir = _root_folder(obstruction, noise_level, folder_name) / "throughput_1d_mono"
    out_dir.mkdir(parents=True, exist_ok=True)

    jobs = [_mask_job(mask_name, charge_val, phase_map_fits, dim, fpm_sam)]
    total_steps = len(jobs) * len(x_coords)
    progress = None
    if tqdm is not None:
        progress = tqdm(total=total_steps, desc="tp_1d_mono", unit="eval")
    else:
        print(
            f"[tp_1d_mono] Starting {total_steps} throughput evaluations "
            f"across {len(jobs)} mask(s). tqdm not available.",
            flush=True,
        )

    out_files = []
    for job in jobs:
        throughput = np.zeros(len(x_coords), dtype=float)
        for ix, ox in enumerate(x_coords, start=1):
            tp = run(
                dim,
                job["runtime_mask_name"],
                wl,
                nsamp,
                fpm_sam,
                charge=job["charge"],
                lyot_stop=lyot_stop,
                binary=False,
                obstruction=obstruction,
                spider=spider,
                get_lyot=False,
                rotate=False,
                greyscale=greyscale,
                offset=(float(ox), 0.0),
                cal_factor=cal_factor,
                phase=None,
                folder_name=folder_name,
                noise_level=noise_level,
                fill_factor=False,
                ghost=ghost,
                phase_shift=phase_shift,
                phase_shift_unit=phase_shift_unit,
                coro_shift=coro_shift,
                mask_opt_params=job["mask_opt_params"],
            )
            throughput[ix - 1] = float(tp)
            if progress is not None:
                progress.set_postfix_str(
                    f"mask={job['output_mask_name']} x={ix}/{len(x_coords)}"
                )
                progress.update(1)
            else:
                print(
                    f"[tp_1d_mono] mask={job['output_mask_name']} x={ix}/{len(x_coords)}",
                    flush=True,
                )

        charge_tag = "none" if job["charge"] is None else str(job["charge"])
        out_file = out_dir / (
            f"{job['output_mask_name']}_{charge_tag}_wl_{wl}_tp_1d_0_to_{max_ld}_step_{step_ld}.txt"
        )
        table = np.column_stack([x_coords, throughput])
        header = "x_lambda_over_D throughput"
        np.savetxt(out_file, table, header=header)

        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(x_coords, throughput, color="tab:blue")
        ax.set_xlabel("x offset [lambda/D]")
        ax.set_ylabel("Throughput")
        ax.set_title(f"Monochromatic 1D throughput: {job['output_mask_name']} wl={wl}")
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig(
            out_dir / f"{job['output_mask_name']}_{charge_tag}_wl_{wl}_tp_1d_0_to_{max_ld}.png",
            dpi=220,
        )
        plt.close(fig)
        out_files.append(out_file)

    if progress is not None:
        progress.close()
    else:
        print("[tp_1d_mono] Completed.", flush=True)

    if len(out_files) == 1:
        return out_files[0]
    return out_files
