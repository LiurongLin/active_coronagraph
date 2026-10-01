"""Minimal coronagraph simulation example.

This example uses a small FQPM coronagraph setup so it runs quickly. It computes
an annular normalized-intensity curve from 0 to 15 lambda/D and can optionally
save a diagnostic figure.

Run from an installed package or from the repository root:

    python examples/minimal_simulation.py --save-figure
"""

from __future__ import annotations

import argparse
import contextlib
import io
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import active_coronagraph
from active_coronagraph.coronagraphs_polished import simple_coro


@contextlib.contextmanager
def working_directory(path: Path):
    old = Path.cwd()
    path.mkdir(parents=True, exist_ok=True)
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(old)


def annular_intensity_curve(
    image: np.ndarray,
    nsamp: int,
    max_radius_ld: float,
    step_ld: float = 1.0,
) -> np.ndarray:
    """Median image value in annuli expressed in lambda/D."""
    y, x = np.indices(image.shape, dtype=float)
    center = (np.array(image.shape) - 1.0) / 2.0
    radius_pixels = np.hypot(x - center[1], y - center[0])
    radius_ld = radius_pixels / float(nsamp)

    edges = np.arange(0.0, max_radius_ld + step_ld, step_ld)
    radii = 0.5 * (edges[:-1] + edges[1:])
    medians = []
    for inner, outer in zip(edges[:-1], edges[1:]):
        mask = (radius_ld >= inner) & (radius_ld < outer)
        medians.append(float(np.median(image[mask])))
    return np.column_stack([radii, medians])


def run_simulation(output_dir: Path, save_figure: bool) -> np.ndarray:
    # Smallest practical optical system for a quick first run:
    # unobstructed entrance pupil, no undersized Lyot stop, and an FQPM mask.
    dim = 32
    nsamp = 2
    fpm_sam = 2

    with working_directory(output_dir):
        # simple_coro writes a FITS product in the current directory. Suppress
        # legacy progress prints so the example output stays focused.
        with contextlib.redirect_stdout(io.StringIO()):
            normalized_psf = simple_coro(
                dim=dim,
                name="FQPM",
                wavelength=1.0,
                nsamp=nsamp,
                fpm_sam=fpm_sam,
                lyot_stop="None_1_1",
                charge=None,
                binary=False,
                obstruction=False,
                get_lyot=False,
                greyscale=None,
            )

        intensity_curve = annular_intensity_curve(normalized_psf, nsamp, 15.0)
        np.savetxt(
            "normalized_intensity_curve_0_to_15_lambda_over_D.txt",
            intensity_curve,
            header="radius_lambda_over_D median_normalized_intensity",
        )

        if save_figure:
            fig, ax = plt.subplots(figsize=(5, 4))
            ax.plot(intensity_curve[:, 0], intensity_curve[:, 1], marker="o")
            ax.set_yscale("log")
            ax.set_title("FQPM normalized intensity curve")
            ax.set_xlabel("radius [lambda/D]")
            ax.set_ylabel("median normalized intensity")
            fig.tight_layout()
            fig.savefig("minimal_fqpm_intensity_curve.png", dpi=150)
            plt.close(fig)

    return intensity_curve


def main() -> None:
    parser = argparse.ArgumentParser(description="Run a minimal active-coronagraph simulation.")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("output/minimal_example"),
        help="Directory for generated FITS and optional figure outputs.",
    )
    parser.add_argument(
        "--save-figure",
        action="store_true",
        help="Save output/minimal_example/minimal_fqpm_intensity_curve.png.",
    )
    args = parser.parse_args()

    intensity_curve = run_simulation(args.output_dir, args.save_figure)

    print(f"active_coronagraph {active_coronagraph.__version__}")
    print("Optical system: unobstructed pupil, FQPM mask, None_1_1 Lyot stop")
    print("Normalized intensity curve from 0 to 15 lambda/D:")
    for radius_ld, median_intensity in intensity_curve:
        print(f"  {radius_ld:4.1f} lambda/D  {median_intensity:.6e}")
    print(f"Generated files are in: {args.output_dir}")


if __name__ == "__main__":
    main()
