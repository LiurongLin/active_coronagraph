# main_polished bundle

This folder contains the project-local code required by `main_polished.py`.

Included Python modules:
- `main_polished.py`
- `basic.py`
- `make_plot.py`
- `coronagraphs.py`
- `coronagraphs_polished.py`
- `main_functions_polished.py`
- `new_mask.py`
- `phase_masks.py`

Included runtime data directories referenced by the code:
- `binary_mask/`
- `phase_screen/`

Notes before running elsewhere:
- `main_polished.py` contains a hardcoded output path:
  - `pixel_noise_mean_output_dir = Path("/media/liurong/My Passport/PLACID")`
  Update that path for the target machine if you need the `pixel_noise` workflow.
- Some workflows also expect generated outputs under directories such as
  `ideal_coro_2rd_mirror_*`. Those are not included here because they are outputs,
  not source dependencies.
- The `add_phase_screen` branch in `main_polished.py` also references
  `2)seeing=1,vmag=8,ZA=30,lag=2/TROIA_phase_screens_new.fits`, which is not part
  of this bundle.

Basic run example:

```bash
active-coronagraph --function main_func --name dual_zone --charge None
```

This uses the default monochromatic wavelength `1.0`, default sampling
`nsamp=10`, default FPM sampling `fpm_sam=10`, and writes FITS, text, and PNG
products under `ideal_coro_2rd_mirror_False_noise_0/`.

Run an FQPM simulation:

```bash
active-coronagraph --function main_func --name FQPM --charge None
```

This changes the focal-plane phase mask to the four-quadrant phase mask. The
output tree and file naming follow the same pattern as the quick-start command,
with `FQPM_None` in the run directory and filenames.

Run a vortex coronagraph simulation:

```bash
active-coronagraph --function main_func --name vortex --charge 2
```

This uses the vortex mask with charge `2`. The final focal-plane FITS product is
written in the generated vortex run directory.

Change the vortex charge:

```bash
active-coronagraph --function main_func --name vortex --charge 4
```

This keeps the same workflow and defaults as the previous vortex example, but
uses charge `4`.

Run a custom monochromatic wavelength grid:

```bash
active-coronagraph --function main_func --name vortex --charge 2 --wavelength 0.95 1.0 1.05
```

This runs one simulation per listed wavelength ratio and writes one run
directory per wavelength.

Run a broadband vortex PSF combination:

```bash
active-coronagraph --function broadband_combine --name vortex --charge 2
```

This runs wavelength ratios `[0.9, 0.95, 1.0, 1.05, 1.1]`, combines the saved
monochromatic PSFs with endpoint-half trapezoid weights, writes a `_bb.fits`
file, and writes a broadband radial text curve ending in `_bb.txt`.

Use the provided phase-screen cubes:

```bash
active-coronagraph --function compare_phase_screens --name FQPM --charge None --phase-screen-folder data/phase_screen
```

This reads every `*.fits` cube in `data/phase_screen/`, runs each realization,
averages each case, and writes a phase-screen comparison plot.

Change the output directory:

```bash
ACTIVE_CORONAGRAPH_OUTPUT_DIR=output/cli_runs \
active-coronagraph --function main_func --name dual_zone --charge None
```

This writes the generated `ideal_coro_2rd_mirror_False_noise_0/` tree under
`output/cli_runs/` instead of the current working directory.

Use all provided binary mask files for a 1D throughput sweep:

```bash
active-coronagraph --function tp_1d_mono --name binary_mask --charge None --tp-1d-max-ld 5 --tp-1d-step-ld 0.5
```

This reads all FITS masks in `data/binary_mask/`, computes throughput along the
positive x direction from `0` to `5 lambda/D`, and writes text and PNG products
under `ideal_coro_2rd_mirror_False_noise_0/throughput_1d_mono/`.

Use one provided binary mask file:

```bash
active-coronagraph --function tp_1d_mono --name binary_mask --charge None --phase-map-fits data/binary_mask/Binary_vortex_sep0.5.fits --tp-1d-max-ld 5 --tp-1d-step-ld 0.5
```

This runs the same 1D throughput sweep, but only for
`Binary_vortex_sep0.5.fits`.

Generate a small 2D throughput map:

```bash
active-coronagraph --function tp_map_mono --name vortex --charge 2 --tp-map-min-ld 0 --tp-map-max-ld 3 --tp-map-step-ld 1
```

This evaluates vortex throughput on an x/y grid in `lambda/D` and writes a text
map plus diagnostic PNG plots.

Sampling and central-obstruction changes are not exposed by CLI flags. Use the
Python API for those settings.

## Using Configuration Files

The current CLI does not support configuration files. There is no `--config`
option.

Path configuration is handled with environment variables and selected CLI path
options:

| Setting | Meaning | Default |
| --- | --- | --- |
| `ACTIVE_CORONAGRAPH_OUTPUT_DIR` | Base directory for generated outputs | current working directory |
| `ACTIVE_CORONAGRAPH_DATA_DIR` | Runtime data root containing `binary_mask/` and `phase_screen/` | repository `data/` directory |
| `--data-dir` | Per-command runtime data root override | unset |
| `--phase-screen-folder` | Per-command phase-screen folder override | `data_dir()/phase_screen` |

Example:

```bash
ACTIVE_CORONAGRAPH_DATA_DIR=data \
ACTIVE_CORONAGRAPH_OUTPUT_DIR=output/cli_runs \
active-coronagraph --function main_func --name dual_zone --charge None
```

CLI arguments override only the specific values they expose. Other simulation
settings remain the defaults defined in `src/active_coronagraph/cli.py`.

## Input Data

Bundled runtime data is stored under `data/`.

`data/binary_mask/` contains FITS phase masks and one PNG overview image. The
CLI uses these FITS files only in the monochromatic throughput workflows
`tp_1d_mono`, `tp_map_mono`, and `tp_directional_mono` when `--name` is
`binary_mask`, `binary_vortex`, or `binary_masks`.

```bash
active-coronagraph --function tp_directional_mono --name binary_mask --charge None --phase-map-fits data/binary_mask/Binary_vortex_sep0.5_45posang_center_axis.fits
```

If `--phase-map-fits` is omitted for these binary-mask throughput workflows, the
code uses all `*.fits` files found in `data/binary_mask/`.

`data/phase_screen/` contains FITS cubes used by `compare_phase_screens`.

```bash
active-coronagraph --function compare_phase_screens --name FQPM --charge None --phase-screen-folder data/phase_screen
```

The supplied phase-screen cubes are FITS files. The code applies each cube slice
as an entrance-pupil phase screen in the phase-screen comparison workflow.

## Output Files

The CLI writes files into a generated directory tree. The root directory is:

```text
ideal_coro_2rd_mirror_<obstruction>[_noise_<noise_level>][<folder_name>]/
```

For normal CLI runs, `obstruction` defaults to `False` and `noise_level`
defaults to `0`, so the default root is:

```text
ideal_coro_2rd_mirror_False_noise_0/
```

Set `ACTIVE_CORONAGRAPH_OUTPUT_DIR` to move this whole tree somewhere else:

```bash
ACTIVE_CORONAGRAPH_OUTPUT_DIR=output/cli_runs \
active-coronagraph --function main_func --name dual_zone --charge None
```

That command writes to:

```text
output/cli_runs/ideal_coro_2rd_mirror_False_noise_0/
```

The next directory level identifies the actual simulation parameters. For the
quick-start command, it is:

```text
lyot_large_cross_dual_zone_None_lambda_1.0_fpm_sam=10_binary_False_obstruction_False/
```

Read that name as:

- `lyot_large_cross`: Lyot stop preset
- `dual_zone_None`: mask name and charge
- `lambda_1.0`: wavelength ratio
- `fpm_sam=10`: focal-plane mask sampling
- `binary_False`: analytic, not binary entrance pupil mode
- `obstruction_False`: no central obstruction

Inside each run directory, common files include:

| File pattern | Meaning |
| --- | --- |
| `m_vor_<name>_<charge>_nsamp=<nsamp>_fpm_sam=<fpm_sam>_rotate=<rotate>_greyscale=<greyscale>_gmethod=<method>_wl_<wavelength>_cal_<cal_factor>.fits` | Saved focal-plane phase mask in radians, written with a leading length-1 axis. |
| `final_focal_plane_coro_mask=<name>_<charge>_sample=<nsamp>_fpm_sam=<fpm_sam>_obstruction_<obstruction>.fits` | Normalized final focal-plane intensity image for a monochromatic run. |
| `final_focal_plane_coro_mask=<name>_<charge>_sample=<nsamp>_fpm_sam=<fpm_sam>_obstruction_<obstruction>_wl_<wl>.txt` | Radial median intensity curve from `annuli_mean()`. |
| `final_focal_plane_coro_mask=<name>_<charge>_sample=<nsamp>_fpm_sam=<fpm_sam>_obstruction_<obstruction>_wl_<wl>_pupil_<dim>_annuli_mean.png` | Plot of the radial curve. |
| `final_focal_plane_coro_mask=<name>_<charge>_sample=<nsamp>_fpm_sam=<fpm_sam>_obstruction_<obstruction>_bb.fits` | Weighted broadband final focal-plane image from `broadband_combine`. |
| `throughput.txt` | Scalar throughput for a single offset evaluation. Used internally by throughput sweeps. |

For most users, the first files to inspect are:

- the final focal-plane FITS file, if you want the 2D coronagraphic image
- the `_wl_<wl>.txt` curve, if you want the radial normalized-intensity values
- the `_annuli_mean.png` plot, if you want a quick visual check

Throughput helper workflows write summary products in subdirectories:

| Workflow | Output subdirectory | Main text output |
| --- | --- | --- |
| `tp_1d_mono` | `throughput_1d_mono/` | `<mask>_<charge>_wl_<wl>_tp_1d_0_to_<max_ld>_step_<step_ld>.txt` |
| `tp_map_mono` | `throughput_map_mono/` | `<mask>_<charge>_wl_<wl>_tp_map_<min_ld>_to_<max_ld>_step_<step_ld>.txt` |
| `tp_directional_mono` | `throughput_directional_mono/` | `<mask>_<charge>_wl_<wl>_tp_0_to_<max_ld>_step_<step_ld>.txt` |

These workflows also write diagnostic PNG plots.

## Interpreting the Main Result

The main FITS product is the coronagraphic final focal-plane intensity,
normalized by the peak intensity of the corresponding no-mask PSF. It is a
contrast-like normalized intensity image, not an absolute calibrated flux.

The radial text curve reports median normalized intensity in annuli. The radial
axis is expressed in `lambda/D`; annuli are sampled out to `15 lambda/D` by
`annuli_mean()`.

Throughput outputs measure the energy transmitted for an off-axis source
relative to a no-mask reference. Throughput offsets use `(x, y)` ordering in
`lambda/D`.

## Troubleshooting

`python3 -m venv .venv` fails on Ubuntu or Debian:

```bash
sudo apt install python3-venv
```

`active-coronagraph: command not found` usually means the package is not
installed in the active environment. Activate the virtual environment and run:

```bash
python3 -m pip install -e .
active-coronagraph --help
```

`ModuleNotFoundError: active_coronagraph` means Python is not using the installed
environment. Re-activate `.venv`, or run from the repository with:

```bash
PYTHONPATH=src python3 -m active_coronagraph.cli --help
```

Wrong data path or missing FITS files: check that the data root contains
`binary_mask/` and `phase_screen/`, then pass `--data-dir` or set
`ACTIVE_CORONAGRAPH_DATA_DIR`.

Invalid coronagraph name: the parser accepts free-form `--name` values, but
unsupported names fail during mask construction with `Unknown focal mask type`.
Use one of the mask names listed in the CLI table above.

Invalid CLI option: run `active-coronagraph --help` and use only listed options.
Sampling, central obstruction, phase greyscale, calibration, and Lyot-stop
selection are not currently CLI flags.

Output directory problems: set `ACTIVE_CORONAGRAPH_OUTPUT_DIR` to a writable
directory. The code creates run directories under that location.

## Python API

For settings that are not exposed by the CLI, use the Python API. This minimal
example writes outputs into `output/python_example` and runs a small FQPM case:

```python
from pathlib import Path
import os

from active_coronagraph.coronagraphs_polished import simple_coro

out = Path("output/python_example")
out.mkdir(parents=True, exist_ok=True)
old = Path.cwd()
os.chdir(out)
try:
    image = simple_coro(
        dim=32,
        name="FQPM",
        wavelength=1.0,
        nsamp=2,
        fpm_sam=2,
        lyot_stop="None_1_1",
        charge=None,
        binary=False,
        obstruction=False,
        get_lyot=False,
        greyscale=None,
    )
finally:
    os.chdir(old)

print(image.shape)
```

For a runnable script, see `examples/minimal_simulation.py`.

## Running Tests

```bash
python3 -m pip install -e ".[dev]"
pytest
```

When running from a source tree without installation:

```bash
PYTHONPATH=src pytest
```

## Citation

If you use this software, cite it using the metadata in `CITATION.cff`.

## License

This project is distributed under the BSD 3-Clause License. See `LICENSE`.

