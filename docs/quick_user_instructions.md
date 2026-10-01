# Quick User Instructions

This page gives the shortest practical path for using Active Coronagraph from
Python. The detailed reference guide is in `docs/user_guide.md`.

## Install

From the repository root:

```bash
python -m pip install -e .
```

For development and tests:

```bash
python -m pip install -e ".[dev]"
pytest
```

Check the install:

```bash
python -c "import active_coronagraph; print(active_coronagraph.__version__)"
```

## Available Masks

Use these strings as the `name` argument in `run()` or `simple_coro()`:

| Mask name | Meaning | `charge` |
| --- | --- | --- |
| `FQPM` | four-quadrant phase mask | `None` |
| `vortex` | vortex phase mask | integer, for example `2` or `4` |
| `ACM` | active coronagraph mask | integer |
| `roddier` | Roddier phase mask | usually `None` |
| `dual_zone` | dual-zone phase mask | `None` |
| `dual_zone_old` | older dual-zone mask | `None` |
| `double_vortex` | double-vortex mask | integer |
| `binary_fits` | phase map loaded from a FITS file | `None` |

The binary-mask throughput helpers also accept `binary_mask`, `binary_masks`,
or `binary_vortex`. Those names run over FITS files in `data/binary_mask/`, or
over one file supplied with `phase_map_fits`.

## Main Functions

Most users should start with these functions:

| Function | Import | Use |
| --- | --- | --- |
| `run()` | `active_coronagraph.main_functions_polished` | Run a simulation and write outputs into the standard output tree |
| `simple_coro()` | `active_coronagraph.coronagraphs_polished` | Run one simulation in the current directory |
| `broadband_combine_final_psf()` | `active_coronagraph.main_functions_polished` | Combine already-saved wavelength-specific PSFs |
| `Phase_masks` | `active_coronagraph.phase_masks` | Generate analytic phase-mask arrays directly |
| `double_vortex_phase_mask()` | `active_coronagraph.new_mask` | Generate a double-vortex phase array directly |
| `throughput_1d_mono()` | `active_coronagraph.new_mask` | Compute 1D monochromatic throughput |
| `throughput_map_mono()` | `active_coronagraph.new_mask` | Compute 2D monochromatic throughput map |
| `directional_throughput_mono()` | `active_coronagraph.new_mask` | Compute throughput along horizontal, vertical, and diagonal directions |

## Run One Simulation

This is the recommended basic pattern:

```python
from active_coronagraph.main_functions_polished import run

image = run(
    dim=100,
    name="vortex",
    wavelength=1.0,
    nsamp=10,
    fpm_sam=10,
    charge=2,
    lyot_stop="large_cross",
    binary=False,
    obstruction=False,
    spider=False,
    get_lyot=False,
    greyscale=8,
    folder_name="_my_vortex_run",
)

print(image.shape)
```

Important arguments:

| Argument | Meaning |
| --- | --- |
| `dim` | base pupil/mask size used by the simulation |
| `name` | mask name, such as `FQPM`, `vortex`, or `dual_zone` |
| `wavelength` | dimensionless wavelength ratio; `1.0` is the reference value |
| `nsamp` | numerical pixels per `lambda/D` |
| `fpm_sam` | focal-plane mask pixels per `lambda/D` |
| `charge` | vortex/ACM/double-vortex charge, or `None` |
| `lyot_stop` | Lyot-stop preset, such as `large_cross` or `None_1_1` |
| `obstruction` | add central obstruction when `True` |
| `spider` | add spider arms when `True` |
| `get_lyot` | use `False` to return normalized final intensity |
| `greyscale` | phase quantization bits; use `None` for no quantization |
| `folder_name` | suffix for the output folder |

Outputs are written below a folder like:

```text
ideal_coro_2rd_mirror_False_my_vortex_run/
```

The final focal-plane FITS file is named like:

```text
final_focal_plane_coro_mask=vortex_2_sample=10_fpm_sam=10_obstruction_False.fits
```

## Run the Minimal Example

The repository includes a small FQPM example:

```bash
python examples/minimal_simulation.py --save-figure
```

It writes to:

```text
output/minimal_example/
```

and creates a normalized intensity curve from 0 to 15 `lambda/D`.

## Change the Mask

Only `name` and `charge` usually need to change.

FQPM:

```python
image = run(
    dim=100,
    name="FQPM",
    wavelength=1.0,
    nsamp=10,
    fpm_sam=10,
    charge=None,
    lyot_stop="large_cross",
    binary=False,
    obstruction=False,
    get_lyot=False,
    greyscale=8,
)
```

Vortex:

```python
image = run(
    dim=100,
    name="vortex",
    wavelength=1.0,
    nsamp=10,
    fpm_sam=10,
    charge=4,
    lyot_stop="large_cross",
    binary=False,
    obstruction=False,
    get_lyot=False,
    greyscale=8,
)
```

Dual-zone:

```python
image = run(
    dim=100,
    name="dual_zone",
    wavelength=1.0,
    nsamp=10,
    fpm_sam=10,
    charge=None,
    lyot_stop="large_cross",
    binary=False,
    obstruction=False,
    get_lyot=False,
    greyscale=8,
)
```

Double vortex:

```python
image = run(
    dim=100,
    name="double_vortex",
    wavelength=1.0,
    nsamp=10,
    fpm_sam=10,
    charge=2,
    lyot_stop="large_cross",
    binary=False,
    obstruction=False,
    get_lyot=False,
    greyscale=None,
    mask_opt_params={
        "dv_separation_ld": 2.0,
        "dv_pa_deg": 0.0,
        "dv_sigma_ld": 1.0,
        "dv_angle_deg": 0.0,
    },
)
```

## Generate a Phase Mask Array

Use `Phase_masks` when you only want the phase mask, not a full coronagraph
simulation:

```python
from active_coronagraph.phase_masks import Phase_masks

masks = Phase_masks(fpm_sam=10, pupil_size=100)

fqpm_phase = masks.FQPM()
vortex_phase = masks.vortex(charge=2)
acm_phase = masks.ACM(charge=2)
roddier_phase = masks.roddier(roddier_d=1.06)
dual_zone_phase = masks.dual_zone()
```

Each returned array has shape:

```text
(pupil_size * fpm_sam, pupil_size * fpm_sam)
```

For a double-vortex phase array:

```python
from active_coronagraph.new_mask import double_vortex_phase_mask

phase = double_vortex_phase_mask(
    dim=100,
    fpm_sam=10,
    charge=2,
    separation_ld=2.0,
    pa_deg=0.0,
    sigma_ld=1.0,
)
```

## Use a FITS Phase Mask

The FITS phase map must have shape:

```text
(dim * fpm_sam, dim * fpm_sam)
```

For the common `dim=100`, `fpm_sam=10` case, that is `(1000, 1000)`.

Example:

```python
from active_coronagraph.main_functions_polished import run

throughput = run(
    dim=100,
    name="binary_fits",
    wavelength=1.0,
    nsamp=10,
    fpm_sam=10,
    charge=None,
    lyot_stop="large_cross",
    binary=False,
    obstruction=False,
    get_lyot=False,
    greyscale=None,
    offset=(1.0, 0.0),
    mask_opt_params={
        "phase_map_fits": "data/binary_mask/Binary_vortex_sep0.5.fits",
    },
)
```

`offset=(x, y)` is in `lambda/D`. When `offset` is supplied, `run()` returns a
scalar throughput instead of an image.

## Broadband PSF

When you want a broadband PSF, run the simulation at several wavelength ratios
across the band and then combine the saved monochromatic FITS files.

The example below uses a 20% bandwidth centered on the reference wavelength:

```text
0.9, 0.95, 1.0, 1.05, 1.1
```

These are dimensionless wavelength ratios `lambda / lambda_0`. The range
`0.9` to `1.1` covers `±10%` around the reference wavelength, or 20% total
bandwidth. The endpoint weights are set to `0.5`, matching the trapezoidal
weighting pattern already used in the code.

```python
import numpy as np

from active_coronagraph.main_functions_polished import (
    broadband_combine_final_psf,
    run,
)

wavelengths = [0.9, 0.95, 1.0, 1.05, 1.1]
weights = np.ones(len(wavelengths))
weights[[0, -1]] = 0.5

for wl in wavelengths:
    run(
        dim=100,
        name="vortex",
        wavelength=wl,
        nsamp=10,
        fpm_sam=10,
        charge=2,
        lyot_stop="large_cross",
        binary=False,
        obstruction=False,
        get_lyot=False,
        greyscale=8,
        folder_name="_broadband_vortex",
    )

broadband_image = broadband_combine_final_psf(
    wavelength=wavelengths,
    fpm_sam=10,
    nsamp=10,
    weights=weights,
    name="vortex",
    charge=2,
    obstruction=False,
    lyot_stop="large_cross",
    folder_name="_broadband_vortex",
)
```

This code first runs five monochromatic simulations, one for each wavelength
ratio. `broadband_combine_final_psf()` then reads those saved FITS files and
writes their weighted average.

The broadband FITS filename ends in:

```text
_bb.fits
```

You can also run the same 20% broadband PSF workflow from the CLI:

```bash
active-coronagraph --function broadband_combine --name vortex --charge 2
```

For this workflow, the CLI defaults to the wavelength ratios
`[0.9, 0.95, 1.0, 1.05, 1.1]`, runs the wavelength-specific simulations, and
then writes the combined broadband product. To use a different wavelength grid,
pass it explicitly:

```bash
active-coronagraph --function broadband_combine \
  --name vortex --charge 2 \
  --wavelength 0.85 0.925 1.0 1.075 1.15
```

## Throughput

1D throughput along `+x`:

```python
from active_coronagraph.new_mask import throughput_1d_mono

out_file = throughput_1d_mono(
    mask_name="vortex",
    charge_val=2,
    wl=1.0,
    max_ld=5.0,
    step_ld=0.5,
)
print(out_file)
```

2D throughput map:

```python
from active_coronagraph.new_mask import throughput_map_mono

out_file = throughput_map_mono(
    mask_name="vortex",
    charge_val=2,
    wl=1.0,
    min_ld=0.0,
    max_ld=5.0,
    step_ld=1.0,
)
print(out_file)
```

Directional throughput:

```python
from active_coronagraph.new_mask import directional_throughput_mono

out_file = directional_throughput_mono(
    mask_name="vortex",
    charge_val=2,
    wl=1.0,
    max_ld=5.0,
    step_ld=0.5,
)
print(out_file)
```

## Phase Screens

Load one 2D phase screen from a FITS cube and pass it as `phase`:

```python
from astropy.io import fits
from active_coronagraph.main_functions_polished import run

with fits.open("data/phase_screen/TROIA_phase_screens_new_jitter5percentLamdaOverD.fits") as hdul:
    phase = hdul[0].data[0]

image = run(
    dim=100,
    name="FQPM",
    wavelength=1.0,
    nsamp=10,
    fpm_sam=10,
    charge=None,
    lyot_stop="None_1_1",
    binary=False,
    obstruction=True,
    get_lyot=False,
    greyscale=8,
    phase=phase,
    folder_name="_phase_screen_example",
)
```

The bundled phase-screen cubes have shape `(100, 100, 100)`, so one slice has
shape `(100, 100)`.

## Command-Line Shortcuts

Basic CLI run:

```bash
active-coronagraph --function main_func --name vortex --charge 2
```

Minimal example:

```bash
python examples/minimal_simulation.py --save-figure
```

CLI help:

```bash
active-coronagraph --help
```

For new work, the Python examples above are usually clearer than the legacy CLI
because they show every parameter being used.
