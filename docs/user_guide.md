# Active Coronagraph User Guide

This guide describes how to install and use the current repository state. It is
based on the package metadata, example script, CLI, tests, bundled data, and the
implementation under `src/active_coronagraph`.

## 1. What This Software Does

Active Coronagraph is a scientific Python package for numerical coronagraph
simulations and focal-plane phase-mask studies. It builds entrance pupils and
Lyot stops, generates focal-plane phase masks, propagates complex optical fields
through a coronagraph model, and writes focal-plane products for later analysis.

The current code supports practical workflows such as:

- running monochromatic simulations and combining wavelength-specific products
  into broadband products;
- comparing phase-mask designs;
- generating and applying analytic focal-plane masks;
- using binary FITS phase masks supplied in `data/binary_mask/`;
- computing radial normalized-intensity or contrast curves from simulated PSFs;
- computing selected monochromatic throughput products;
- averaging simulations over supplied phase-screen FITS cubes in legacy
  workflows.

Implemented mask names available through the current polished propagation path
include:

- `FQPM`: four-quadrant phase mask;
- `vortex`: vortex phase mask with user-supplied charge;
- `ACM`: active coronagraph mask based on a cosine phase pattern;
- `roddier`: Roddier phase mask;
- `dual_zone`: dual-zone phase mask;
- `dual_zone_old`: older dual-zone parameterization;
- `double_vortex`: testable double-vortex phase mask;
- `binary_fits`: phase map read from a FITS file, used internally by binary-mask
  throughput helpers.

The software produces FITS files for phase masks and final focal-plane
intensity images, text tables for radial intensity/contrast or throughput
products, and PNG diagnostic plots in selected workflows.

## 2. Requirements

The package requires Python 3.10 or newer, as declared in `pyproject.toml`.

Required Python dependencies are:

| Dependency | Used for |
| --- | --- |
| `numpy` | array operations and numerical calculations |
| `matplotlib` | plotting and diagnostic figures |
| `astropy` | FITS input/output |
| `scipy` | numerical/image utilities used by legacy plotting helpers |
| `tqdm` | progress bars for longer throughput and CLI workflows |
| `hcipy` | pupil grids, apertures, spiders, and optical utilities |

Development-only dependency:

| Dependency | Used for |
| --- | --- |
| `pytest` | test suite |

The code is filesystem-oriented and writes generated products to directories
under the current working directory unless `ACTIVE_CORONAGRAPH_OUTPUT_DIR` is
set. Matplotlib also needs a writable cache directory. On restricted systems,
set `MPLCONFIGDIR` to a writable directory such as `/tmp/matplotlib-cache`.

Computational cost scales quickly with `dim`, `nsamp`, `fpm_sam`, throughput
grid size, and number of phase-screen realizations. The minimal example uses
`dim=32`, `nsamp=2`, and `fpm_sam=2` specifically so it runs quickly. Several
legacy CLI workflows use defaults such as `dim=100`, `nsamp=10`, and
`fpm_sam=10`, which create substantially larger arrays and many output files.

## 3. Installation

Clone the repository and enter it:

```bash
git clone <repository-url>
cd active_coronagraph
```

Create and activate a virtual environment:

```bash
python -m venv .venv
source .venv/bin/activate
```

On systems where `python` is not Python 3, use `python3` instead:

```bash
python3 -m venv .venv
source .venv/bin/activate
```

Install the package for normal runtime use:

```bash
python -m pip install .
```

For an editable development install:

```bash
python -m pip install -e .
```

To include the development test dependency:

```bash
python -m pip install -e ".[dev]"
```

The repository metadata defines a console script named `active-coronagraph`.
There is no PyPI installation command documented in the repository, so install
from the source tree unless a package release is provided separately.

## 4. Verify the Installation

Check that the package imports:

```bash
python -c "import active_coronagraph; print(active_coronagraph.__version__)"
```

The current package version is exposed as `active_coronagraph.__version__`.

Check the CLI help:

```bash
active-coronagraph --help
```

From a source tree, this equivalent module command also works:

```bash
python -m active_coronagraph.cli --help
```

Run the test suite if development dependencies are installed:

```bash
pytest
```

If running directly from a source tree without installation:

```bash
PYTHONPATH=src pytest
```

## 5. First Simulation

The intended beginner example is `examples/minimal_simulation.py`. Run it from
the repository root or from an installed package environment:

```bash
python examples/minimal_simulation.py
```

To also save a PNG plot:

```bash
python examples/minimal_simulation.py --save-figure
```

By default, the example writes to:

```text
output/minimal_example/
```

You can choose another output directory:

```bash
python examples/minimal_simulation.py --output-dir output/my_first_run
```

The script simulates a small, quick FQPM coronagraph case:

| Quantity | Value in the example |
| --- | --- |
| `dim` | `32` |
| `nsamp` | `2` |
| `fpm_sam` | `2` |
| `wavelength` | `1.0` |
| `name` | `"FQPM"` |
| `charge` | `None` |
| `lyot_stop` | `"None_1_1"` |
| `binary` pupil mode | `False` |
| central obstruction | `False` |
| spiders | default `False` |
| `get_lyot` | `False` |
| greyscale quantization | `None` |

The wavelength value is a dimensionless scaling factor used by the simulation
code. The default `1.0` corresponds to the reference wavelength in the code's
normalization.

The example calls `active_coronagraph.coronagraphs_polished.simple_coro()`.
Because `get_lyot=False` and `offset=None`, `simple_coro()` returns a normalized
2D focal-plane intensity image and writes a FITS file in the current output
directory.

The example then computes the median normalized intensity in annuli from 0 to
15 `lambda/D` and writes:

```text
normalized_intensity_curve_0_to_15_lambda_over_D.txt
```

The first column is annulus radius in `lambda/D`; the second column is median
normalized intensity in that annulus.

Additional files written by the underlying simulation include:

```text
m_vor_FQPM_None_nsamp=2_fpm_sam=2_rotate=False_greyscale=None_gmethod=legacy_wl_1.0_cal_1.0.fits
final_focal_plane_coro_mask=FQPM_None_sample=2_fpm_sam=2_obstruction_False.fits
```

If `--save-figure` is used, the example also writes:

```text
minimal_fqpm_intensity_curve.png
```

## 6. Understanding the Minimal Example

The example has three main pieces:

- `working_directory(path)`: creates the output directory and temporarily runs
  the simulation there because `simple_coro()` writes files to the current
  working directory.
- `annular_intensity_curve(image, nsamp, max_radius_ld, step_ld=1.0)`: converts
  pixel radius to `lambda/D` using `radius_ld = radius_pixels / nsamp` and
  computes median image values in radial annuli.
- `run_simulation(output_dir, save_figure)`: configures and runs the FQPM
  simulation, writes the text curve, and optionally writes a plot.

Important parameters in `run_simulation()`:

| Variable | Meaning | Unit/convention | Effect |
| --- | --- | --- | --- |
| `dim` | Linear pupil/grid scale used by the simulation setup | pixels and `lambda/D` field scale in code conventions | Larger values increase array sizes and field extent |
| `nsamp` | Numerical focal-plane sampling | pixels per `lambda/D` | Used to convert image-pixel radius to `lambda/D`; larger values increase output array size |
| `fpm_sam` | Focal-plane mask sampling | mask pixels per `lambda/D` | Sets the resolution of the phase mask before resampling to the simulation grid |
| `wavelength` | Wavelength scaling factor | dimensionless ratio to the reference wavelength | Used in pupil sizing and final mask phase scaling |
| `name` | Mask name passed to `simple_coro()` | string | Selects the implemented FQPM mask |
| `charge` | Vortex/ACM/double-vortex charge | integer or `None` | Not used by FQPM, so the example passes `None` |
| `lyot_stop` | Lyot-stop preset name | string | `"None_1_1"` uses the preset with fraction `1.0` and secondary fraction `1.0` |
| `binary` | Entrance-pupil mode | boolean | `False` uses `make_aperture()` rather than `make_entrance_pupil()` |
| `obstruction` | Central obstruction flag | boolean | `False` makes an unobstructed circular aperture |
| `get_lyot` | Whether to return/save Lyot-normalized complex-field products | boolean | `False` makes `simple_coro()` return normalized final intensity |
| `greyscale` | Phase quantization setting | integer, `True`, `False`, or `None` | `None` disables quantization in the example |
| `max_radius_ld` | Maximum annular-profile radius | `lambda/D` | The example uses `15.0` |
| `step_ld` | Annulus width | `lambda/D` | Default is `1.0` |
| `output_dir` | Destination directory | filesystem path | The example changes into this directory before running |

## 7. Basic Python API

The package does not currently expose a compact high-level simulation API from
`active_coronagraph.__init__`; it only exposes `__version__`. Users should import
the specific functions they need from implementation modules.

For a direct script similar to the minimal example:

```python
from pathlib import Path
import os
import numpy as np

from active_coronagraph.coronagraphs_polished import simple_coro

output_dir = Path("output/my_api_run")
output_dir.mkdir(parents=True, exist_ok=True)

old_cwd = Path.cwd()
try:
    os.chdir(output_dir)
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
    os.chdir(old_cwd)

print(image.shape)
print(np.nanmin(image), np.nanmax(image))
```

For the standard output directory structure used by the CLI, use
`main_functions_polished.run()`:

```python
from active_coronagraph.main_functions_polished import run

image = run(
    dim=32,
    name="FQPM",
    wavelength=1.0,
    nsamp=2,
    fpm_sam=2,
    charge=None,
    lyot_stop="None_1_1",
    binary=False,
    obstruction=False,
    get_lyot=False,
    greyscale=None,
    folder_name="_my_run",
)
```

This writes under a root directory whose name is built from the obstruction,
noise, and folder settings, for example:

```text
ideal_coro_2rd_mirror_False_my_run/
```

Important return behavior:

| Function call pattern | Return value | Files written |
| --- | --- | --- |
| `simple_coro(..., get_lyot=False, offset=None)` | normalized 2D intensity array | phase-mask FITS and final focal-plane FITS |
| `simple_coro(..., get_lyot=False, offset=(x, y))` | scalar throughput | phase-mask FITS and `throughput.txt` |
| `simple_coro(..., get_lyot=True)` | complex final focal-plane field | phase-mask FITS; optional Lyot-plane files are controlled internally |
| `run(..., do_res=True)` | residual-energy scalar from `res_ene()` | standard run directory products for that workflow |

The current API is filesystem-oriented. Several functions write products as a
side effect in the current working directory or a generated run directory.

## 8. Configuration

Configuration is currently path-based and uses environment variables plus a few
CLI options. There is no separate configuration file or configuration object.

The path helpers are in `src/active_coronagraph/config.py`:

| Parameter / function | Meaning | Unit | Default | Allowed values |
| --- | --- | --- | --- | --- |
| `ACTIVE_CORONAGRAPH_DATA_DIR` / `data_dir()` | Runtime data root containing `binary_mask/` and `phase_screen/` | path | `<repo>/data` | Any existing data directory with expected subdirectories |
| `ACTIVE_CORONAGRAPH_OUTPUT_DIR` / `output_dir()` | Base directory for generated simulation outputs | path | current working directory | Any writable directory |
| `phase_screen_dir()` | Directory containing phase-screen FITS cubes | path | `data_dir() / "phase_screen"` | Derived from data root |
| `binary_mask_dir()` | Directory containing binary-mask FITS files | path | `data_dir() / "binary_mask"` | Derived from data root |

Example:

```bash
ACTIVE_CORONAGRAPH_OUTPUT_DIR=/path/to/output \
active-coronagraph --function main_func --name dual_zone --charge None
```

The CLI also accepts selected path overrides:

```bash
active-coronagraph --function compare_phase_screens --data-dir data
```

CLI defaults are module-level variables in `cli.py`. Important defaults include:

| Variable | Meaning | Default |
| --- | --- | --- |
| `dim` | simulation dimension | `100` |
| `nsamp` | numerical sampling | `10` |
| `fpm_sam` | focal-plane mask sampling | `10` |
| `wavelength` | resolved wavelength list | `[1.0]` for monochromatic workflows; `[0.9, 0.95, 1.0, 1.05, 1.1]` for broadband workflows |
| `lyot_stop` | selected Lyot stop preset | first entry of `lyot_stop_default`, currently `"large_cross"` |
| `obstruction` | central obstruction flag | `False` |
| `spider` | spider flag | `False` |
| `greyscale` / `gs` | greyscale passed by CLI workflows | `None` in `main_func()` |
| `noise_level` | Gaussian phase-noise level | `0` |
| `phase_shift` | focal-plane shift as `(DY, DX)` | `(0.0, 0.0)` |
| `phase_shift_unit` | shift units | `"array"` |
| `coro_shift` | half-pixel coronagraph centering shift | `True` |

Only some CLI options are plumbed into every workflow. For example,
`--phase-shift`, `--phase-shift-unit`, `--coro-shift`, data path options, and
throughput options are used by their corresponding branches. `--wavelength`
sets the wavelength-ratio list used by workflows that iterate over the global
CLI wavelength list.

## 9. Coronagraph Types

The following mask names are accepted by
`active_coronagraph.coronagraphs_polished.focal_mask_new()` and therefore by
`simple_coro()` / `run()` when routed through the polished path.

| Code name | Scientific name / description | Required parameters | Optional parameters | Example |
| --- | --- | --- | --- | --- |
| `FQPM` | Four-quadrant phase mask | `charge=None` is customary | `vortex_center` through `simple_coro()` | `simple_coro(..., name="FQPM", charge=None)` |
| `vortex` | Vortex phase mask | integer `charge` | `rotate`, `vortex_center` | `simple_coro(..., name="vortex", charge=2)` |
| `ACM` | Active coronagraph mask using `2.4048*cos(charge*theta)` before phase wrapping | integer `charge` | none exposed beyond common options | `simple_coro(..., name="ACM", charge=2)` |
| `roddier` | Roddier phase mask | `charge=None` is customary | `roddier_d`, default `1.06` | `simple_coro(..., name="roddier", charge=None, roddier_d=1.06)` |
| `dual_zone` | Dual-zone phase mask | none beyond common arguments | fixed defaults `a1=0.515`, `a2=0.705`, `z1=0.47`, `z2=0.92` inside `focal_mask_new()` | `simple_coro(..., name="dual_zone", charge=None)` |
| `dual_zone_old` | Older dual-zone phase mask | none beyond common arguments | fixed defaults in `Phase_masks.dual_zone_old()` | `simple_coro(..., name="dual_zone_old", charge=None)` |
| `double_vortex` | Double-vortex phase mask | non-`None` integer `charge` | `dv_*` values through `mask_opt_params` or throughput CLI options | `run(..., name="double_vortex", charge=2, mask_opt_params={"dv_separation_ld": 2.0})` |
| `binary_fits` | Phase map read from a FITS file | `opt_params["phase_map_fits"]` | common options | usually used by binary-mask throughput helpers |

Double-vortex example:

```python
from active_coronagraph.main_functions_polished import run

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
        "dv_sigma_ld": 1.0,
        "dv_angle_deg": 0.0,
        "dv_separation_ld": 2.0,
        "dv_pa_deg": 0.0,
    },
)
```

Binary FITS masks are exposed most conveniently by the throughput helpers using
`mask_name="binary_mask"`, `mask_name="binary_masks"`, or
`mask_name="binary_vortex"`. Those names cause the helpers to run one job per
FITS file in `binary_mask_dir()`, or one job for `phase_map_fits` if provided.

## 10. Phase Masks

Analytic phase-mask generation is implemented in `phase_masks.py` by the
`Phase_masks` class:

```python
from active_coronagraph.phase_masks import Phase_masks

masks = Phase_masks(fpm_sam=10, pupil_size=100)
phase = masks.FQPM()
```

The constructor arguments are:

| Argument | Meaning |
| --- | --- |
| `fpm_sam` | focal-plane mask sampling in pixels per `lambda/D` |
| `pupil_size` | linear field size used to build the mask |

Generated analytic phase arrays have shape:

```text
(pupil_size * fpm_sam, pupil_size * fpm_sam)
```

User-facing methods:

| Function | Purpose | Arguments | Output |
| --- | --- | --- | --- |
| `dual_zone(show=False, a1=0.515, a2=0.705, z1=0.47, z2=0.92)` | Build a two-zone radial phase mask | radii `a1`, `a2` in `lambda/D` scaled by `fpm_sam`; phase fractions `z1`, `z2` multiplied by `2*pi` | 2D float phase array in radians |
| `dual_zone_old(show=False, a1=0.475, a2=0.8, z1=0.34, z2=0.74)` | Older two-zone implementation | same style as `dual_zone()` | 2D float phase array in radians |
| `FQPM(center=0)` | Build a four-quadrant phase mask | optional center offset | 2D float array containing `0` and `pi` regions |
| `vortex(charge, rotate=False, center=0)` | Build a wrapped vortex phase | integer charge, optional 45-degree rotation, optional center | 2D float phase array wrapped through complex angle arithmetic into a nonnegative range |
| `ACM(charge, rotate=False)` | Build active coronagraph mask phase | integer charge | 2D float phase array |
| `roddier(roddier_d)` | Build a circular Roddier mask | diameter parameter scaled by `fpm_sam` | 2D float phase array, aperture times `pi` |

`coronagraphs_polished.focal_mask_new()` converts these phase arrays into a
complex focal-plane mask:

```python
from active_coronagraph.coronagraphs_polished import focal_mask_new

complex_mask = focal_mask_new(
    name="vortex",
    dim=100,
    nsamp=10,
    fpm_sam=10,
    wavelength=1.0,
    charge=2,
    greyscale=None,
)
```

Important conventions implemented in `focal_mask_new()`:

- the initial analytic mask is generated on a grid of `dim * fpm_sam` pixels;
- the mask is resampled to the numerical grid with
  `np.repeat(..., scale)` where `scale = nsamp / fpm_sam`;
- this repeat operation expects `nsamp / fpm_sam` to be usable as an integer
  repeat count;
- phase values are transformed through wrapping operations involving `pi` and
  scaled by `1 / wavelength`;
- `greyscale` quantizes phase to `2**greyscale` levels when it is not `None`;
- `greyscale_method` accepts `"legacy"` or `"new"`;
- `fill_factor=True` zeros every 10th row and column in the amplitude mask;
- `noise_level` adds Gaussian phase noise with standard deviation
  `noise_level * 2*pi`;
- the function writes a FITS file named `m_vor_...fits` containing
  `np.angle(m_vor) + pi`;
- the return value is a complex-valued mask array.

The double-vortex generator in `new_mask.py` is also directly usable:

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

It returns a phase array in `[0, 2*pi)` with shape
`(dim * fpm_sam, dim * fpm_sam)`. The tests verify this shape and range.

## 11. Pupil and Optical Parameters

Entrance pupils and Lyot stops are built in `coronagraphs_polished.py`.

`make_aperture()` builds an analytic HCIPy aperture sampled onto an array:

```python
from active_coronagraph.coronagraphs_polished import make_aperture

pupil = make_aperture(
    dim=100,
    oversample=8,
    obstruction=True,
    spider=False,
    obstruction_ratio=0.25,
    phase=None,
)
```

Parameters:

| Parameter | Meaning | Default |
| --- | --- | --- |
| `dim` | linear aperture dimension used for the pupil grid | required |
| `oversample` | supersampling factor passed to HCIPy aperture evaluation | required |
| `obstruction` | include central obstruction | `True` |
| `spider` | include four spider arms | `False` |
| `obstruction_ratio` | central obstruction ratio | `0.25` |
| `lyot` | build a Lyot stop instead of an entrance pupil | `False` |
| `lyot_fraction` | Lyot-stop outer diameter fraction when `lyot=True` | `None`, treated as `1.0` |
| `nsamp` | Lyot-stop sampling multiplier when `lyot=True` | `10` |
| `nor` | normalize field intensity with `Isum()` | `True` |
| `phase` | optional phase array applied as `exp(1j * phase)` | `None` |

`make_entrance_pupil(dim, d=None, nor=True)` builds a binary circular pupil
using project helper `circle_mask()`. If `d` is `None`, it uses `dim` as the
diameter. With `nor=True`, it normalizes total intensity to 100.

`make_lyot_stop()` wraps `make_aperture(..., lyot=True)`:

```python
from active_coronagraph.coronagraphs_polished import make_lyot_stop

lyot = make_lyot_stop(
    dim=100,
    fraction=0.95,
    sec_fraction=1.1,
    obstruction=True,
    nsamp=10,
)
```

The Lyot stop uses:

```text
obstruction_ratio = 0.25 * sec_fraction / fraction
```

`simple_coro()` provides named Lyot-stop presets:

| Preset | `fraction` | `sec_fraction` |
| --- | ---: | ---: |
| `small_cross` | `0.98` | `1.2` |
| `large_cross` | `0.95` | `1.3` |
| `zernike_cross` | `0.98` | `1.05` |
| `None` | `1.0` | `1.2` |
| `None_1_1` | `1.0` | `1.0` |
| `None_1_1.1` | `1.0` | `1.1` |
| `None_1_1.3` | `1.0` | `1.3` |
| `None_1_1.4` | `1.0` | `1.4` |
| `None_1_1.5` | `1.0` | `1.5` |

The spider geometry in the polished module uses a pupil diameter of `4.0` m, a
spider width of `1.958e-2` m, and a spider offset of `[0.0, 0.37251]` m before
normalization to the sampled grid.

## 12. Sampling

Sampling is represented by several parameters:

| Parameter | Meaning in code |
| --- | --- |
| `dim` | base linear dimension used for pupil size and field/mask construction |
| `nsamp` | numerical sampling in focal-plane pixels per `lambda/D`; final FFT padding uses `nsamp * dim` |
| `fpm_sam` | focal-plane mask sampling in mask pixels per `lambda/D`; analytic masks are generated with size `dim * fpm_sam` |
| `oversample` | HCIPy supersampling factor for evaluating apertures |

For a normalized intensity image returned by `simple_coro()`, the minimal
example interprets radius as:

```python
radius_ld = radius_pixels / nsamp
```

This convention is also consistent with `make_plot.annuli_mean()`, which builds
annuli using `nsamp` and reports radial bins in `lambda/D`.

The current `focal_mask_new()` implementation rescales the phase mask with
`np.repeat()` using:

```python
scale = nsamp / fpm_sam
```

In normal documented examples, use integer-compatible values such as
`nsamp=10`, `fpm_sam=10` or `nsamp=2`, `fpm_sam=2`.

Increasing sampling generally increases array sizes and runtime. The repository
does not document quantitative convergence or accuracy thresholds, so users
should validate sampling choices for their own scientific use case.

## 13. Wavelength

Wavelength is represented as a dimensionless scaling factor, not a physical
length with units. The CLI help describes it as a wavelength ratio
`lambda/lambda_0`, and default workflows use `1.0`.

In `simple_coro()` and related functions, `wavelength` affects:

- the aperture size passed to `make_aperture(round_to_even(dim / wavelength))`;
- the final focal-plane mask phase scaling in `focal_mask_new()` through
  multiplication by `1 / wavelength`.

The Lyot stop and final output grid use the base `dim` and `nsamp` values so
all wavelength-specific PSFs in a broadband run have the same array shape and
can be combined.

Monochromatic operation is the clearest beginner workflow, but the code also
supports broadband products by combining simulations that have already been run
at multiple wavelength ratios. The broadband workflow is not a separate optical
propagator; it is a weighted average of saved wavelength-specific outputs.

Broadband helpers:

- `main_functions_polished.broadband_combine_final_psf()`;
- `main_functions_polished.broadband_combine_tp()`;
- CLI workflows such as `broadband_combine`, `broadband_kit`, and `tp_bb`.

`broadband_combine_final_psf()` expects that matching monochromatic final-PSF
FITS files already exist for every wavelength in the wavelength list. It reads
those files, applies the supplied weights, and writes a broadband FITS product
with `_bb.fits` in the filename.

Example Python workflow for a 20% bandwidth centered on the reference
wavelength:

```text
0.9, 0.95, 1.0, 1.05, 1.1
```

These dimensionless wavelength ratios cover `±10%` around `lambda_0`, or 20%
total bandwidth. The first and last weights are `0.5`, matching the
trapezoidal weighting pattern used by the existing workflow.

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
        folder_name="_broadband_example",
    )

broadband_psf = broadband_combine_final_psf(
    wavelength=wavelengths,
    fpm_sam=10,
    nsamp=10,
    weights=weights,
    name="vortex",
    charge=2,
    obstruction=False,
    lyot_stop="large_cross",
    folder_name="_broadband_example",
)
```

This example writes one monochromatic final-PSF FITS file per wavelength and
then writes:

```text
final_focal_plane_coro_mask=vortex_2_sample=10_fpm_sam=10_obstruction_False_bb.fits
```

inside the `lambda_1.0` run directory for the selected output folder.

The CLI `broadband_combine` workflow uses the same 20% wavelength grid by
default, runs each wavelength-specific simulation, and then combines the saved
FITS products. Passing `--wavelength` overrides the default grid.

## 14. Using Provided Data

The `data/` directory contains runtime data for selected workflows:

```text
data/binary_mask/
data/phase_screen/
```

`data/binary_mask/` contains three binary focal-plane phase-mask FITS files and
a preview PNG:

```text
Binary_vortex_sep0.5.fits
Binary_vortex_sep0.5_45posang_SLM_axis.fits
Binary_vortex_sep0.5_45posang_center_axis.fits
binary_masks_cyclic_cmap.png
```

The tests verify that bundled binary mask FITS files have:

- shape `(1000, 1000)`;
- floating-point data;
- FITS header `BUNIT = "rad"`.

`data/phase_screen/` contains TROIA phase-screen FITS cubes:

```text
TROIA_phase_screens_new_jitter0percentLamdaOverD.fits
TROIA_phase_screens_new_jitter5percentLamdaOverD.fits
TROIA_phase_screens_new_jitter10percentLamdaOverD.fits
TROIA_phase_screens_new_jitter20percentLamdaOverD.fits
```

The tests verify that these phase-screen files have shape `(100, 100, 100)` and
floating-point data. Their FITS headers do not currently define `BUNIT`.

The repository `data/README.md` states that provenance, generation method,
license, and expected citation still need to be documented before publication.

The software locates runtime data through:

```python
from active_coronagraph.config import binary_mask_dir, phase_screen_dir
```

or with `ACTIVE_CORONAGRAPH_DATA_DIR`:

```bash
ACTIVE_CORONAGRAPH_DATA_DIR=/path/to/data active-coronagraph --function compare_phase_screens
```

If users provide their own binary FITS phase masks for `binary_fits`, the code
requires the phase-map array shape to equal:

```text
(dim * fpm_sam, dim * fpm_sam)
```

For the default throughput helpers, `dim=100` and `fpm_sam=10`, so the expected
shape is `(1000, 1000)`. The code treats values as phase in radians.

## 15. Using a Phase Screen

Phase-screen support exists in two forms.

First, the low-level simulation functions accept a phase array:

```python
from astropy.io import fits
from active_coronagraph.main_functions_polished import run

with fits.open("data/phase_screen/TROIA_phase_screens_new_jitter5percentLamdaOverD.fits") as hdul:
    phase_cube = hdul[0].data

phase = phase_cube[0]

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
    folder_name="_phase_example",
)
```

The phase array is passed to `make_aperture()` and applied to the entrance pupil
as:

```python
A = A * np.exp(1j * phase)
```

The supplied phase-screen cubes are `(100, 100, 100)`. A selected slice has
shape `(100, 100)`, matching `dim=100` at `wavelength=1.0` in the default
phase-screen workflows.

The CLI also includes a phase-screen comparison workflow:

```bash
active-coronagraph --function compare_phase_screens --name FQPM --charge None
```

This reads every `*.fits` file from `phase_screen_dir()`, runs each realization
in each cube, averages the resulting PSFs with `combine_phase_screen_psf()`, and
writes a comparison plot.

To point at a different phase-screen folder:

```bash
active-coronagraph --function compare_phase_screens \
  --name FQPM --charge None \
  --phase-screen-folder /path/to/phase_screen
```

The code does not currently document physical units for the phase-screen values.
Use the supplied files and any custom replacements with care, and preserve array
shape compatibility with the aperture being simulated.

## 16. Outputs

The exact outputs depend on which entry point is used.

For a direct `simple_coro()` call with `get_lyot=False` and no `offset`, the
return value is a normalized 2D intensity array:

```python
norm_final = (abs(final_field) ** 2) / wo_mask_peak
```

The corresponding FITS file is:

```text
final_focal_plane_coro_mask=<name>_<charge>_sample=<nsamp>_fpm_sam=<fpm_sam>_obstruction_<obstruction>.fits
```

The FITS data is written as `[norm_final]`, so readers should expect the primary
data array to have a leading length-1 axis.

Every `focal_mask_new()` call writes a phase-mask FITS file named like:

```text
m_vor_<name>_<charge>_nsamp=<nsamp>_fpm_sam=<fpm_sam>_rotate=<rotate>_greyscale=<greyscale>_gmethod=<method>_wl_<wavelength>_cal_<cal_factor>.fits
```

For throughput mode, `simple_coro(..., get_lyot=False, offset=(x, y))` returns a
scalar throughput and writes:

```text
throughput.txt
```

The `offset` convention for throughput is `(x, y)` in `lambda/D`.

`main_functions_polished.run()` writes into a generated output tree. The root is
computed from obstruction, noise level, and optional folder name:

```text
ideal_coro_2rd_mirror_<obstruction>[_noise_<noise_level>][<folder_name>]/
```

The per-run directory includes the Lyot stop, mask name, charge, wavelength,
`fpm_sam`, binary flag, obstruction flag, greyscale setting, and optional phase
shift.

`make_plot.annuli_mean()` reads final focal-plane FITS files and writes radial
median curves. For monochromatic products, the filename is:

```text
final_focal_plane_coro_mask=<name>_<charge>_sample=<nsamp>_fpm_sam=<fpm_sam>_obstruction_<obstruction>_wl_<wl>.txt
```

For broadband products, the filename ends in:

```text
_bb.txt
```

`broadband_combine_final_psf()` writes a weighted-average broadband FITS file
next to the wavelength-`1.0` product, with the filename:

```text
final_focal_plane_coro_mask=<name>_<charge>_sample=<nsamp>_fpm_sam=<fpm_sam>_obstruction_<obstruction>_bb.fits
```

The annular bins use rings of width `nsamp / 10` pixels and report radii in
`lambda/D`.

Throughput helper outputs:

| Helper / CLI function | Output directory | Main text output |
| --- | --- | --- |
| `throughput_1d_mono()` / `tp_1d_mono` | `throughput_1d_mono/` | `<mask>_<charge>_wl_<wl>_tp_1d_0_to_<max_ld>_step_<step_ld>.txt` |
| `throughput_map_mono()` / `tp_map_mono` | `throughput_map_mono/` | `<mask>_<charge>_wl_<wl>_tp_map_<min_ld>_to_<max_ld>_step_<step_ld>.txt` |
| `directional_throughput_mono()` / `tp_directional_mono` | `throughput_directional_mono/` | `<mask>_<charge>_wl_<wl>_tp_0_to_<max_ld>_step_<step_ld>.txt` |

These helpers also write PNG diagnostic plots.

## 17. Command-Line Workflows

The installed CLI entry point is:

```bash
active-coronagraph --function <workflow>
```

The available workflow choices are defined in `cli.py`:

```text
main_func
res_ene
annuli_mean
summary_plot
broadband_combine
tp_bb
plot_tp
fpm_plot
broadband_kit
add_fill_factor
ghost_im
compare_phase_screens
tp_directional_mono
tp_map_mono
tp_1d_mono
```

Basic run:

```bash
active-coronagraph --function main_func --name FQPM --charge None
```

Vortex run:

```bash
active-coronagraph --function main_func --name vortex --charge 2
```

Broadband final-PSF simulation and combination is available through:

```bash
active-coronagraph --function broadband_combine --name vortex --charge 2
```

By default, this runs the 20% bandwidth wavelength grid
`[0.9, 0.95, 1.0, 1.05, 1.1]`, combines the saved monochromatic FITS files with
endpoint-half trapezoid weights, and writes the `_bb.fits` product.

To override the wavelength grid:

```bash
active-coronagraph --function broadband_combine \
  --name vortex --charge 2 \
  --wavelength 0.85 0.925 1.0 1.075 1.15
```

Monochromatic 1D throughput:

```bash
active-coronagraph --function tp_1d_mono \
  --name vortex --charge 2 \
  --tp-1d-max-ld 5 \
  --tp-1d-step-ld 0.5
```

Binary-mask throughput using all FITS files in `data/binary_mask/`:

```bash
active-coronagraph --function tp_1d_mono \
  --name binary_mask --charge None \
  --tp-1d-max-ld 5 \
  --tp-1d-step-ld 0.5
```

Binary-mask throughput using one specific FITS file:

```bash
active-coronagraph --function tp_1d_mono \
  --name binary_mask --charge None \
  --phase-map-fits data/binary_mask/Binary_vortex_sep0.5.fits \
  --tp-1d-max-ld 5 \
  --tp-1d-step-ld 0.5
```

`plot_tp` is currently a placeholder branch that prints that the original
function is commented out.

## 18. Running Your Own Simulations

For new simulations, prefer starting from `examples/minimal_simulation.py` or a
small direct call to `run()`. A practical workflow is:

1. Choose a mask name supported by `focal_mask_new()`.
2. Choose `dim`, `nsamp`, and `fpm_sam`.
3. Choose `wavelength=1.0` unless you are intentionally studying wavelength
   scaling.
4. Choose a Lyot-stop preset.
5. Decide whether to include a central obstruction and spiders.
6. Decide whether phase quantization (`greyscale`) should be used.
7. Set `ACTIVE_CORONAGRAPH_OUTPUT_DIR` or pass `folder_name` to keep outputs
   organized.
8. Run a small case first, then increase sampling after validating file shapes
   and output locations.

Example:

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
    folder_name="_vortex_charge2",
)
```

For a custom FITS phase mask, use throughput helpers or route through
`binary_fits` with `mask_opt_params`:

```python
from active_coronagraph.main_functions_polished import run

tp = run(
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

For custom phase screens, load a FITS cube with Astropy, select one 2D slice,
and pass it as `phase=...` as shown in the phase-screen section.

## 19. Limitations and Current Caveats

- The public API is not yet concentrated in `active_coronagraph.__init__`; users
  import from implementation modules.
- Many routines write files as side effects.
- The CLI is a workflow wrapper with module-level defaults. Some options
  are workflow-specific and are not applied globally.
- Several plotting helpers assume specific generated directory names.
- The bundled data README does not yet document provenance, generation method,
  license, or citation information.
- The repository tests are currently smoke and shape tests; scientific users
  should add or run numerical validation appropriate to their study.
