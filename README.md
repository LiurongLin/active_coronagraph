# Active Coronagraph

Scientific Python code for coronagraph simulations and focal-plane mask studies.

## Repository Layout

- `src/active_coronagraph/`: installable Python package.
- `data/`: bundled runtime FITS/PNG data used by selected workflows.
- `examples/`: small usage examples.
- `tests/`: smoke tests for imports, configuration, and the CLI.
- `docs/`: publication documentation scaffold.
- `scripts/`: reserved for maintenance or one-off helper scripts.

## Installation

```bash
python -m pip install -e ".[dev]"
```

For runtime-only use:

```bash
python -m pip install -e .
```

## Basic CLI Usage

```bash
active-coronagraph --function main_func --name dual_zone --charge None
```

Run a 20% broadband PSF workflow centered on the reference wavelength:

```bash
active-coronagraph --function broadband_combine --name vortex --charge 2
```

The broadband workflow uses wavelength ratios `[0.9, 0.95, 1.0, 1.05, 1.1]`
by default and combines the wavelength-specific PSFs with endpoint-half
trapezoid weights.

Generated simulation outputs are written under the current working directory by
default.

## Examples

Run the minimal FQPM example:

```bash
python examples/minimal_simulation.py --save-figure
```

Run the broadband vortex example and save a diagnostic image:

```bash
cd examples
python test_example.py
```

This writes `test_example_broadband_vortex_image.png` and simulation products
under `examples/ideal_coro_2rd_mirror_False_broadband_example/`.

## User Documentation

For a short practical guide, see
[`docs/quick_user_instructions.md`](docs/quick_user_instructions.md). It lists
the available masks, main functions, and copy-paste Python examples.

For the longer reference guide, see [`docs/user_guide.md`](docs/user_guide.md).

## Runtime Configuration

Paths can be configured with environment variables:

- `ACTIVE_CORONAGRAPH_DATA_DIR`: directory containing `binary_mask/` and `phase_screen/`.
- `ACTIVE_CORONAGRAPH_OUTPUT_DIR`: base directory for generated simulation outputs.

The CLI also accepts:

```bash
active-coronagraph --function compare_phase_screens --data-dir data
```

## Tests

```bash
pytest
```

## License

This project is distributed under the BSD 3-Clause License. See `LICENSE`.

## Before Public Release

Expand `data/README.md` with provenance, license, and citation information for
the bundled FITS data before public release.
