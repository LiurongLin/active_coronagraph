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

The legacy source-tree entry point is still available:

```bash
python main_polished.py --function main_func --name dual_zone --charge None
```

Generated simulation outputs are written under the current working directory by
default.

## User Documentation

For a short practical guide, see
[`docs/quick_user_instructions.md`](docs/quick_user_instructions.md). It lists
the available masks, main functions, and copy-paste Python examples.

For the longer reference guide, see [`docs/user_guide.md`](docs/user_guide.md).

## Runtime Configuration

Paths can be configured with environment variables:

- `ACTIVE_CORONAGRAPH_DATA_DIR`: directory containing `binary_mask/` and `phase_screen/`.
- `ACTIVE_CORONAGRAPH_OUTPUT_DIR`: base directory for generated simulation outputs.
- `ACTIVE_CORONAGRAPH_PIXEL_NOISE_OUTPUT_DIR`: base directory for pixel-noise mean outputs.
- `ACTIVE_CORONAGRAPH_PHASE_SCREEN_FILE`: optional single phase-screen FITS cube for the legacy `add_phase_screen` workflow.

The CLI also accepts:

```bash
active-coronagraph --function compare_phase_screens --data-dir data
active-coronagraph --function add_phase_screen --phase-screen-file path/to/screens.fits
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
