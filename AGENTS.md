# AGENTS.md

Guidance for coding agents and maintainers working in this repository.

## Project Purpose

`active-coronagraph` is scientific Python software for coronagraph simulations
and focal-plane mask studies. The code models pupil and Lyot-stop geometry,
generates focal-plane phase masks, runs coronagraph propagation workflows,
computes contrast/throughput products, and reads/writes FITS outputs used by
the simulation workflows.

Preserve scientific behavior unless the change is explicitly intended,
documented, and covered by tests.

## High-Level Architecture

- `src/active_coronagraph/`: installable Python package.
- `src/active_coronagraph/cli.py`: command-line workflows and legacy workflow
  orchestration.
- `src/active_coronagraph/config.py`: configurable runtime data and output
  paths.
- `src/active_coronagraph/main_functions_polished.py`: run-directory
  construction, simulation orchestration, broadband combination, and
  post-processing helpers.
- `src/active_coronagraph/coronagraphs_polished.py`: main coronagraph model,
  pupil/stop generation, focal-plane masks, ghost image, simple coronagraph
  runs, and residual-energy calculation.
- `src/active_coronagraph/coronagraphs.py`: older/legacy coronagraph routines
  retained for compatibility and reference.
- `src/active_coronagraph/phase_masks.py`: analytic phase-mask definitions.
- `src/active_coronagraph/new_mask.py`: double-vortex and throughput-map
  workflows, including binary FITS mask support.
- `src/active_coronagraph/basic.py`: numerical utilities such as FFT/IFFT,
  padding, radial masks, FITS writing, and coordinate helpers.
- `src/active_coronagraph/make_plot.py`: plotting and analysis helpers for
  generated simulation products.
- `data/`: bundled FITS/PNG runtime data.
- `examples/`: small examples.
- `tests/`: smoke tests for imports, configuration, and CLI availability.
- `docs/`: documentation scaffold.
- `main_polished.py`: backward-compatible wrapper around `active_coronagraph.cli`.

## Core Scientific Algorithm Areas

Treat these files as scientific core code:

- `src/active_coronagraph/coronagraphs_polished.py`
- `src/active_coronagraph/coronagraphs.py`
- `src/active_coronagraph/phase_masks.py`
- `src/active_coronagraph/new_mask.py`
- `src/active_coronagraph/basic.py`
- `src/active_coronagraph/main_functions_polished.py`

Changes here can alter numerical results. Do not refactor algorithms, array
indexing, normalization, phase wrapping, propagation order, sampling, or file
naming casually.

## Files and Data Not to Modify Casually

- `data/binary_mask/*.fits`
- `data/phase_screen/*.fits`
- `data/binary_mask/*.png`
- `src/active_coronagraph/coronagraphs_polished.py`
- `src/active_coronagraph/phase_masks.py`
- `src/active_coronagraph/new_mask.py`
- `src/active_coronagraph/basic.py`

Bundled FITS files are runtime scientific data. Do not regenerate, compress,
rename, or remove them without documenting provenance, expected shape/dtype, and
result changes.

`LICENSE` is currently a placeholder and must be replaced with the selected
open-source license before public release.

## Coding Conventions

- Use package-relative imports inside `src/active_coronagraph`, for example
  `from .config import data_dir`.
- Keep public CLI behavior backward compatible where practical.
- Keep `main_polished.py` as a lightweight compatibility wrapper only.
- Prefer `pathlib.Path` for filesystem paths.
- Do not add hard-coded local machine paths.
- Use `config.py` or CLI arguments for runtime data/output paths.
- Keep generated outputs outside package source. By default, outputs are rooted
  at `ACTIVE_CORONAGRAPH_OUTPUT_DIR` or the current working directory.
- Avoid broad style-only rewrites in scientific modules unless separately
  requested.
- If touching legacy code with wildcard imports, keep changes minimal and verify
  imports/tests.
- Do not change units, coordinate ordering, or sign conventions without explicit
  tests and documentation.

## Installation

Development install:

```bash
python -m pip install -e ".[dev]"
```

Runtime-only install:

```bash
python -m pip install -e .
```

The package requires Python 3.10 or newer. Runtime dependencies are listed in
`pyproject.toml`: `numpy`, `matplotlib`, `astropy`, `scipy`, `tqdm`, and
`hcipy`.

## Tests

Run the current test suite with:

```bash
pytest
```

When working from a source tree without installation, use:

```bash
PYTHONPATH=src pytest
```

Current tests are smoke tests. For scientific changes, add focused numerical
tests that compare shapes, units, normalization, saved FITS products, and
selected reference values.

## Basic Example

After installation, run:

```bash
active-coronagraph --function main_func --name dual_zone --charge None
```

The legacy wrapper also works from the source tree:

```bash
python main_polished.py --function main_func --name dual_zone --charge None
```

Generated simulation outputs are written under the current working directory by
default. To keep outputs elsewhere:

```bash
ACTIVE_CORONAGRAPH_OUTPUT_DIR=/path/to/output active-coronagraph --function main_func --name dual_zone --charge None
```

Runtime data can be redirected with:

```bash
ACTIVE_CORONAGRAPH_DATA_DIR=/path/to/data active-coronagraph --function compare_phase_screens
```

## Numerical Validation Expectations

For changes that may affect numerical behavior:

- Run the existing test suite.
- Add tests that would fail if the intended scientific behavior changed.
- Validate representative FITS output shapes and dtype.
- Validate normalization behavior for final focal-plane products.
- Validate throughput products for expected dimensions and coordinate order.
- Compare selected numerical outputs against a known reference before and after
  the change.
- Document intentional result changes in the PR/commit notes and, when
  appropriate, `CHANGELOG.md`.

Do not accept numerical changes based only on successful imports or CLI `--help`.

## Units and Coordinate Conventions

Preserve the existing conventions unless explicitly changing them with tests:

- Separation and throughput map coordinates are expressed in `lambda/D`.
- Focal-plane mask sampling uses `fpm_sam`.
- Numerical pupil/focal sampling uses `nsamp`.
- `phase_shift` is documented as `(DY, DX)` and can be specified in `array` or
  `fpm` units.
- Throughput offsets passed to `simple_coro(offset=...)` use `(x, y)` ordering;
  `new_mask.py` calls this out explicitly.
- FITS masks in `data/binary_mask` are expected by binary-mask workflows.
- Phase-screen workflows expect FITS cubes under `data/phase_screen` unless
  configured otherwise.

If a task requires changing any convention above, update code, tests, docs, and
examples together.
