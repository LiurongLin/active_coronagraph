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
