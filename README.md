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
python main_polished.py --function main_func
```
