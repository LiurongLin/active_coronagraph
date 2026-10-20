"""User-supplied masks use synthetic data and the shared propagation pipeline."""

import numpy as np
import pytest
from astropy.io import fits

from active_coronagraph.coronagraphs_polished import focal_mask_new, simple_coro
from active_coronagraph.new_mask import load_phase_map
from active_coronagraph.phase_masks import Phase_masks


def test_user_phase_map_preserves_values_and_coordinate_order(tmp_path):
    phase = np.arange(16, dtype=np.float32).reshape(4, 4) / 10
    path = tmp_path / "user_phase.fits"
    fits.writeto(path, phase)
    loaded = load_phase_map(path, dim=2, fpm_sam=2)
    np.testing.assert_array_equal(loaded, phase)
    assert loaded.dtype == np.dtype(float)


@pytest.mark.parametrize("phase", [np.zeros((3, 4)), np.zeros((1, 4, 4)), np.full((4, 4), np.nan)])
def test_user_phase_map_rejects_invalid_data(tmp_path, phase):
    path = tmp_path / "invalid_phase.fits"
    fits.writeto(path, phase)
    with pytest.raises(ValueError, match="new_mask"):
        load_phase_map(path, dim=2, fpm_sam=2)


def test_user_phase_map_requires_explicit_file():
    with pytest.raises(ValueError, match="requires"):
        load_phase_map(None, dim=2, fpm_sam=2)


@pytest.mark.parametrize("wavelength", [0.9, 1.0])
def test_user_mask_matches_analytic_reference_and_saved_products(monkeypatch, tmp_path, wavelength):
    monkeypatch.chdir(tmp_path)
    phase = Phase_masks(fpm_sam=2, pupil_size=10).FQPM()
    path = tmp_path / "user_phase.fits"
    fits.writeto(path, phase)
    options = {"phase_map_fits": str(path)}
    mask_args = dict(dim=10, nsamp=4, fpm_sam=2, wavelength=wavelength, charge=None)
    loaded_mask = focal_mask_new("new_mask", opt_params=options, **mask_args)
    analytic_mask = focal_mask_new("FQPM", **mask_args)
    np.testing.assert_array_equal(loaded_mask, analytic_mask)
    assert loaded_mask.shape == (40, 40)
    np.testing.assert_allclose(np.abs(loaded_mask), 1.0)

    run_args = dict(**mask_args, lyot_stop="None_1_1", binary=False,
                    obstruction=False, get_lyot=False, greyscale=None)
    for offset in [None, (1.0, 0.0), (0.0, 1.0)]:
        loaded = simple_coro(name="new_mask", mask_opt_params=options, offset=offset, **run_args)
        reference = simple_coro(name="FQPM", offset=offset, **run_args)
        np.testing.assert_array_equal(loaded, reference)
        if offset is None:
            saved = fits.getdata(next(tmp_path.glob("final_focal_plane_coro_mask=new_mask*.fits")))
            assert saved.shape == (1, 40, 40)
            assert saved.dtype.kind == "f"
            np.testing.assert_array_equal(saved[0], loaded)
