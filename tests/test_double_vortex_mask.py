import numpy as np
import pytest

from active_coronagraph.new_mask import double_vortex_phase_mask


def test_double_vortex_phase_mask_has_expected_shape_and_range():
    phase = double_vortex_phase_mask(dim=6, fpm_sam=4, charge=2)

    assert phase.shape == (24, 24)
    assert np.all(np.isfinite(phase))
    assert np.min(phase) >= 0.0
    assert np.max(phase) < 2 * np.pi


def test_double_vortex_flip_lr_is_numpy_left_right_flip():
    kwargs = dict(
        dim=5,
        fpm_sam=3,
        charge=4,
        separation_ld=1.5,
        pa_deg=30.0,
        angle_deg=10.0,
        flux_ratio=0.7,
        sigma_ld=0.9,
    )
    unflipped = double_vortex_phase_mask(**kwargs, flip_lr=False)
    flipped = double_vortex_phase_mask(**kwargs, flip_lr=True)

    np.testing.assert_allclose(flipped, np.fliplr(unflipped))


def test_double_vortex_requires_charge():
    with pytest.raises(ValueError, match="requires a non-None charge"):
        double_vortex_phase_mask(dim=4, fpm_sam=4, charge=None)
