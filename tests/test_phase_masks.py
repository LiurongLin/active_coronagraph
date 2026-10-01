import numpy as np

from active_coronagraph.phase_masks import Phase_masks


def test_dual_zone_assigns_documented_radial_phase_regions():
    masks = Phase_masks(fpm_sam=10, pupil_size=4)
    phase = masks.dual_zone(a1=0.5, a2=1.0, z1=0.25, z2=0.75)

    assert phase.shape == (40, 40)
    np.testing.assert_allclose(phase[20, 20], 0.5 * np.pi)
    np.testing.assert_allclose(phase[20, 27], 1.5 * np.pi)
    np.testing.assert_allclose(phase[0, 0], 0.0)


def test_fqpm_has_opposite_quadrants_at_pi_and_zero():
    masks = Phase_masks(fpm_sam=2, pupil_size=2)
    phase = masks.FQPM()

    expected = np.array(
        [
            [np.pi, np.pi, 0.0, 0.0],
            [np.pi, np.pi, 0.0, 0.0],
            [0.0, 0.0, np.pi, np.pi],
            [0.0, 0.0, np.pi, np.pi],
        ]
    )
    np.testing.assert_allclose(phase, expected)
