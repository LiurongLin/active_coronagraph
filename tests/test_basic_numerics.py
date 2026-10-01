import numpy as np
import pytest

from active_coronagraph.basic import Isum, pad_array, phase_2pi, phi_ramp, r_theta


def test_phase_2pi_wraps_complex_phase_to_nonnegative_interval():
    angles = np.array([-np.pi / 2, 0.0, np.pi / 2, np.pi])
    wrapped = phase_2pi(np.exp(1j * angles))

    np.testing.assert_allclose(
        wrapped,
        np.array([3 * np.pi / 2, 0.0, np.pi / 2, np.pi]),
        atol=1e-15,
    )
    assert np.all(wrapped >= 0)
    assert np.all(wrapped <= 2 * np.pi)


def test_pad_array_centers_input_and_preserves_values():
    array = np.array([[1, 2], [3, 4]])
    padded = pad_array(array, pad_size=6)

    assert padded.shape == (6, 6)
    np.testing.assert_array_equal(padded[2:4, 2:4], array)
    assert np.count_nonzero(padded) == np.count_nonzero(array)
    assert padded.sum() == array.sum()


def test_pad_array_rejects_arrays_larger_than_target():
    with pytest.raises(ValueError):
        pad_array(np.ones((3, 2)), pad_size=2)


def test_phi_ramp_matches_documented_linear_formula():
    image = np.zeros((3, 5))
    ramp = phi_ramp(image, npx=2, npy=3)

    lx = np.linspace(-0.5, 0.5, 5) * np.pi * 2 * 2
    ly = np.linspace(-0.5, 0.5, 3) * np.pi * 3 * 2
    expected_x, expected_y = np.meshgrid(lx, ly)
    np.testing.assert_allclose(ramp, expected_x + expected_y)


def test_r_theta_uses_y_x_grid_relative_to_center():
    image = np.zeros((3, 3))
    radius, theta = r_theta(image, xc=1, yc=1)

    np.testing.assert_allclose(radius[1, 1], 0.0)
    np.testing.assert_allclose(radius[1, 2], 1.0)
    np.testing.assert_allclose(theta[1, 2], 0.0)
    np.testing.assert_allclose(theta[2, 1], np.pi / 2)


def test_isum_returns_real_total_intensity():
    field = np.array([1 + 1j, 2 - 1j])
    assert Isum(field) == pytest.approx(7.0)
