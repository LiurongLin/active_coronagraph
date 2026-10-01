import subprocess
import sys

import numpy as np

from active_coronagraph.cli import resolve_wavelengths, trapezoid_weights


def test_cli_help_runs():
    result = subprocess.run(
        [sys.executable, "-m", "active_coronagraph.cli", "--help"],
        check=True,
        text=True,
        capture_output=True,
    )
    assert "--function" in result.stdout


def test_broadband_workflows_default_to_20_percent_band():
    assert resolve_wavelengths("broadband_combine", None) == [0.9, 0.95, 1.0, 1.05, 1.1]


def test_monochromatic_workflows_default_to_reference_wavelength():
    assert resolve_wavelengths("main_func", None) == [1.0]


def test_user_wavelengths_override_workflow_defaults():
    assert resolve_wavelengths("broadband_combine", [0.8, 1.0, 1.2]) == [0.8, 1.0, 1.2]


def test_trapezoid_weights_half_weight_endpoints():
    np.testing.assert_allclose(
        trapezoid_weights([0.9, 0.95, 1.0, 1.05, 1.1]),
        np.array([0.5, 1.0, 1.0, 1.0, 0.5]),
    )
