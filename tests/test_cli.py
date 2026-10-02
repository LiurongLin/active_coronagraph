import subprocess
import sys
from pathlib import Path

import numpy as np

from active_coronagraph import cli
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


def test_main_func_restores_cwd_after_nested_annuli_mean(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(cli, "wavelength", [1.0])

    def fake_run(*args, **kwargs):
        return None

    def fake_annuli_mean(*args, **kwargs):
        nested = (
            tmp_path
            / "ideal_coro_2rd_mirror_False_noise_0_phase_screen_compare"
            / "no_screen"
            / "lyot_large_cross_FQPM_None_lambda_1.0_fpm_sam=10_binary_False_obstruction_False"
        )
        nested.mkdir(parents=True)
        # Reproduce the legacy helper behavior that caused the reported bug:
        # it enters a nested output path, then only climbs back two levels.
        monkeypatch.chdir(nested)
        monkeypatch.chdir(Path("../.."))
        return np.ones(3)

    monkeypatch.setattr(cli, "run", fake_run)
    monkeypatch.setattr(cli, "annuli_mean", fake_annuli_mean)

    cli.main_func("FQPM", None, folder="_phase_screen_compare/no_screen")

    assert Path.cwd() == tmp_path
    expected_plot = (
        tmp_path
        / "ideal_coro_2rd_mirror_False_noise_0_phase_screen_compare"
        / "no_screen"
        / "lyot_large_cross_FQPM_None_lambda_1.0_fpm_sam=10_binary_False_obstruction_False"
        / "final_focal_plane_coro_mask=FQPM_None_sample=10_fpm_sam=10_obstruction_False_wl_1.0_pupil_100_annuli_mean.png"
    )
    assert expected_plot.is_file()
