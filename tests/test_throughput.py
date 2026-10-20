"""Keep throughput sampling and output conventions stable for supported masks."""

import numpy as np

from active_coronagraph import main_functions_polished as workflows
from active_coronagraph.new_mask import (
    directional_throughput_mono,
    throughput_1d_mono,
    throughput_map_mono,
)


def test_retained_throughput_sampling_and_saved_tables(monkeypatch, tmp_path):
    monkeypatch.setenv("ACTIVE_CORONAGRAPH_OUTPUT_DIR", str(tmp_path))
    calls = []

    def measure(*args, **kwargs):
        assert args[1] == "vortex"
        assert kwargs["charge"] == 2
        x, y = kwargs["offset"]
        calls.append((x, y))
        return 1.0 + x + 2.0 * y

    monkeypatch.setattr(workflows, "run", measure)
    options = dict(dim=4, fpm_sam=2, max_ld=1.0, step_ld=1.0)
    line = throughput_1d_mono("vortex", 2, **options)
    np.testing.assert_array_equal(np.loadtxt(line), [[0, 1], [1, 2]])
    assert calls == [(0, 0), (1, 0)]

    calls.clear()
    grid = throughput_map_mono("vortex", 2, **options)
    np.testing.assert_array_equal(np.loadtxt(grid), [[1, 2], [3, 4]])
    assert calls == [(0, 0), (1, 0), (0, 1), (1, 1)]
    assert "rows=y_lambda_over_D cols=x_lambda_over_D" in grid.read_text()

    calls.clear()
    directional = directional_throughput_mono("vortex", 2, **options)
    radii = np.array([-1.0, 0.0, 1.0])
    np.testing.assert_allclose(
        np.loadtxt(directional),
        np.column_stack([radii, 1 + radii, 1 + 2 * radii, 1 + 3 * radii / np.sqrt(2)]),
    )
    assert len(calls) == 9
    assert len(list(tmp_path.rglob("*.png"))) == 5


def test_user_mask_forwarded_to_all_throughput_workflows(monkeypatch, tmp_path):
    from astropy.io import fits
    import matplotlib.axes

    monkeypatch.setenv("ACTIVE_CORONAGRAPH_OUTPUT_DIR", str(tmp_path))
    phase = np.arange(16, dtype=float).reshape(4, 4) / 10
    path = tmp_path / "my_phase.fits"
    fits.writeto(path, phase)
    calls = []
    previews = []
    original_imshow = matplotlib.axes.Axes.imshow

    def imshow(ax, image, *args, **kwargs):
        previews.append(np.array(image))
        return original_imshow(ax, image, *args, **kwargs)

    def measure(*args, **kwargs):
        assert args[1] == "new_mask"
        assert kwargs["charge"] is None
        assert kwargs["mask_opt_params"] == {"phase_map_fits": str(path)}
        calls.append(kwargs["offset"])
        return 0.25

    monkeypatch.setattr(workflows, "run", measure)
    monkeypatch.setattr(matplotlib.axes.Axes, "imshow", imshow)
    for workflow in [throughput_1d_mono, throughput_map_mono, directional_throughput_mono]:
        output = workflow("new_mask", None, dim=2, fpm_sam=2,
                          max_ld=1, step_ld=1, phase_map_fits=path)
        assert output.is_file()
        assert output.name.startswith("new_mask_my_phase_")
    assert len(calls) == 15
    assert sum(np.array_equal(preview, phase) for preview in previews) == 2
