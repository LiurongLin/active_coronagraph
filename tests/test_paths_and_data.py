from pathlib import Path

import numpy as np
from astropy.io import fits

from active_coronagraph.main_functions_polished import _root_folder, _run_dir


def test_root_folder_respects_output_dir_environment(monkeypatch, tmp_path):
    monkeypatch.setenv("ACTIVE_CORONAGRAPH_OUTPUT_DIR", str(tmp_path))

    root = _root_folder(obstruction=True, noise_level=0.25, folder_name="_case")

    assert root == tmp_path / "ideal_coro_2rd_mirror_True_noise_0.25_case"


def test_run_dir_preserves_phase_shift_suffix_and_offset_mode():
    run_dir = _run_dir(
        lyot_stop="large_cross",
        name="vortex",
        charge=2,
        wavelength=1.0,
        nsamp=10,
        fpm_sam=10,
        binary=False,
        obstruction=True,
        greyscale=8,
        offset=(1.0, 0.0),
        phase_shift=(0.5, -0.25),
        phase_shift_unit="fpm",
    )

    assert run_dir == Path(
        "lyot_large_cross_vortex_2_"
        "sam_10_fpm_sam=10_binary_False_obstruction_True_greyscale_8_offset_True"
        "_phase_shift_0.5_-0.25_fpm"
    )


def test_bundled_binary_mask_fits_files_have_expected_shape_and_units():
    for path in sorted(Path("data/binary_mask").glob("*.fits")):
        with fits.open(path, memmap=False) as hdul:
            data = hdul[0].data
            header = hdul[0].header

        assert data.shape == (1000, 1000)
        assert np.issubdtype(data.dtype, np.floating)
        assert header["BUNIT"] == "rad"


def test_bundled_phase_screen_cubes_have_expected_shape():
    for path in sorted(Path("data/phase_screen").glob("*.fits")):
        with fits.open(path, memmap=False) as hdul:
            data = hdul[0].data

        assert data.shape == (100, 100, 100)
        assert np.issubdtype(data.dtype, np.floating)


def test_source_tree_has_no_machine_specific_runtime_paths():
    scanned_suffixes = {".py", ".md", ".toml", ".cff"}
    forbidden = [
        "/" + "media" + "/" + "liurong",
        "My " + "Passport",
        "pixel_noise" + "_mean_output_dir",
        "2)" + "seeing=1,vmag=8,ZA=30,lag=2",
        "TROIA_phase_screens_new" + ".fits",
    ]

    for path in Path(".").rglob("*"):
        if any(part in {".git", ".venv", "__pycache__", ".pytest_cache"} for part in path.parts):
            continue
        if not path.is_file() or path.suffix not in scanned_suffixes:
            continue

        text = path.read_text(encoding="utf-8")
        for value in forbidden:
            assert value not in text, f"{value!r} found in {path}"
