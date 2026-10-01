from pathlib import Path

from active_coronagraph import config


def test_default_data_dirs_exist():
    assert (config.data_dir() / "binary_mask").is_dir()
    assert (config.data_dir() / "phase_screen").is_dir()


def test_output_dir_env_override(monkeypatch, tmp_path):
    monkeypatch.setenv("ACTIVE_CORONAGRAPH_OUTPUT_DIR", str(tmp_path))
    assert config.output_dir() == Path(tmp_path).resolve()
