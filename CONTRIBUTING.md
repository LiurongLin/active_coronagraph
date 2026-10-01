# Contributing

Please open an issue before large scientific or algorithmic changes.

For development:

```bash
python -m pip install -e ".[dev]"
pytest
```

Keep numerical algorithm changes separate from packaging, documentation, and
style-only changes so they can be reviewed independently.
