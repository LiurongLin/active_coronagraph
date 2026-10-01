"""Backward-compatible entry point for the active_coronagraph CLI."""

import sys
from pathlib import Path

src_dir = Path(__file__).resolve().parent / "src"
if src_dir.exists():
    sys.path.insert(0, str(src_dir))

from active_coronagraph.cli import main


if __name__ == "__main__":
    main()
