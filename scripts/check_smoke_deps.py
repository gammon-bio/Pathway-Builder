#!/usr/bin/env python3
import importlib.util
import sys

MISSING = [m for m in ("numpy", "pandas") if importlib.util.find_spec(m) is None]
if MISSING:
    mods = ", ".join(MISSING)
    sys.exit(
        f"Missing Python dependency(ies): {mods}. Install extras via pip install -e '.[dev]' before running make toy-bulk."
    )
