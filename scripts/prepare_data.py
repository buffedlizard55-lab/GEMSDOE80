#!/usr/bin/env python3
"""Prepare hash-verified inputs, registry index, and 52-channel zero-catalogue feature matrix."""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from vendor.gems57.grid import load_grid
from scripts.restore_data import DATA_DIR, EVIDENCE_DIR, EXT_DIR, restore_all
from src.gems80.features import STATIC_FEATURE_NAMES, build_or_load_static_features


def main() -> None:
    t0 = time.time()
    receipt = restore_all()
    reg_path = ROOT / ".cache" / "registry_index.json"
    if not reg_path.is_file():
        import subprocess
        subprocess.run([sys.executable, str(ROOT / "scripts" / "build_registry_index.py")], check=True)

    grid = load_grid(DATA_DIR)
    cache_path = Path("/tmp/gems80_cache/x_static.npy")
    mm, _ = build_or_load_static_features(grid.footprint, DATA_DIR, EXT_DIR, cache_path)
    summary = {
        "evidence_class": "FEATURE-MATRIX-RECEIPT",
        "shape": list(mm.shape),
        "dtype": str(mm.dtype),
        "n_valid_cells": int(grid.footprint.sum()),
        "n_static_features": len(STATIC_FEATURE_NAMES),
        "feature_names": list(STATIC_FEATURE_NAMES),
        "elapsed_sec": round(time.time() - t0, 2),
        "core_sha256_verified": all(v["verified"] for v in receipt["core_files"].values()),
    }
    (EVIDENCE_DIR / "feature_matrix_receipt.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({"status": "OK", "shape": summary["shape"], "elapsed_sec": summary["elapsed_sec"]}, indent=2))


if __name__ == "__main__":
    main()
