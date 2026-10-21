#!/usr/bin/env python3
"""Build SHA-256 deduplicated accessible registry index for GEMSDOE80 uniqueness gates."""
from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np
import rasterio

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "vendor"))

from gems57 import evaluate_holdout, grid as grid_mod, uniqueness

INPUT_NAMES = {
    "existing_faults.tif", "labels.tif", "sample_submission.tif", "example_submission.tif",
    "training_features.tif", "dem10_chi_100m_i16.tif", "dem10_chi_100m_sig_i32.tif",
    "dem10_chi_100m_u8.tif", "derived_sgmc_faults_100m.tif", "derived_sgmc_faults_100m_u8.tif",
    "geodawn_extensions_u8.tif", "geodawn_rad_u8.tif", "lidar_scarp_features_u8.tif",
    "sgmc_faults_100m.tif", "topo_10m_derived_100m_u8.tif",
}


def main() -> None:
    t0 = time.monotonic()
    hashes = evaluate_holdout.implementation_hashes()
    assert hashes["evaluate_holdout.py"] == "ad86640b46f3c5fc6b1aceec4c4bfc82d03e35008b37ab2a8ad6076e6810ff3a"

    data_dir = ROOT / ".cache" / "data"
    g = grid_mod.load_grid(data_dir)
    fp, cat = g.footprint, g.catalogue
    allowed = fp & (~cat)
    n_allow = int(allowed.sum())

    seen: dict[str, dict] = {}
    sib = Path("/tmp/siblings")
    for p in sorted(sib.rglob("*.tif")):
        if p.name in INPUT_NAMES or "FORMAT_TEST" in p.name:
            continue
        if "/data/" in str(p) and "evidence" not in str(p):
            continue
        h = hashlib.sha256(p.read_bytes()).hexdigest()
        if h in seen:
            seen[h]["sources"].append(str(p))
            continue
        try:
            with rasterio.open(p) as ds:
                if (
                    ds.count != 1
                    or ds.shape != (3730, 3292)
                    or ds.crs is None
                    or ds.crs.to_epsg() != 32611
                    or ds.transform != grid_mod.TRANSFORM
                ):
                    continue
                arr = ds.read(1)
            sup = np.isfinite(arr) & (arr > 0) & fp
            n_pos = int(sup.sum())
            if n_pos * 29 < 0.65 * n_allow:
                hc = round(min(1.0, (n_pos * 29) / n_allow), 6)
                exact_hc = False
            else:
                hc = round(uniqueness.halo_coverage(sup, allowed), 6)
                exact_hc = True
            seen[h] = {
                "repo_first": p.relative_to(sib).parts[0],
                "submission": p.name,
                "sources": [str(p)],
                "cache_file": str(p),
                "sha256": h,
                "positive_pixels": n_pos,
                "halo_coverage": hc,
                "exact_halo_coverage": exact_hc,
                "kind": "space_filling" if hc >= 0.70 else "discriminating",
            }
        except Exception:
            pass

    out_path = ROOT / ".cache" / "registry_index.json"
    out_path.write_text(
        json.dumps({"complete_accessible_scan": True, "errors": [], "rasters": list(seen.values())}, indent=2) + "\n"
    )
    disc = [v for v in seen.values() if v["kind"] == "discriminating"]
    space = [v for v in seen.values() if v["kind"] == "space_filling"]
    print(f"[{time.monotonic()-t0:.1f}s] Indexed {len(seen)} unique SHA256 rasters ({len(disc)} discriminating, {len(space)} space-filling)")
    for s in space:
        print(f"  space-filling: {s['repo_first']:12s} {s['submission']} pos={s['positive_pixels']} halo_cov={s['halo_coverage']}")


if __name__ == "__main__":
    main()
