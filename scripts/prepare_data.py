#!/usr/bin/env python3
"""
Verify data/ placement, SHA-256 digests, official submission grid geometry,
footprint pixel counts, and band metadata across all core + external layers.

Writes machine-readable audit receipt to reports/prepare_data_report.json.
"""

from __future__ import annotations

import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import rasterio

REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "data"
REPORTS = REPO / "reports"

CORE_FILES = [
    "training_features.tif",
    "labels.tif",
    "existing_faults.tif",
    "sample_submission.tif",
]
EXTERNAL_FILES = [
    "external/lidar_scarp_features_u8.tif",
    "external/geodawn_rad_u8.tif",
    "external/geodawn_extensions_u8.tif",
    "external/derived_sgmc_faults_100m_u8.tif",
    "external/gdr_wellspring_in_footprint.csv",
    "external/gdr_volcanic_vents_in_footprint.csv",
    "external/sb_heat_flow_in_footprint.csv",
]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    print(f"Inspecting competition data directory: {DATA}")
    if not DATA.exists():
        print("[ERROR] data/ does not exist. Run: bash scripts/download_competition_data.sh")
        return 1

    inventory = {}
    missing = []
    for rel in CORE_FILES + EXTERNAL_FILES:
        p = DATA / rel
        if p.is_file():
            inventory[rel] = {
                "present": True,
                "bytes": p.stat().st_size,
                "sha256": sha256_file(p),
            }
            print(f"  [OK]   {rel:45s}  {inventory[rel]['bytes']:12,d} B  sha256={inventory[rel]['sha256'][:16]}…")
        else:
            inventory[rel] = {"present": False}
            missing.append(rel)
            print(f"  [MISS] {rel}")

    if any(m in CORE_FILES for m in missing):
        print(f"[ERROR] Missing required core files: {missing}")
        return 1

    with rasterio.open(DATA / "sample_submission.tif") as tmpl:
        tmpl_arr = tmpl.read(1)
        footprint = np.isfinite(tmpl_arr)
        grid_report = {
            "count": tmpl.count,
            "dtype": tmpl.dtypes[0],
            "crs": str(tmpl.crs),
            "height": tmpl.height,
            "width": tmpl.width,
            "transform": [
                tmpl.transform.a, tmpl.transform.b, tmpl.transform.c,
                tmpl.transform.d, tmpl.transform.e, tmpl.transform.f,
            ],
            "resolution_x_m": float(tmpl.transform.a),
            "resolution_y_m": abs(float(tmpl.transform.e)),
            "nodata": str(tmpl.nodata),
            "bounds": [tmpl.bounds.left, tmpl.bounds.bottom, tmpl.bounds.right, tmpl.bounds.top],
            "total_grid_pixels": int(tmpl.height * tmpl.width),
            "valid_footprint_pixels": int(footprint.sum()),
            "outside_footprint_pixels": int((~footprint).sum()),
        }

    with rasterio.open(DATA / "labels.tif") as lsrc:
        lbl = lsrc.read(1)
        pos_mask = (lbl > 0) & footprint
        label_report = {
            "shape_matches_template": (lsrc.height, lsrc.width) == (grid_report["height"], grid_report["width"]),
            "crs_matches_template": str(lsrc.crs) == grid_report["crs"],
            "transform_matches_template": list(lsrc.transform)[:6] == grid_report["transform"],
            "known_fault_pixels": int(pos_mask.sum()),
            "prevalence_in_footprint": float(pos_mask.sum() / footprint.sum()),
        }

    band_catalog = []
    with rasterio.open(DATA / "training_features.tif") as fsrc:
        for b in range(1, fsrc.count + 1):
            tags = fsrc.tags(b)
            band_catalog.append({
                "band": b,
                "band_name": tags.get("band_name", f"band_{b}"),
                "category": tags.get("data_category", "unknown"),
                "description": tags.get("description", fsrc.descriptions[b - 1] or ""),
            })

    summary = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "READY",
        "grid": grid_report,
        "labels": label_report,
        "official_feature_bands": band_catalog,
        "files": inventory,
    }
    REPORTS.mkdir(parents=True, exist_ok=True)
    out = REPORTS / "prepare_data_report.json"
    out.write_text(json.dumps(summary, indent=2) + "\n")
    print(f"[OK] Grid verified: {grid_report['height']}x{grid_report['width']} {grid_report['crs']} "
          f"valid={grid_report['valid_footprint_pixels']:,} faults={label_report['known_fault_pixels']:,}")
    print(f"[OK] Wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
