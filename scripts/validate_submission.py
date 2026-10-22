#!/usr/bin/env python3
"""Independent programmatic submission GeoTIFF validator for DrivenData GEMS.

Checks every requirement without human intervention:
  1. Single-band GeoTIFF (`count == 1`), `dtype == float32`
  2. Exact grid shape `(3730, 3292)`, CRS `EPSG:32611`, and exact `Affine` transform
     matching `sample_submission.tif`
  3. Zero `NaN` or `Inf` inside valid footprint (and finite across full raster when `-zeros.tif`)
  4. Strict `[0.0, 1.0]` range (`min >= 0.0`, `max <= 1.0` — never fails `[0, 1]` check)
  5. Zero positive predictions on `existing_faults.tif` known catalogue pixels
  6. Non-trivial sparse emission (`5,000 <= n_dots <= 120,000`), exact SHA-256 receipt
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import rasterio

ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / ".cache" / "data"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def validate_geotiff(
    sub_path: Path,
    template_path: Path = DATA_DIR / "sample_submission.tif",
    catalogue_path: Path = DATA_DIR / "existing_faults.tif",
) -> dict:
    checks: dict[str, bool] = {}
    notes: list[str] = []

    if not sub_path.is_file():
        return {"passed": False, "error": f"File not found: {sub_path}"}

    with rasterio.open(template_path) as tpl, rasterio.open(catalogue_path) as cat_ds, rasterio.open(sub_path) as ds:
        tpl_arr = tpl.read(1)
        valid = np.isfinite(tpl_arr) if tpl.nodata is None or np.isnan(tpl.nodata) else (tpl_arr != tpl.nodata)
        cat = (cat_ds.read(1) > 0) & valid

        checks["single_band"] = ds.count == 1
        checks["dtype_float32"] = ds.dtypes[0] == "float32"
        checks["shape_match"] = ds.shape == tpl.shape
        checks["crs_epsg_32611"] = (ds.crs == tpl.crs) and (ds.crs.to_epsg() == 32611)
        checks["transform_exact"] = ds.transform == tpl.transform

        arr = ds.read(1)
        in_foot = arr[valid]
        out_foot = arr[~valid]

        checks["no_nan_inside_footprint"] = bool(np.isfinite(in_foot).all())
        checks["all_pixels_finite"] = bool(np.isfinite(arr).all())
        min_in = float(np.nanmin(in_foot))
        max_in = float(np.nanmax(in_foot))
        min_all = float(np.nanmin(arr))
        max_all = float(np.nanmax(arr))

        checks["range_0_1_inside_footprint"] = (min_in >= 0.0) and (max_in <= 1.0)
        checks["range_0_1_entire_raster"] = (min_all >= 0.0) and (max_all <= 1.0)
        checks["outside_footprint_zero"] = bool(np.all(out_foot == 0.0))

        cat_leak_dots = int((arr[cat] > 0.0).sum())
        checks["zero_catalogue_dots"] = cat_leak_dots == 0

        n_pos = int((in_foot > 0.0).sum())
        n_ones = int((in_foot == 1.0).sum())
        checks["non_degenerate_emission"] = 5_000 <= n_pos <= 120_000

        if not checks["range_0_1_entire_raster"]:
            notes.append(f"Out-of-range values found: [{min_all}, {max_all}]")

    passed = all(checks.values())
    return {
        "evidence_class": "VALIDATOR-RECEIPT",
        "file_name": sub_path.name,
        "file_path": str(sub_path),
        "file_size_bytes": sub_path.stat().st_size,
        "file_size_mib": round(sub_path.stat().st_size / (1024 * 1024), 3),
        "sha256": sha256_file(sub_path),
        "shape": [int(arr.shape[0]), int(arr.shape[1])],
        "band_count": 1,
        "dtype": "float32",
        "crs": "EPSG:32611",
        "min_value": min_all,
        "max_value": max_all,
        "nan_count_inside_footprint": int((~np.isfinite(in_foot)).sum()),
        "nan_count_total": int((~np.isfinite(arr)).sum()),
        "positive_pixels": n_pos,
        "binary_unit_pixels": n_ones,
        "catalogue_overlap_pixels": cat_leak_dots,
        "checks": checks,
        "notes": notes,
        "passed": passed,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("submission", type=Path, help="Path to candidate .tif submission")
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON output receipt path")
    args = parser.parse_args()

    res = validate_geotiff(args.submission)
    if args.out is not None:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(res, indent=2) + "\n")
    print(json.dumps(res, indent=2))
    if not res["passed"]:
        sys.exit(1)


if __name__ == "__main__":
    main()
