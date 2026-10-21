#!/usr/bin/env python3
"""
Inspect data/ and report whether the competition inputs are present.

The competition inputs are auth-gated on DrivenData (login required) and are
NOT downloaded automatically by this repo. See scripts/download_competition_data.sh
for the official source URLs and manual placement instructions.

Expected layout (after the user places the files manually):
  data/
    training_features.tif   multi-band features, EPSG:32611, 100 m
    labels.tif              fault label raster (training ground truth)
    sample_submission.tif   official submission template (grid to match exactly)
    1m_DEM_links.csv        links to 1 m DEM downloads

This script prints a grid report for sample_submission.tif if present so that
generate_submission.py can be re-run against the exact competition grid.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import rasterio

REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "data"

EXPECTED = ["training_features.tif", "labels.tif", "sample_submission.tif", "1m_DEM_links.csv"]


def main():
    print(f"data/ directory: {DATA}")
    if not DATA.exists():
        print("  (missing - create data/ and place competition inputs there)")
        return 1
    present = {p.name: p for p in DATA.iterdir() if p.is_file()}
    for name in EXPECTED:
        print(f"  {'OK ' if name in present else 'MISS'} {name}")
    if "sample_submission.tif" in present:
        p = present["sample_submission.tif"]
        with rasterio.open(p) as src:
            report = {
                "file": str(p),
                "count": src.count,
                "dtype": src.dtypes[0],
                "crs": str(src.crs),
                "height": src.height,
                "width": src.width,
                "transform": [src.transform.a, src.transform.b, src.transform.c,
                              src.transform.d, src.transform.e, src.transform.f],
                "resolution_x": float(src.transform.a),
                "resolution_y": abs(float(src.transform.e)),
                "nodata": src.nodata,
                "bounds": [src.bounds.left, src.bounds.bottom, src.bounds.right, src.bounds.top],
            }
        print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
