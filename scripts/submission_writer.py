#!/usr/bin/env python3
"""Write a legal GEMS Prize GeoTIFF pair (zeros-outside and NaN-outside).

Fixes DrivenData 'Predicted values must be in range [0, 1]':
  * zeros twin: nodata=None, every cell finite and in [0, 1]
  * nan twin: NaN strictly outside the official footprint, [0, 1] inside
  * CRS / shape / transform copied from sample_submission.tif
  * catalogue pixels forced to 0 (test labels are off-catalogue)
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import zipfile
from pathlib import Path

import numpy as np
import rasterio

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gems80.io_grid import load_labels, load_template, sha256_file  # noqa: E402


def write_pair(mask: np.ndarray, template: Path, labels_path: Path, out_dir: Path, stem: str) -> dict:
    tpl = load_template(template)
    foot = tpl["foot"]
    labels = load_labels(labels_path, foot)
    meta = tpl["meta"]
    clean = (np.asarray(mask, dtype=bool) & foot & ~labels).astype(np.float32)
    assert clean.shape == (meta["height"], meta["width"])

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    zeros_path = out_dir / f"{stem}-zeros.tif"
    nan_path = out_dir / f"{stem}-nan.tif"

    arr_zeros = np.where(foot, np.clip(clean, 0.0, 1.0), np.float32(0.0)).astype(np.float32)
    meta_z = meta.copy()
    meta_z["nodata"] = None
    with rasterio.open(zeros_path, "w", **meta_z) as dst:
        dst.write(arr_zeros, 1)
        dst.set_band_description(1, "fault_probability")

    arr_nan = np.where(foot, np.clip(clean, 0.0, 1.0), np.float32(np.nan)).astype(np.float32)
    meta_n = meta.copy()
    meta_n["nodata"] = np.nan
    with rasterio.open(nan_path, "w", **meta_n) as dst:
        dst.write(arr_nan, 1)
        dst.set_band_description(1, "fault_probability")

    zip_path = zeros_path.with_suffix(".zip")
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        zf.write(zeros_path, arcname=zeros_path.name)

    def audit(path: Path, mode: str) -> dict:
        with rasterio.open(path) as src:
            a = src.read(1)
            nodata = src.nodata
            shape = [src.height, src.width]
            transform = [float(src.transform.a), float(src.transform.b), float(src.transform.c),
                         float(src.transform.d), float(src.transform.e), float(src.transform.f)]
            crs = str(src.crs)
            dtype = src.dtypes[0]
            count = src.count
        in_f = a[foot]
        out_f = a[~foot]
        finite_in = bool(np.isfinite(in_f).all())
        in_range = bool(finite_in and ((in_f >= 0) & (in_f <= 1)).all())
        if mode == "zeros":
            outside_ok = bool(np.isfinite(a).all() and (out_f == 0).all())
            full_01 = bool(np.isfinite(a).all() and ((a >= 0) & (a <= 1)).all())
        else:
            outside_ok = bool(np.isnan(out_f).all())
            full_01 = in_range
        return {
            "path": str(path),
            "sha256": sha256_file(path),
            "bytes": path.stat().st_size,
            "mode": mode,
            "crs": crs,
            "shape": shape,
            "transform": transform,
            "dtype": dtype,
            "count": count,
            "nodata": None if nodata is None or (isinstance(nodata, float) and np.isnan(nodata)) else nodata,
            "positive_px": int((np.isfinite(a) & (a > 0.5) & foot).sum()),
            "on_catalogue_px": int((np.isfinite(a) & (a > 0.5) & labels).sum()),
            "nan_inside_footprint": int((~np.isfinite(in_f)).sum()),
            "values_in_0_1_inside": in_range,
            "outside_compliant": outside_ok,
            "full_grid_in_0_1": full_01,
            "format_valid": bool(count == 1 and dtype == "float32" and "32611" in crs
                                 and shape == [meta["height"], meta["width"]]
                                 and in_range and finite_in and outside_ok
                                 and int((np.isfinite(a) & (a > 0.5) & labels).sum()) == 0),
        }

    bundle = {
        "stem": stem,
        "zeros": audit(zeros_path, "zeros"),
        "nan": audit(nan_path, "nan"),
        "zip_sha256": sha256_file(zip_path),
        "zip_path": str(zip_path),
        "template_report": tpl["report"],
    }
    sidecar = out_dir / f"{stem}.audit.json"
    sidecar.write_text(json.dumps(bundle, indent=2) + "\n")
    return bundle


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mask-tif", type=Path, required=True, help="Any single-band raster; >0.5 treated as emit")
    ap.add_argument("--template", type=Path, default=ROOT / "data" / "sample_submission.tif")
    ap.add_argument("--labels", type=Path, default=ROOT / "data" / "labels.tif")
    ap.add_argument("--out-dir", type=Path, default=ROOT / "docs" / "downloads")
    ap.add_argument("--stem", required=True)
    args = ap.parse_args(argv)
    with rasterio.open(args.mask_tif) as src:
        mask = src.read(1) > 0.5
    bundle = write_pair(mask, args.template, args.labels, args.out_dir, args.stem)
    print(json.dumps(bundle, indent=2))
    return 0 if bundle["zeros"]["format_valid"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
