#!/usr/bin/env python3
"""Validate a single-band GeoTIFF against the competition's sample raster.

This checks file-format safety only. It does not train a model, score a raster,
prove novelty, or certify that its predictions have scientific merit.
Requires rasterio and numpy; install from the project's approved environment.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("prediction", type=Path, help="Candidate single-band GeoTIFF")
    parser.add_argument("--template", required=True, type=Path,
                        help="Official sample-submission GeoTIFF")
    args = parser.parse_args()

    try:
        import numpy as np
        import rasterio
    except ImportError as exc:
        print(json.dumps({"format_valid": False, "download_ready": False,
                          "error": f"Missing dependency: {exc}",
                          "install": "Install rasterio and numpy in the project environment."},
                         indent=2))
        return 2

    report: dict[str, object] = {
        "format_valid": False,
        "download_ready": False,
        "prediction": str(args.prediction),
        "template": str(args.template),
        "checks": {},
        "errors": [],
    }
    errors: list[str] = report["errors"]  # type: ignore[assignment]
    checks: dict[str, object] = report["checks"]  # type: ignore[assignment]

    for path, label in ((args.prediction, "prediction"), (args.template, "template")):
        if not path.is_file():
            errors.append(f"{label} file does not exist: {path}")
    if errors:
        print(json.dumps(report, indent=2))
        return 1

    try:
        with rasterio.open(args.prediction) as pred, rasterio.open(args.template) as tmpl:
            checks["driver"] = pred.driver
            checks["band_count"] = pred.count
            checks["template_band_count"] = tmpl.count
            checks["shape_matches"] = (pred.height, pred.width) == (tmpl.height, tmpl.width)
            checks["crs_matches"] = pred.crs == tmpl.crs
            checks["transform_matches"] = pred.transform == tmpl.transform
            checks["template_has_crs"] = tmpl.crs is not None
            checks["prediction_dtype"] = list(pred.dtypes)

            if pred.driver != "GTiff":
                errors.append(f"Prediction driver must be GTiff, got {pred.driver!r}")
            if pred.count != 1:
                errors.append(f"Prediction must have exactly one band, got {pred.count}")
            if tmpl.count < 1:
                errors.append("Template has no raster band")
            if not checks["shape_matches"]:
                errors.append("Prediction height/width do not match template")
            if not checks["crs_matches"]:
                errors.append("Prediction CRS does not exactly match template CRS")
            if not checks["transform_matches"]:
                errors.append("Prediction affine transform does not exactly match template")
            if not checks["template_has_crs"]:
                errors.append("Template has no CRS; cannot certify submission alignment")

            if pred.count == 1 and (pred.height, pred.width) == (tmpl.height, tmpl.width):
                values = pred.read(1)
                pred_mask = pred.read_masks(1) > 0
                template_mask = tmpl.read_masks(1) > 0
                checks["footprint_matches"] = bool(np.array_equal(pred_mask, template_mask))
                if not checks["footprint_matches"]:
                    errors.append("Prediction valid-data footprint does not match template")
                footprint = template_mask
                if not np.any(footprint):
                    errors.append("Template footprint is empty")
                else:
                    inside = values[footprint]
                    finite = bool(np.isfinite(inside).all())
                    in_range = bool(finite and (inside >= 0).all() and (inside <= 1).all())
                    checks["no_nan_or_infinite_inside_footprint"] = finite
                    checks["values_in_0_1_inside_footprint"] = in_range
                    checks["footprint_pixel_count"] = int(footprint.sum())
                    if not finite:
                        errors.append("Prediction contains NaN or infinite values inside template footprint")
                    elif not in_range:
                        errors.append("Prediction values inside template footprint must be in [0, 1]")
    except Exception as exc:  # rasterio errors include invalid/corrupt raster formats
        errors.append(f"Could not read or validate GeoTIFF: {type(exc).__name__}: {exc}")

    report["format_valid"] = not errors
    report["download_ready"] = False
    report["interpretation"] = (
        "FORMAT CHECK PASSED ONLY; scientific validation and holdout are still required."
        if not errors else "FORMAT CHECK FAILED; NOT READY TO DOWNLOAD OR SUBMIT."
    )
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if not errors else 1


if __name__ == "__main__":
    sys.exit(main())
