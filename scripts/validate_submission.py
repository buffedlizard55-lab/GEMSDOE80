#!/usr/bin/env python3
"""
Validate a single-band GeoTIFF against the competition's sample raster.

Format-safety check only. It does not train a model, score a raster against
labels, prove novelty, or certify scientific merit.

Usage:
    python scripts/validate_submission.py PREDICTION.tif --template OFFICIAL_SAMPLE.tif [--labels LABELS.tif] [--features FEATURES.tif]

Exit codes: 0 = format OK, 1 = format problems, 2 = missing dependencies/files.

JSON report is printed to stdout and written next to the prediction as
PREDICTION.validate.json. Report fields follow the project run-card convention:
format_valid, download_ready, values_in_0_1, nan_inside_footprint, No scientific evaluation is performed here; use the shared evaluator.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]

def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def build_report(prediction: Path, template: Path):
    report = {
        "format_valid": False,
        "download_ready": False,
        "organizer_confirmed": False,
        "prediction": str(prediction),
        "template": str(template),
        "checks": {},
        "errors": [],
        "value_range": None,
        "nan_inside_footprint": None,
        "crs": None,
        "shape": None,
        "transform": None,
        "sha256": None,
        "holdout": {"available": False},
        "leakage_canary": {"available": False},
        "lane_check": {
            "available": False,
            "note": "no registry raster provided with public competition data",
        },
    }
    errors = report["errors"]
    checks = report["checks"]

    for path, label in ((prediction, "prediction"), (template, "template")):
        if not path.is_file():
            errors.append(f"{label} file does not exist: {path}")
    if errors:
        return report

    try:
        import numpy as np  # noqa: F401  (used below via module-level import)
        import rasterio
    except ImportError as exc:
        errors.append(f"Missing dependency: {exc}")
        errors.append("Install rasterio and numpy in the project environment.")
        return report

    try:
        with rasterio.open(prediction) as pred, rasterio.open(template) as tmpl:
            report["crs"] = str(pred.crs)
            report["shape"] = [pred.height, pred.width]
            report["transform"] = [
                pred.transform.a,
                pred.transform.b,
                pred.transform.c,
                pred.transform.d,
                pred.transform.e,
                pred.transform.f,
            ]
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

            values = pred.read(1).astype(np.float64)
            pred_mask = pred.read_masks(1) > 0
            template_mask = tmpl.read_masks(1) > 0
            valid = pred_mask & np.isfinite(values)

            checks["finite_values"] = bool(np.isfinite(values[valid]).all()) if valid.any() else True
            checks["no_nan_or_infinite_inside_footprint"] = bool(
                np.isfinite(values[pred_mask]).all()
            ) if pred_mask.any() else True
            in_range = bool(((values[valid] >= 0.0) & (values[valid] <= 1.0)).all()) if valid.any() else False
            report["values_in_0_1"] = in_range
            report["value_range"] = (
                [float(values[valid].min()), float(values[valid].max())] if valid.any() else None
            )
            if not in_range:
                errors.append("Prediction values must be in range [0, 1]")

            nan_inside = int((~valid & pred_mask).sum())  # NaN where prediction claims data
            report["nan_inside_footprint"] = nan_inside
            if nan_inside:
                errors.append(f"{nan_inside} NaN pixels inside the prediction footprint")

            # The template, never the candidate's nodata mask, defines required coverage.
            if checks["shape_matches"]:
                required = template_mask
                holes = int((required & (~pred_mask | ~np.isfinite(values))).sum())
                checks["template_footprint_nonempty"] = bool(required.any())
                checks["template_footprint_covered"] = holes == 0
                report["nan_inside_footprint"] = int((required & ~np.isfinite(values)).sum())
                report["missing_inside_footprint"] = holes
                if not required.any():
                    errors.append("Template footprint is empty")
                if holes:
                    errors.append(f"{holes} missing/nonfinite pixels inside template footprint")
            else:
                checks["template_footprint_covered"] = False
                report["nan_inside_footprint"] = None
            # Masked finite sentinels must not hide the original range error.
            finite = np.isfinite(values)
            if np.any(finite & ((values < 0) | (values > 1))):
                errors.append("Finite values outside [0, 1], including masked sentinels")
                report["values_in_0_1"] = False
            if np.isinf(values).any():
                errors.append("Infinite values are not permitted")

            report["sha256"] = sha256_file(prediction)

    except Exception as exc:  # raster read failure
        errors.append(f"Failed to read rasters: {exc}")
        return report

    report["format_valid"] = not errors
    report["format_ready"] = report["format_valid"]
    report["download_ready"] = False  # Promotion is a separate, evidence-gated selector.
    report["promotion_blockers"] = [
        "Shared hide-and-recover HOLDOUT-DTI and CI unavailable",
        "Per-feature leakage canary unavailable",
        "Surface and final-dot registry comparisons unavailable",
        "Current holdout best and weekly cap unverified",
    ]
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("prediction", type=Path, help="Candidate single-band GeoTIFF")
    parser.add_argument("--template", required=True, type=Path, help="Official sample-submission GeoTIFF")
    parser.add_argument("--labels", type=Path, default=None, help="Optional label raster for holdout DTI")
    parser.add_argument("--features", type=Path, default=None, help="Optional feature raster for leakage canary")
    args = parser.parse_args(argv)

    try:
        import rasterio  # noqa: F401
    except ImportError:
        rasterio = None  # type: ignore

    report = build_report(args.prediction, args.template)

    if args.labels or args.features:
        report["holdout"]["note"] = (
            "Disabled: the former private evaluator was incorrect. Restore the shared "
            "evaluate_holdout.py and cached stack; refit each fold using visible faults only."
        )
        report["leakage_canary"]["note"] = "Not evaluated; no valid shared holdout available."

    out_path = args.prediction.with_name(args.prediction.name + ".validate.json")
    out_path.write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))

    return 0 if report["format_valid"] else 1


if __name__ == "__main__":
    sys.exit(main())
