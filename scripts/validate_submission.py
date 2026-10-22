#!/usr/bin/env python3
"""
Validate a single-band GeoTIFF against the competition's sample raster.

Format-safety check only. It does not train a model, score a raster against
labels, prove novelty, or certify scientific merit.

Usage:
    python scripts/validate_submission.py PREDICTION.tif --template SAMPLE_TEMPLATE.tif

Exit codes: 0 = format OK, 1 = format problems, 2 = missing dependencies/files.

JSON report is printed to stdout and written next to the prediction as
PREDICTION.validate.json. This tool checks format only. It deliberately does not
calculate HOLDOUT-DTI or leakage canaries; those require a hide-and-rebuild
experiment and are run by scripts/run_h80_experiments.py.
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
        "readiness_scope": (
            "Format checks apply only against the caller-supplied template. This report does not "
            "authenticate the template or guarantee organizer acceptance."
        ),
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
            if tmpl.count != 1:
                errors.append(f"Template must have exactly one band, got {tmpl.count}")
            if pred.dtypes[0] != "float32":
                errors.append(f"Prediction must use float32, got {pred.dtypes[0]!r}")
            if not checks["shape_matches"]:
                errors.append("Prediction height/width do not match template")
            if not checks["crs_matches"]:
                errors.append("Prediction CRS does not exactly match template CRS")
            if not checks["transform_matches"]:
                errors.append("Prediction affine transform does not exactly match template")
            if not checks["template_has_crs"]:
                errors.append("Template has no CRS; cannot certify submission alignment")

            # Array comparisons are only meaningful when dimensions agree. Fail closed
            # rather than broadcasting masks or returning a partial report on mismatch.
            if not checks["shape_matches"] or pred.count != 1 or tmpl.count != 1:
                report["sha256"] = sha256_file(prediction)
                report["format_valid"] = False
                report["submit_format_ok"] = False
                report["download_ready"] = False
                report["format_scope"] = (
                    "Checks use the caller-supplied template; they do not establish its official provenance "
                    "or guarantee organizer acceptance."
                )
                return report

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

            nan_inside = int((~np.isfinite(values) & template_mask).sum())
            report["nan_inside_footprint"] = nan_inside
            report["missing_inside_footprint"] = nan_inside
            if nan_inside:
                errors.append(f"{nan_inside} NaN pixels inside the template footprint")
            checks["no_nan_or_infinite_inside_footprint"] = nan_inside == 0
            # Finite values anywhere (including a numeric nodata sentinel) must be in [0, 1].
            finite_all = np.isfinite(values)
            if finite_all.any() and not bool(((values[finite_all] >= 0.0) & (values[finite_all] <= 1.0)).all()):
                in_range = False
                report["values_in_0_1"] = False
                if "Prediction values must be in range [0, 1]" not in errors:
                    errors.append("Prediction values must be in range [0, 1]")

            footprint_overlap = bool((pred_mask & template_mask).any()) if pred_mask.any() else False
            template_footprint_covered = bool(np.all(pred_mask[template_mask])) if template_mask.any() else False
            checks["footprint_overlaps_template"] = footprint_overlap
            checks["template_footprint_covered"] = template_footprint_covered
            if not template_footprint_covered:
                errors.append("Prediction mask does not cover the complete template footprint")

            report["sha256"] = sha256_file(prediction)

    except Exception as exc:  # raster read failure
        errors.append(f"Failed to read rasters: {exc}")
        return report

    report["format_valid"] = not errors
    # Format-valid files are OK to upload (will not trip "Predicted values must be in range [0, 1]").
    # Organizer confirmation is a separate field and is never self-asserted.
    report["submit_format_ok"] = report["format_valid"]
    report["download_ready"] = report["format_valid"]
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("prediction", type=Path, help="Candidate single-band GeoTIFF")
    parser.add_argument("--template", required=True, type=Path, help="Official sample-submission GeoTIFF")
    args = parser.parse_args(argv)

    report = build_report(args.prediction, args.template)

    out_path = args.prediction.with_name(args.prediction.name + ".validate.json")
    out_path.write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))

    return 0 if report["format_valid"] else 1


if __name__ == "__main__":
    sys.exit(main())
