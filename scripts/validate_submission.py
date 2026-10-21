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
format_valid, download_ready, values_in_0_1, nan_inside_footprint, plus
holdout DTI / leakage canary / lane check when labels/features/registry rasters
are available.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]

ALPHA, BETA = 0.2, 0.8
KERNEL_R_M = 300.0
RESOLUTION_M = 100.0
BOOTSTRAP_SEED = 20261021
N_BOOTSTRAP = 2000


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def auc_from_ranks(scores, labels):
    """Mann-Whitney U AUC with average ranks for ties (no sklearn required)."""
    scores = np.asarray(scores, dtype=np.float64)
    labels = np.asarray(labels, dtype=np.int64)
    pos = labels == 1
    neg = labels == 0
    n_pos = int(pos.sum())
    n_neg = int(neg.sum())
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    order = np.argsort(scores, kind="mergesort")
    ranks = np.empty(scores.size, dtype=np.float64)
    sorted_scores = scores[order]
    i = 0
    while i < scores.size:
        j = i
        while j + 1 < scores.size and sorted_scores[j + 1] == sorted_scores[i]:
            j += 1
        avg_rank = 0.5 * ((i + 1) + (j + 1))
        ranks[order[i : j + 1]] = avg_rank
        i = j + 1
    u = float(ranks[pos].sum()) - n_pos * (n_pos + 1) / 2.0
    return u / (n_pos * n_neg)


def triangular_kernel(d, r):
    return np.clip(1.0 - d / r, 0.0, None)


def dti_weighted(pred, gt, resolution_m=RESOLUTION_M, alpha=ALPHA, beta=BETA):
    """Distance-weighted Tversky index per the competition definition."""
    R = int(math.ceil(KERNEL_R_M / resolution_m))
    p = np.where(np.isfinite(pred), pred, 0.0).astype(np.float64)
    g = (np.asarray(gt) > 0).astype(np.float64)

    offs = []
    for dy in range(-R, R + 1):
        for dx in range(-R, R + 1):
            d = math.hypot(dy, dx)
            if 0 < d <= R:
                offs.append((dy, dx, float(d)))
    k_by_off = [(dy, dx, triangular_kernel(np.array([d]), R)[0]) for dy, dx, d in offs]

    H, W = p.shape
    padded = np.pad(p, R, mode="constant", constant_values=0.0)
    local_max = np.full((H, W), -1.0, dtype=np.float64)
    for dy, dx, k in k_by_off:
        shifted = padded[R + dy : R + dy + H, R + dx : R + dx + W]
        local_max = np.maximum(local_max, shifted * k)

    try:
        from scipy.ndimage import distance_transform_edt

        dist_to_gt = distance_transform_edt((g == 0).astype(np.uint8))
    except Exception:
        dist_to_gt = np.full((H, W), np.inf)
    kmap = triangular_kernel(dist_to_gt, R)

    gt_mask = g > 0
    pred_pos = p > 0

    tp_w = float(np.where(gt_mask, local_max, 0.0).sum())
    fn_w = float((1.0 - np.where(gt_mask, local_max, 0.0)).clip(min=0.0).sum())
    fp_w = float((p[pred_pos] * (1.0 - kmap[pred_pos])).sum())

    dti = tp_w / (tp_w + alpha * fp_w + beta * fn_w + 1e-12)
    return {"TP_w": tp_w, "FP_w": fp_w, "FN_w": fn_w, "DTI": float(dti)}


def bootstrap_ci(scores, n_bootstrap=N_BOOTSTRAP, seed=BOOTSTRAP_SEED):
    rng = np.random.default_rng(seed)
    scores = np.asarray([s for s in scores if np.isfinite(s)], dtype=np.float64)
    if scores.size == 0:
        return (float("nan"), float("nan"))
    boots = np.empty(n_bootstrap)
    for i in range(n_bootstrap):
        boots[i] = rng.choice(scores, size=scores.size, replace=True).mean()
    return float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))


def holdout_segment_mask(gt, holdout_frac=1.0 / 3.0):
    """Whole fault segments withheld (connected components of the label raster)."""
    try:
        from scipy import ndimage

        labeled, n = ndimage.label(gt > 0)
    except Exception:
        return (gt > 0).copy()
    sizes = np.bincount(labeled.ravel())[1:]
    order = np.argsort(-sizes)
    n_hold = max(1, int(round(n * holdout_frac)))
    hold_labels = set(int(v) for v in order[:n_hold])
    return np.isin(labeled, list(hold_labels)) if hold_labels else np.zeros_like(gt, dtype=bool)


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

            footprint_overlap = bool((pred_mask & template_mask).any()) if pred_mask.any() else False
            checks["footprint_overlaps_template"] = footprint_overlap

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


def holdout_dti(report, prediction, template, labels_path, features_path):
    import rasterio

    with rasterio.open(prediction) as pred, rasterio.open(labels_path) as ls, rasterio.open(
        template
    ) as tmpl:
        values = pred.read(1).astype(np.float64)
        gt = ls.read(1).astype(np.float64)
    if gt.shape != values.shape:
        report["holdout"]["note"] = "labels grid does not match submission grid"
        return
    holdout = holdout_segment_mask(gt)
    gt_held = np.where(holdout, gt, 0.0)
    res = dti_weighted(values, gt_held)
    lo, hi = bootstrap_ci(np.array([res["DTI"]]))
    report["holdout"] = {
        "available": True,
        "evaluator_version": "GEMS distance-weighted Tversky alpha=0.2 beta=0.8 R=300m",
        "withheld_positive_count": int((gt_held > 0).sum()),
        "score": res["DTI"],
        "ci95": [lo, hi],
        "TP_w": res["TP_w"],
        "FP_w": res["FP_w"],
        "FN_w": res["FN_w"],
        "status": "COMPUTED (single-image scalar; CI approximate)",
    }

    # Leakage canary: each feature band vs held-out labels.
    with rasterio.open(features_path) as fs:
        aucs = []
        for b in range(1, fs.count + 1):
            band = fs.read(b).astype(np.float64)
            if band.shape != gt.shape:
                continue
            valid_b = np.isfinite(band) & holdout
            if not valid_b.any():
                continue
            desc = fs.tags(b).get("description", f"band {b}")
            aucs.append({"band": b, "description": desc, "auc": auc_from_ranks(band[valid_b], gt[valid_b].astype(int))})
    max_auc = max((a["auc"] for a in aucs if np.isfinite(a["auc"])), default=float("nan"))
    report["leakage_canary"] = {
        "available": True,
        "per_band_auc": aucs,
        "max_auc": float(max_auc),
        "leakage_flag": bool(max_auc > 0.90),
        "status": "COMPUTED",
        "note": "AUC > 0.90 on a single feature means leakage until proven otherwise",
    }


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

    if report["format_valid"] and args.labels and args.features and np is not None:
        if args.labels.is_file() and args.features.is_file():
            holdout_dti(report, args.prediction, args.template, args.labels, args.features)
        else:
            report["holdout"]["note"] = "labels or features file not found"
            report["leakage_canary"]["note"] = "labels or features file not found"

    out_path = args.prediction.with_name(args.prediction.name + ".validate.json")
    out_path.write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))

    return 0 if report["format_valid"] else 1


if __name__ == "__main__":
    sys.exit(main())
