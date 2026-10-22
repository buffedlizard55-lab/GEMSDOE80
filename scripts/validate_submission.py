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
sys.path.insert(0, str(REPO / "src"))

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
            if d <= R:  # include the pixel itself (k(0)=1)
                offs.append((dy, dx, float(d)))
    k_by_off = [(dy, dx, triangular_kernel(np.array([d]), R)[0]) for dy, dx, d in offs]

    H, W = p.shape
    padded = np.pad(p, R, mode="constant", constant_values=0.0)
    local_max = np.zeros((H, W), dtype=np.float64)
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

    tp_w = float(local_max[gt_mask].sum()) if gt_mask.any() else 0.0
    fn_w = float((1.0 - local_max[gt_mask]).clip(min=0.0).sum()) if gt_mask.any() else 0.0
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
    """Whole fault segments withheld (1-based connected component IDs of the label raster)."""
    try:
        from scipy import ndimage

        labeled, n = ndimage.label(gt > 0)
    except Exception:
        return (gt > 0).copy()
    if n == 0:
        return np.zeros_like(gt, dtype=bool)
    sizes = np.bincount(labeled.ravel())[1:]
    order = np.argsort(-sizes) + 1  # 1-based component labels (0 is background!)
    n_hold = max(1, int(round(n * holdout_frac)))
    hold_labels = set(int(v) for v in order[:n_hold])
    return np.isin(labeled, list(hold_labels))


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
    from gems80.holdout import score_holdout, withhold_segments
    from gems80.uniqueness import check_registry, load_compact_registry

    with rasterio.open(prediction) as pred, rasterio.open(labels_path) as ls, rasterio.open(
        template
    ) as tmpl:
        values = pred.read(1).astype(np.float64)
        gt = ls.read(1).astype(np.float64)
        foot = tmpl.read_masks(1) > 0
    if gt.shape != values.shape:
        report["holdout"]["note"] = "labels grid does not match submission grid"
        return

    cat_bool = (gt > 0) & foot
    if cat_bool.any():
        split = withhold_segments(cat_bool, foot)
        if values.max() > 0.5 and np.isin(np.unique(values[np.isfinite(values)]), [0.0, 1.0]).all():
            h_res = score_holdout((values > 0.5) & foot, split, foot, bootstrap=True)
            report["holdout"] = {
                "available": True,
                "label": h_res["label"],
                "evaluator_version": h_res["evaluator_version"],
                "withheld_positive_count": h_res["withheld_positive_count"],
                "score": h_res["score"],
                "ci95": h_res["ci95"],
                "TP_w": h_res["TP_w"],
                "FP_w": h_res["FP_w"],
                "FN_w": h_res["FN_w"],
                "spatial_quadrant_dti": h_res["spatial_quadrant_dti"],
                "status": "COMPUTED",
                "note": (
                    "Direct evaluation of this GeoTIFF against labels.tif. Note: live submission "
                    "rasters explicitly mask existing_faults.tif (= labels.tif) plus a 200 m collar "
                    "to avoid wasting false positives on known faults, so only 200-300 m fringe "
                    "kernels touch labels.tif; see research/h80_experiments_summary.json and "
                    "data/cache_h80v2/holdout_exp3_choquet.tif for blind hide-and-recover HOLDOUT-DTI."
                ),
            }
        else:
            res = dti_weighted(np.where(foot & ~split["visible"], values, 0.0), split["hidden"] & foot)
            report["holdout"] = {
                "available": True,
                "label": "HOLDOUT-DTI",
                "evaluator_version": "GEMS distance-weighted Tversky alpha=0.2 beta=0.8 R=300m (continuous)",
                "withheld_positive_count": int((split["hidden"] & foot).sum()),
                "score": res["DTI"],
                "ci95": [res["DTI"], res["DTI"]],
                "TP_w": res["TP_w"],
                "FP_w": res["FP_w"],
                "FN_w": res["FN_w"],
                "status": "COMPUTED",
            }
        eval_domain = foot & ~split["visible"]
        pos_mask = split["hidden"] & eval_domain
        neg_mask = (~split["hidden"]) & eval_domain
    else:
        pos_mask = (gt > 0) & foot
        neg_mask = (gt <= 0) & foot

    # Leakage canary: each feature band evaluated on held-out positives vs off-visible negatives
    nodata_floor = -1e30
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    with rasterio.open(features_path) as fs:
        aucs = []
        for b in range(1, fs.count + 1):
            band = fs.read(b).astype(np.float64)
            if band.shape != gt.shape:
                continue
            ok = np.isfinite(band) & (band > nodata_floor)
            pos_idx = np.flatnonzero((pos_mask & ok).ravel())
            neg_pool = np.flatnonzero((neg_mask & ok).ravel())
            if pos_idx.size == 0 or neg_pool.size == 0:
                continue
            if neg_pool.size > 200_000:
                neg_idx = rng.choice(neg_pool, size=200_000, replace=False)
            else:
                neg_idx = neg_pool
            sample_idx = np.concatenate([pos_idx, neg_idx])
            y_bin = np.concatenate([np.ones(pos_idx.size, dtype=int), np.zeros(neg_idx.size, dtype=int)])
            auc_val = float(auc_from_ranks(band.ravel()[sample_idx], y_bin))
            desc = fs.tags(b).get("description") or (fs.descriptions[b - 1] if fs.descriptions and fs.descriptions[b - 1] else f"band {b}")
            aucs.append({"band": b, "description": desc, "auc": round(auc_val, 6)})
    max_auc = max((max(a["auc"], 1.0 - a["auc"]) for a in aucs if np.isfinite(a["auc"])), default=0.0)
    report["leakage_canary"] = {
        "available": True,
        "per_band_auc": aucs,
        "max_auc": round(float(max_auc), 6),
        "leakage_flag": bool(max_auc > 0.90),
        "status": "COMPUTED",
        "note": "AUC > 0.90 on a single feature means leakage until proven otherwise",
    }

    # Lane uniqueness check if registry rasters are available AND grid matches competition dimensions
    reg_dir = REPO / "data" / "registry"
    if reg_dir.is_dir() and values.shape == (3730, 3292):
        reg_paths = [p for p in sorted(reg_dir.glob("*.tif")) if p.resolve() != prediction.resolve()]
        if reg_paths:
            reg = load_compact_registry(reg_paths, foot)
            u = check_registry((values > 0.5) & foot, values.astype(np.float32), foot, reg)
            report["lane_check"] = {
                "available": True,
                "registry_count": len(reg_paths),
                "duplicate_lane": u["duplicate_lane"],
                "max_abs_spearman_surface": u["max_abs_spearman_surface"],
                "max_abs_spearman_dots": u["max_abs_spearman_dots"],
                "max_dot_fraction_within_3px": u["max_dot_fraction_within_3px"],
                "thresholds": u["thresholds"],
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
