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
# 8-connectivity keeps diagonal fault rasters in one segment, so pieces of one fault
# cannot be split between the visible and withheld sets (which would leak).
EIGHT = np.ones((3, 3), dtype=bool)


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
    """k(d) = max(1 - d/R, 0), as in the official rules (problem description, metric section)."""
    return np.clip(1.0 - np.asarray(d, dtype=np.float64) / r, 0.0, None)


def _window_offsets(radius_px: float):
    """All integer offsets (dy, dx) with Euclidean distance <= radius, including the centre."""
    rc = int(math.ceil(radius_px))
    out = []
    for dy in range(-rc, rc + 1):
        for dx in range(-rc, rc + 1):
            d = math.hypot(dy, dx)
            if d <= radius_px + 1e-12:
                out.append((dy, dx, d))
    return rc, out


def dti_components(pred, gt, domain=None, resolution_m=RESOLUTION_M, alpha=ALPHA, beta=BETA):
    """Distance-weighted Tversky components, following the competition definition.

    TP_w = sum_{g in G} max_{x: d(x,g) <= R} p(x) k(d(x,g))            (centre x = g included)
    FN_w = sum_{g in G} [1 - max_{x: d(x,g) <= R} p(x) k(d(x,g))]
    FP_w = sum_{x: p(x) > 0} p(x) [1 - max_{g in G} k(d(x,g))]
    DTI  = TP_w / (TP_w + alpha FP_w + beta FN_w + eps)

    `domain` (bool array) restricts every term to the scored pixels; pixels outside it
    contribute nothing. Also returns per-fault-segment sums for a segment-level bootstrap.
    """
    from scipy import ndimage

    R = KERNEL_R_M / resolution_m  # 3.0 px at 100 m
    H, W = pred.shape
    dom = np.ones((H, W), dtype=bool) if domain is None else np.asarray(domain, dtype=bool)
    p = np.where(np.isfinite(pred) & dom, np.clip(np.nan_to_num(pred, nan=0.0), 0.0, 1.0), 0.0).astype(np.float64)
    g = (np.asarray(gt) > 0) & dom

    rc, offs = _window_offsets(R)
    padded = np.pad(p, rc, mode="constant", constant_values=0.0)
    local_max = np.zeros((H, W), dtype=np.float64)
    for dy, dx, d in offs:
        shifted = padded[rc + dy: rc + dy + H, rc + dx: rc + dx + W]
        local_max = np.maximum(local_max, shifted * float(triangular_kernel(d, R)))

    n_seg = 0
    if g.any():
        seg_labels, n_seg = ndimage.label(g, structure=EIGHT)  # 8-connected fault segments
        dist, (iy, ix) = ndimage.distance_transform_edt(~g, return_indices=True)
        near_seg = seg_labels[iy, ix]          # segment of the nearest GT pixel
        k_near = triangular_kernel(dist, R)
    else:
        near_seg = np.zeros((H, W), dtype=np.int64)
        k_near = np.zeros((H, W), dtype=np.float64)

    tp_pix = np.where(g, local_max, 0.0)
    fn_pix = np.where(g, np.clip(1.0 - local_max, 0.0, None), 0.0)
    fp_pix = np.where(p > 0, p * (1.0 - k_near), 0.0)

    tp_w, fn_w, fp_w = float(tp_pix.sum()), float(fn_pix.sum()), float(fp_pix.sum())
    dti = tp_w / (tp_w + alpha * fp_w + beta * fn_w + 1e-12)

    if n_seg:
        seg_idx = seg_labels[g]
        tp_seg = np.bincount(seg_idx, weights=tp_pix[g], minlength=n_seg + 1)[1:]
        fn_seg = np.bincount(seg_idx, weights=fn_pix[g], minlength=n_seg + 1)[1:]
        fp_seg = np.bincount(near_seg[fp_pix > 0].ravel(), weights=fp_pix[fp_pix > 0], minlength=n_seg + 1)[1:]
        fp_bg = float(np.bincount(near_seg[fp_pix > 0].ravel(), weights=fp_pix[fp_pix > 0], minlength=n_seg + 1)[0])
    else:
        tp_seg = fn_seg = fp_seg = np.zeros(0)
        fp_bg = fp_w
    return {
        "TP_w": tp_w, "FP_w": fp_w, "FN_w": fn_w, "DTI": float(dti),
        "n_gt_pixels": int(g.sum()), "n_segments": int(n_seg),
        "segments": {"tp": tp_seg, "fn": fn_seg, "fp": fp_seg, "fp_background": fp_bg},
    }


def dti_weighted(pred, gt, resolution_m=RESOLUTION_M, alpha=ALPHA, beta=BETA, domain=None):
    """Distance-weighted Tversky index (scalar), see dti_components."""
    res = dti_components(pred, gt, domain=domain, resolution_m=resolution_m, alpha=alpha, beta=beta)
    return {k: res[k] for k in ("TP_w", "FP_w", "FN_w", "DTI")}


def segment_bootstrap_ci(seg, alpha=ALPHA, beta=BETA, n_bootstrap=N_BOOTSTRAP, seed=BOOTSTRAP_SEED):
    """95% CI by resampling whole fault segments with replacement (cluster bootstrap).

    Each replicate draws segments with replacement and recomputes DTI from the summed
    per-segment TP/FN/FP. FP mass not attributable to any segment is held fixed.
    """
    tp, fn, fp = seg["tp"], seg["fn"], seg["fp"]
    n = tp.size
    if n == 0:
        return (float("nan"), float("nan"))
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, n, size=(n_bootstrap, n))
    tp_b = tp[idx].sum(axis=1)
    fn_b = fn[idx].sum(axis=1)
    fp_b = fp[idx].sum(axis=1) + seg["fp_background"]
    dti_b = tp_b / (tp_b + alpha * fp_b + beta * fn_b + 1e-12)
    return float(np.percentile(dti_b, 2.5)), float(np.percentile(dti_b, 97.5))


def holdout_segment_mask(gt, holdout_frac=1.0 / 3.0, seed=BOOTSTRAP_SEED):
    """Withhold a seeded random 1/3 of whole fault segments (8-connected components of GT).

    Random (not largest-first) so that the held-out set is not biased toward the longest faults.
    Spatially-blocked holdout and buffered visible masking are TODO (see README limitations).
    """
    from scipy import ndimage

    labeled, n = ndimage.label(gt > 0, structure=EIGHT)
    if n == 0:
        return np.zeros(gt.shape, dtype=bool)
    rng = np.random.default_rng(seed)
    n_hold = max(1, int(round(n * holdout_frac)))
    hold_labels = rng.choice(np.arange(1, n + 1), size=n_hold, replace=False)
    return np.isin(labeled, hold_labels)


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
            if nan_inside:
                errors.append(f"{nan_inside} NaN pixels inside the prediction footprint")

            footprint_overlap = bool((pred_mask & template_mask).any()) if pred_mask.any() else False
            checks["footprint_overlaps_template"] = footprint_overlap

            # The TEMPLATE, never the candidate's own nodata mask, defines required coverage.
            # A candidate cannot hide holes by declaring them nodata.
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
                report["missing_inside_footprint"] = None

            # Finite sentinels (e.g. -9999 masked as nodata) must still fail the range rule.
            finite_raw = np.isfinite(values)
            if np.any(finite_raw & ((values < 0) | (values > 1))):
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
    # Promotion is a separate, evidence-gated selector step; format validity never sets this.
    report["download_ready"] = False
    report["promotion_blockers"] = [
        "Shared hide-and-recover HOLDOUT-DTI and CI unavailable unless --labels/--features supplied",
        "Per-feature leakage canary unavailable unless --labels/--features supplied",
        "Registry duplicate check (rank correlation > 0.90 or > 70% dots within 3 px) not run",
        "Current holdout best and weekly cap unverified",
    ]
    return report


def holdout_dti(report, prediction, template, labels_path, features_path):
    """Internal hide-and-recover holdout (NOT the organizer's private test set).

    * Withheld positives = a seeded random 1/3 of whole fault segments (holdout_segment_mask).
    * Visible faults are not scored: pixels within 300 m (R) of a visible fault are removed
      from the scoring domain, so predictions there count neither as TP nor FP. This is an
      internal protocol choice; the official rules do not state a visible-fault mask.
    * Scoring domain also requires the template footprint (valid template pixels).
    * Canary: each feature band's AUC, held-out positives vs in-domain negatives.
    """
    import rasterio
    from scipy import ndimage

    with rasterio.open(prediction) as pred, rasterio.open(labels_path) as ls, rasterio.open(template) as tmpl:
        values = pred.read(1).astype(np.float64)
        gt = ls.read(1).astype(np.float64)
        footprint = tmpl.read_masks(1) > 0
    if gt.shape != values.shape:
        report["holdout"]["note"] = "labels grid does not match submission grid"
        return

    gt_pos = gt > 0
    held = holdout_segment_mask(gt)
    gt_held = gt_pos & held
    visible = gt_pos & ~held
    # visible-fault buffer of R = 3 px (300 m) removed from scoring domain
    yy, xx = np.mgrid[-3:4, -3:4]
    disk = (yy ** 2 + xx ** 2) <= 9  # Euclidean radius 3 px = 300 m
    buf = ndimage.binary_dilation(visible, structure=disk)
    domain = footprint & ~buf

    res = dti_components(values, gt_held.astype(np.float64), domain=domain)
    lo, hi = segment_bootstrap_ci(res["segments"])
    report["holdout"] = {
        "available": True,
        "evaluator_version": "GEMS distance-weighted Tversky alpha=0.2 beta=0.8 R=300m (project re-implementation, see scripts/validate_submission.py)",
        "protocol": "random 1/3 of fault segments withheld; visible-fault 300 m buffer excluded from domain (internal)",
        "withheld_positive_count": int(gt_held.sum()),
        "withheld_segments": int(ndimage.label(gt_held)[1]),
        "score": res["DTI"],
        "ci95": [lo, hi],
        "ci_method": "cluster bootstrap over fault segments, 2000 replicates, seed 20261021",
        "TP_w": res["TP_w"],
        "FP_w": res["FP_w"],
        "FN_w": res["FN_w"],
        "status": "HOLDOUT-DTI (internal evaluator; not organizer-confirmed)",
    }

    with rasterio.open(features_path) as fs:
        aucs = []
        for b in range(1, fs.count + 1):
            band = fs.read(b).astype(np.float64)
            if band.shape != gt.shape:
                continue
            ok = np.isfinite(band) & domain
            pos = ok & gt_held
            neg = ok & ~gt_pos & ~held  # background pixels never used as a withheld positive
            if not pos.any() or not neg.any():
                continue
            rng = np.random.default_rng(BOOTSTRAP_SEED + b)
            neg_idx = np.flatnonzero(neg)
            if neg_idx.size > 200_000:
                neg_idx = rng.choice(neg_idx, size=200_000, replace=False)
            pos_idx = np.flatnonzero(pos)
            if pos_idx.size > 200_000:
                pos_idx = rng.choice(pos_idx, size=200_000, replace=False)
            scores = np.concatenate([band.ravel()[pos_idx], band.ravel()[neg_idx]])
            labs = np.concatenate([np.ones(pos_idx.size, int), np.zeros(neg_idx.size, int)])
            desc = fs.tags(b).get("description", f"band {b}")
            aucs.append({"band": b, "description": desc, "auc": auc_from_ranks(scores, labs)})
    finite = [a["auc"] for a in aucs if np.isfinite(a["auc"])]
    max_auc = max(finite) if finite else float("nan")
    report["leakage_canary"] = {
        "available": True,
        "per_band_auc": aucs,
        "max_auc": float(max_auc),
        "leakage_flag": bool(max_auc > 0.90) if finite else None,
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
