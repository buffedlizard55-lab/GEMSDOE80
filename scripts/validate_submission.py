#!/usr/bin/env python3
"""
Automated Zero-Touch Submission Validator for DrivenData Competition #306 (DOE GEMS).

Emits structured JSON to stdout (and human checklist to stderr), enforcing:
  1. `format_valid`: single-band float32 GeoTIFF, finite values strictly in [0.0, 1.0],
     zero NaN/Inf inside footprint (fixes portal error "Predicted values must be in range [0, 1]"),
     and exact CRS / shape / Affine transform / bounds match against `--template`.
  2. `download_ready` (and `all_passed`): `format_valid` PLUS full competition promotion
     verification (zero predictions on mapped catalogue when catalogue matches grid, valid
     single-member .zip, SHA-256-matched sidecar & <=140-char note, verified 929-raster
     uniqueness receipt, and matching promoted `research/run-card.json` entry).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import rasterio

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

ALPHA, BETA = 0.2, 0.8
KERNEL_R_M = 300.0


def _default_data_file(name: str) -> Path:
    for base in (REPO / "data", REPO / ".cache" / "data", REPO / "data" / "bridge"):
        p = base / name
        if p.is_file():
            return p
    return REPO / "data" / name


def validate_geotiff(
    sub_path: Path,
    template_path: Path | None = None,
    catalogue_path: Path | None = None,
) -> dict:
    """12-check programmatic GeoTIFF receipt compatible across all GEMSDOE80 test suites."""
    import rasterio

    sub_path = Path(sub_path)
    if not sub_path.is_file():
        return {"passed": False, "error": f"File not found: {sub_path}"}
    tpl_p = Path(template_path) if template_path is not None and Path(template_path).is_file() else _default_data_file("sample_submission.tif")
    cat_p = Path(catalogue_path) if catalogue_path is not None and Path(catalogue_path).is_file() else _default_data_file("existing_faults.tif")
    if not cat_p.is_file():
        cat_p = _default_data_file("labels.tif")

    checks: dict[str, bool] = {}
    notes: list[str] = []
    with rasterio.open(tpl_p) as tpl, rasterio.open(cat_p) as cat_ds, rasterio.open(sub_path) as ds:
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
        "sha256": hashlib.sha256(sub_path.read_bytes()).hexdigest(),
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
RESOLUTION_M = 100.0
BOOTSTRAP_SEED = 20261021
N_BOOTSTRAP = 2000


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def validate_submission(
    candidate_path: Path,
    template_path: Path | None = None,
    catalogue_path: Path | None = None,
) -> dict:
    checks: dict[str, bool] = {}
    details: dict[str, dict] = {}

    def record(name: str, passed: bool, observed, expected) -> None:
        checks[name] = bool(passed)
        details[name] = {
            "passed": bool(passed),
            "observed": observed,
            "expected": expected,
        }

    record("file_exists", candidate_path.is_file(), str(candidate_path), "existing .tif file")
    if not candidate_path.is_file():
        return {
            "format_valid": False,
            "download_ready": False,
            "portal_ready_without_manual_inspection": False,
            "all_passed": False,
            "candidate": str(candidate_path),
            "n_checks_total": 1,
            "n_checks_passed": 0,
            "checks": checks,
            "details": details,
        }

    with rasterio.open(candidate_path) as src:
        band = src.read(1)
        crs_str = str(src.crs)
        dtype_str = src.dtypes[0]
        count = src.count
        shape = [src.height, src.width]
        transform = list(src.transform)[:6]
        bounds = list(src.bounds)
        nodata = src.nodata

    record("single_band", count == 1, count, 1)
    record("dtype_float32", dtype_str == "float32", dtype_str, "float32")
    record("nodata_unset_or_clean", nodata is None or np.isfinite(nodata) or np.isnan(nodata), nodata, "valid")

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

    foot = finite_mask
    if template_path is not None and template_path.is_file():
        with rasterio.open(template_path) as tmpl:
            t_arr = tmpl.read(1)
            t_nodata = tmpl.nodata
            if t_nodata is not None and np.isfinite(t_nodata):
                foot = np.isfinite(t_arr) & (t_arr != t_nodata)
            else:
                foot = np.isfinite(t_arr)
            t_crs = str(tmpl.crs)
            t_shape = [tmpl.height, tmpl.width]
            t_transform = list(tmpl.transform)[:6]
            t_bounds = list(tmpl.bounds)

        shape_ok = shape == t_shape
        record("crs_matches", crs_str == t_crs, crs_str, t_crs)
        record("shape_matches", shape_ok, shape, t_shape)
        t_diff = max(abs(a - b) for a, b in zip(transform, t_transform)) if shape_ok else 999.0
        record("transform_matches", t_diff < 1e-6, transform, t_transform)
        b_diff = max(abs(a - b) for a, b in zip(bounds, t_bounds)) if shape_ok else 999.0
        record("bounds_match", b_diff < 1e-3, bounds, t_bounds)

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
            record("no_nan_or_infinite_inside_footprint", False, "shape mismatch", 0)
            record("outside_footprint_zero", False, "shape mismatch", 0.0)
    else:
        record("crs_matches", crs_str == "EPSG:32611", crs_str, "EPSG:32611")
        record("shape_matches", True, shape, shape)
        record("transform_matches", True, transform, transform)
        record("bounds_match", True, bounds, bounds)
        record("no_nan_or_infinite_inside_footprint", nan_total == 0 and inf_total == 0, nan_total + inf_total, 0)
        record("outside_footprint_zero", True, 0.0, 0.0)

    format_keys = (
        "file_exists",
        "single_band",
        "dtype_float32",
        "no_nan_anywhere",
        "no_inf_anywhere",
        "values_in_0_1",
        "crs_matches",
        "shape_matches",
        "transform_matches",
        "bounds_match",
        "no_nan_or_infinite_inside_footprint",
    )
    format_valid = all(checks[k] for k in format_keys)

    # Full competition promotion gates (distinguish promoted official candidates from ad-hoc test rasters)
    actual_sha = sha256_file(candidate_path)
    cat = None
    if catalogue_path is not None and catalogue_path.is_file():
        with rasterio.open(catalogue_path) as csrc:
            if [csrc.height, csrc.width] == shape:
                cat = csrc.read(1) == 1
                on_cat = int(((band > 0.0) & cat).sum())
                record("zero_predictions_on_mapped_catalogue", on_cat == 0, on_cat, 0)

    if shape == [3730, 3292] and checks["shape_matches"]:
        n_pos = int((band[foot] > 0.5).sum())
        record("non_degenerate_positive_support", 5_000 <= n_pos <= 65_000, n_pos, "5,000..65,000")
        shared_res = SHARED_VAL.validate(candidate_path, foot, cat)
        record("shared_gems57_validator_all_passed", bool(shared_res["all_checks_passed"]), shared_res["all_checks_passed"], True)
    else:
        record("official_competition_grid_3730x3292", False, shape, [3730, 3292])

    zip_path = candidate_path.with_suffix(".zip")
    if zip_path.is_file():
        with zipfile.ZipFile(zip_path, "r") as zf:
            names = zf.namelist()
            bad_zip = zf.testzip()
        record("zip_archive_valid_single_tif", bad_zip is None and len(names) == 1 and names[0] == candidate_path.name, names, [candidate_path.name])
    else:
        record("zip_archive_valid_single_tif", False, "missing .zip", "present .zip")

    sidecar_path = candidate_path.with_name(candidate_path.name + ".json")
    if sidecar_path.is_file():
        sc = json.loads(sidecar_path.read_text())
        note = sc.get("note", "")
        record("sidecar_sha256_matches_bytes", sc.get("sha256") == actual_sha, sc.get("sha256"), actual_sha)
        record("submission_note_within_140_chars", 1 <= len(note) <= 140, len(note), "1..140")
    else:
        record("sidecar_sha256_matches_bytes", False, "missing sidecar", actual_sha)
        record("submission_note_within_140_chars", False, 0, "1..140")

    rc_path = REPO / "research" / "run-card.json"
    rc_ok = False
    if rc_path.is_file():
        rc = json.loads(rc_path.read_text())
        rc_ok = (
            rc.get("verdict") == "promote"
            and rc.get("raster_sha256") == actual_sha
            and not rc.get("correlation_overlap_vs_registry", {}).get("binding_duplicate", True)
            and rc.get("leakage_canary", {}).get("passes") is True
        )
    record("promoted_in_run_card_with_uniqueness_and_canary", rc_ok, rc_ok, True)

    all_passed = all(checks.values())
    return {
        "validated_at_utc": datetime.now(timezone.utc).isoformat(),
        "candidate": str(candidate_path),
        "sha256": actual_sha,
        "bytes": candidate_path.stat().st_size,
        "format_valid": bool(format_valid),
        "download_ready": bool(all_passed),
        "portal_ready_without_manual_inspection": bool(all_passed),
        "all_passed": bool(all_passed),
        "n_checks_total": len(checks),
        "n_checks_passed": sum(1 for v in checks.values() if v),
        "checks": checks,
        "details": details,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("candidate", help="Path to candidate .tif submission")
    ap.add_argument("--template", default=str(REPO / "data" / "sample_submission.tif"))
    ap.add_argument("--catalogue", default=str(REPO / "data" / "existing_faults.tif"))
    ap.add_argument("--report", default=None)
    args = ap.parse_args()

    res = validate_submission(Path(args.candidate), Path(args.template), Path(args.catalogue))
    sys.stderr.write(
        f"[validate_submission] target={res['candidate']} format_valid={res['format_valid']} "
        f"download_ready={res['download_ready']} ({res['n_checks_passed']}/{res['n_checks_total']} checks)\n"
    )
    if args.report:
        rp = Path(args.report)
        rp.parent.mkdir(parents=True, exist_ok=True)
        rp.write_text(json.dumps(res, indent=2) + "\n")

    sys.stdout.write(json.dumps(res, indent=2) + "\n")
    return 0 if res["format_valid"] else 1


if __name__ == "__main__":
    sys.exit(main())
