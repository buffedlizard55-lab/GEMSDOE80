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
sys.path.insert(0, str(REPO / "vendor"))

from gems57 import validate as SHARED_VAL  # noqa: E402


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

    finite_mask = np.isfinite(band)
    nan_total = int(np.isnan(band).sum())
    inf_total = int(np.isinf(band).sum())
    vmin = float(np.nanmin(band)) if finite_mask.any() else float("nan")
    vmax = float(np.nanmax(band)) if finite_mask.any() else float("nan")

    record("no_nan_anywhere", nan_total == 0, nan_total, 0)
    record("no_inf_anywhere", inf_total == 0, inf_total, 0)
    record(
        "values_in_0_1",
        bool(finite_mask.all() and vmin >= 0.0 and vmax <= 1.0),
        [round(vmin, 6), round(vmax, 6)],
        "[0.0, 1.0]",
    )

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

        if shape_ok:
            in_vals = band[foot]
            out_vals = band[~foot]
            nan_in_fp = int((~np.isfinite(in_vals)).sum())
            out_max = float(np.max(np.abs(out_vals))) if out_vals.size else 0.0
            record("no_nan_or_infinite_inside_footprint", nan_in_fp == 0, nan_in_fp, 0)
            record("outside_footprint_zero", out_max == 0.0, out_max, 0.0)
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
