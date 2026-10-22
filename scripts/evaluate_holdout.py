#!/usr/bin/env python3
"""Hide-and-rebuild spatial holdout for the current H80-C detector.

Unlike scoring a prebuilt submission raster, this script recomputes every
catalogue-dependent term from visible faults in each holdout split, masks those
visible pixels, and scores only withheld whole fault segments. It writes a JSON
report, not a submission TIFF.

Usage:
  python scripts/evaluate_holdout.py \
      --template data/sample_submission.tif \
      --labels data/labels.tif \
      --features data/training_features.tif [--report reports/holdout.json]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import rasterio

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from gems80.emission import gradient_magnitude  # noqa: E402
from gems80.holdout import score_holdout, withhold_segments  # noqa: E402
from gems80.io_grid import NODATA_F32, load_template  # noqa: E402
from gems80.metric import EVALUATOR_VERSION  # noqa: E402
from generate_submission import build_h80c  # noqa: E402
from run_h80_experiments import leakage_canary  # noqa: E402

# Band ids verified in the existing H80 run card; descriptions are retained in the report.
BAND_IDS = {"dilatation": 8, "strain_2ndinv": 4, "tmi_hg": 3, "tilt_tc": 6,
            "conductivity": 17, "rtp": 2, "shear": 7, "earthquake_density": 16}


def assert_same_grid(dataset, reference, name: str) -> None:
    mismatches = []
    if dataset.height != reference.height or dataset.width != reference.width:
        mismatches.append("shape")
    if dataset.crs != reference.crs:
        mismatches.append("CRS")
    if dataset.transform != reference.transform:
        mismatches.append("transform")
    if mismatches:
        raise ValueError(f"{name} grid does not match template ({', '.join(mismatches)})")


def evaluate(template_path: Path, labels_path: Path, features_path: Path) -> dict:
    tpl = load_template(template_path)
    foot = tpl["foot"]

    with rasterio.open(template_path) as template, rasterio.open(labels_path) as label_src, rasterio.open(features_path) as fs:
        assert_same_grid(label_src, template, "labels")
        assert_same_grid(fs, template, "features")
        if fs.count < max(BAND_IDS.values()):
            raise ValueError(f"features raster has {fs.count} bands; H80-C canary requires at least {max(BAND_IDS.values())}")
        labels_raw = label_src.read(1)
        labels = (labels_raw >= 1) & foot & np.isfinite(labels_raw)

        bands = {}
        descriptions = {}
        arrays = {}
        for name, band_id in BAND_IDS.items():
            arr = fs.read(band_id).astype(np.float32)
            desc = fs.descriptions[band_id - 1] or f"band {band_id}"
            valid = np.isfinite(arr) & (arr > NODATA_F32 * np.float32(0.5))
            arrays[name] = (arr, valid)
            bands[name] = (arr, valid)
            descriptions[name] = {"band": band_id, "description": desc}

        dil, v8 = arrays["dilatation"]
        inv, v4 = arrays["strain_2ndinv"]
        hg, v3 = arrays["tmi_hg"]
        tc, v6 = arrays["tilt_tc"]
        cond, v17 = arrays["conductivity"]
        rtp, v2 = arrays["rtp"]
        valid_model = foot & v8 & v2 & v6
        if int(valid_model.sum()) == 0:
            raise ValueError("no valid H80-C feature pixels inside the template footprint")

        split = withhold_segments(labels, foot)
        # This is the essential hide-and-rebuild step: only split['visible'] enters
        # the catalogue-derived feature (distance residual and exclusion collar).
        dots, _ = build_h80c(dil, rtp, tc, valid_model, foot, split["visible"])
        holdout = score_holdout(dots, split, foot, bootstrap=True)

        bands["dilatation_grad"] = (gradient_magnitude(np.where(v8, dil, 0.0)), v8)
        bands["inv_grad"] = (gradient_magnitude(np.where(v4, inv, 0.0)), v4)
        canary = leakage_canary(bands, split["hidden"], foot, split["visible"])

    return {
        "label": "HOLDOUT-DTI",
        "evaluator_version": EVALUATOR_VERSION,
        "method": "hide whole connected catalogue segments with 3 px buffer; derive catalogue-distance residual and exclusion collar from visible faults only; mask visible fault pixels exactly; pooled DTI on withheld positives",
        "withheld_positive_count": int(split["n_hidden_px"]),
        "held_components": int(split["held_ids"].size),
        "total_components": int(split["n_components"]),
        "holdout_dti": holdout,
        "leakage_canary": canary,
        "feature_bands": descriptions,
        "submission_tif_written": False,
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--template", type=Path, default=ROOT / "data" / "sample_submission.tif")
    parser.add_argument("--labels", type=Path, default=ROOT / "data" / "labels.tif")
    parser.add_argument("--features", type=Path, default=ROOT / "data" / "training_features.tif")
    parser.add_argument("--report", type=Path, default=None)
    args = parser.parse_args(argv)
    missing = [str(p) for p in (args.template, args.labels, args.features) if not p.is_file()]
    if missing:
        print(json.dumps({"status": "BLOCKED", "missing_inputs": missing,
                          "message": "No holdout score computed; required rasters are unavailable."}, indent=2))
        return 2
    try:
        result = evaluate(args.template, args.labels, args.features)
    except Exception as exc:
        print(json.dumps({"status": "ERROR", "error": str(exc), "score_computed": False}, indent=2))
        return 1
    payload = json.dumps(result, indent=2)
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(payload + "\n")
    print(payload)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
