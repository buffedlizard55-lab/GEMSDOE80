#!/usr/bin/env python3
"""Rebuild the GEMSDOE80 H80-C unique submission GeoTIFF.

Selected lane (experiment 3 of 3): catalogue-residual along-fabric
directional derivative of dilatation, corroborated by magnetic tilt/curvature,
4 px catalogue collar, 3 px NMS.

Requires data/sample_submission.tif, data/labels.tif, data/training_features.tif.
Writes zeros + nan twins under artifacts/candidates/ (git-ignored); it does not publish or promote them.
"""
from __future__ import annotations

import argparse
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.ndimage import binary_dilation, distance_transform_edt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from gems80.emission import (  # noqa: E402
    along_fabric_derivative,
    distance_residual,
    nms_dots,
    percentile_norm,
)
from gems80.io_grid import load_band, load_labels, load_template  # noqa: E402
from submission_writer import write_pair  # noqa: E402


def build_h80c(dil, rtp, tc, valid, foot, labels, buffer_px=4, sep=3, max_k=34033):
    dist = distance_transform_edt(~labels)
    near = binary_dilation(labels, iterations=buffer_px)
    emit = foot & valid & ~near
    af = along_fabric_derivative(np.where(valid, dil, 0.0), np.where(valid, rtp, 0.0))
    resid = distance_residual(af, dist, valid)
    score = percentile_norm(resid, emit) * (0.50 + 0.50 * percentile_norm(tc, valid))
    score = np.where(dil > 0, score, score * np.float32(0.20))
    score = np.where(emit, score, 0.0).astype(np.float32)
    dots = nms_dots(score, emit, sep=sep, max_k=max_k)
    return dots, score


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--template", default=str(ROOT / "data" / "sample_submission.tif"))
    ap.add_argument("--features", default=str(ROOT / "data" / "training_features.tif"))
    ap.add_argument("--labels", default=str(ROOT / "data" / "labels.tif"))
    ap.add_argument("--out-dir", default=str(ROOT / "artifacts" / "candidates"))
    ap.add_argument("--out", default=None, help="Ignored compatibility flag (H80-C writes a zeros/nan pair).")
    args = ap.parse_args()

    template, features, labels_p = Path(args.template), Path(args.features), Path(args.labels)
    if not template.is_file():
        print("[STATUS] BLOCKED: NOT SUBMISSION-READY: missing sample_submission.tif")
        print("[ACTION] Place the official template in data/ and re-run.")
        return 2
    if not features.is_file() or not labels_p.is_file():
        print("[STATUS] BLOCKED: NOT SUBMISSION-READY: missing training_features.tif or labels.tif")
        return 2

    tpl = load_template(template)
    foot = tpl["foot"]
    labels = load_labels(labels_p, foot)
    dil, v8, *_ = load_band(features, 8)
    rtp, v2, *_ = load_band(features, 2)
    tc, v6, *_ = load_band(features, 6)
    valid = foot & v8 & v2 & v6
    dots, _ = build_h80c(dil, rtp, tc, valid, foot, labels)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    stem = f"gemsdoe80-h80c-strain-accom-resid-{stamp}"
    bundle = write_pair(dots, template, labels_p, Path(args.out_dir), stem)
    print("[OK] zeros", bundle["zeros"]["path"], "sha256", bundle["zeros"]["sha256"])
    print("[OK] format_valid", bundle["zeros"]["format_valid"], "positive_px", bundle["zeros"]["positive_px"])
    return 0 if bundle["zeros"]["format_valid"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
