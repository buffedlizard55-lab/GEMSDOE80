#!/usr/bin/env python3
"""Hide-and-recover holdout DTI for a candidate GeoTIFF.

Usage:
  python scripts/evaluate_holdout.py PREDICTION.tif \\
      --template data/sample_submission.tif --labels data/labels.tif

Prints a JSON object labelled HOLDOUT-DTI (never a projection / organizer score).
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

from gems80.holdout import score_holdout, withhold_segments  # noqa: E402
from gems80.io_grid import load_labels, load_template  # noqa: E402


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("prediction", type=Path)
    ap.add_argument("--template", type=Path, default=ROOT / "data" / "sample_submission.tif")
    ap.add_argument("--labels", type=Path, default=ROOT / "data" / "labels.tif")
    args = ap.parse_args(argv)

    tpl = load_template(args.template)
    foot = tpl["foot"]
    labels = load_labels(args.labels, foot)
    with rasterio.open(args.prediction) as src:
        pred = src.read(1)
    pred_bool = np.isfinite(pred) & (pred > 0.5)
    split = withhold_segments(labels, foot)
    report = score_holdout(pred_bool, split, foot)
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
