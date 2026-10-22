#!/usr/bin/env python3
"""Verify official GeoDAWN inputs and materialize the 34-channel static memmap cache."""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gems80.features_multiphysics import build_or_load_static_stack
from gems80.io_grid import feature_descriptions, load_labels, load_template

REQUIRED = [
    "training_features.tif",
    "labels.tif",
    "existing_faults.tif",
    "sample_submission.tif",
]


def main() -> int:
    data = ROOT / "data"
    missing = [f for f in REQUIRED if not (data / f).is_file()]
    if missing:
        print("MISSING inputs in data/:", ", ".join(missing), file=sys.stderr)
        print("Run `python scripts/restore_mirrors.py` first.", file=sys.stderr)
        return 2

    tmpl = load_template(data / "sample_submission.tif")
    foot = tmpl["foot"]
    lab = load_labels(data / "labels.tif", foot)
    cat = load_labels(data / "existing_faults.tif", foot)
    descs = feature_descriptions(data / "training_features.tif")

    mm, aux, static_names = build_or_load_static_stack(data, foot, tmpl["meta"]["transform"])

    summary = {
        "template": tmpl["report"],
        "labels_positive_px": int(lab.sum()),
        "existing_faults_positive_px": int(cat.sum()),
        "labels_equal_existing_faults": bool((lab == cat).all()),
        "feature_bands": len(descs),
        "static_multiphysics_channels": len(static_names),
        "static_cache_path": "data/cache_h80v2/static_f16.mmap",
        "static_cache_shape": list(mm.shape),
        "aux_fields": sorted(list(aux.keys())),
    }
    out = ROOT / "research" / "data-summary.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
