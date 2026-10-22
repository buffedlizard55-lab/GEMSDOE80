#!/usr/bin/env python3
"""Prepare and cache the 50-channel static multi-physics + updip-advected + geothermal
+ deconfounded scarp feature stack in `data/processed/`, and run the mandatory
59-feature single-feature leakage canary audit (`AUC <= 0.90`, Protocol Rule 4)
on the spatially-blocked buffered holdout (`gems57-pooled-hide-v2`)."""
from __future__ import annotations

import json
import resource
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "vendor"))

from gems57 import grid as GRID, holdout as HO  # noqa: E402
from gems80.features import build_or_load_static_stack, run_leakage_canary_audit  # noqa: E402


def rss_mb() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


def main() -> None:
    t0 = time.time()
    raw_dir = ROOT / "data" / "raw"
    ext_dir = raw_dir / "external"
    cache_dir = ROOT / "data" / "processed"
    evid_dir = ROOT / "evidence"
    evid_dir.mkdir(parents=True, exist_ok=True)

    g = GRID.load_grid(ROOT / "data" / "bridge")
    print(f"[prepare_data] Loaded GEMS grid {g.shape}, footprint={int(g.footprint.sum()):,}, catalogue={int(g.catalogue.sum()):,}")

    stack = build_or_load_static_stack(raw_dir, ext_dir, cache_dir, g.footprint, force_rebuild=False)
    print(f"[prepare_data] Static stack ready: shape={stack.X_mmap.shape}, RSS={rss_mb():.1f} MB ({time.time() - t0:.1f}s)")

    folds20, _ = HO.buffered_component_draw(g, seed=20, hide_frac=0.20, buffer_px=3)
    canary = run_leakage_canary_audit(stack, folds20, g.footprint, max_allowed_auc=0.90, seed=80)
    canary_path = evid_dir / "h80_leakage_canary.json"
    canary_path.write_text(json.dumps(canary, indent=2) + "\n")
    print(
        f"[prepare_data] Leakage canary PASS: {canary['n_features_tested']} features tested, "
        f"max_symmetric_auc={canary['max_symmetric_auc']:.5f} <= 0.90 "
        f"(top={canary['top_feature']}), RSS={rss_mb():.1f} MB, total={time.time() - t0:.1f}s"
    )


if __name__ == "__main__":
    main()
