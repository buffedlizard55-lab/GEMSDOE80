#!/usr/bin/env python3
"""Execute the Strictly Budgeted 3-Experiment GEMSDOE80 v2 Ablation Suite.

Peak RAM strictly bounded < 1.1 GiB on 3.8 GiB / 0-swap CPU sandboxes via float16 memmaps,
compact registry pre-indexing, and streaming 47-feature leakage-canary evaluation.
"""

from __future__ import annotations

import gc
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gems80.choquet_submodular import dinkelbach_choquet_allocate
from gems80.emission import auc_mannwhitney, nms_dots
from gems80.features_multiphysics import build_or_load_static_stack
from gems80.holdout import score_holdout, withhold_segments
from gems80.io_grid import load_labels, load_template
from gems80.two_regime_model import fit_predict_phys_oof, fit_predict_struct_oof, fuse_two_regime
from gems80.uniqueness import check_registry, check_surface_preplacement, load_compact_registry


def streaming_leakage_canary(
    mm_static: np.memmap,
    static_names: list[str],
    topo_eval: dict[str, np.ndarray],
    split: dict,
    foot: np.ndarray,
    seed: int = 77,
    max_neg: int = 200_000,
    auc_limit: float = 0.90,
) -> dict:
    """Audit all 47 static + fold-visible topological features in < 5 MB peak RAM."""
    fp_idx = np.flatnonzero(foot.ravel())
    vis_fp = split["visible"].ravel()[fp_idx]
    hid_fp = split["hidden"].ravel()[fp_idx]
    eval_fp = ~vis_fp

    pos_rows = np.flatnonzero(eval_fp & hid_fp)
    neg_pool = np.flatnonzero(eval_fp & (~hid_fp))
    rng = np.random.default_rng(seed)
    neg_rows = rng.choice(neg_pool, size=min(max_neg, int(neg_pool.size)), replace=False)
    sample_fp_rows = np.concatenate([pos_rows, neg_rows])
    sample_flat = fp_idx[sample_fp_rows]
    y = np.concatenate([np.ones(pos_rows.size, dtype=np.int64), np.zeros(neg_rows.size, dtype=np.int64)])

    rows = []
    max_auc = 0.0
    flagged = False
    for j, sname in enumerate(static_names):
        vec = np.nan_to_num(mm_static[sample_fp_rows, j].astype(np.float64), nan=0.0)
        auc = float(auc_mannwhitney(vec, y))
        dev = max(auc, 1.0 - auc)
        max_auc = max(max_auc, dev)
        hit = dev > auc_limit
        flagged = flagged or hit
        rows.append({"feature": sname, "domain": "static_multiphysics", "auc": round(auc, 5), "effective_auc": round(dev, 5), "leakage": bool(hit)})

    for tname, tarr in topo_eval.items():
        vec = np.nan_to_num(tarr.ravel()[sample_flat].astype(np.float64), nan=0.0)
        auc = float(auc_mannwhitney(vec, y))
        dev = max(auc, 1.0 - auc)
        max_auc = max(max_auc, dev)
        hit = dev > auc_limit
        flagged = flagged or hit
        rows.append({"feature": tname, "domain": "fold_visible_topology", "auc": round(auc, 5), "effective_auc": round(dev, 5), "leakage": bool(hit)})

    rows.sort(key=lambda r: -r["effective_auc"])
    return {
        "leakage_flag": bool(flagged),
        "threshold_auc": auc_limit,
        "max_auc": round(max_auc, 5),
        "n_features_audited": len(rows),
        "top_15_features": rows[:15],
        "all_features": rows,
    }


def main() -> int:
    t0 = time.time()
    data = ROOT / "data"
    tmpl = load_template(data / "sample_submission.tif")
    foot = tmpl["foot"]
    cat = load_labels(data / "existing_faults.tif", foot)

    split = withhold_segments(cat, foot)
    visible = split["visible"]
    eval_mask = foot & ~visible
    print(
        f"[setup] Holdout split ready: visible={split['n_visible_px']} px, "
        f"hidden={split['n_hidden_px']} px ({len(split['held_ids'])}/{split['n_components']} comps)",
        flush=True,
    )

    mm_static, aux_fields, static_names = build_or_load_static_stack(data, foot, tmpl["meta"]["transform"])
    reg_paths = sorted((data / "registry").glob("*.tif"))
    registry = load_compact_registry(reg_paths, foot)
    print(
        f"[setup] Loaded {len(static_names)}-ch memmap stack + {len(registry)} compact registry TIFs in {time.time() - t0:.1f}s",
        flush=True,
    )

    # Zero-cost reference baseline: matched-support random dots (37,654 dots off-visible)
    rng = np.random.default_rng(20261022)
    rand_surf = np.zeros(foot.shape, dtype=np.float32)
    rand_surf[eval_mask] = rng.random(int(eval_mask.sum()), dtype=np.float32)
    rand_dots = nms_dots(rand_surf, eval_mask, sep=3, max_k=37_654)
    rand_hold = score_holdout(rand_dots, split, foot, bootstrap=True)
    print(f"[ref] Random control ({int(rand_dots.sum())} NMS dots): HOLDOUT-DTI = {rand_hold['score']:.4f} CI95={rand_hold['ci95']}", flush=True)
    del rand_surf, rand_dots
    gc.collect()

    # =========================================================================
    # EXPERIMENT 1 OF 3: H80-Exp1-MPhys-NMS (Pure Label-Free 34-ch Multi-Physics OOF)
    # =========================================================================
    t_e1 = time.time()
    p_phys, phys_info = fit_predict_phys_oof(mm_static, visible, foot, seed=20261022)
    e1_pre_u = check_surface_preplacement(p_phys, foot, registry)
    assert not e1_pre_u["duplicate_surface"], f"Exp 1 pre-placement duplicate: {e1_pre_u['max_abs_spearman']}"

    dots_e1 = nms_dots(p_phys, eval_mask, sep=3, max_k=37_654)
    hold_e1 = score_holdout(dots_e1, split, foot, bootstrap=True)
    uniq_e1 = check_registry(dots_e1, p_phys, foot, registry)
    dt_e1 = time.time() - t_e1
    print(
        f"[Exp 1/3] H80-Exp1-MPhys-NMS ({int(dots_e1.sum())} dots, {dt_e1:.1f}s): "
        f"HOLDOUT-DTI = {hold_e1['score']:.4f} CI95={hold_e1['ci95']} | "
        f"max|rho|={uniq_e1['max_abs_spearman_surface']:.4f}, max_ov3px={uniq_e1['max_dot_fraction_within_3px']:.4f}, dup={uniq_e1['duplicate_lane']}",
        flush=True,
    )

    # =========================================================================
    # EXPERIMENT 2 OF 3: H80-Exp2-TwoRegime-NMS (+ Nested Self-Supervised Crack-Tip & Bridges)
    # =========================================================================
    t_e2 = time.time()
    p_struct, allowed_eval, topo_eval, struct_info = fit_predict_struct_oof(
        mm_static, aux_fields, visible, foot, nested_seeds=(101, 202, 303), seed=20261022
    )

    canary_report = streaming_leakage_canary(mm_static, static_names, topo_eval, split, foot, seed=77)
    del topo_eval
    gc.collect()
    print(
        f"[canary] Audited all {canary_report['n_features_audited']} features: "
        f"max holdout AUC = {canary_report['max_auc']:.4f} (threshold 0.90, flag={canary_report['leakage_flag']})",
        flush=True,
    )
    assert not canary_report["leakage_flag"], f"LEAKAGE CANARY TRIPPED: {canary_report}"

    p_fused = fuse_two_regime(p_struct, p_phys, allowed_eval, greenfield_weight=0.22)
    e2_pre_u = check_surface_preplacement(p_fused, foot, registry)
    assert not e2_pre_u["duplicate_surface"], f"Exp 2 pre-placement duplicate: {e2_pre_u['max_abs_spearman']}"

    dots_e2 = nms_dots(p_fused, allowed_eval, sep=3, max_k=37_654)
    hold_e2 = score_holdout(dots_e2, split, foot, bootstrap=True)
    uniq_e2 = check_registry(dots_e2, p_fused, foot, registry)
    dt_e2 = time.time() - t_e2
    print(
        f"[Exp 2/3] H80-Exp2-TwoRegime-NMS ({int(dots_e2.sum())} dots, {dt_e2:.1f}s): "
        f"HOLDOUT-DTI = {hold_e2['score']:.4f} CI95={hold_e2['ci95']} | "
        f"max|rho|={uniq_e2['max_abs_spearman_surface']:.4f}, max_ov3px={uniq_e2['max_dot_fraction_within_3px']:.4f}, dup={uniq_e2['duplicate_lane']}",
        flush=True,
    )

    # =========================================================================
    # EXPERIMENT 3 OF 3: H80-Exp3-Choquet-Submodular (Exact 6-Shell Dinkelbach Allocator)
    # =========================================================================
    t_e3 = time.time()
    e3_pre_u = check_surface_preplacement(p_fused, foot, registry)
    assert not e3_pre_u["duplicate_surface"], f"Exp 3 pre-placement duplicate: {e3_pre_u['max_abs_spearman']}"

    dots_e3, alloc_stats = dinkelbach_choquet_allocate(
        p_fused,
        allowed_eval,
        eval_mask=eval_mask,
        min_dots=0,
        max_dots=70_000,
        min_sep_px=2.2,
        n_outer_iters=6,
        max_candidates=420_000,
    )
    hold_e3 = score_holdout(dots_e3, split, foot, bootstrap=True)
    uniq_e3 = check_registry(dots_e3, p_fused, foot, registry)
    assert not uniq_e3["duplicate_lane"], f"Exp 3 post-placement duplicate: {uniq_e3}"
    dt_e3 = time.time() - t_e3
    print(
        f"[Exp 3/3] H80-Exp3-Choquet-Submodular ({int(dots_e3.sum())} dots, {dt_e3:.1f}s): "
        f"HOLDOUT-DTI = {hold_e3['score']:.4f} CI95={hold_e3['ci95']} | "
        f"expected_DTI={alloc_stats['expected_dti']:.4f} | "
        f"max|rho|={uniq_e3['max_abs_spearman_surface']:.4f}, max_ov3px={uniq_e3['max_dot_fraction_within_3px']:.4f}, dup={uniq_e3['duplicate_lane']}",
        flush=True,
    )

    # Also write Holdout Exp 3 GeoTIFF via shared submission_writer internals so CLI
    # `scripts/evaluate_holdout.py` can independently verify the exact score on disk!
    holdout_tif = data / "cache_h80v2" / "holdout_exp3_choquet.tif"
    meta_h = tmpl["meta"].copy()
    meta_h["nodata"] = None
    import rasterio
    with rasterio.open(holdout_tif, "w", **meta_h) as dst:
        dst.write((dots_e3 & eval_mask).astype(np.float32), 1)

    np.savez_compressed(
        data / "cache_h80v2" / "holdout_experiment_maps.npz",
        p_phys=p_phys.astype(np.float16),
        p_struct=p_struct.astype(np.float16),
        p_fused=p_fused.astype(np.float16),
        dots_e1=np.packbits(dots_e1.ravel()),
        dots_e2=np.packbits(dots_e2.ravel()),
        dots_e3=np.packbits(dots_e3.ravel()),
    )

    out_report = {
        "budget_spent": "3 of 3 pre-registered experiments",
        "elapsed_seconds_total": round(time.time() - t0, 2),
        "registry_rasters_checked": len(registry),
        "leakage_canary_all_47_features": canary_report,
        "reference_random_control": rand_hold,
        "reference_prior_repo_h80c_v1": {
            "label": "HOLDOUT-DTI",
            "score": 0.0507,
            "note": "Prior GEMSDOE80 v1 unsupervised strain-transfer heuristic",
        },
        "experiments": [
            {
                "experiment_index": 1,
                "id": "H80-Exp1-MPhys-NMS",
                "hypothesis": "H80-G2 + H80-G3 pure 34-channel label-free multi-physics 4-quadrant OOF GBDT + NMS(sep=3)",
                "runtime_sec": round(dt_e1, 2),
                "n_dots": int(dots_e1.sum()),
                "holdout": hold_e1,
                "pre_placement_uniqueness": {
                    "duplicate_surface": e1_pre_u["duplicate_surface"],
                    "max_abs_spearman": e1_pre_u["max_abs_spearman"],
                },
                "uniqueness": uniq_e1,
                "model_diagnostics": phys_info,
            },
            {
                "experiment_index": 2,
                "id": "H80-Exp2-TwoRegime-NMS",
                "hypothesis": "H80-G1 + H80-G2 + H80-G3 Two-Regime OOF (nested self-supervised crack-tip cones + fracture-arrest + collinear gap bridges fused with M_phys) + NMS(sep=3)",
                "runtime_sec": round(dt_e2, 2),
                "n_dots": int(dots_e2.sum()),
                "holdout": hold_e2,
                "pre_placement_uniqueness": {
                    "duplicate_surface": e2_pre_u["duplicate_surface"],
                    "max_abs_spearman": e2_pre_u["max_abs_spearman"],
                },
                "uniqueness": uniq_e2,
                "model_diagnostics": struct_info,
            },
            {
                "experiment_index": 3,
                "id": "H80-Exp3-Choquet-Submodular",
                "hypothesis": "H80-G1..G4 Full Two-Regime OOF posterior + Exact 6-Shell Discrete Choquet-Survival Dinkelbach Fractional Submodular Allocator",
                "runtime_sec": round(dt_e3, 2),
                "n_dots": int(dots_e3.sum()),
                "holdout": hold_e3,
                "allocator_stats": alloc_stats,
                "pre_placement_uniqueness": {
                    "duplicate_surface": e3_pre_u["duplicate_surface"],
                    "max_abs_spearman": e3_pre_u["max_abs_spearman"],
                },
                "uniqueness": uniq_e3,
            },
        ],
    }
    report_path = ROOT / "research" / "h80_experiments_summary.json"
    report_path.write_text(json.dumps(out_report, indent=2) + "\n")
    print(f"[done] Saved full 3-experiment ablation ledger -> {report_path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
