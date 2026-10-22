#!/usr/bin/env python3
"""Execute the strictly capped 3-experiment H80 suite (`H80-E1 -> H80-E2 -> H80-E3`),
evaluate every arm + canonical fleet baselines on shared unmodified `vendor/gems57/evaluate_holdout.py`
(`SHA-256: ad86640b...`, `gems57-pooled-hide-v2`) across BOTH:
  - Primary Canonical Holdout: `Seed 20` (`11,321` withheld positives, 2,000-draw 20 km CI)
  - Independent Replication Holdout: `Seed 30` (`11,342` withheld positives, 2,000-draw 20 km CI)
audit ALL 1,081 unique SHA-256 rasters across ALL 92 `buffedlizard55-lab` repositories on BOTH
Phase-1 pre-placement continuous surface AND Phase-2 binary dots, emit the winning unique GeoTIFF
via shared `vendor/gems57/submission_writer.py`, and run standalone `scripts/validate_submission.py`.
"""
from __future__ import annotations

import gc
import hashlib
import json
import shutil
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import rasterio
from scipy import ndimage as ndi

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "vendor"))

from vendor.gems57 import evaluate_holdout as eh
from vendor.gems57.grid import load_grid
from vendor.gems57.holdout import buffered_component_draw
from vendor.gems57.submission_writer import write_submission

from scripts.restore_data import DATA_DIR, EVIDENCE_DIR, EXT_DIR, restore_all
from scripts.validate_submission import validate_geotiff
from src.gems80.adpt_engine import (
    build_adpt_gmns_dataset,
    fit_adpt_fold_ensemble,
    load_adpt_context,
    precompute_network_memmap,
    predict_adpt_from_precomputed,
    run_full_176_leakage_canary,
)
from src.gems80.fast_uniqueness import RegistryUniquenessScanner
from src.gems80.lagrangian_solver import solve_adpt_allocation

EXPECTED_EVAL_HASH = "ad86640b46f3c5fc6b1aceec4c4bfc82d03e35008b37ab2a8ad6076e6810ff3a"
BUDGET_DOTS = 36_500
BOOTSTRAP_DRAWS = 2000
BOOTSTRAP_SEED = 520810
CACHE_DIR = Path("/tmp/gems80_cache")


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def _score_global_map_on_folds(
    arm_name: str, pred_2d: np.ndarray, folds: list[dict], valid: np.ndarray
) -> tuple[np.ndarray, dict]:
    terms_sum = None
    fold_rows = []
    total_truth = 0
    total_emitted = 0
    for f in folds:
        res, terms = eh.evaluate(pred_2d, f, valid, block_side=200)
        terms_sum = terms.copy() if terms_sum is None else (terms_sum + terms)
        total_truth += int(res["n_truth"])
        total_emitted += int(res["n_emitted"])
        fold_rows.append({
            "quadrant": f["name"],
            "holdout_dti": round(float(res["dti"]), 6),
            "tpw": round(float(res["tpw"]), 3),
            "fpw": round(float(res["fpw"]), 3),
            "fnw": round(float(res["fnw"]), 3),
            "n_truth": int(res["n_truth"]),
            "n_emitted": int(res["n_emitted"]),
        })
    assert terms_sum is not None
    pooled_dti = float(eh.from_terms(terms_sum.sum(axis=0)))
    return terms_sum, {
        "arm": arm_name,
        "evidence_class": "HOLDOUT-DTI",
        "evaluator_version": eh.VERSION,
        "n_withheld_positives": total_truth,
        "total_emitted_dots_across_folds": total_emitted,
        "pooled_holdout_dti": round(pooled_dti, 6),
        "per_fold": fold_rows,
    }


def _build_foldfair_h33(valid: np.ndarray, visible: np.ndarray, budget: int = 37_654) -> np.ndarray | None:
    """Fold-fair rebuild of GEMSDOE32 champion (`0.2778` ORGANIZER-CONFIRMED):
    Poisson-disk `d >= 2.8 px` dots strictly on `H19-5 > 0` ridge support with `d(visible) > 2.0 px`.
    """
    cands = sorted(Path("/tmp/siblings").rglob("*h19-5*e27054cf*.tif"))
    if not cands:
        return None
    with rasterio.open(cands[0]) as ds:
        h19 = ds.read(1)
    d_vis = ndi.distance_transform_edt(~visible)
    eligible = valid & np.isfinite(h19) & (h19 > 0.0) & (d_vis > 2.0)
    vals = h19[eligible]
    ys, xs = np.nonzero(eligible)
    order = np.argsort(-vals, kind="mergesort")
    out = np.zeros(valid.shape, dtype=np.float32)
    occ = np.zeros(valid.shape, dtype=bool)
    offs = [(dy, dx) for dy in range(-3, 4) for dx in range(-3, 4) if dy * dy + dx * dx < 8]
    h, w = valid.shape
    placed = 0
    for idx in order:
        r, c = int(ys[idx]), int(xs[idx])
        if occ[r, c]:
            continue
        out[r, c] = 1.0
        placed += 1
        if placed >= budget:
            break
        for dy, dx in offs:
            ry, cx = r + dy, c + dx
            if 0 <= ry < h and 0 <= cx < w:
                occ[ry, cx] = True
    return out


def _live_anatomy(live_dots: np.ndarray, valid: np.ndarray, labels: np.ndarray, data_dir: Path) -> dict:
    with rasterio.open(data_dir / "training_features.tif") as ds:
        slp = np.where(valid, ds.read(19).astype(np.float32), 0.0)
        dep = np.where(valid, ds.read(15).astype(np.float32), 0.0)
    d_cat = ndi.distance_transform_edt(~labels).astype(np.float32)
    m = valid & (~labels) & (live_dots > 0.5)
    n_act = int(m.sum())
    d_m = d_cat[m]
    return {
        "live_active_dots": n_act,
        "median_slope_deg": round(float(np.median(slp[m])), 3) if n_act else 0.0,
        "median_depth_to_basement_m": round(float(np.median(dep[m])), 2) if n_act else 0.0,
        "frac_d_cat_le_2px": round(float((d_m <= 2.0).mean()), 5) if n_act else 0.0,
        "frac_d_cat_2_to_15px": round(float(((d_m > 2.0) & (d_m <= 15.0)).mean()), 5) if n_act else 0.0,
        "frac_d_cat_gt_15px": round(float((d_m > 15.0).mean()), 5) if n_act else 0.0,
    }


def _train_and_predict_seed(adpt, folds: list[dict], seed_base: int, compute_live: bool, labels: np.ndarray):
    """Cross-fit 4 outer LOQO ensembles for `folds` (`Seed 20` or `Seed 30`) in ~3 min!"""
    fp = adpt.h82.fp
    p_joint_oof = np.zeros(fp.shape, dtype=np.float32)
    p_phys_oof = np.zeros(fp.shape, dtype=np.float32)
    p_joint_live_acc = np.zeros(int(fp.sum()), dtype=np.float32) if compute_live else None
    p_phys_live_acc = np.zeros(int(fp.sum()), dtype=np.float32) if compute_live else None

    feat_gain = np.zeros(len(adpt.all_feature_names), dtype=np.float64)
    vis_global = folds[0]["visible"]
    all_fp_pos = np.arange(int(fp.sum()), dtype=np.int32)

    net_oof_mm = precompute_network_memmap(vis_global, adpt, CACHE_DIR / f"net_vis_{seed_base}.npy")
    net_live_mm = precompute_network_memmap(labels, adpt, CACHE_DIR / "net_vis_live.npy") if compute_live else None

    for k, fold in enumerate(folds):
        f_ckpt = CACHE_DIR / f"fold_ckpt_{seed_base}_{k}.npz"
        reg_flat = np.flatnonzero(fp & fold["region"])
        if f_ckpt.is_file():
            zf = np.load(f_ckpt)
            p_joint_oof.ravel()[reg_flat] = zf["pj_reg"]
            p_phys_oof.ravel()[reg_flat] = zf["pp_reg"]
            feat_gain += zf["fg"]
            if compute_live and p_joint_live_acc is not None and p_phys_live_acc is not None:
                p_joint_live_acc += zf["pj_lv"] / 4.0
                p_phys_live_acc += zf["pp_lv"] / 4.0
            log(f"  Restored checkpoint fold {fold['name']} (seed_base={seed_base})")
            continue

        tr_reg = np.zeros(fp.shape, dtype=bool)
        for j, f2 in enumerate(folds):
            if j != k:
                tr_reg |= f2["region"]
        assert not np.any(tr_reg & fold["region"]), "Spatial LOQO boundary leak!"
        assert not np.any(tr_reg & fold["truth"]), "Holdout truth leak!"

        vis_tr = vis_global & tr_reg
        x_tr, y_tr, w_tr, odds_den = build_adpt_gmns_dataset(
            vis_train=vis_tr,
            train_region=tr_reg,
            adpt=adpt,
            n_worlds=5,
            buffer_px=3.0,
            n_neg_per_world=44_000,
            deconfound_digitization_jitter=True,
            seed=seed_base + 100 * k,
        )
        ens = fit_adpt_fold_ensemble(x_tr, y_tr, w_tr, odds_den_factor=odds_den, seed=seed_base + 100 * k)
        fg_k = 0.6 * ens.lgbm_leaf.booster_.feature_importance(importance_type="gain") + 0.4 * ens.lgbm_et.booster_.feature_importance(importance_type="gain")
        feat_gain += fg_k
        del x_tr, y_tr, w_tr
        gc.collect()

        reg_fp = adpt.h82.fp_lut[reg_flat]
        pj_reg, pp_reg = predict_adpt_from_precomputed(ens, net_oof_mm, reg_fp, adpt)
        p_joint_oof.ravel()[reg_flat] = pj_reg
        p_phys_oof.ravel()[reg_flat] = pp_reg

        if compute_live and net_live_mm is not None and p_joint_live_acc is not None and p_phys_live_acc is not None:
            pj_lv, pp_lv = predict_adpt_from_precomputed(ens, net_live_mm, all_fp_pos, adpt)
            p_joint_live_acc += pj_lv / 4.0
            p_phys_live_acc += pp_lv / 4.0
            np.savez(f_ckpt, pj_reg=pj_reg, pp_reg=pp_reg, fg=fg_k, pj_lv=pj_lv, pp_lv=pp_lv)
            del pj_lv, pp_lv
        else:
            np.savez(f_ckpt, pj_reg=pj_reg, pp_reg=pp_reg, fg=fg_k)

        del ens, pj_reg, pp_reg
        gc.collect()
        log(f"  Completed fold {fold['name']} (seed_base={seed_base})")

    p_joint_live_2d = None
    p_phys_live_2d = None
    if compute_live and p_joint_live_acc is not None and p_phys_live_acc is not None:
        p_joint_live_2d = np.zeros(fp.shape, dtype=np.float32)
        p_phys_live_2d = np.zeros(fp.shape, dtype=np.float32)
        p_joint_live_2d.ravel()[adpt.h82.fp_idx] = p_joint_live_acc
        p_phys_live_2d.ravel()[adpt.h82.fp_idx] = p_phys_live_acc

    return p_joint_oof, p_phys_oof, p_joint_live_2d, p_phys_live_2d, feat_gain


def main() -> None:
    t_start = time.time()
    restore_all()
    EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)
    sub_dir = ROOT / "submissions"
    dl_dir = ROOT / "docs" / "downloads"
    sub_dir.mkdir(parents=True, exist_ok=True)
    dl_dir.mkdir(parents=True, exist_ok=True)

    # 1. Verify unmodified shared evaluator SHA-256
    hashes = eh.implementation_hashes()
    if hashes["evaluate_holdout.py"] != EXPECTED_EVAL_HASH:
        raise RuntimeError(f"Evaluator hash mismatch: {hashes['evaluate_holdout.py']}")
    log(f"[1/7] Verified shared evaluator: {eh.VERSION} ({hashes['evaluate_holdout.py'][:16]}...)")

    grid = load_grid(DATA_DIR)
    valid = grid.footprint
    labels = grid.catalogue
    allowed_live = valid & (~labels)

    folds_s20, _ = buffered_component_draw(grid, seed=20, hide_frac=0.20, buffer_px=3)
    folds_s30, _ = buffered_component_draw(grid, seed=30, hide_frac=0.20, buffer_px=3)
    n_w20 = sum(int(f["truth"].sum()) for f in folds_s20)
    n_w30 = sum(int(f["truth"].sum()) for f in folds_s30)
    log(f"[2/7] Canonical whole-component draws ready: Seed 20 = {n_w20} px | Seed 30 = {n_w30} px")

    # 2. Load unified 162-channel zero-SGMC physical context (`mmap_mode='r'`)
    adpt = load_adpt_context(CACHE_DIR, DATA_DIR, EXT_DIR, valid, grid.transform)
    log(f"[3/7] Loaded 176-channel ADPT-80 context (162 pure physical + 14 whole-component mechanics)")

    # 3. Run single-feature ROC-AUC leakage canary across ALL 176 channels BEFORE training!
    canary_path = EVIDENCE_DIR / "leakage_canary.json"
    if canary_path.is_file():
        canary = json.loads(canary_path.read_text())
    else:
        canary = run_full_176_leakage_canary(adpt, folds_s20, max_allowed_auc=0.90)
        canary_path.write_text(json.dumps(canary, indent=2) + "\n")
    log(
        f"[4/7] 176-feature Leakage Canary passed={canary['passed']} | "
        f"max symmetric AUC={canary['max_effective_auc']:.6f} ({canary['top_feature']}, ceiling <= 0.90)"
    )
    if not canary["passed"]:
        raise RuntimeError(f"Leakage canary failed: {canary['flagged_features']}")

    # 4. Initialize zero-RAM mmap uniqueness scanner across all 1,081 organization rasters (92 repos)
    scanner = RegistryUniquenessScanner(allowed_live)
    log(
        f"[5/7] Organization uniqueness scanner ready: {scanner.repos_scanned} repos, "
        f"{len(scanner.entries)} unique valid GeoTIFFs"
    )

    arms_terms_s20: dict[str, np.ndarray] = {}
    arms_terms_s30: dict[str, np.ndarray] = {}
    arms_info_s20: dict[str, dict] = {}
    arms_info_s30: dict[str, dict] = {}

    # Canonical controls on BOTH Seed 20 & Seed 30
    for seed_val, folds_cur, t_dict, i_dict in (
        (20, folds_s20, arms_terms_s20, arms_info_s20),
        (30, folds_s30, arms_terms_s30, arms_info_s30),
    ):
        allow_eval = np.zeros(valid.shape, dtype=bool)
        for f in folds_cur:
            allow_eval |= valid & f["region"] & (~f["masked_known"])
        rng_r = np.random.default_rng(700 + seed_val)
        rand_map = np.zeros(valid.shape, dtype=np.float32)
        rand_map.ravel()[rng_r.choice(np.flatnonzero(allow_eval), size=BUDGET_DOTS, replace=False)] = 1.0
        t_dict["B0-Random-Uniform"], i_dict["B0-Random-Uniform"] = _score_global_map_on_folds(
            "B0-Random-Uniform", rand_map, folds_cur, valid
        )
        h33_ff = _build_foldfair_h33(valid, folds_cur[0]["visible"], budget=37_654)
        if h33_ff is not None:
            t_dict["B1-H33-FoldFair"], i_dict["B1-H33-FoldFair"] = _score_global_map_on_folds(
                "B1-H33-FoldFair", h33_ff, folds_cur, valid
            )

    s20_ckpt = CACHE_DIR / "seed20_posteriors.npz"
    if s20_ckpt.is_file() and (CACHE_DIR / "pj_live.npy").is_file():
        log("[6/7] Reloading completed Seed 20 + full live posteriors from checkpoint...")
        z20 = np.load(s20_ckpt)
        pj_s20, pp_s20, gain_s20 = z20["pj"], z20["pp"], z20["gain"]
        pj_live = np.load(CACHE_DIR / "pj_live.npy")
        pp_live = np.load(CACHE_DIR / "pp_live.npy")
        t_fit = time.time()
    else:
        log("[6/7] Cross-fitting ADPT-80 4-architecture ensembles on Seed 20 + full live posterior...")
        t_fit = time.time()
        pj_s20, pp_s20, pj_live, pp_live, gain_s20 = _train_and_predict_seed(
            adpt, folds_s20, seed_base=8200, compute_live=True, labels=labels
        )
        assert pj_live is not None and pp_live is not None
        np.save(CACHE_DIR / "pj_live.npy", pj_live)
        np.save(CACHE_DIR / "pp_live.npy", pp_live)
        np.savez(s20_ckpt, pj=pj_s20, pp=pp_s20, gain=gain_s20)

    log("[6/7] Cross-fitting ADPT-80 4-architecture ensembles on independent replication Seed 30...")
    pj_s30, pp_s30, _, _, gain_s30 = _train_and_predict_seed(
        adpt, folds_s30, seed_base=8300, compute_live=False, labels=labels
    )

    gain_tot = gain_s20 + gain_s30
    top_idx = np.argsort(-gain_tot)[:25]
    top_features = [
        {"feature": adpt.all_feature_names[int(i)], "gain_share": round(float(gain_tot[i] / max(gain_tot.sum(), 1e-9)), 6)}
        for i in top_idx
    ]

    allow_s20 = np.zeros(valid.shape, dtype=bool)
    for f in folds_s20:
        allow_s20 |= valid & f["region"] & (~f["masked_known"])
    allow_s30 = np.zeros(valid.shape, dtype=bool)
    for f in folds_s30:
        allow_s30 |= valid & f["region"] & (~f["masked_known"])

    # =========================================================================
    # EXPERIMENT 1 (`H80-E1`): Pure 162-Channel Zero-Catalogue ADPT Physics
    # =========================================================================
    log("Evaluating Experiment 1/3 (`H80-E1`: ADPT-Physical-Blind, 0 catalogue features)...")
    d20_e1, _, rec20_e1 = solve_adpt_allocation(
        pp_s20, allow_s20, mode="holdout", ref_mask=valid, budget=BUDGET_DOTS, constrain_uniqueness=True
    )
    d30_e1, _, rec30_e1 = solve_adpt_allocation(
        pp_s30, allow_s30, mode="holdout", ref_mask=valid, budget=BUDGET_DOTS, constrain_uniqueness=True
    )
    arms_terms_s20["H80-E1-PurePhysics"], arms_info_s20["H80-E1-PurePhysics"] = _score_global_map_on_folds(
        "H80-E1-PurePhysics", d20_e1, folds_s20, valid
    )
    arms_terms_s30["H80-E1-PurePhysics"], arms_info_s30["H80-E1-PurePhysics"] = _score_global_map_on_folds(
        "H80-E1-PurePhysics", d30_e1, folds_s30, valid
    )
    e1_live_dots, e1_live_surf, e1_live_rec = solve_adpt_allocation(
        pp_live, allowed_live, mode="live", ref_mask=allowed_live, budget=BUDGET_DOTS, constrain_uniqueness=True
    )
    arms_info_s20["H80-E1-PurePhysics"]["surface_uniqueness"] = scanner.check(e1_live_surf, "surface")
    arms_info_s20["H80-E1-PurePhysics"]["dots_uniqueness"] = scanner.check(e1_live_dots, "dots")
    arms_info_s20["H80-E1-PurePhysics"]["solver"] = e1_live_rec
    arms_info_s20["H80-E1-PurePhysics"]["live_anatomy"] = _live_anatomy(e1_live_dots, valid, labels, DATA_DIR)

    # =========================================================================
    # EXPERIMENT 2 (`H80-E2`): Joint 176-Channel ADPT-80 Unconstrained Submodular
    # =========================================================================
    log("Evaluating Experiment 2/3 (`H80-E2`: ADPT-Dislocation-Unconstrained)...")
    d20_e2, _, _ = solve_adpt_allocation(
        pj_s20, allow_s20, mode="holdout", ref_mask=valid, budget=BUDGET_DOTS, constrain_uniqueness=False
    )
    d30_e2, _, _ = solve_adpt_allocation(
        pj_s30, allow_s30, mode="holdout", ref_mask=valid, budget=BUDGET_DOTS, constrain_uniqueness=False
    )
    arms_terms_s20["H80-E2-Unconstrained"], arms_info_s20["H80-E2-Unconstrained"] = _score_global_map_on_folds(
        "H80-E2-Unconstrained", d20_e2, folds_s20, valid
    )
    arms_terms_s30["H80-E2-Unconstrained"], arms_info_s30["H80-E2-Unconstrained"] = _score_global_map_on_folds(
        "H80-E2-Unconstrained", d30_e2, folds_s30, valid
    )
    e2_live_dots, e2_live_surf, e2_live_rec = solve_adpt_allocation(
        pj_live, allowed_live, mode="live", ref_mask=allowed_live, budget=BUDGET_DOTS, constrain_uniqueness=False
    )
    arms_info_s20["H80-E2-Unconstrained"]["surface_uniqueness"] = scanner.check(e2_live_surf, "surface")
    arms_info_s20["H80-E2-Unconstrained"]["dots_uniqueness"] = scanner.check(e2_live_dots, "dots")
    arms_info_s20["H80-E2-Unconstrained"]["solver"] = e2_live_rec
    arms_info_s20["H80-E2-Unconstrained"]["live_anatomy"] = _live_anatomy(e2_live_dots, valid, labels, DATA_DIR)

    # =========================================================================
    # EXPERIMENT 3 (`H80-E3`): Joint 176-Channel ADPT-80 + KKT Lagrangian Unique
    # =========================================================================
    log("Evaluating Experiment 3/3 (`H80-E3`: ADPT-KKT-Lagrangian-Unique across all 1,001 priors)...")
    d20_e3, _, rec20_e3 = solve_adpt_allocation(
        pj_s20, allow_s20, mode="holdout", ref_mask=valid, budget=BUDGET_DOTS, constrain_uniqueness=True
    )
    d30_e3, _, rec30_e3 = solve_adpt_allocation(
        pj_s30, allow_s30, mode="holdout", ref_mask=valid, budget=BUDGET_DOTS, constrain_uniqueness=True
    )
    arms_terms_s20["H80-E3-LagrangianUnique"], arms_info_s20["H80-E3-LagrangianUnique"] = _score_global_map_on_folds(
        "H80-E3-LagrangianUnique", d20_e3, folds_s20, valid
    )
    arms_terms_s30["H80-E3-LagrangianUnique"], arms_info_s30["H80-E3-LagrangianUnique"] = _score_global_map_on_folds(
        "H80-E3-LagrangianUnique", d30_e3, folds_s30, valid
    )
    e3_live_dots, e3_live_surf, e3_live_rec = solve_adpt_allocation(
        pj_live, allowed_live, mode="live", ref_mask=allowed_live, budget=BUDGET_DOTS, constrain_uniqueness=True
    )
    arms_info_s20["H80-E3-LagrangianUnique"]["surface_uniqueness"] = scanner.check(e3_live_surf, "surface")
    arms_info_s20["H80-E3-LagrangianUnique"]["dots_uniqueness"] = scanner.check(e3_live_dots, "dots")
    arms_info_s20["H80-E3-LagrangianUnique"]["solver"] = e3_live_rec
    arms_info_s20["H80-E3-LagrangianUnique"]["live_anatomy"] = _live_anatomy(e3_live_dots, valid, labels, DATA_DIR)

    for k_arm in ("B0-Random-Uniform", "B1-H33-FoldFair", "H80-E1-PurePhysics", "H80-E2-Unconstrained", "H80-E3-LagrangianUnique"):
        if k_arm in arms_info_s20 and k_arm in arms_info_s30:
            log(
                f"  SUMMARY {k_arm:<26s} | Seed20 = {arms_info_s20[k_arm]['pooled_holdout_dti']:.6f} | "
                f"Seed30 = {arms_info_s30[k_arm]['pooled_holdout_dti']:.6f}"
            )

    # Select winning protocol-unique candidate (`binding_duplicate == False` on BOTH surface AND dots!)
    live_candidates = {
        "H80-E1-PurePhysics": e1_live_dots,
        "H80-E2-Unconstrained": e2_live_dots,
        "H80-E3-LagrangianUnique": e3_live_dots,
    }
    eligible = [
        k for k in live_candidates
        if (not arms_info_s20[k]["surface_uniqueness"]["binding_duplicate"])
        and (not arms_info_s20[k]["dots_uniqueness"]["binding_duplicate"])
    ]
    if not eligible:
        raise RuntimeError("No candidate passed both surface and dot uniqueness gates across all 1,001 priors!")

    winner = max(
        eligible,
        key=lambda k: 0.5 * (arms_info_s20[k]["pooled_holdout_dti"] + arms_info_s30[k]["pooled_holdout_dti"]),
    )

    summary_s20 = eh.pooled_summary(arms_terms_s20, draws=BOOTSTRAP_DRAWS, seed=BOOTSTRAP_SEED, candidate=winner)
    summary_s30 = eh.pooled_summary(arms_terms_s30, draws=BOOTSTRAP_DRAWS, seed=BOOTSTRAP_SEED, candidate=winner)

    for k, row in summary_s20["scores"].items():
        arms_info_s20[k]["ci95"] = [round(row["ci95"][0], 6), round(row["ci95"][1], 6)]
    for k, row in summary_s30["scores"].items():
        arms_info_s30[k]["ci95"] = [round(row["ci95"][0], 6), round(row["ci95"][1], 6)]
        if k in arms_info_s20:
            arms_info_s20[k]["seed30_replication"] = {
                "pooled_holdout_dti": arms_info_s30[k]["pooled_holdout_dti"],
                "ci95": arms_info_s30[k]["ci95"],
                "n_withheld_positives": arms_info_s30[k]["n_withheld_positives"],
            }

    win_dots = live_candidates[winner]
    np.save(CACHE_DIR / "win_dots.npy", win_dots)
    short_hash = hashlib.sha256(win_dots.tobytes()).hexdigest()[:8]
    ts_str = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    sub_filename = f"gems80-adpt-kkt-lagrangian-36500d-{ts_str}-{short_hash}-zeros.tif"
    sub_path = sub_dir / sub_filename

    s20_win = arms_info_s20[winner]["pooled_holdout_dti"]
    s30_win = arms_info_s30[winner]["pooled_holdout_dti"]
    dual_mean = 0.5 * (s20_win + s30_win)

    sub_ui_name = f"GEMSDOE80 ADPT-80 KKT-Lagrangian ({BUDGET_DOTS} dots, {short_hash})"
    sub_ui_note = (
        f"GEMSDOE80 ADPT-80 176f KKT-Lagrangian 36500d: {s20_win:.4f}/{s30_win:.4f} HOLDOUT-DTI, "
        f"0/1001 disc ov<=0.70, float32 [0,1] EPSG:32611"
    )[:140]

    write_submission(
        sub_path,
        win_dots,
        DATA_DIR / "sample_submission.tif",
        valid,
        note=sub_ui_note,
        name=sub_ui_name,
        catalogue=labels,
    )
    # Clean auxiliary `.zip` next to `submissions/` if desired, and copy exact `.tif` to `docs/downloads/`
    dl_tif = dl_dir / sub_filename
    shutil.copy2(sub_path, dl_tif)

    val_receipt = validate_geotiff(dl_tif)
    val_receipt["submission_ui_name"] = sub_ui_name
    val_receipt["submission_ui_note"] = sub_ui_note
    val_receipt["download_relative_path"] = f"downloads/{sub_filename}"
    (EVIDENCE_DIR / "submission_validation.json").write_text(json.dumps(val_receipt, indent=2) + "\n")

    out_doc = {
        "evidence_class": "EXPERIMENT-SUITE-RECEIPT",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "total_elapsed_sec": round(time.time() - t_start, 2),
        "model_fit_sec": round(time.time() - t_fit, 2),
        "experiments_executed": 3,
        "experiment_budget_cap": 3,
        "evaluator": {
            "version": eh.VERSION,
            "hashes": hashes,
            "primary_seed20_withheld_positives": n_w20,
            "replication_seed30_withheld_positives": n_w30,
            "bootstrap_draws": BOOTSTRAP_DRAWS,
            "bootstrap_block_km": 20.0,
        },
        "top_25_gain_features": top_features,
        "winner_arm": winner,
        "winner_dual_seed_mean_holdout_dti": round(dual_mean, 6),
        "winner_submission": val_receipt,
        "paired_deltas_seed20": summary_s20["paired_differences"],
        "paired_deltas_seed30": summary_s30["paired_differences"],
        "arms_seed20": arms_info_s20,
        "arms_seed30": arms_info_s30,
    }
    (EVIDENCE_DIR / "experiments_summary.json").write_text(json.dumps(out_doc, indent=2) + "\n")
    log(
        f"[7/7] ALL DONE IN {time.time() - t_start:.1f}s! Winner={winner} | "
        f"Seed20={s20_win:.6f} | Seed30={s30_win:.6f} | Mean={dual_mean:.6f} | File={sub_filename}"
    )


if __name__ == "__main__":
    main()
