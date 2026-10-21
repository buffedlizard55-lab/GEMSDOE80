#!/usr/bin/env python3
"""Execute the strictly capped 3-experiment H80 suite (`E1 -> E2 -> E3`),
evaluate every arm + canonical fleet baselines on shared `vendor/gems57/evaluate_holdout.py`
(`SHA-256: ad86640b...`, `gems57-pooled-hide-v2`) across BOTH:
  - Primary Canonical Holdout: `Seed 20` (`11,321` withheld positives, 1,000-draw 20 km CI)
  - Independent Replication Holdout: `Seed 30` (`11,342` withheld positives, 1,000-draw 20 km CI)
audit all 454 registry rasters on BOTH surface and dots, emit the winning unique GeoTIFF via
shared `vendor/gems57/submission_writer.py`, and validate it without human intervention.
"""
from __future__ import annotations

import gc
import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import rasterio
from scipy import ndimage as ndi

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from vendor.gems57 import evaluate_holdout as eh
from vendor.gems57.grid import Grid, load_grid
from vendor.gems57.holdout import buffered_component_draw
from vendor.gems57.submission_writer import write_submission

from scripts.restore_data import DATA_DIR, EVIDENCE_DIR, EXT_DIR, restore_all
from scripts.validate_submission import validate_geotiff
from src.gems80.allocator import (
    dinkelbach_submodular_sequence,
    priority_disk_thin,
    select_threshold_on_nested_subfold,
)
from src.gems80.fast_uniqueness import RegistryUniquenessScanner
from src.gems80.features import build_or_load_static_features
from src.gems80.model import (
    _quadrant_partition_masks,
    fit_stage1_physical_oof,
    fit_stage2_fracture_meta_and_calibrators,
    predict_stage1_ensemble_all_valid,
    predict_stage2_for_visible_network,
    run_leakage_canary,
)

EXPECTED_EVAL_HASH = "ad86640b46f3c5fc6b1aceec4c4bfc82d03e35008b37ab2a8ad6076e6810ff3a"
LIVE_REF_DENSITY = 37654.0 / 5106385.0  # Exact active density of H33-2-B2 (`0.0073739`)


def _sigmoid01(x: np.ndarray) -> np.ndarray:
    c = np.clip(x.astype(np.float64), -25.0, 25.0)
    return (1.0 / (1.0 + np.exp(-c))).astype(np.float32)


def _eval_arm(arm_name: str, fold_preds: list[np.ndarray], folds: list[dict], valid: np.ndarray) -> tuple[np.ndarray, dict]:
    terms_sum = None
    fold_rows = []
    total_truth = 0
    total_emitted = 0
    for q_idx, (pred, f) in enumerate(zip(fold_preds, folds)):
        res, terms = eh.evaluate(pred, f, valid, block_side=200)
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


def _live_anatomy(
    live_dots: np.ndarray, valid: np.ndarray, labels: np.ndarray, data_dir: Path
) -> dict:
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


def main() -> None:
    t_start = time.time()
    restore_all()
    EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)
    sub_dir = ROOT / "submissions"
    sub_dir.mkdir(parents=True, exist_ok=True)

    # 1. Verify shared evaluator SHA-256 bit-for-bit
    hashes = eh.implementation_hashes()
    if hashes["evaluate_holdout.py"] != EXPECTED_EVAL_HASH:
        raise RuntimeError(f"Evaluator hash mismatch: {hashes['evaluate_holdout.py']}")
    print(f"[1/7] Shared evaluator verified: {eh.VERSION} ({hashes['evaluate_holdout.py'][:16]}...)", flush=True)

    grid = load_grid(DATA_DIR)
    valid = grid.footprint
    labels = grid.catalogue
    allowed_live = valid & (~labels)
    n_valid = int(valid.sum())
    quad_masks = _quadrant_partition_masks(valid)

    folds_s20, _ = buffered_component_draw(grid, seed=20, hide_frac=0.20, buffer_px=3)
    folds_s30, _ = buffered_component_draw(grid, seed=30, hide_frac=0.20, buffer_px=3)
    n_w20 = sum(int(f["truth"].sum()) for f in folds_s20)
    n_w30 = sum(int(f["truth"].sum()) for f in folds_s30)
    print(f"[2/7] Canonical whole-component draws ready: Seed 20 = {n_w20} px | Seed 30 = {n_w30} px", flush=True)

    # 2. Build/load 52-channel zero-catalogue static feature matrix
    cache_path = Path("/tmp/gems80_cache/x_static.npy")
    X_static, _ = build_or_load_static_features(valid, DATA_DIR, EXT_DIR, cache_path)
    print(f"[3/7] Static zero-catalogue matrix ready: {X_static.shape} ({time.time() - t_start:.1f}s)", flush=True)

    # 3. Single-feature Leakage Canary BEFORE training any model!
    canary = run_leakage_canary(X_static, folds_s20, valid, max_allowed_auc=0.90)
    (EVIDENCE_DIR / "leakage_canary.json").write_text(json.dumps(canary, indent=2) + "\n")
    print(
        f"[4/7] Leakage canary passed={canary['passed']} | max single-feature AUC={canary['max_effective_auc']} "
        f"({canary['top_feature']}, threshold <= 0.90)",
        flush=True,
    )
    if not canary["passed"]:
        raise RuntimeError(f"Leakage canary failed: {canary['flagged_features']}")

    # Pre-index all 454 unique SHA-256 registry rasters for bit-exact surface + dot gates
    print("[5/7] Indexing 454 unique SHA-256 registry rasters for surface + dot uniqueness gates...", flush=True)
    scanner = RegistryUniquenessScanner(allowed_live, max_cells=250_000, seed=20)

    arms_terms_s20: dict[str, np.ndarray] = {}
    arms_terms_s30: dict[str, np.ndarray] = {}
    arms_info: dict[str, dict] = {}
    arms_info_s30: dict[str, dict] = {}

    # Baseline B0: Deterministic Random Matched-Density Control
    for seed_num, f_list, t_dict, i_dict in (
        (20, folds_s20, arms_terms_s20, arms_info),
        (30, folds_s30, arms_terms_s30, arms_info_s30),
    ):
        rng_ctrl = np.random.default_rng(seed_num)
        ctrl_preds = []
        for f in f_list:
            act_f = valid & f["region"] & (~f["masked_known"])
            k_f = int(round(LIVE_REF_DENSITY * act_f.sum()))
            p0 = np.zeros(valid.shape, dtype=np.float32)
            idx0 = rng_ctrl.choice(np.flatnonzero(act_f), size=k_f, replace=False)
            p0.ravel()[idx0] = 1.0
            ctrl_preds.append(p0)
        t_dict["B0-Random-Matched"], i_dict["B0-Random-Matched"] = _eval_arm(
            "B0-Random-Matched", ctrl_preds, f_list, valid
        )
        del ctrl_preds

    # Baseline B1: H33-2-B2 (`0.2778` ORGANIZER-CONFIRMED live leader, masked on each fold's active cells)
    h33_path = Path("/tmp/siblings/GEMSDOE32/submissions/h33-h33-2-b2-20261004T220000Z-e5eb6e7e-zeros.tif")
    if h33_path.is_file():
        with rasterio.open(h33_path) as ds:
            h33_arr = np.clip(np.nan_to_num(ds.read(1).astype(np.float32), nan=0.0), 0.0, 1.0)
        arms_terms_s20["B1-H33-Live02778"], arms_info["B1-H33-Live02778"] = _eval_arm(
            "B1-H33-Live02778", [h33_arr] * 4, folds_s20, valid
        )
        arms_terms_s30["B1-H33-Live02778"], arms_info_s30["B1-H33-Live02778"] = _eval_arm(
            "B1-H33-Live02778", [h33_arr] * 4, folds_s30, valid
        )
        del h33_arr

    # Baseline B2: GEMSDOE70 P8-OOF38k (`0.2478` ORGANIZER-CONFIRMED live, `0.1113` OOF benchmark)
    g70_path = Path("/tmp/siblings/GEMSDOE70/submissions/p8-oof-submod-38k-20261015T062535Z-7045d5ea-zeros.tif")
    if g70_path.is_file():
        with rasterio.open(g70_path) as ds:
            g70_arr = np.clip(np.nan_to_num(ds.read(1).astype(np.float32), nan=0.0), 0.0, 1.0)
        arms_terms_s20["B2-G70-LiveSub"], arms_info["B2-G70-LiveSub"] = _eval_arm(
            "B2-G70-LiveSub", [g70_arr] * 4, folds_s20, valid
        )
        arms_terms_s30["B2-G70-LiveSub"], arms_info_s30["B2-G70-LiveSub"] = _eval_arm(
            "B2-G70-LiveSub", [g70_arr] * 4, folds_s30, valid
        )
        del g70_arr

    # =========================================================================
    # EXPERIMENT 1 (`H80-E1`): ADPT-Physical-OOF (Zero-Catalogue 52-Channel GBDT+LGBM)
    # =========================================================================
    print("[6/7] Running Experiment 1/3 (`H80-E1`: ADPT-Physical-OOF)...", flush=True)
    t_e1 = time.time()
    s1_oof_v, s1_models, s1_feat_imp = fit_stage1_physical_oof(X_static, labels, valid, quad_masks)
    s1_oof_2d = np.full(valid.shape, -30.0, dtype=np.float32)
    s1_oof_2d[valid] = s1_oof_v

    for f_list, t_dict, i_dict in (
        (folds_s20, arms_terms_s20, arms_info),
        (folds_s30, arms_terms_s30, arms_info_s30),
    ):
        preds_e1 = []
        for f in f_list:
            act_f = valid & f["region"] & (~f["masked_known"])
            k_f = int(round(LIVE_REF_DENSITY * act_f.sum()))
            preds_e1.append(priority_disk_thin(s1_oof_2d, act_f, n_budget=k_f, min_spacing_px=2.85))
        t_dict["H80-E1"], i_dict["H80-E1"] = _eval_arm("H80-E1", preds_e1, f_list, valid)
        del preds_e1

    arms_info["H80-E1"]["elapsed_sec"] = round(time.time() - t_e1, 2)

    s1_live_v = predict_stage1_ensemble_all_valid(X_static, s1_models, n_valid)
    s1_live_surf = np.zeros(valid.shape, dtype=np.float32)
    s1_live_surf[allowed_live] = _sigmoid01(s1_live_v)[~labels[valid]]
    e1_live_dots = priority_disk_thin(s1_live_surf, allowed_live, n_budget=37_654, min_spacing_px=2.85)
    arms_info["H80-E1"]["surface_uniqueness"] = scanner.check(s1_live_surf, "surface")
    arms_info["H80-E1"]["dots_uniqueness"] = scanner.check(e1_live_dots, "dots")
    arms_info["H80-E1"]["live_anatomy"] = _live_anatomy(e1_live_dots, valid, labels, DATA_DIR)
    print(
        f"      -> H80-E1 Seed20={arms_info['H80-E1']['pooled_holdout_dti']:.6f} | "
        f"Seed30={arms_info_s30['H80-E1']['pooled_holdout_dti']:.6f} | "
        f"disc_ov={arms_info['H80-E1']['dots_uniqueness']['max_discriminating_overlap_3px']:.4f}",
        flush=True,
    )

    # =========================================================================
    # EXPERIMENT 2 (`H80-E2`): ADPT-Curvilinear-FractureMeta
    # =========================================================================
    print("[6/7] Running Experiment 2/3 (`H80-E2`: ADPT-Curvilinear-FractureMeta)...", flush=True)
    t_e2 = time.time()
    (
        meta_oof_s20_v,
        p_oof_s20_v,
        q_oof_s20_v,
        meta_models,
        p_cals,
        q_cals,
    ) = fit_stage2_fracture_meta_and_calibrators(s1_oof_v, X_static, grid, folds_s20, quad_masks)

    # Independent Seed-30 OOF fracture evaluation using exact same quadrant models & Seed-30's `visible`!
    meta_oof_s30_v, p_oof_s30_v, q_oof_s30_v = predict_stage2_for_visible_network(
        s1_oof_v, X_static, folds_s30[0]["visible"], valid, meta_models, p_cals, q_cals, quad_masks=quad_masks
    )

    for lg_v, f_list, t_dict, i_dict in (
        (meta_oof_s20_v, folds_s20, arms_terms_s20, arms_info),
        (meta_oof_s30_v, folds_s30, arms_terms_s30, arms_info_s30),
    ):
        m2d = np.full(valid.shape, -30.0, dtype=np.float32)
        m2d[valid] = lg_v
        preds_e2 = []
        for f in f_list:
            act_f = valid & f["region"] & (~f["masked_known"])
            k_f = int(round(LIVE_REF_DENSITY * act_f.sum()))
            preds_e2.append(priority_disk_thin(m2d, act_f, n_budget=k_f, min_spacing_px=2.85))
        t_dict["H80-E2"], i_dict["H80-E2"] = _eval_arm("H80-E2", preds_e2, f_list, valid)
        del preds_e2, m2d

    arms_info["H80-E2"]["elapsed_sec"] = round(time.time() - t_e2, 2)

    meta_live_logit_v, p_live_v, q_live_v = predict_stage2_for_visible_network(
        s1_live_v, X_static, labels, valid, meta_models, p_cals, q_cals, quad_masks=None
    )
    meta_live_surf = np.zeros(valid.shape, dtype=np.float32)
    meta_live_surf[allowed_live] = p_live_v[~labels[valid]]
    e2_live_dots = priority_disk_thin(meta_live_surf, allowed_live, n_budget=37_654, min_spacing_px=2.85)
    arms_info["H80-E2"]["surface_uniqueness"] = scanner.check(meta_live_surf, "surface")
    arms_info["H80-E2"]["dots_uniqueness"] = scanner.check(e2_live_dots, "dots")
    arms_info["H80-E2"]["live_anatomy"] = _live_anatomy(e2_live_dots, valid, labels, DATA_DIR)
    print(
        f"      -> H80-E2 Seed20={arms_info['H80-E2']['pooled_holdout_dti']:.6f} | "
        f"Seed30={arms_info_s30['H80-E2']['pooled_holdout_dti']:.6f} | "
        f"disc_ov={arms_info['H80-E2']['dots_uniqueness']['max_discriminating_overlap_3px']:.4f}",
        flush=True,
    )

    # =========================================================================
    # EXPERIMENT 3 (`H80-E3`): ADPT-Dinkelbach-Submodular (Theorem 3 Exact Allocator)
    # =========================================================================
    print("[6/7] Running Experiment 3/3 (`H80-E3`: ADPT-Dinkelbach-Submodular)...", flush=True)
    t_e3 = time.time()
    p_s20_2d = np.zeros(valid.shape, dtype=np.float32)
    q_s20_2d = np.zeros(valid.shape, dtype=np.float32)
    p_s20_2d[valid] = p_oof_s20_v
    q_s20_2d[valid] = q_oof_s20_v

    p_s30_2d = np.zeros(valid.shape, dtype=np.float32)
    q_s30_2d = np.zeros(valid.shape, dtype=np.float32)
    p_s30_2d[valid] = p_oof_s30_v
    q_s30_2d[valid] = q_oof_s30_v

    subfold_calibrations = []
    preds_e3_s20 = []
    preds_e3_s30 = []
    for q_idx, q_mask in enumerate(quad_masks):
        # Strictly inside quadrant `q_idx`'s training complement (`d_q >= 12`), tune `(gamma*, frac*)` on whole-component sub-draw!
        d_q = ndi.distance_transform_edt(~q_mask).astype(np.float32)
        train_dom = valid & (d_q >= 12.0)
        sub_grid = Grid(
            footprint=train_dom,
            catalogue=folds_s20[q_idx]["visible"] & train_dom,
            crs=grid.crs,
            transform=grid.transform,
            height=grid.height,
            width=grid.width,
        )
        sf_list, _ = buffered_component_draw(sub_grid, seed=95 + q_idx, hide_frac=0.20, buffer_px=3, boundary_px=6)
        sub_truth = np.zeros(valid.shape, dtype=bool)
        sub_reg_u = np.zeros(valid.shape, dtype=bool)
        for sf in sf_list:
            sub_truth |= sf["truth"]
            sub_reg_u |= sf["region"]
        sub_act = train_dom & sub_reg_u & (~sf_list[0]["masked_known"])

        best_cfg = {"best_dti": -1.0, "gamma": 1.35, "best_frac": LIVE_REF_DENSITY}
        for gam in (1.15, 1.35, 1.55):
            cal_res = select_threshold_on_nested_subfold(p_s20_2d, q_s20_2d, sub_act, sub_truth, gamma=gam)
            if cal_res["best_dti"] > best_cfg["best_dti"]:
                best_cfg = {"gamma": gam, **cal_res}
        subfold_calibrations.append({"fold": folds_s20[q_idx]["name"], **best_cfg})

        act20 = valid & folds_s20[q_idx]["region"] & (~folds_s20[q_idx]["masked_known"])
        k20 = int(round(best_cfg["best_frac"] * act20.sum()))
        sy20, sx20, _ = dinkelbach_submodular_sequence(
            p_s20_2d, q_s20_2d, act20, gamma=best_cfg["gamma"], max_select=k20, candidate_cap=k20 * 5
        )
        pr20 = np.zeros(valid.shape, dtype=np.float32)
        pr20[sy20, sx20] = 1.0
        preds_e3_s20.append(pr20)

        act30 = valid & folds_s30[q_idx]["region"] & (~folds_s30[q_idx]["masked_known"])
        k30 = int(round(best_cfg["best_frac"] * act30.sum()))
        sy30, sx30, _ = dinkelbach_submodular_sequence(
            p_s30_2d, q_s30_2d, act30, gamma=best_cfg["gamma"], max_select=k30, candidate_cap=k30 * 5
        )
        pr30 = np.zeros(valid.shape, dtype=np.float32)
        pr30[sy30, sx30] = 1.0
        preds_e3_s30.append(pr30)

    arms_terms_s20["H80-E3"], arms_info["H80-E3"] = _eval_arm("H80-E3", preds_e3_s20, folds_s20, valid)
    arms_terms_s30["H80-E3"], arms_info_s30["H80-E3"] = _eval_arm("H80-E3", preds_e3_s30, folds_s30, valid)
    arms_info["H80-E3"]["elapsed_sec"] = round(time.time() - t_e3, 2)
    arms_info["H80-E3"]["nested_subfold_calibrations"] = subfold_calibrations
    del preds_e3_s20, preds_e3_s30

    live_gamma = float(np.median([c["gamma"] for c in subfold_calibrations]))
    live_frac = float(np.mean([c["best_frac"] for c in subfold_calibrations]))
    n_live_e3 = int(round(live_frac * allowed_live.sum()))
    p_live_2d = np.zeros(valid.shape, dtype=np.float32)
    q_live_2d = np.zeros(valid.shape, dtype=np.float32)
    p_live_2d[valid] = p_live_v
    q_live_2d[valid] = q_live_v
    sy_l, sx_l, _ = dinkelbach_submodular_sequence(
        p_live_2d, q_live_2d, allowed_live, gamma=live_gamma, max_select=n_live_e3, candidate_cap=n_live_e3 * 5
    )
    e3_live_dots = np.zeros(valid.shape, dtype=np.float32)
    e3_live_dots[sy_l, sx_l] = 1.0

    arms_info["H80-E3"]["surface_uniqueness"] = scanner.check(meta_live_surf, "surface")
    arms_info["H80-E3"]["dots_uniqueness"] = scanner.check(e3_live_dots, "dots")
    arms_info["H80-E3"]["live_anatomy"] = _live_anatomy(e3_live_dots, valid, labels, DATA_DIR)
    print(
        f"      -> H80-E3 Seed20={arms_info['H80-E3']['pooled_holdout_dti']:.6f} | "
        f"Seed30={arms_info_s30['H80-E3']['pooled_holdout_dti']:.6f} | "
        f"live_dots={n_live_e3} | disc_ov={arms_info['H80-E3']['dots_uniqueness']['max_discriminating_overlap_3px']:.4f}",
        flush=True,
    )

    # =========================================================================
    # 7. Exact Paired 20 km Physical Block Bootstrap CIs (`eh.pooled_summary`)
    # =========================================================================
    candidate_arms = {"H80-E1": e1_live_dots, "H80-E2": e2_live_dots, "H80-E3": e3_live_dots}
    eligible = [
        k for k in ("H80-E1", "H80-E2", "H80-E3")
        if (not arms_info[k]["surface_uniqueness"]["binding_duplicate"])
        and (not arms_info[k]["dots_uniqueness"]["binding_duplicate"])
    ]
    if not eligible:
        raise RuntimeError("No candidate passed both surface and dot uniqueness gates!")

    winner = max(eligible, key=lambda k: arms_info[k]["pooled_holdout_dti"])
    summary_s20 = eh.pooled_summary(arms_terms_s20, draws=1000, seed=520810, candidate=winner)
    summary_s30 = eh.pooled_summary(arms_terms_s30, draws=1000, seed=520810, candidate=winner)

    for k, ci_row in summary_s20["scores"].items():
        arms_info[k]["ci95"] = [round(ci_row["ci95"][0], 6), round(ci_row["ci95"][1], 6)]
    for k, ci_row in summary_s30["scores"].items():
        arms_info_s30[k]["ci95"] = [round(ci_row["ci95"][0], 6), round(ci_row["ci95"][1], 6)]
        if k in arms_info:
            arms_info[k]["seed30_replication"] = {
                "pooled_holdout_dti": arms_info_s30[k]["pooled_holdout_dti"],
                "ci95": arms_info_s30[k]["ci95"],
                "n_withheld_positives": arms_info_s30[k]["n_withheld_positives"],
            }

    # Emit winning unique GeoTIFF via shared `vendor/gems57/submission_writer.py`
    ts_str = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    arm_slug = winner.lower()
    sub_name = f"adpt80-{arm_slug}-{ts_str}-zeros.tif"
    sub_path = sub_dir / sub_name
    win_dots = candidate_arms[winner]
    n_win_dots = int((win_dots > 0.5).sum())

    sub_ui_name = f"GEMSDOE80 ADPT-80 {winner} ({n_win_dots:,} dots)"
    sub_ui_note = (
        f"ADPT-80 {winner}: 360-deg dip-polarity tomography + Radon fracture meta + "
        f"submodular DTI ({arms_info[winner]['pooled_holdout_dti']:.4f} HOLDOUT-DTI)."
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

    val_receipt = validate_geotiff(sub_path)
    val_receipt["submission_ui_name"] = sub_ui_name
    val_receipt["submission_ui_note"] = sub_ui_note
    (EVIDENCE_DIR / "submission_validation.json").write_text(json.dumps(val_receipt, indent=2) + "\n")

    results_doc = {
        "evidence_class": "EXPERIMENT-SUITE-RECEIPT",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "total_elapsed_sec": round(time.time() - t_start, 2),
        "experiments_executed": 3,
        "experiment_budget_cap": 3,
        "evaluator": {
            "version": eh.VERSION,
            "hashes": hashes,
            "primary_seed20_withheld_positives": n_w20,
            "replication_seed30_withheld_positives": n_w30,
            "bootstrap_draws": 1000,
            "bootstrap_block_km": 20.0,
        },
        "stage1_top_features": s1_feat_imp["top_20_gain_features"],
        "winner_arm": winner,
        "winner_submission": val_receipt,
        "paired_deltas_seed20": summary_s20["paired_differences"],
        "paired_deltas_seed30": summary_s30["paired_differences"],
        "arms_seed20": arms_info,
        "arms_seed30": arms_info_s30,
    }
    (EVIDENCE_DIR / "experiments_summary.json").write_text(json.dumps(results_doc, indent=2) + "\n")
    print(f"[7/7] COMPLETE in {time.time() - t_start:.1f}s! Winner={winner} -> {sub_path.name}", flush=True)


if __name__ == "__main__":
    main()
