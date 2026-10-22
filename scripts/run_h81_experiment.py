#!/usr/bin/env python3
"""GEMSDOE80 session 2026-10-22 — ONE new experiment (H81), preregistered here.

Candidates (full ranking in research/h81_candidates_20261022.md):
  H81-1 (TESTED HERE): clay-cap / magnetite-destruction conduit gate.
    Surface = OOF self-masked meta surface x (0.60 + 0.40 * conduit_rank),
    conduit_rank = mean rank01 of static-stack channels
    [mt_clay_cap_feeder(18), mt_cond_grad_hp(19), mt_cond_worm_1500m(20),
     hydro_demag_trough(22), hydro_shallow_demag_u150(23),
     strain_transtension(48)].  Emission identical to the E80-2 control
    (submodular_pack_300m, skeleton prefilter, collar d_vis>2.40, N=37,654
    per full grid) so the ONLY difference vs control is the conduit gate.
    Non-fault mimic: playa evaporite / lacustrine clay conductivity
    (suppressed because the gate only re-weights meta-surface dots, which
    already require fracture-transfer support, and collar excludes catalogue).
  Controls (not new experiments): E80-2 pure meta transfer and the shipped
    E80-3 hero hybrid allocator, re-evaluated on identical folds/terms.

Promotion rule (fixed before running): ship H81-1 only if its pooled
HOLDOUT-DTI (gems57-pooled-hide-v2, K=22,663) exceeds BOTH controls.
Otherwise ship the E80-2 pure-transfer emission (itself a new, unique dot
set vs the shipped hero) only if it beats the hero control pooled; if
neither beats the hero, verdict=negative and no new file is promoted.
Every number is HOLDOUT-DTI; no live/organizer score is claimed.
"""
from __future__ import annotations

import gc
import hashlib
import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np
import rasterio
from scipy import ndimage as ndi
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "vendor"))
sys.path.insert(0, str(ROOT / "scripts"))

from gems57 import evaluate_holdout as EH, grid as GRID, holdout as HO, submission_writer as SW  # noqa: E402
from gems80.features import build_or_load_static_stack  # noqa: E402
from gems80.forensics import compute_spatial_profile, run_full_uniqueness_audit  # noqa: E402
from gems80.models import build_full_catalogue_surfaces, run_holdout_oof_for_seed, submodular_pack_300m, _rank01  # noqa: E402

spec = importlib.util.spec_from_file_location("run_experiments", ROOT / "scripts" / "run_experiments.py")
RE = importlib.util.module_from_spec(spec)
spec.loader.exec_module(RE)  # provides pack_h80_hybrid_dots, compute_sub_vent_conduit_map

RAW_DIR = ROOT / "data" / "raw"
EXT_DIR = RAW_DIR / "external"
SCORED_DIR = RAW_DIR / "scored"
CACHE_DIR = ROOT / "data" / "processed"
EVID_DIR = ROOT / "evidence"
SUB_DIR = ROOT / "submissions"
DL_DIR = ROOT / "docs" / "downloads"
CONDUIT_COLS = [18, 19, 20, 22, 23, 48]


def conduit_map(stack, fp) -> np.ndarray:
    X = stack.X_mmap
    ranks = np.zeros(stack.valid_idx.size, dtype=np.float32)
    for c in CONDUIT_COLS:
        ranks += _rank01(np.asarray(X[:, c], dtype=np.float32))
    out = np.zeros(fp.size, dtype=np.float32)
    out[stack.valid_idx] = ranks / len(CONDUIT_COLS)
    return out.reshape(fp.shape)


def conduit_canary(stack, folds, fp, conduit) -> dict:
    rng = np.random.default_rng(81)
    pos, neg = [], []
    for f in folds:
        reg = f["region"]
        truth = f["truth"] & reg
        neg_pool = reg & (~f["visible"]) & (~f["hidden_all"])
        p_idx = np.flatnonzero(truth)
        n_pool = np.flatnonzero(neg_pool)
        n_idx = rng.choice(n_pool, size=min(len(n_pool), 4 * len(p_idx)), replace=False)
        pos.append(conduit.ravel()[p_idx]); neg.append(conduit.ravel()[n_idx])
    pos = np.concatenate(pos); neg = np.concatenate(neg)
    y = np.concatenate([np.ones(len(pos)), np.zeros(len(neg))])
    s = np.concatenate([pos, neg])
    auc = float(roc_auc_score(y, s))
    return {"feature": "h81_conduit_composite(cols18,19,20,22,23,48)",
            "symmetric_auc": round(max(auc, 1.0 - auc), 5), "raw_auc": round(auc, 5),
            "n_pos": int(len(pos)), "pass_le_0.90": bool(max(auc, 1 - auc) <= 0.90)}


def main() -> None:
    t0 = time.time()
    g = GRID.load_grid(ROOT / "data" / "bridge")
    fp, cat = g.footprint, g.catalogue
    print(f"[h81] grid {g.shape} footprint={int(fp.sum()):,} catalogue={int(cat.sum()):,}", flush=True)

    stack = build_or_load_static_stack(RAW_DIR, EXT_DIR, CACHE_DIR, fp, force_rebuild=False)
    conduit = conduit_map(stack, fp)
    print(f"[h81] conduit map built ({time.time()-t0:.0f}s)", flush=True)

    folds20, quad20 = HO.buffered_component_draw(g, seed=20, hide_frac=0.20, buffer_px=3)
    folds30, quad30 = HO.buffered_component_draw(g, seed=30, hide_frac=0.20, buffer_px=3)
    canary = conduit_canary(stack, folds20, fp, conduit)
    print(f"[h81] conduit canary: {canary}", flush=True)
    assert canary["pass_le_0.90"], "conduit composite leakage canary failed"

    with rasterio.open(SCORED_DIR / "gems24-h25-1-dotted-h19-5-d2-8-20261002-e56ea318af89-nan.tif") as s:
        d28_bool = np.isfinite(s.read(1)) & (s.read(1) > 0) & fp
    with rasterio.open(SCORED_DIR / "GEMSDOE32__gemsdoe32-h33-h33-2-b2-20261004T220000Z-e5eb6e7e-zeros.tif") as s:
        h33_bool = np.isfinite(s.read(1)) & (s.read(1) > 0) & fp

    oofs, foldsets = {}, {}
    for seed, folds, quad in ((20, folds20, quad20), (30, folds30, quad30)):
        cache = CACHE_DIR / f"h80_oof_seed{seed}.npz"
        if cache.exists():
            z = np.load(cache); oofs[seed] = {k: z[k].astype(np.float32) for k in ("oof_e1_phys", "oof_e2_selfmask", "oof_e3_hero")}
            print(f"[h81] OOF seed{seed} from cache", flush=True)
        else:
            oofs[seed] = run_holdout_oof_for_seed(stack, folds, quad, fp, seed=seed)
            np.savez_compressed(cache, **{k: v.astype(np.float16) for k, v in oofs[seed].items()})
        foldsets[seed] = folds

    sub_vent = RE.compute_sub_vent_conduit_map(stack, fp)
    arms_terms_all, summaries = {}, {}
    for seed in (20, 30):
        oof = oofs[seed]; folds = foldsets[seed]
        gated = (oof["oof_e2_selfmask"] * (0.60 + 0.40 * conduit)).astype(np.float32)
        total_active = float(sum(int((f["region"] & fp & (~f["masked_known"])).sum()) for f in folds))
        terms = {}
        def accum(name, pred, f):
            _, t = EH.evaluate(pred, f, valid=fp, block_side=200)
            terms[name] = t if name not in terms else terms[name] + t
        for f in folds:
            reg = f["region"] & fp; vis = f["visible"] & fp
            act = reg & (~f["masked_known"])
            d_vis = ndi.distance_transform_edt(~vis).astype(np.float32)
            q_share = float(act.sum()) / total_active
            q37 = max(1, int(round(37654 * q_share))); q41 = max(1, int(round(41600 * q_share)))
            d_e2, _ = submodular_pack_300m(oof["oof_e2_selfmask"], act & (d_vis > 2.40), target_n=q37, support_quantile=0.86, gamma=2.1, skeleton_prefilter=True)
            accum("CTRL_E80_2_puremeta_37654", d_e2.astype(np.float32), f)
            d_h81, _ = submodular_pack_300m(gated, act & (d_vis > 2.40), target_n=q37, support_quantile=0.86, gamma=2.1, skeleton_prefilter=True)
            accum("H81_1_conduit_gated_37654", d_h81.astype(np.float32), f)
            d_hero, _ = RE.pack_h80_hybrid_dots(oof["oof_e2_selfmask"], oof["oof_e1_phys"], sub_vent, stack.updip_energy_full, stack.slope10_full, vis, d_vis, act, fp, target_n=q41)
            accum("CTRL_E80_3_hero_hybrid_41600", d_hero.astype(np.float32), f)
            accum("GEMSDOE32_H33_fold_fair", (d28_bool & (d_vis > 2.0) & act).astype(np.float32), f)
        summaries[seed] = EH.pooled_summary(terms, draws=1000, seed=520810 + seed, candidate="H81_1_conduit_gated_37654")
        arms_terms_all[seed] = terms
        for k, v in summaries[seed]["scores"].items():
            print(f"[h81] seed{seed} {k}: HOLDOUT-DTI={v['dti']:.5f} CI={v['ci95']}", flush=True)
        del gated
        gc.collect()

    pooled_terms = {k: arms_terms_all[20][k] + arms_terms_all[30][k] for k in arms_terms_all[20]}
    pooled = EH.pooled_summary(pooled_terms, draws=1000, seed=520899, candidate="H81_1_conduit_gated_37654")
    for k, v in pooled["scores"].items():
        print(f"[h81] POOLED {k}: HOLDOUT-DTI={v['dti']:.5f} CI={v['ci95']} K={v['withheld_positive_pixels']}", flush=True)
    print(f"[h81] paired vs H81-1: {json.dumps(pooled['paired_differences'], indent=1)}", flush=True)

    h81_dti = pooled["scores"]["H81_1_conduit_gated_37654"]["dti"]
    e2_dti = pooled["scores"]["CTRL_E80_2_puremeta_37654"]["dti"]
    hero_dti = pooled["scores"]["CTRL_E80_3_hero_hybrid_41600"]["dti"]

    # Full-catalogue surfaces (cached by run_experiments convention)
    full_cache = CACHE_DIR / "h80_full_surfaces.npz"
    if full_cache.exists():
        zf = np.load(full_cache); p_meta_full = zf["p_meta_full"].astype(np.float32); p_phys_full = zf["p_phys_full"].astype(np.float32)
    else:
        p_meta_full, p_phys_full = build_full_catalogue_surfaces(stack, cat, fp, seed=8080)
        np.savez_compressed(full_cache, p_meta_full=p_meta_full.astype(np.float16), p_phys_full=p_phys_full.astype(np.float16))

    d_cat = ndi.distance_transform_edt(~cat).astype(np.float32)
    if h81_dti > e2_dti and h81_dti > hero_dti:
        lane, surf_full = "H81-1 conduit-gated meta", (p_meta_full * (0.60 + 0.40 * conduit)).astype(np.float32)
        stem = "gems80-h81-conduit-gated-meta-37654-20261022-zeros"
        name = "gems80-h81-conduit-gated-meta-37654-zeros"
        note = "H81-1 clay-cap/demag conduit-gated self-masked meta, skeleton-submod 300m, collar2.4, 37654 dots [0,1]"
        verdict = "promote"
    elif e2_dti > hero_dti:
        lane, surf_full = "H81-2 pure self-masked meta transfer (no hybrid dilution)", p_meta_full
        stem = "gems80-h81-puremeta-transfer-37654-20261022-zeros"
        name = "gems80-h81-puremeta-transfer-37654-zeros"
        note = "H81-2 pure self-masked meta transfer, skeleton-submod 300m, collar2.4, 37654 dots, no vent dilution [0,1]"
        verdict = "promote"
    else:
        lane, surf_full, verdict = "none", None, "negative"
    print(f"[h81] decision lane={lane} verdict={verdict} note_len={len(note) if verdict=='promote' else 0}", flush=True)

    result = {"schema": "gems80.h81-experiment.v1", "evaluator_version": EH.VERSION,
              "evidence_class": "HOLDOUT-DTI", "conduit_canary": canary,
              "seed_20": summaries[20], "seed_30": summaries[30], "pooled": pooled,
              "decision_lane": lane, "verdict": verdict}

    if verdict == "promote":
        assert len(note) <= 140
        dots, picked = submodular_pack_300m(surf_full, fp & (~cat) & (d_cat > 2.40), target_n=37654, support_quantile=0.86, gamma=2.1, skeleton_prefilter=True)
        print(f"[h81] final dots={picked} min_dcat={float(d_cat[dots].min()):.3f}", flush=True)
        # overlap vs shipped hero (our own previous file, reported separately from registry)
        with rasterio.open(SUB_DIR / "gems80-h80-dpvc-updip-vent-submod-41600-20261021-zeros.tif") as s:
            hero_dots_prev = (s.read(1) > 0) & fp
        from gems57 import uniqueness as UQ
        halo_hero = UQ._near(hero_dots_prev, radius_px=3.0)
        hero_overlap = float(halo_hero[dots].mean()) if dots.any() else 0.0
        uniq = run_full_uniqueness_audit(surf_full, dots, SCORED_DIR, fp, cat, stem)
        # Pass-2 fix (2026-10-22): uniqueness is a promotion GATE, checked
        # before writing. The first run wrote H81-2 (overlap 0.8593 > 0.70)
        # and it had to be withdrawn as a duplicate.
        if uniq["post_placement_binding_like_for_like_gate"]["duplicate"] or not uniq["pre_placement_surface_gate"]["passed"]:
            result.update({"verdict": "negative-duplicate", "submission_stem": stem,
                           "uniqueness": uniq, "dots": int(dots.sum())})
            (EVID_DIR / "h81_experiment.json").write_text(json.dumps(result, indent=2, default=str) + "\n")
            print("[h81] uniqueness gate FAILED — duplicate, not writing submission", flush=True)
            return
        profile = compute_spatial_profile(dots, stack.slope10_full, stack.lidar_valid_full, cat, fp)
        pred = np.zeros(fp.shape, dtype=np.float32); pred[dots] = 1.0
        sub_tif = SUB_DIR / f"{stem}.tif"
        receipt = SW.write_submission(sub_tif, pred, ROOT / "data" / "bridge" / "sample_submission.tif", fp,
                                      note=note, name=name, catalogue=cat,
                                      metadata={"lane": lane, "holdout_dti_pooled": pooled["scores"]})
        import shutil
        for ext in (".tif", ".zip", ".json"):
            shutil.copy2(SUB_DIR / f"{stem}{ext}", DL_DIR / f"{stem}{ext}")
        sha = hashlib.sha256((SUB_DIR / f"{stem}.tif").read_bytes()).hexdigest()
        result.update({"submission_stem": stem, "submission_name": name, "submission_note": note,
                       "submission_note_chars": len(note), "raster_sha256": sha,
                       "raster_bytes": (SUB_DIR / f"{stem}.tif").stat().st_size,
                       "dots": int(dots.sum()), "overlap_3px_vs_shipped_hero": round(hero_overlap, 5),
                       "uniqueness": uniq, "spatial_profile": profile, "validator_receipt": receipt})
        print(f"[h81] wrote {stem} sha256={sha} validator_pass={receipt['validator']['all_checks_passed']}", flush=True)

    (EVID_DIR / "h81_experiment.json").write_text(json.dumps(result, indent=2, default=str) + "\n")
    print(f"[h81] done in {time.time()-t0:.0f}s -> evidence/h81_experiment.json", flush=True)


if __name__ == "__main__":
    main()
