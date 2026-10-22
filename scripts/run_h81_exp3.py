#!/usr/bin/env python3
"""H81 Experiment 2 (Pass 2 remedy): core-weighted hybrid allocator frontier.

Experiment 1 result: H81-1 conduit gate negative (pooled 0.18403 < pure-meta
0.18897); H81-2 pure-meta transfer beat hero holdout (0.18897 > 0.16128) but
FAILED the binding uniqueness gate (0.8593 > 0.70 vs GEMSDOE77 H85-EDTI) and
is logged DUPLICATE — withdrawn, do not submit.

This experiment tests ONE allocator change on cached OOF surfaces:
shares (core, collar, blind) = (0.78, 0.07, 0.15) at N=41,600 in
pack_h80_hybrid_dots, vs the shipped hero shares (0.635, 0.085, 0.280).
Promotion requires BOTH: pooled HOLDOUT-DTI > hero control on identical
folds AND final-dot uniqueness gates pass (pre |rho|<=0.90, discriminating
3-px overlap <=0.70). Otherwise verdict=negative, shipped hero stays primary.
"""
from __future__ import annotations
import gc, hashlib, importlib.util, json, shutil, sys, time
from pathlib import Path
import numpy as np, rasterio
from scipy import ndimage as ndi

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "vendor")); sys.path.insert(0, str(ROOT / "scripts"))
from gems57 import evaluate_holdout as EH, grid as GRID, holdout as HO, submission_writer as SW  # noqa: E402
from gems80.features import build_or_load_static_stack  # noqa: E402
from gems80.forensics import compute_spatial_profile, run_full_uniqueness_audit  # noqa: E402
from gems80.models import build_full_catalogue_surfaces  # noqa: E402
spec = importlib.util.spec_from_file_location("run_experiments", ROOT / "scripts" / "run_experiments.py")
RE = importlib.util.module_from_spec(spec); spec.loader.exec_module(RE)

RAW_DIR = ROOT / "data" / "raw"; EXT_DIR = RAW_DIR / "external"; SCORED_DIR = RAW_DIR / "scored"
CACHE_DIR = ROOT / "data" / "processed"; EVID_DIR = ROOT / "evidence"; SUB_DIR = ROOT / "submissions"; DL_DIR = ROOT / "docs" / "downloads"
SHARES_NEW = (0.67, 0.08, 0.25)

def main():
    t0 = time.time()
    g = GRID.load_grid(ROOT / "data" / "bridge"); fp, cat = g.footprint, g.catalogue
    stack = build_or_load_static_stack(RAW_DIR, EXT_DIR, CACHE_DIR, fp, force_rebuild=False)
    sub_vent = RE.compute_sub_vent_conduit_map(stack, fp)
    pooled_terms_by_seed = []
    summaries = {}
    for seed in (20, 30):
        folds, _ = HO.buffered_component_draw(g, seed=seed, hide_frac=0.20, buffer_px=3)
        z = np.load(CACHE_DIR / f"h80_oof_seed{seed}.npz")
        oof = {k: z[k].astype(np.float32) for k in ("oof_e1_phys", "oof_e2_selfmask")}
        total_active = float(sum(int((f["region"] & fp & (~f["masked_known"])).sum()) for f in folds))
        terms = {}
        def accum(name, pred, f):
            _, t = EH.evaluate(pred, f, valid=fp, block_side=200)
            terms[name] = t if name not in terms else terms[name] + t
        for f in folds:
            reg = f["region"] & fp; vis = f["visible"] & fp; act = reg & (~f["masked_known"])
            d_vis = ndi.distance_transform_edt(~vis).astype(np.float32)
            q41 = max(1, int(round(41600 * float(act.sum()) / total_active)))
            d_hero, _ = RE.pack_h80_hybrid_dots(oof["oof_e2_selfmask"], oof["oof_e1_phys"], sub_vent, stack.updip_energy_full, stack.slope10_full, vis, d_vis, act, fp, target_n=q41)
            accum("CTRL_hero_41600", d_hero.astype(np.float32), f)
            d_new, _ = RE.pack_h80_hybrid_dots(oof["oof_e2_selfmask"], oof["oof_e1_phys"], sub_vent, stack.updip_energy_full, stack.slope10_full, vis, d_vis, act, fp, target_n=q41, shares=SHARES_NEW)
            accum("EXP3_core67_41600", d_new.astype(np.float32), f)
        summaries[seed] = EH.pooled_summary(terms, draws=1000, seed=520810 + seed, candidate="EXP3_core67_41600")
        pooled_terms_by_seed.append(terms)
        for k, v in summaries[seed]["scores"].items():
            print(f"[exp3] seed{seed} {k}: HOLDOUT-DTI={v['dti']:.5f} CI={v['ci95']}", flush=True)
        del oof; gc.collect()
    pooled_terms = {k: pooled_terms_by_seed[0][k] + pooled_terms_by_seed[1][k] for k in pooled_terms_by_seed[0]}
    pooled = EH.pooled_summary(pooled_terms, draws=1000, seed=520899, candidate="EXP3_core67_41600")
    for k, v in pooled["scores"].items():
        print(f"[exp3] POOLED {k}: HOLDOUT-DTI={v['dti']:.5f} CI={v['ci95']} K={v['withheld_positive_pixels']}", flush=True)
    print(f"[exp3] paired: {json.dumps(pooled['paired_differences'])}", flush=True)

    result = {"schema": "gems80.h81-exp3.v1", "evidence_class": "HOLDOUT-DTI", "evaluator_version": EH.VERSION,
              "shares_new": SHARES_NEW, "seed_20": summaries[20], "seed_30": summaries[30], "pooled": pooled}
    new_dti = pooled["scores"]["EXP3_core67_41600"]["dti"]; hero_dti = pooled["scores"]["CTRL_hero_41600"]["dti"]

    if new_dti > hero_dti:
        zf = np.load(CACHE_DIR / "h80_full_surfaces.npz")
        p_meta_full = zf["p_meta_full"].astype(np.float32); p_phys_full = zf["p_phys_full"].astype(np.float32)
        d_cat = ndi.distance_transform_edt(~cat).astype(np.float32)
        surf_full = np.where(fp & (~cat), 0.55 * p_meta_full + 0.45 * sub_vent, 0.0).astype(np.float32)
        dots, picked = RE.pack_h80_hybrid_dots(p_meta_full, p_phys_full, sub_vent, stack.updip_energy_full, stack.slope10_full, cat, d_cat, fp & (~cat), fp, target_n=41600, shares=SHARES_NEW)
        print(f"[exp3] final dots={picked} min_dcat={float(d_cat[dots].min()):.3f}", flush=True)
        stem = "gems80-h81-core67-hybrid-41600-20261022-zeros"
        uniq = run_full_uniqueness_audit(surf_full, dots, SCORED_DIR, fp, cat, stem)
        print(f"[exp3] uniqueness: pre_rho={uniq['pre_placement_surface_gate']['max_abs_spearman_rho']} worst_overlap={uniq['post_placement_binding_like_for_like_gate']['worst_forward_overlap']} duplicate={uniq['post_placement_binding_like_for_like_gate']['duplicate']}", flush=True)
        result["uniqueness"] = uniq
        gates_pass = bool(uniq["pre_placement_surface_gate"]["passed"]) and not bool(uniq["post_placement_binding_like_for_like_gate"]["duplicate"])
        if gates_pass:
            name = "gems80-h81-core67-hybrid-41600-zeros"
            note = "H81-E3 core-weighted hybrid 67/8/25 skeleton-submod 300m, 41600 dots, updip+vent+self-masked meta [0,1]"
            assert len(note) <= 140, len(note)
            profile = compute_spatial_profile(dots, stack.slope10_full, stack.lidar_valid_full, cat, fp)
            pred = np.zeros(fp.shape, dtype=np.float32); pred[dots] = 1.0
            receipt = SW.write_submission(SUB_DIR / f"{stem}.tif", pred, ROOT / "data" / "bridge" / "sample_submission.tif", fp, note=note, name=name, catalogue=cat, metadata={"lane": "H81-E3 core67 hybrid", "holdout_dti_pooled": pooled["scores"]})
            for ext in (".tif", ".zip", ".json"):
                shutil.copy2(SUB_DIR / f"{stem}{ext}", DL_DIR / f"{stem}{ext}")
            result.update({"verdict": "promote", "submission_stem": stem, "submission_name": name, "submission_note": note,
                           "submission_note_chars": len(note), "dots": int(dots.sum()),
                           "raster_sha256": hashlib.sha256((SUB_DIR / f"{stem}.tif").read_bytes()).hexdigest(),
                           "raster_bytes": (SUB_DIR / f"{stem}.tif").stat().st_size,
                           "spatial_profile": profile, "validator_receipt": receipt})
        else:
            result["verdict"] = "negative-duplicate"
    else:
        result["verdict"] = "negative-holdout"
    (EVID_DIR / "h81_exp3.json").write_text(json.dumps(result, indent=2, default=str) + "\n")
    print(f"[exp3] verdict={result['verdict']} done {time.time()-t0:.0f}s", flush=True)

if __name__ == "__main__":
    main()
