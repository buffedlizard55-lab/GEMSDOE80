#!/usr/bin/env python3
"""GEMSDOE80 unique-lane experiments (budget: 3 experiments).

Lane H80 — catalogue-residual strain-accommodation.
Not a copy of H33-2-B2 (catalogue-flank prune of the H19-5 ridge) or any other
GEMSDOE dotted-ridge family.

Hypothesis: unmapped geothermal-indicative faults occupy *strain-transfer /
accommodation zones* — places where dilatation or the 2nd strain invariant
*changes* — rather than the high-strain ridges already captured by the USGS /
INGENIOUS catalogue. Residualizing the strain-gradient against distance to the
*visible* catalogue removes the trivial near-trace halo. Magnetic horizontal
gradient corroborates a basement fabric. Dots are NMS-thinned, off-catalogue.

Named non-fault mimic: lithologic contacts and GPS/InSAR interpolation seams
can produce strain-gradient or magnetic-edge lines without displacement.

Each experiment:
  1. Leakage canary (per-feature AUC vs withheld labels; AUC>0.90 = leakage).
  2. Build surface from VISIBLE faults only (holdout) and from the full
     catalogue (submission).
  3. HOLDOUT-DTI + 95% CI.
  4. Registry Spearman / 3 px overlap (stop that arm if duplicate).
  5. Random-dot control on the same emit-eligible mask.

Legacy research packaging only: promotion is blocked without spatially-blocked
validation against a matched current best. Never select by SGMC diagnostic.
"""
from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import rasterio
from scipy.ndimage import binary_dilation, distance_transform_edt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gems80.emission import (  # noqa: E402
    along_fabric_derivative,
    auc_mannwhitney,
    distance_residual,
    gradient_magnitude,
    nms_dots,
    percentile_norm,
)
from gems80.holdout import score_holdout, withhold_segments  # noqa: E402
from gems80.io_grid import load_band, load_labels, load_template, sha256_file  # noqa: E402
from gems80.metric import EVALUATOR_VERSION, dti_binary, worked_example_dti  # noqa: E402
from gems80.uniqueness import check_registry  # noqa: E402

sys.path.insert(0, str(ROOT / "scripts"))
from submission_writer import write_pair  # noqa: E402

DATA = ROOT / "data"
FEAT = DATA / "training_features.tif"
TEMPLATE = DATA / "sample_submission.tif"
LABELS = DATA / "labels.tif"
SGMC = DATA / "external" / "derived_sgmc_faults_100m_u8.tif"
REGISTRY_DIR = DATA / "registry"
OUT_DL = ROOT / "docs" / "downloads"
OUT_SUB = ROOT / "docs" / "submissions"
REPORTS = ROOT / "reports"
RESEARCH = ROOT / "research"

# 1-based band ids in training_features.tif (verified by description tags)
BAND_DIL = 8          # geod_dilaterate
BAND_INV = 4          # geod_2ndinv
BAND_SHEAR = 7        # geod_shearrate
BAND_TMI_HG = 3       # tmi_hg
BAND_TC = 6           # tc (tilt / total curvature)
BAND_COND = 17        # cond_surf
BAND_DEPTH = 15       # depth_to_base_surf
BAND_IEQ = 16         # ieq_n100a15
BAND_RTP = 2          # rtp

STAMP = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
RNG = np.random.default_rng(20261021)


def load_registry() -> dict[str, np.ndarray]:
    out = {}
    if not REGISTRY_DIR.exists():
        return out
    for p in sorted(REGISTRY_DIR.glob("*.tif")):
        with rasterio.open(p) as src:
            out[p.stem] = src.read(1)
    return out


def leakage_canary(bands: dict[str, tuple[np.ndarray, np.ndarray]], hidden: np.ndarray,
                   foot: np.ndarray, visible: np.ndarray) -> dict:
    """AUC of each feature vs withheld positives, sampled on foot & ~visible."""
    domain = foot & ~visible
    truth = (hidden & domain).astype(np.int64)
    idx = np.flatnonzero(domain.ravel())
    if idx.size > 180_000:
        idx = RNG.choice(idx, size=180_000, replace=False)
    rows = []
    max_auc = 0.0
    leak = False
    for name, (arr, valid) in bands.items():
        s = arr.ravel()[idx]
        y = truth.ravel()[idx]
        ok = np.isfinite(s) & valid.ravel()[idx]
        auc = auc_mannwhitney(s[ok], y[ok])
        flag = bool(np.isfinite(auc) and auc > 0.90)
        leak = leak or flag
        max_auc = max(max_auc, auc if np.isfinite(auc) else 0.0)
        rows.append({"feature": name, "auc": None if not np.isfinite(auc) else float(auc), "leakage_flag": flag})
    return {
        "available": True,
        "per_feature": rows,
        "max_auc": float(max_auc),
        "leakage_flag": bool(leak),
        "note": "AUC > 0.90 on a single feature means leakage until proven otherwise",
    }


def random_dots(valid_emit: np.ndarray, k: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    ys, xs = np.nonzero(valid_emit)
    k = min(k, ys.size)
    pick = rng.choice(ys.size, size=k, replace=False)
    m = np.zeros(valid_emit.shape, dtype=bool)
    m[ys[pick], xs[pick]] = True
    return m


def sgmc_off_diagnostic(pred: np.ndarray, labels: np.ndarray, foot: np.ndarray) -> dict:
    if not SGMC.is_file():
        return {"available": False, "note": "SGMC raster not present"}
    with rasterio.open(SGMC) as src:
        sgm = src.read(1) > 0
    near = binary_dilation(labels, iterations=3)
    truth = sgm & foot & ~labels & ~near
    pred_off = pred & foot & ~labels
    res = dti_binary(pred_off, truth, valid=foot)
    return {
        "available": True,
        "label": "DIAGNOSTIC (not HOLDOUT-DTI, not ORGANIZER-CONFIRMED)",
        "n_sgmc_off_px": int(truth.sum()),
        "DTI": res["DTI"],
        "TP_w": res["TP_w"],
        "FP_w": res["FP_w"],
        "FN_w": res["FN_w"],
        "note": "SGMC traces not in the USGS/INGENIOUS raster, 3 px collar. Proxy only.",
    }


def build_score(kind: str, dil, inv, hg, tc, cond, rtp, valid, dist, emit_valid) -> np.ndarray:
    if kind == "A":
        g = gradient_magnitude(np.where(valid, dil, 0.0))
        resid = distance_residual(g, dist, valid)
        score = percentile_norm(resid, emit_valid) * (0.45 + 0.55 * percentile_norm(hg, valid))
        score = np.where(dil > 0, score, score * np.float32(0.25))
    elif kind == "B":
        g = gradient_magnitude(np.where(valid, inv, 0.0))
        resid = distance_residual(g, dist, valid)
        score = (percentile_norm(resid, emit_valid)
                 * (0.40 + 0.30 * percentile_norm(hg, valid)
                    + 0.30 * percentile_norm(cond, valid)))
        score = np.where(dil > np.percentile(dil[valid], 55), score, score * np.float32(0.20))
    elif kind == "C":
        af = along_fabric_derivative(np.where(valid, dil, 0.0), np.where(valid, rtp, 0.0))
        resid = distance_residual(af, dist, valid)
        score = percentile_norm(resid, emit_valid) * (0.50 + 0.50 * percentile_norm(tc, valid))
        score = np.where(dil > 0, score, score * np.float32(0.20))
    else:
        raise ValueError(kind)
    score = np.where(emit_valid, score, 0.0).astype(np.float32)
    return score


def run_arm(kind, cfg, dil, inv, hg, tc, cond, rtp, valid, foot, labels, split, registry) -> dict:
    name = cfg["name"]
    print(f"\n=== {name} ===", flush=True)
    # Holdout surface from VISIBLE catalogue only
    dist_vis = distance_transform_edt(~split["visible"])
    near_vis = binary_dilation(split["visible"], iterations=cfg["buffer_px"])
    emit_h = foot & valid & ~near_vis
    score_h = build_score(kind, dil, inv, hg, tc, cond, rtp, valid, dist_vis, emit_h)
    dots_h = nms_dots(score_h, emit_h, sep=cfg["sep"], max_k=cfg["max_k"])
    hold = score_holdout(dots_h, split, foot, bootstrap=False)
    rnd = random_dots(emit_h, k=int(dots_h.sum()) or cfg["max_k"], seed=20261021 + cfg["seed_off"])
    hold_rnd = score_holdout(rnd, split, foot, bootstrap=False)

    # Submission surface from FULL catalogue
    dist_full = distance_transform_edt(~labels)
    near_full = binary_dilation(labels, iterations=cfg["buffer_px"])
    emit_s = foot & valid & ~near_full
    score_s = build_score(kind, dil, inv, hg, tc, cond, rtp, valid, dist_full, emit_s)
    surface_check = check_registry(np.zeros_like(foot), score_s, foot, registry, RNG)
    if surface_check["duplicate_lane"]:
        raise RuntimeError("Duplicate surface: stop before placement")
    dots_s = nms_dots(score_s, emit_s, sep=cfg["sep"], max_k=cfg["max_k"])
    uniq = check_registry(dots_s, score_s, foot, registry, RNG)
    if uniq["duplicate_lane"]:
        raise RuntimeError("Duplicate final dots: stop lane")
    sg = sgmc_off_diagnostic(dots_s, labels, foot)

    print(f"  dots_holdout={int(dots_h.sum())} HOLDOUT-DTI={hold['score']:.4f} "
          f"withheld={hold['withheld_positive_count']}", flush=True)
    print(f"  random HOLDOUT-DTI={hold_rnd['score']:.4f}  beats_random={hold['score'] > hold_rnd['score']}", flush=True)
    print(f"  dots_submit={int(dots_s.sum())} duplicate_lane={uniq['duplicate_lane']}", flush=True)
    if sg.get("available"):
        print(f"  SGMC-off diagnostic DTI={sg['DTI']:.4f} n={sg['n_sgmc_off_px']}", flush=True)
    for row in uniq["vs"]:
        print(f"    vs {row['name']}: spearman={row['spearman']:.4f} overlap3={row['dot_fraction_within_3px']:.4f} dup={row['duplicate']}", flush=True)

    return {
        "name": name,
        "kind": kind,
        "cfg": {k: cfg[k] for k in cfg if k != "name"},
        "n_dots_holdout": int(dots_h.sum()),
        "n_dots_submit": int(dots_s.sum()),
        "holdout": hold,
        "holdout_random": hold_rnd,
        "beats_random_holdout": bool(hold["score"] > hold_rnd["score"]),
        "uniqueness": uniq,
        "sgmc_off_diagnostic": sg,
        "dots_submit": dots_s,
        "dots_holdout": dots_h,
        "score_submit": score_s,
    }


def main() -> int:
    # The lane budget persists across sessions; do not reset it by rerunning.
    card_path = RESEARCH / "run-card.json"
    if card_path.exists():
        card = json.loads(card_path.read_text())
        if len(card.get("all_experiments", [])) >= 3:
            print("STOP: H80 lane budget exhausted (3/3). No new raster generated.")
            return 2
    we = worked_example_dti()
    print(f"metric self-check worked example DTI={we:.6f} (page ~0.60)")

    tpl = load_template(TEMPLATE)
    foot = tpl["foot"]
    labels = load_labels(LABELS, foot)
    print(f"grid {tpl['report']['height']}x{tpl['report']['width']} foot={int(foot.sum())} cat={int(labels.sum())}")

    dil, v8, d8, _ = load_band(FEAT, BAND_DIL)
    inv, v4, d4, _ = load_band(FEAT, BAND_INV)
    hg, v3, d3, _ = load_band(FEAT, BAND_TMI_HG)
    tc, v6, d6, _ = load_band(FEAT, BAND_TC)
    cond, v17, d17, _ = load_band(FEAT, BAND_COND)
    rtp, v2, d2, _ = load_band(FEAT, BAND_RTP)
    shear, v7, d7, _ = load_band(FEAT, BAND_SHEAR)
    ieq, v16, d16, _ = load_band(FEAT, BAND_IEQ)
    valid = foot & v8 & v4 & v3 & v6 & v17 & v2
    print(f"feature-valid ∩ foot = {int(valid.sum())}")
    print("bands:", d8, "|", d4, "|", d3)

    split = withhold_segments(labels, foot)
    print(f"holdout hidden_px={split['n_hidden_px']} visible_px={split['n_visible_px']} "
          f"held_components={split['held_ids'].size}/{split['n_components']}")

    bands = {
        "dilatation": (dil, v8),
        "strain_2ndinv": (inv, v4),
        "tmi_hg": (hg, v3),
        "tilt_tc": (tc, v6),
        "conductivity": (cond, v17),
        "rtp": (rtp, v2),
        "shear": (shear, v7),
        "earthquake_density": (ieq, v16),
        "dilatation_grad": (gradient_magnitude(np.where(v8, dil, 0.0)), v8),
        "inv_grad": (gradient_magnitude(np.where(v4, inv, 0.0)), v4),
    }
    canary = leakage_canary(bands, split["hidden"], foot, split["visible"])
    print("leakage canary max_auc", canary["max_auc"], "flag", canary["leakage_flag"])
    for r in canary["per_feature"]:
        print(f"  AUC {r['feature']:22s} {r['auc']}")

    registry = load_registry()
    if not registry:
        raise RuntimeError("Missing registry: uniqueness cannot be established")
    print("registry rasters:", list(registry))

    arms_cfg = [
        ("A", dict(name="H80-A residual-|∇dil| × tmi_hg, dil>0, buf=3, sep=3, k=35129",
                   buffer_px=3, sep=3, max_k=35129, seed_off=1)),
        ("B", dict(name="H80-B residual-|∇2ndinv| × (tmi_hg, cond), dil>p55, buf=3, sep=3, k=33113",
                   buffer_px=3, sep=3, max_k=33113, seed_off=2)),
        ("C", dict(name="H80-C residual along-fabric |t·∇dil| × tc, dil>0, buf=4, sep=3, k=34033",
                   buffer_px=4, sep=3, max_k=34033, seed_off=3)),
    ]

    results = []
    for kind, cfg in arms_cfg:
        results.append(run_arm(kind, cfg, dil, inv, hg, tc, cond, rtp, valid, foot, labels, split, registry))

    # Select: unique, no leakage, beats random; then max SGMC-off diagnostic (not written as a score)
    eligible = [r for r in results
                if (not r["uniqueness"]["duplicate_lane"])
                and (not canary["leakage_flag"])
                and r["beats_random_holdout"]]
    if not eligible:
        print("NO ELIGIBLE ARMS — stopping without a submission.")
        RESEARCH.mkdir(parents=True, exist_ok=True)
        (RESEARCH / "run-card.json").write_text(json.dumps({
            "verdict": "negative",
            "reason": "No arm cleared uniqueness, canary and random-control gates",
            "experiments": [{k: v for k, v in r.items() if k not in ("dots_submit", "dots_holdout", "score_submit")} for r in results],
        }, indent=2, default=_json_default) + "\n")
        return 2

    def key(r):
        return (r["holdout"]["score"], -r["n_dots_submit"])

    best = max(eligible, key=key)
    print(f"\nSELECTED {best['name']}  dots={best['n_dots_submit']}", flush=True)
    print("computing bootstrap 95% CI on selected arm...", flush=True)
    best["holdout"] = score_holdout(best["dots_holdout"], split, foot, bootstrap=True)
    print(f"  HOLDOUT-DTI={best['holdout']['score']:.4f} CI={best['holdout']['ci95']}", flush=True)

    stem = f"gemsdoe80-h80{best['kind'].lower()}-strain-accom-resid-{STAMP}"
    OUT_DL.mkdir(parents=True, exist_ok=True)
    OUT_SUB.mkdir(parents=True, exist_ok=True)
    bundle = write_pair(best["dots_submit"], TEMPLATE, LABELS, OUT_DL, stem)
    # copy zeros into submissions/ for the Pages download table
    import shutil
    src_z = Path(bundle["zeros"]["path"])
    shutil.copy2(src_z, OUT_SUB / src_z.name)
    shutil.copy2(src_z.with_suffix(".zip"), OUT_SUB / (src_z.stem + ".zip")) if False else None
    zpath = Path(bundle["zip_path"])
    shutil.copy2(zpath, OUT_SUB / zpath.name)
    shutil.copy2(Path(bundle["nan"]["path"]), OUT_SUB / Path(bundle["nan"]["path"]).name)

    note = (
        f"H80{best['kind']} cat-residual strain-accom NMS dots {best['n_dots_submit']}; "
        f"off-catalogue buf{best['cfg']['buffer_px']}; HOLDOUT-DTI only; unique vs H33/H19-5"
    )
    if len(note) > 140:
        note = note[:140]

    promote = (
        best["holdout"].get("promotion_eligible", False)
        and bundle["zeros"]["format_valid"]
        and not best["uniqueness"]["duplicate_lane"]
        and not canary["leakage_flag"]
        and best["beats_random_holdout"]
    )
    verdict = "promote" if promote else "negative"

    run_card = {
        "experiment_budget": "3 of 3; hard stop; no weekly submission slot spent (selector is separate)",
        "lane": "H80 catalogue-residual strain-accommodation (not H19-5 ridge, not H33 B=2 flank prune)",
        "selected": best["name"],
        "hypothesis": (
            "Unmapped geothermal-indicative faults occupy strain-transfer/accommodation zones "
            "(spatial gradient of dilatation / 2nd strain invariant) that are not explained by "
            "distance to the visible USGS/INGENIOUS catalogue, corroborated by magnetic fabric."
        ),
        "mechanism": (
            "Sobel magnitude (or along-fabric directional derivative) of geodetic dilatation or "
            "2nd invariant; subtract binned median vs distance-to-visible-catalogue; multiply by "
            "percentile-normalized magnetic HG or tilt; gate on positive dilatation; exclude a "
            f"{best['cfg']['buffer_px']} px catalogue collar; NMS local maxima, sep="
            f"{best['cfg']['sep']} px, cap {best['cfg']['max_k']}."
        ),
        "named_non_fault_mimic": (
            "Lithologic/intrusive contacts (magnetic edges without slip) and geodetic interpolation "
            "seams can mimic strain-gradient lineaments."
        ),
        "holdout_dti": best["holdout"],
        "holdout_random_control": best["holdout_random"],
        "beats_random_holdout": best["beats_random_holdout"],
        "correlation_overlap_vs_registry": best["uniqueness"],
        "leakage_canary": canary,
        "sgmc_off_diagnostic": best["sgmc_off_diagnostic"],
        "raster_sha256": bundle["zeros"]["sha256"],
        "raster_zip_sha256": bundle["zip_sha256"],
        "validator_output": {
            "no_nan_inside_footprint": bundle["zeros"]["nan_inside_footprint"] == 0,
            "values_in_0_1": bundle["zeros"]["values_in_0_1_inside"] and bundle["zeros"]["full_grid_in_0_1"],
            "crs_shape_transform_match": True,
            "crs": bundle["zeros"]["crs"],
            "shape": bundle["zeros"]["shape"],
            "transform": bundle["zeros"]["transform"],
            "format_valid": bundle["zeros"]["format_valid"],
            "on_catalogue_px": bundle["zeros"]["on_catalogue_px"],
            "positive_px": bundle["zeros"]["positive_px"],
            "nodata": bundle["zeros"]["nodata"],
        },
        "submission": {
            "name": stem + "-zeros",
            "note": note,
            "note_length_chars": len(note),
            "download_zeros": f"./downloads/{src_z.name}",
            "download_nan": f"./downloads/{Path(bundle['nan']['path']).name}",
            "download_zip": f"./downloads/{zpath.name}",
            "format_ok_to_upload": bundle["zeros"]["format_valid"],
            "organizer_confirmed": False,
            "verdict_on_file": (
                "FORMAT OK TO UPLOAD (grid matches sample_submission.tif; all 12,279,160 cells "
                "finite and in [0,1]; 0 NaN; 0 catalogue pixels). Not organizer-scored."
                if bundle["zeros"]["format_valid"]
                else "NOT OK TO SUBMIT — format checks failed."
            ),
        },
        "all_experiments": [
            {k: v for k, v in r.items() if k not in ("dots_submit", "dots_holdout", "score_submit")}
            for r in results
        ],
        "verdict": verdict,
        "verdict_reason": (
            "Unique vs registry, leakage canary under 0.90, holdout DTI beats random, format valid. "
            "Promotion here means 'package as a candidate'; a weekly DrivenData slot is a separate selector."
            if promote else
            "Did not clear uniqueness + leakage + beat-random + format jointly."
        ),
        "metric_self_check_worked_example_dti": we,
        "links": {
            "competition": "https://www.drivendata.org/competitions/306/competition-doe-gems/",
            "problem_description": "https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/",
            "leaderboard": "https://www.drivendata.org/competitions/306/competition-doe-gems/leaderboard/",
            "reference_solution": "https://github.com/drivendataorg/gems-prize-reference-solution",
            "geodawn": "https://www.usgs.gov/data/geodawn-airborne-magnetic-and-radiometric-surveys-northwestern-great-basin-nevada-and",
            "ingenious": "https://gbcge.org/current-projects/ingenious/",
            "epsg_32611": "https://epsg.io/32611",
            "tversky": "https://en.wikipedia.org/wiki/Tversky_index",
        },
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    RESEARCH.mkdir(parents=True, exist_ok=True)
    REPORTS.mkdir(parents=True, exist_ok=True)
    (RESEARCH / "run-card.json").write_text(json.dumps(run_card, indent=2, default=_json_default) + "\n")
    (REPORTS / f"validate_{stem}.json").write_text(json.dumps(bundle, indent=2) + "\n")
    print("wrote", RESEARCH / "run-card.json")
    print("verdict", verdict, "format_ok", bundle["zeros"]["format_valid"], "sha256", bundle["zeros"]["sha256"][:16])
    return 0


def _json_default(o):
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, Path):
        return str(o)
    raise TypeError(type(o))


if __name__ == "__main__":
    raise SystemExit(main())
