"""Spatially blocked, whole-segment hide-and-recover holdout (parallel-run protocol).

- Withhold whole connected components of the visible catalogue, plus a buffer.
- Catalogue-based features must be derived only from the remaining visible faults.
- Visible faults are masked pixel-exactly in the prediction before scoring.
- Score pooled DTI (alpha 0.2, beta=0.8, 300 m triangular kernel).
- 95% CI uses an exact Voronoi multinomial cluster bootstrap across the M withheld
  components (fixing v1's biased `np.unique(rng.choice(...))` which dropped 36.8%
  of withheld components per draw and undervalued TP_w while inflating FP_w).
"""

from __future__ import annotations

import numpy as np
from scipy.ndimage import binary_dilation, distance_transform_edt, label

from .emission import auc_mannwhitney
from .metric import ALPHA, BETA, EPS, EVALUATOR_VERSION, dti_binary, dti_from_components, kernel

HOLD_FRAC = 0.30
BUFFER_PX = 3
BOOTSTRAP_N = 1000
BOOTSTRAP_SEED = 20261021


def connected_components(mask: np.ndarray) -> tuple[np.ndarray, int]:
    structure = np.ones((3, 3), dtype=int)
    comp, n = label(mask, structure=structure)
    return comp, int(n)


def withhold_segments(catalogue: np.ndarray, foot: np.ndarray, hold_frac: float = HOLD_FRAC,
                      buffer_px: int = BUFFER_PX, seed: int = BOOTSTRAP_SEED) -> dict:
    cat = catalogue & foot
    comp, n = connected_components(cat)
    if n == 0:
        raise RuntimeError("no catalogue components to withhold")
    sizes = np.bincount(comp.ravel(), minlength=n + 1)
    rng = np.random.default_rng(seed)
    order = rng.permutation(np.arange(1, n + 1))
    target = hold_frac * float(cat.sum())
    acc = 0
    held_ids = []
    for cid in order:
        held_ids.append(int(cid))
        acc += int(sizes[cid])
        if acc >= target and len(held_ids) >= max(1, int(0.15 * n)):
            break
    held_ids = np.array(held_ids, dtype=int)
    hidden = np.isin(comp, held_ids)
    hidden_buf = binary_dilation(hidden, iterations=buffer_px) if buffer_px else hidden
    visible = cat & ~hidden_buf
    return {
        "n_components": n,
        "held_ids": held_ids,
        "comp": comp,
        "hidden": hidden & foot,
        "hidden_buffer": hidden_buf & foot,
        "visible": visible,
        "n_hidden_px": int((hidden & foot).sum()),
        "n_visible_px": int(visible.sum()),
        "hold_frac_pixels": float((hidden & foot).sum()) / float(cat.sum()),
    }


def bootstrap_dti_ci(pred: np.ndarray, hidden: np.ndarray, comp: np.ndarray, held_ids: np.ndarray,
                     valid: np.ndarray, n_boot: int = BOOTSTRAP_N, seed: int = BOOTSTRAP_SEED) -> tuple[float, float]:
    """Unbiased Voronoi-partition multinomial cluster bootstrap across withheld components.

    Every pixel in `valid` is partitioned into the Voronoi cell V_m of its nearest
    withheld fault component C_m (m = 0 .. M-1), giving exact additive components:
        sum_m tp_m == TP_w,   sum_m fp_m == FP_w,   sum_m fn_m == FN_w.
    Resampling component multiplicities u ~ Multinomial(M, [1/M]*M) preserves exact
    first moments E[u_m] = 1 so the bootstrap distribution is centered cleanly on the
    pooled ratio estimator DTI(1) without component-collision loss.
    """
    p = np.asarray(pred, dtype=bool) & valid
    g = np.asarray(hidden, dtype=bool) & valid
    m_count = int(held_ids.size)
    if m_count == 0 or not g.any():
        return 0.0, 0.0

    # Map component IDs 1..N -> 0..M-1
    max_cid = int(max(int(comp.max()), int(held_ids.max())))
    lut = np.full(max_cid + 1, -1, dtype=np.int32)
    lut[held_ids] = np.arange(m_count, dtype=np.int32)

    # Distance from every pixel to nearest withheld truth pixel + Voronoi owner index
    dg, nearest_rc = distance_transform_edt(~g, return_indices=True)
    owner_cid = comp[nearest_rc[0], nearest_rc[1]]
    owner_idx = lut[owner_cid]

    # Exact per-component TP_w and FN_w
    if p.any():
        dp = distance_transform_edt(~p)
        kw_g = kernel(dp[g])
    else:
        kw_g = np.zeros(int(g.sum()), dtype=np.float64)
    g_owner = lut[comp[g]]
    tp_vec = np.bincount(g_owner, weights=kw_g, minlength=m_count).astype(np.float64)
    sz_vec = np.bincount(g_owner, minlength=m_count).astype(np.float64)
    fn_vec = sz_vec - tp_vec

    # Exact per-Voronoi-cell FP_w
    if p.any():
        fp_weights = 1.0 - kernel(dg[p])
        p_owner = np.clip(owner_idx[p], 0, m_count - 1)
        fp_vec = np.bincount(p_owner, weights=fp_weights, minlength=m_count).astype(np.float64)
    else:
        fp_vec = np.zeros(m_count, dtype=np.float64)

    rng = np.random.default_rng(seed + 7)
    probs = np.full(m_count, 1.0 / m_count, dtype=np.float64)
    mult = rng.multinomial(m_count, probs, size=n_boot).astype(np.float64)
    tp_b = mult @ tp_vec
    fp_b = mult @ fp_vec
    fn_b = mult @ fn_vec
    denom = tp_b + ALPHA * fp_b + BETA * fn_b + EPS
    scores = np.where(denom > 0.0, tp_b / denom, 0.0)
    return float(np.percentile(scores, 2.5)), float(np.percentile(scores, 97.5))


def quadrant_dti_breakdown(pred: np.ndarray, truth: np.ndarray, foot: np.ndarray) -> dict[str, float]:
    """Per-quadrant DTI (NW, NE, SW, SE split at footprint median row/col)."""
    ys, xs = np.nonzero(foot)
    rmid = int(np.median(ys))
    cmid = int(np.median(xs))
    rr, cc = np.indices(foot.shape)
    quads = {
        "NW": foot & (rr < rmid) & (cc < cmid),
        "NE": foot & (rr < rmid) & (cc >= cmid),
        "SW": foot & (rr >= rmid) & (cc < cmid),
        "SE": foot & (rr >= rmid) & (cc >= cmid),
    }
    return {name: float(dti_binary(pred, truth, valid=qmask)["DTI"]) for name, qmask in quads.items()}


def score_holdout(pred_bool: np.ndarray, split: dict, foot: np.ndarray, bootstrap: bool = True) -> dict:
    """Mask visible faults pixel-exactly, score against withheld positives."""
    pred = np.asarray(pred_bool, dtype=bool) & foot & ~split["visible"]
    truth = split["hidden"] & foot
    res = dti_binary(pred, truth, valid=foot)
    if bootstrap:
        lo, hi = bootstrap_dti_ci(pred, split["hidden"], split["comp"], split["held_ids"], foot)
        ci = [lo, hi]
    else:
        ci = None
    return {
        "label": "HOLDOUT-DTI",
        "evaluator_version": EVALUATOR_VERSION,
        "withheld_positive_count": int(truth.sum()),
        "n_components_total": int(split["n_components"]),
        "n_components_held": int(split["held_ids"].size),
        "hold_frac_pixels": split["hold_frac_pixels"],
        "score": res["DTI"],
        "ci95": ci,
        "TP_w": res["TP_w"],
        "FP_w": res["FP_w"],
        "FN_w": res["FN_w"],
        "n_pred_after_visible_mask": int(pred.sum()),
        "spatial_quadrant_dti": quadrant_dti_breakdown(pred, truth, foot),
    }


def leakage_canary(
    features: dict[str, np.ndarray],
    split: dict,
    foot: np.ndarray,
    rng: np.random.Generator,
    max_neg: int = 200_000,
    auc_limit: float = 0.90,
) -> dict:
    """Protocol Rule 4: per-feature holdout AUC > 0.90 means leakage until proven otherwise."""
    eval_domain = foot & ~split["visible"]
    pos_idx = np.flatnonzero((split["hidden"] & eval_domain).ravel())
    neg_pool = np.flatnonzero((~split["hidden"] & eval_domain).ravel())
    neg_idx = rng.choice(neg_pool, size=min(max_neg, int(neg_pool.size)), replace=False)
    idx = np.concatenate([pos_idx, neg_idx])
    y = np.concatenate([np.ones(pos_idx.size, dtype=np.int64), np.zeros(neg_idx.size, dtype=np.int64)])

    rows = []
    max_auc = 0.0
    flagged = False
    for name, arr in features.items():
        s = np.nan_to_num(arr.ravel()[idx], nan=0.0).astype(np.float64)
        auc = float(auc_mannwhitney(s, y))
        dev_auc = max(auc, 1.0 - auc)
        max_auc = max(max_auc, dev_auc)
        hit = dev_auc > auc_limit
        flagged = flagged or hit
        rows.append({"feature": name, "auc": round(auc, 5), "effective_auc": round(dev_auc, 5), "leakage": bool(hit)})
    rows.sort(key=lambda r: -r["effective_auc"])
    return {
        "leakage_flag": bool(flagged),
        "threshold_auc": auc_limit,
        "max_auc": round(max_auc, 5),
        "top_15_features": rows[:15],
        "all_features": rows,
    }
