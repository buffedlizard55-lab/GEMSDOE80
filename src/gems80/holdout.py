"""Whole-component hide-and-recover proxy (NOT a spatially blocked split).

- Withhold whole connected components of the visible catalogue, plus a buffer.
- Catalogue-based features must be derived only from the remaining visible faults.
- Visible faults are masked pixel-exactly in the prediction before scoring.
- Score pooled DTI (alpha 0.2, beta 0.8, 300 m triangular kernel).
- Uncertainty is a conditional spatial-block bootstrap, not model generalization.
- Segment sampling is NOT a spatially-blocked train/test split.
"""

from __future__ import annotations

import numpy as np
from scipy.ndimage import binary_dilation, label

from .metric import EVALUATOR_VERSION, dti_binary

HOLD_FRAC = 0.30
BUFFER_PX = 3
BOOTSTRAP_N = 80
BOOTSTRAP_SEED = 20261021


def connected_components(mask: np.ndarray) -> tuple[np.ndarray, int]:
    structure = np.ones((3, 3), dtype=int)
    comp, n = label(mask, structure=structure)
    return comp, int(n)


def withhold_segments(catalogue: np.ndarray, foot: np.ndarray, hold_frac: float = HOLD_FRAC,
                      buffer_px: int = BUFFER_PX, seed: int = BOOTSTRAP_SEED) -> dict:
    if not 0 < hold_frac < 1 or buffer_px < 0:
        raise ValueError("hold_frac must be in (0,1); buffer_px must be nonnegative")
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
                     valid: np.ndarray, n_boot: int = 1000, seed: int = BOOTSTRAP_SEED):
    """Conditional 95% percentile interval from additive spatial-block contributions.

    Compute nearest-neighbour credit on the entire domain BEFORE partitioning,
    preserving cross-block kernel credit. Resample blocks with multiplicity and
    recompute the pooled ratio (never deduplicate sampled components).
    This does not include training/split uncertainty or guarantee coverage.
    """
    from scipy.ndimage import distance_transform_edt
    from .metric import kernel, dti_from_components

    p = np.asarray(pred, dtype=bool) & valid
    g = np.asarray(hidden, dtype=bool) & valid
    if n_boot < 2:
        raise ValueError("n_boot must be at least two")
    credit = kernel(distance_transform_edt(~p)) if p.any() else np.zeros(p.shape)
    fp_credit = 1 - kernel(distance_transform_edt(~g)) if g.any() else np.ones(g.shape)
    blocks = []
    ys = np.linspace(0, p.shape[0], min(4, p.shape[0]) + 1, dtype=int)
    xs = np.linspace(0, p.shape[1], min(4, p.shape[1]) + 1, dtype=int)
    for y0, y1 in zip(ys[:-1], ys[1:]):
        for x0, x1 in zip(xs[:-1], xs[1:]):
            sl = np.s_[y0:y1, x0:x1]
            if not valid[sl].any():
                continue
            tp = float(credit[sl][g[sl]].sum())
            blocks.append([tp, float(fp_credit[sl][p[sl]].sum()), float(g[sl].sum()) - tp])
    if len(blocks) < 4:
        return None
    rng = np.random.default_rng(seed + 7)
    parts = np.asarray(blocks)
    draws = rng.integers(len(parts), size=(n_boot, len(parts)))
    totals = parts[draws].sum(axis=1)
    scores = [dti_from_components(*row) for row in totals]
    return tuple(float(x) for x in np.percentile(scores, [2.5, 97.5]))


def score_holdout(pred_bool: np.ndarray, split: dict, foot: np.ndarray, bootstrap: bool = True) -> dict:
    """Mask visible faults pixel-exactly, score against withheld positives."""
    pred = np.asarray(pred_bool, dtype=bool) & foot & ~split["visible"]
    truth = split["hidden"] & foot
    res = dti_binary(pred, truth, valid=foot)
    if bootstrap:
        interval = bootstrap_dti_ci(pred, split["hidden"], split["comp"], split["held_ids"], foot)
        ci = list(interval) if interval is not None else None
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
        "ci95_method": "conditional 4x4 spatial-block percentile bootstrap; 1000 draws; fixed predictions/split; coverage unverified",
        "split_method": "random whole connected components plus buffer; NOT spatially-blocked",
        "promotion_eligible": False,
        "TP_w": res["TP_w"],
        "FP_w": res["FP_w"],
        "FN_w": res["FN_w"],
        "n_pred_after_visible_mask": int(pred.sum()),
    }
