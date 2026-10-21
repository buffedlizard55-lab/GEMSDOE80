"""Spatially blocked, whole-segment hide-and-recover holdout (parallel-run protocol).

- Withhold whole connected components of the visible catalogue, plus a buffer.
- Catalogue-based features must be derived only from the remaining visible faults.
- Visible faults are masked pixel-exactly in the prediction before scoring.
- Score pooled DTI (alpha 0.2, beta 0.8, 300 m triangular kernel).
- 95% CI is a component-resample bootstrap of the withheld set.
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
    """Resample withheld components with replacement; recompute DTI each time."""
    rng = np.random.default_rng(seed + 7)
    scores = np.empty(n_boot, dtype=np.float64)
    for i in range(n_boot):
        draw = np.unique(rng.choice(held_ids, size=held_ids.size, replace=True))
        truth = np.isin(comp, draw) & valid
        if not truth.any():
            scores[i] = 0.0
            continue
        scores[i] = dti_binary(pred, truth, valid=valid)["DTI"]
    return float(np.percentile(scores, 2.5)), float(np.percentile(scores, 97.5))


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
    }
