"""Spatially blocked, whole-segment hide-and-recover holdout (parallel-run protocol).

- Withhold whole connected components of the visible catalogue, plus a buffer.
- Catalogue-based features must be derived only from the remaining visible faults.
- Visible faults are masked pixel-exactly in the prediction before scoring.
- Score pooled DTI (alpha 0.2, beta 0.8, 300 m triangular kernel).
- Optional uncertainty is a conditional spatial-block bootstrap of per-pixel DTI contributions;
  it is reported only with at least 20 nonempty blocks and is not a refit-across-splits interval.
"""

from __future__ import annotations

import numpy as np
from scipy.ndimage import binary_dilation, label

from .metric import EVALUATOR_VERSION, dti_binary, dti_contribution_maps, dti_from_components

HOLD_FRAC = 0.30
BUFFER_PX = 3
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


def spatial_block_dti_ci(pred: np.ndarray, truth: np.ndarray, valid: np.ndarray,
                         block_px: int = 50, n_boot: int = 1000,
                         seed: int = BOOTSTRAP_SEED) -> tuple[float | None, float | None, int]:
    """Conditional spatial-block bootstrap for the pooled DTI ratio.

    Contributions are first calculated on the full fixed prediction/truth field,
    then aggregated over non-overlapping 50-pixel (5 km at 100 m) blocks. Blocks
    are resampled with replacement as units. This estimates spatial sampling
    variation conditional on one hide-and-recover split; it does not include
    uncertainty from refitting on different visible-catalogue splits.
    """
    if block_px <= 0 or n_boot < 1:
        raise ValueError("block_px and n_boot must be positive")
    pred = np.asarray(pred, dtype=bool)
    truth = np.asarray(truth, dtype=bool)
    valid = np.asarray(valid, dtype=bool)
    if pred.shape != truth.shape or pred.shape != valid.shape:
        raise ValueError("prediction, truth, and valid mask must have the same shape")

    tp_map, fp_map, fn_map = dti_contribution_maps(pred, truth, valid)
    block_sums = []
    h, w = valid.shape
    for y0 in range(0, h, block_px):
        y1 = min(y0 + block_px, h)
        for x0 in range(0, w, block_px):
            x1 = min(x0 + block_px, w)
            if not valid[y0:y1, x0:x1].any():
                continue
            block_sums.append((float(tp_map[y0:y1, x0:x1].sum()),
                               float(fp_map[y0:y1, x0:x1].sum()),
                               float(fn_map[y0:y1, x0:x1].sum())))
    n_blocks = len(block_sums)
    if n_blocks < 20:
        return None, None, n_blocks

    sums = np.asarray(block_sums, dtype=np.float64)
    rng = np.random.default_rng(seed)
    scores = np.empty(n_boot, dtype=np.float64)
    for i in range(n_boot):
        counts = rng.multinomial(n_blocks, np.full(n_blocks, 1.0 / n_blocks))
        tp, fp, fn = counts @ sums
        scores[i] = dti_from_components(tp, fp, fn)
    lo, hi = np.percentile(scores, [2.5, 97.5])
    return float(lo), float(hi), n_blocks


def score_holdout(pred_bool: np.ndarray, split: dict, foot: np.ndarray, bootstrap: bool = True) -> dict:
    """Mask visible faults pixel-exactly, score against withheld positives."""
    pred = np.asarray(pred_bool, dtype=bool) & foot & ~split["visible"]
    truth = split["hidden"] & foot
    res = dti_binary(pred, truth, valid=foot)
    if bootstrap:
        lo, hi, n_blocks = spatial_block_dti_ci(pred, truth, foot)
        ci = [lo, hi] if lo is not None and hi is not None else None
        ci_status = ("conditional spatial-block bootstrap; does not include refit/split uncertainty"
                     if ci is not None else f"unavailable: only {n_blocks} nonempty spatial blocks (minimum 20)")
    else:
        ci = None
        n_blocks = None
        ci_status = "not computed"
    return {
        "label": "HOLDOUT-DTI",
        "evaluator_version": EVALUATOR_VERSION,
        "withheld_positive_count": int(truth.sum()),
        "n_components_total": int(split["n_components"]),
        "n_components_held": int(split["held_ids"].size),
        "hold_frac_pixels": split["hold_frac_pixels"],
        "score": res["DTI"],
        "ci95": ci,
        "ci95_method": ci_status,
        "spatial_block_count": n_blocks,
        "TP_w": res["TP_w"],
        "FP_w": res["FP_w"],
        "FN_w": res["FN_w"],
        "n_pred_after_visible_mask": int(pred.sum()),
    }
