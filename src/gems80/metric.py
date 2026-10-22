"""Distance-weighted Tversky index (official GEMS Prize metric).

Source (verbatim equations): DrivenData competition #306 problem description
https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/#performance-metric

    k(d) = max(1 - d/R, 0),  R = 300 m = 3 px at 100 m
    TP_w = sum_{g in G} max_{x: d(x,g) <= R} p(x) k(d(x,g))
    FP_w = sum_{x: p(x)>0} p(x) [1 - max_{g in G} k(d(x,g))]
    FN_w = |G| - TP_w
    DTI(alpha=0.2, beta=0.8) = TP_w / (TP_w + alpha FP_w + beta FN_w + eps)

Worked example on that page: TP_w=3.00, FP_w=1.89, FN_w=2.00 -> 3/(3+0.378+1.6)=0.60265 (~0.60).
"""

from __future__ import annotations

import math

import numpy as np
from scipy.ndimage import distance_transform_edt

ALPHA = 0.2
BETA = 0.8
RADIUS_PX = 3.0
RADIUS_M = 300.0
EPS = 1e-12
EVALUATOR_VERSION = (
    "GEMSDOE80-evaluate_holdout.py-v2 "
    "(DTI alpha=0.2 beta=0.8 R=300m triangular kernel, Voronoi multinomial cluster CI)"
)


def kernel(d, radius: float = RADIUS_PX):
    return np.maximum(1.0 - np.asarray(d, dtype=np.float64) / radius, 0.0)


def kernel_offsets(radius: float = RADIUS_PX) -> list[tuple[int, int, float]]:
    """All integer grid offsets (dy, dx, k) with strictly positive kernel weight k > 0."""
    r_int = int(math.ceil(radius))
    out: list[tuple[int, int, float]] = []
    for dy in range(-r_int, r_int + 1):
        for dx in range(-r_int, r_int + 1):
            d = math.hypot(dy, dx)
            kv = max(1.0 - d / radius, 0.0)
            if kv > 0.0:
                out.append((dy, dx, float(kv)))
    out.sort(key=lambda t: (-t[2], abs(t[0]) + abs(t[1])))
    return out


KERNEL_OFFSETS = kernel_offsets(RADIUS_PX)


def dti_from_components(tp: float, fp: float, fn: float, alpha: float = ALPHA, beta: float = BETA) -> float:
    return float(tp / (tp + alpha * fp + beta * fn + EPS))


def dti_binary(pred_bool: np.ndarray, truth_bool: np.ndarray, valid: np.ndarray | None = None) -> dict:
    """Exact DTI for binary {0,1} predictions via Euclidean distance transforms."""
    p = np.asarray(pred_bool, dtype=bool)
    g = np.asarray(truth_bool, dtype=bool)
    if valid is not None:
        valid = np.asarray(valid, dtype=bool)
        p = p & valid
        g = g & valid
    n = int(g.sum())
    n_pred = int(p.sum())
    if n == 0:
        return {"TP_w": 0.0, "FP_w": float(n_pred), "FN_w": 0.0, "n_truth": 0, "n_pred": n_pred, "DTI": 0.0}
    if n_pred == 0:
        return {"TP_w": 0.0, "FP_w": 0.0, "FN_w": float(n), "n_truth": n, "n_pred": 0, "DTI": 0.0}
    dp = distance_transform_edt(~p)
    tp = float(kernel(dp[g]).sum())
    fn = float(n) - tp
    dg = distance_transform_edt(~g)
    fp = float((1.0 - kernel(dg[p])).sum())
    return {
        "TP_w": tp,
        "FP_w": fp,
        "FN_w": fn,
        "n_truth": n,
        "n_pred": n_pred,
        "DTI": dti_from_components(tp, fp, fn),
    }


def dti_continuous(pred_prob: np.ndarray, truth_bool: np.ndarray, valid: np.ndarray | None = None) -> dict:
    """Exact DTI for general continuous predictions p(x) in [0, 1]."""
    p = np.where(np.isfinite(pred_prob), np.clip(pred_prob, 0.0, 1.0), 0.0).astype(np.float64)
    g = np.asarray(truth_bool, dtype=bool)
    if valid is not None:
        valid = np.asarray(valid, dtype=bool)
        p = np.where(valid, p, 0.0)
        g = g & valid
    n = int(g.sum())
    pos = p > 0.0
    if n == 0:
        fp = float(p[pos].sum())
        return {"TP_w": 0.0, "FP_w": fp, "FN_w": 0.0, "n_truth": 0, "mass": fp, "DTI": 0.0}
    if not pos.any():
        return {"TP_w": 0.0, "FP_w": 0.0, "FN_w": float(n), "n_truth": n, "mass": 0.0, "DTI": 0.0}
    h, w = p.shape
    local_max = np.zeros((h, w), dtype=np.float64)
    for dy, dx, kv in KERNEL_OFFSETS:
        y0, y1 = max(0, dy), min(h, h + dy)
        x0, x1 = max(0, dx), min(w, w + dx)
        local_max[y0:y1, x0:x1] = np.maximum(
            local_max[y0:y1, x0:x1],
            p[y0 - dy : y1 - dy, x0 - dx : x1 - dx] * kv,
        )
    tp = float(local_max[g].sum())
    fn = float(n) - tp
    dg = distance_transform_edt(~g)
    fp = float((p[pos] * (1.0 - kernel(dg[pos]))).sum())
    return {
        "TP_w": tp,
        "FP_w": fp,
        "FN_w": fn,
        "n_truth": n,
        "mass": float(p[pos].sum()),
        "DTI": dti_from_components(tp, fp, fn),
    }


def worked_example_dti() -> float:
    """Published worked example: 3 / (3 + 0.2*1.89 + 0.8*2) = 0.602651..."""
    return dti_from_components(3.00, 1.89, 2.00)
