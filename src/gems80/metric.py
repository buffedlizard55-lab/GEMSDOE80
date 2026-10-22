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
EVALUATOR_VERSION = "GEMSDOE80-evaluate_holdout.py-v2 (DTI alpha=0.2 beta=0.8 R=300m triangular kernel; conditional spatial-block bootstrap CI)"


def kernel(d, radius: float = RADIUS_PX):
    return np.maximum(1.0 - np.asarray(d, dtype=np.float64) / radius, 0.0)


def dti_from_components(tp: float, fp: float, fn: float, alpha: float = ALPHA, beta: float = BETA) -> float:
    return float(tp / (tp + alpha * fp + beta * fn + EPS))


def dti_contribution_maps(pred_bool: np.ndarray, truth_bool: np.ndarray,
                          valid: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return per-pixel TP, FP and FN contributions for binary predictions.

    Each truth pixel receives nearest-prediction triangular credit; each predicted
    pixel pays false-positive cost according to distance to the nearest truth.
    These maps also support spatial-block uncertainty estimates.
    """
    p = np.asarray(pred_bool, dtype=bool)
    g = np.asarray(truth_bool, dtype=bool)
    if p.shape != g.shape:
        raise ValueError("prediction and truth arrays must have the same shape")
    if valid is not None:
        valid = np.asarray(valid, dtype=bool)
        if valid.shape != p.shape:
            raise ValueError("valid mask must have the same shape as prediction")
        p = p & valid
        g = g & valid
    n = int(g.sum())
    n_pred = int(p.sum())
    if n == 0:
        return (np.zeros(p.shape, dtype=np.float64),
                p.astype(np.float64), np.zeros(p.shape, dtype=np.float64))
    if n_pred == 0:
        return (np.zeros(p.shape, dtype=np.float64),
                np.zeros(p.shape, dtype=np.float64), g.astype(np.float64))
    dp = distance_transform_edt(~p)
    dg = distance_transform_edt(~g)
    tp_map = g * kernel(dp)
    fn_map = g.astype(np.float64) - tp_map
    fp_map = p * (1.0 - kernel(dg))
    return tp_map, fp_map, fn_map


def dti_binary(pred_bool: np.ndarray, truth_bool: np.ndarray, valid: np.ndarray | None = None) -> dict:
    """Exact DTI for binary {0,1} predictions via Euclidean distance transforms."""
    p = np.asarray(pred_bool, dtype=bool)
    g = np.asarray(truth_bool, dtype=bool)
    if valid is not None:
        valid = np.asarray(valid, dtype=bool)
        if valid.shape != p.shape:
            raise ValueError("valid mask must have the same shape as prediction")
        p = p & valid
        g = g & valid
    tp_map, fp_map, fn_map = dti_contribution_maps(p, g)
    tp, fp, fn = float(tp_map.sum()), float(fp_map.sum()), float(fn_map.sum())
    return {
        "TP_w": tp,
        "FP_w": fp,
        "FN_w": fn,
        "n_truth": int(g.sum()),
        "n_pred": int(p.sum()),
        "DTI": dti_from_components(tp, fp, fn),
    }


def worked_example_dti() -> float:
    """Published worked example: 3 / (3 + 0.2*1.89 + 0.8*2) = 0.602651..."""
    return dti_from_components(3.00, 1.89, 2.00)
