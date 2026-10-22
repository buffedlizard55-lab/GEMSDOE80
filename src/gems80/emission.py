"""NMS dotted emission and catalogue-residual strain-accommodation surface."""

from __future__ import annotations

import numpy as np
from scipy.ndimage import maximum_filter, sobel


def percentile_norm(arr: np.ndarray, valid: np.ndarray, lo: float = 2.0, hi: float = 98.0) -> np.ndarray:
    out = np.zeros(arr.shape, dtype=np.float32)
    if not valid.any():
        return out
    vals = arr[valid]
    p_lo, p_hi = np.percentile(vals, [lo, hi])
    if not np.isfinite(p_lo) or not np.isfinite(p_hi) or p_hi <= p_lo:
        return out
    n = np.clip((arr - p_lo) / (p_hi - p_lo), 0.0, 1.0)
    out[valid] = n[valid].astype(np.float32)
    return out


def gradient_magnitude(arr: np.ndarray) -> np.ndarray:
    gx = sobel(arr.astype(np.float64), axis=1)
    gy = sobel(arr.astype(np.float64), axis=0)
    return np.hypot(gx, gy).astype(np.float32)


def along_fabric_derivative(field: np.ndarray, fabric_src: np.ndarray) -> np.ndarray:
    """Absolute directional derivative of `field` along the strike of `fabric_src` edges.

    Magnetic-contact strike is perpendicular to the fabric gradient.
    """
    fx = sobel(fabric_src.astype(np.float64), axis=1)
    fy = sobel(fabric_src.astype(np.float64), axis=0)
    nrm = np.hypot(fx, fy) + 1e-12
    # unit strike: rotate gradient 90 degrees
    tx = -fy / nrm
    ty = fx / nrm
    gx = sobel(field.astype(np.float64), axis=1)
    gy = sobel(field.astype(np.float64), axis=0)
    return np.abs(tx * gx + ty * gy).astype(np.float32)


def distance_residual(signal: np.ndarray, dist_px: np.ndarray, valid: np.ndarray, max_bin: int = 40) -> np.ndarray:
    """Remove the radial median of `signal` vs distance-to-visible-catalogue.

    High residual = structure not explained by merely sitting near a mapped trace.
    """
    bins = np.clip(dist_px.astype(np.int32), 0, max_bin)
    expected_bin = np.zeros(max_bin + 1, dtype=np.float64)
    sig = signal.astype(np.float64)
    vr = valid.ravel()
    br = bins.ravel()[vr]
    sr = sig.ravel()[vr]
    for b in range(max_bin + 1):
        sel = sr[br == b]
        if sel.size:
            expected_bin[b] = float(np.median(sel))
    expected = expected_bin[bins]
    resid = np.maximum(sig - expected, 0.0)
    resid[~valid] = 0.0
    return resid.astype(np.float32)


def nms_dots(score: np.ndarray, valid: np.ndarray, sep: int, max_k: int, min_k: int = 8000) -> np.ndarray:
    """Keep spatially separated local maxima of `score` (binary mask)."""
    fill = np.where(valid, score, -1.0).astype(np.float32)
    # Break residual plateaus so maximum_filter returns isolated peaks, not whole flats.
    H, W = score.shape
    yy = np.arange(H, dtype=np.float32)[:, None]
    xx = np.arange(W, dtype=np.float32)[None, :]
    fill = fill + np.float32(1e-7) * np.mod(yy * np.float32(12.9898) + xx * np.float32(78.233), 1.0)
    fill = np.where(valid, fill, np.float32(-1.0))
    size = 2 * sep + 1
    local = maximum_filter(fill, size=size)
    peaks = valid & (fill >= local) & (score > 0)
    ys, xs = np.nonzero(peaks)
    if ys.size == 0:
        return np.zeros(score.shape, dtype=bool)
    vals = fill[ys, xs]
    if ys.size > 250_000:
        # keep the strongest candidates only before the greedy pass
        top = np.argpartition(-vals, 250_000)[:250_000]
        ys, xs, vals = ys[top], xs[top], vals[top]
    order = np.argsort(-vals)
    keep = np.zeros(score.shape, dtype=bool)
    suppressed = np.zeros(score.shape, dtype=bool)
    taken = 0
    for i in order:
        y = int(ys[i])
        x = int(xs[i])
        if suppressed[y, x]:
            continue
        keep[y, x] = True
        taken += 1
        y0, y1 = max(0, y - sep), min(H, y + sep + 1)
        x0, x1 = max(0, x - sep), min(W, x + sep + 1)
        suppressed[y0:y1, x0:x1] = True
        if taken >= max_k:
            break
    return keep


def auc_mannwhitney(scores: np.ndarray, labels: np.ndarray) -> float:
    scores = np.asarray(scores, dtype=np.float64)
    labels = np.asarray(labels, dtype=np.int64)
    pos = labels == 1
    neg = labels == 0
    n_pos = int(pos.sum())
    n_neg = int(neg.sum())
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    order = np.argsort(scores, kind="mergesort")
    ranks = np.empty(scores.size, dtype=np.float64)
    sorted_scores = scores[order]
    i = 0
    while i < scores.size:
        j = i
        while j + 1 < scores.size and sorted_scores[j + 1] == sorted_scores[i]:
            j += 1
        avg_rank = 0.5 * ((i + 1) + (j + 1))
        ranks[order[i : j + 1]] = avg_rank
        i = j + 1
    u = float(ranks[pos].sum()) - n_pos * (n_pos + 1) / 2.0
    return u / (n_pos * n_neg)
