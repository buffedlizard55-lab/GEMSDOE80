"""Discrete emission operators for the GEMS Distance-Weighted Tversky Index
(alpha = 0.2, beta = 0.8, triangular kernel radius R = 300m = 3.0 px).

Implements:
  1. `spaced_nms_ranked`: fast non-maximum-suppressed blue-noise ranking with
     minimum Euclidean separation `min_sep_px` (e.g. 2.8 px = 280m, the H33/D2.8 standard).
  2. `skeleton_ranked`: 1-px 3x3 ridge local-maximum skeleton ranking.
  3. `celf_submodular_ranked`: Cost-Effective Lazy Forward (CELF, Leskovec et al. 2007)
     min-heap submodular maximization of expected 300m triangular-kernel coverage
     sum_g p_cal(g) * max_{x in S} k(d(x, g)), returning the exact ordered pick list
     so both prefix budgets and Dinkelbach fractional-programming stopping rules run in < 2s.
"""

from __future__ import annotations

import heapq
import math
import numpy as np
from scipy import ndimage as ndi

from gems57 import metric as M


def disk_offsets(radius_px: float) -> tuple[np.ndarray, np.ndarray]:
    r = int(math.floor(radius_px))
    dy, dx = [], []
    for j in range(-r, r + 1):
        for i in range(-r, r + 1):
            if j * j + i * i < radius_px * radius_px - 1e-9:
                dy.append(j)
                dx.append(i)
    return np.array(dy, dtype=np.int32), np.array(dx, dtype=np.int32)


def spaced_nms_ranked(surface: np.ndarray, allowed: np.ndarray, min_sep_px: float = 2.8, max_picks: int = 80_000):
    """Greedy non-maximum suppression with exact Euclidean exclusion disk `min_sep_px`.

    Returns `(rank_grid, picked_scores)` where `rank_grid[y, x] in {1..K}` gives exact 1-based pick rank.
    """
    h, w = surface.shape
    cand = allowed & (surface > 0)
    ys, xs = np.nonzero(cand)
    if ys.size == 0:
        return np.zeros((h, w), dtype=np.int32), np.empty(0, dtype=np.float32)

    vals = surface[ys, xs]
    if vals.size > max_picks * 8:
        top_idx = np.argpartition(vals, -(max_picks * 8))[-(max_picks * 8):]
        ys, xs, vals = ys[top_idx], xs[top_idx], vals[top_idx]

    order = np.argsort(-vals, kind="mergesort")
    ys, xs, vals = ys[order], xs[order], vals[order]

    ody, odx = disk_offsets(min_sep_px)
    blocked = np.zeros((h, w), dtype=bool)
    rank_grid = np.zeros((h, w), dtype=np.int32)
    picked_vals = []
    k = 0
    for y, x, v in zip(ys.tolist(), xs.tolist(), vals.tolist()):
        if blocked[y, x]:
            continue
        k += 1
        rank_grid[y, x] = k
        picked_vals.append(v)
        if k >= max_picks:
            break
        ny = y + ody
        nx = x + odx
        ok = (ny >= 0) & (ny < h) & (nx >= 0) & (nx < w)
        blocked[ny[ok], nx[ok]] = True

    return rank_grid, np.array(picked_vals, dtype=np.float32)


def skeleton_ranked(surface: np.ndarray, allowed: np.ndarray, smooth_sigma: float = 0.7):
    """3x3 local-maximum ridge skeleton ranking on `allowed` pixels."""
    s = ndi.gaussian_filter(surface, sigma=smooth_sigma, mode="nearest") if smooth_sigma > 0 else surface
    mx = ndi.maximum_filter(s, size=3, mode="nearest")
    cand = allowed & (s >= mx - 1e-7) & (s > 0)
    ys, xs = np.nonzero(cand)
    if ys.size == 0:
        return np.zeros(surface.shape, dtype=np.int32)
    order = np.argsort(-s[ys, xs], kind="mergesort")
    rank_grid = np.zeros(surface.shape, dtype=np.int32)
    rank_grid[ys[order], xs[order]] = np.arange(1, order.size + 1, dtype=np.int32)
    return rank_grid


def celf_submodular_ranked(
    posterior: np.ndarray,
    allowed: np.ndarray,
    max_picks: int = 80_000,
    candidate_pool: int = 320_000,
    min_sep_px: float = 2.85,
):
    """Exact CELF lazy-greedy submodular maximizer of expected 300m triangular-kernel
    coverage `F(S) = sum_g posterior(g) * max_{x in S} k(d(x, g))` subject to minimum
    inter-dot spacing `min_sep_px` (prevents 100m carpet double-packing while preferring
    elongated multi-pixel fault crests over isolated 1-pixel spikes).
    """
    p = np.where(np.isfinite(posterior), np.maximum(posterior, 0.0), 0.0).astype(np.float32)
    h, w = p.shape
    off_dy = M.OFF_DY.astype(np.int32)
    off_dx = M.OFF_DX.astype(np.int32)
    off_k = M.OFF_K.astype(np.float32)
    sep_dy, sep_dx = disk_offsets(min_sep_px) if min_sep_px > 0 else (np.empty(0, np.int32), np.empty(0, np.int32))
    blocked = np.zeros((h, w), dtype=bool)

    # Initial upper bound = exact convolution sum_off p(x+off) * k(off)
    init_gain = np.zeros((h, w), dtype=np.float32)
    for dy, dx, kw in zip(off_dy.tolist(), off_dx.tolist(), off_k.tolist()):
        sy0, sy1 = max(0, dy), min(h, h + dy)
        sx0, sx1 = max(0, dx), min(w, w + dx)
        ty0, ty1 = max(0, -dy), min(h, h - dy)
        tx0, tx1 = max(0, -dx), min(w, w - dx)
        init_gain[ty0:ty1, tx0:tx1] += p[sy0:sy1, sx0:sx1] * kw

    cand = allowed & (init_gain > 0)
    ys, xs = np.nonzero(cand)
    if ys.size == 0:
        return np.zeros((h, w), dtype=np.int32), np.empty(0, dtype=np.float32)

    gains = init_gain[ys, xs]
    if gains.size > candidate_pool:
        idx = np.argpartition(gains, -candidate_pool)[-candidate_pool:]
        ys, xs, gains = ys[idx], xs[idx], gains[idx]

    # Build min-heap of (-gain, tie_id, y, x)
    heap = [(-float(g), idx, int(y), int(x)) for idx, (y, x, g) in enumerate(zip(ys, xs, gains))]
    heapq.heapify(heap)

    cover = np.zeros((h, w), dtype=np.float32)
    rank_grid = np.zeros((h, w), dtype=np.int32)
    marginal_gains = []
    picked = 0

    while heap and picked < max_picks:
        neg_g, idx, y, x = heapq.heappop(heap)
        if blocked[y, x]:
            continue
        # Recompute exact marginal gain on the 29-pixel 300m kernel stencil
        ny = y + off_dy
        nx = x + off_dx
        ok = (ny >= 0) & (ny < h) & (nx >= 0) & (nx < w)
        nyv, nxv, kw = ny[ok], nx[ok], off_k[ok]
        delta_k = kw - cover[nyv, nxv]
        pos = delta_k > 0.0
        new_g = float((p[nyv[pos], nxv[pos]] * delta_k[pos]).sum()) if pos.any() else 0.0

        # If heap is empty or `new_g` still beats the next best upper bound `-heap[0][0]`, accept!
        if not heap or new_g >= -heap[0][0] - 1e-9:
            if new_g <= 0.0:
                break
            picked += 1
            rank_grid[y, x] = picked
            marginal_gains.append(new_g)
            cover[nyv, nxv] = np.maximum(cover[nyv, nxv], kw)
            if sep_dy.size:
                by = y + sep_dy
                bx = x + sep_dx
                bok = (by >= 0) & (by < h) & (bx >= 0) & (bx < w)
                blocked[by[bok], bx[bok]] = True
        else:
            heapq.heappush(heap, (-new_g, idx, y, x))

    return rank_grid, np.array(marginal_gains, dtype=np.float32)
