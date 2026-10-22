"""Exact 6-Shell Discrete Choquet-Survival Dinkelbach Submodular Allocator (H80-G4).

Mathematical derivation:
  Let Y_y ~ Bernoulli(p(y)) be calibrated pixel fault probabilities and k(d) = max(1 - d/3, 0).
  On Z^2, k(||y - x||_2) > 0 holds on exactly 25 integer offsets partitioned into 6 shells:
      v_0 = 1.0000 (d=0),  v_1 = 2/3 (d=1),  v_2 = 1 - sqrt(2)/3 (d=sqrt(2)),
      v_3 = 1/3 (d=2),     v_4 = 1 - sqrt(5)/3 (d=sqrt(5)),  v_5 = 1 - sqrt(8)/3 (d=sqrt(8)),  v_6 = 0.

  1. Exact Expected Self-Credit (Choquet survival integral — zero clipping bias):
       E[max_{y in B_3(x)} k(||y - x||_2) Y_y]
         = sum_{m=0}^5 (v_m - v_{m+1}) * [1 - exp( sum_{||y-x||_2 <= r_m} ln(1 - p(y)) )]
     giving exact expected false-positive increment E[FP_w(x)] = 1 - E[max k Y].

  2. Exact Submodular Expected True-Positive Marginal Coverage:
       Delta E[TP_w](x | S) = sum_{y in B_3(x)} p(y) * max(0, k(||y - x||_2) - C_S(y))
     where C_S(y) = max_{s in S} k(||y - s||_2) is non-decreasing as dots are added, making
     E[TP_w](S) strictly monotone submodular.

  3. Exact Dinkelbach Fractional Programming (Dinkelbach 1967, Management Science 13:492-498):
     Maximizing E[DTI](S) = E[TP_w](S) / (E[TP_w](S) + alpha * E[FP_w](S) + beta * N_truth)
     is solved iteratively via linear-fractional parameter lambda_t = E[DTI](S_{t-1}) with
     Minoux (1978) lazy-greedy priority-queue evaluation:
         gain_lambda(x | S) = (1 - alpha * lambda) * Delta E[TP_w](x | S) - alpha * lambda * E[FP_w(x)].
"""

from __future__ import annotations

import heapq
import math

import numpy as np
from scipy.ndimage import convolve

from .metric import ALPHA, BETA, EPS, KERNEL_OFFSETS

SHELL_SPECS: list[tuple[float, float]] = [
    (0.0, 1.0),
    (1.0, 2.0 / 3.0),
    (math.sqrt(2.0), 1.0 - math.sqrt(2.0) / 3.0),
    (2.0, 1.0 / 3.0),
    (math.sqrt(5.0), 1.0 - math.sqrt(5.0) / 3.0),
    (math.sqrt(8.0), 1.0 - math.sqrt(8.0) / 3.0),
]


def choquet_expected_self_credit(prob: np.ndarray) -> np.ndarray:
    """Exact closed-form E[max_{y in B_3(x)} k(||y - x||_2) Y_y] on the integer lattice."""
    p = np.clip(np.nan_to_num(prob, nan=0.0), 0.0, 0.999999).astype(np.float64)
    log_surv = np.log1p(-p)
    yy, xx = np.mgrid[-3:4, -3:4]
    dist_mat = np.hypot(yy, xx)

    total = np.zeros(p.shape, dtype=np.float64)
    for m, (r_m, v_m) in enumerate(SHELL_SPECS):
        v_next = SHELL_SPECS[m + 1][1] if m + 1 < len(SHELL_SPECS) else 0.0
        disk = (dist_mat <= r_m + 1e-9).astype(np.float64)
        sum_log_q = convolve(log_surv, disk, mode="constant", cval=0.0)
        hit_at_least_vm = -np.expm1(sum_log_q)
        total += (v_m - v_next) * hit_at_least_vm
    return np.clip(total, 0.0, 1.0).astype(np.float32)


def dinkelbach_choquet_allocate(
    prob: np.ndarray,
    allowed: np.ndarray,
    *,
    eval_mask: np.ndarray | None = None,
    min_dots: int = 26000,
    max_dots: int = 40000,
    min_sep_px: float = 2.0,
    n_outer_iters: int = 3,
    max_candidates: int = 260000,
) -> tuple[np.ndarray, dict]:
    """Allocate binary dots maximizing exact Choquet expected DTI via Dinkelbach + Minoux heap."""
    cov_mask = allowed if eval_mask is None else np.asarray(eval_mask, dtype=bool)
    # Strictly zero out visible/outside-footprint cells so candidate dots near known faults never
    # receive phantom coverage credit from visible fault pixels excluded by the official scorer!
    p_cov = np.where(cov_mask & np.isfinite(prob), np.clip(prob, 0.0, 1.0), 0.0).astype(np.float64)
    n_truth_est = float(p_cov.sum())
    if n_truth_est <= 0.0 or not allowed.any():
        return np.zeros(prob.shape, dtype=bool), {"n_dots": 0, "expected_dti": 0.0, "lambda_trace": []}

    e_self = choquet_expected_self_credit(p_cov).astype(np.float64)
    e_fp = 1.0 - e_self

    # Initial upper bound on Delta E[TP_w](x | empty) = sum_y p_cov(y) k(x, y)
    init_tp = np.zeros(p_cov.shape, dtype=np.float64)
    h, w = p_cov.shape
    for dy, dx, kv in KERNEL_OFFSETS:
        y0, y1 = max(0, dy), min(h, h + dy)
        x0, x1 = max(0, dx), min(w, w + dx)
        init_tp[y0:y1, x0:x1] += p_cov[y0 - dy : y1 - dy, x0 - dx : x1 - dx] * kv

    cand_mask = allowed & (init_tp > 0.0)
    cand_r, cand_c = np.nonzero(cand_mask)
    if cand_r.size > max_candidates:
        top_idx = np.argpartition(init_tp[cand_r, cand_c], -max_candidates)[-max_candidates:]
        cand_r, cand_c = cand_r[top_idx], cand_c[top_idx]

    n_cand = int(cand_r.size)
    offs = [(int(dy), int(dx), float(kv)) for dy, dx, kv in KERNEL_OFFSETS]

    # Precompute exclusion offsets for minimum Euclidean dot spacing `min_sep_px`
    sep_r = int(math.ceil(min_sep_px))
    excl_offs = [
        (dy, dx)
        for dy in range(-sep_r, sep_r + 1)
        for dx in range(-sep_r, sep_r + 1)
        if math.hypot(dy, dx) < min_sep_px - 1e-6
    ]

    lam = 0.24
    lam_trace: list[float] = []
    best_emitted = np.zeros((h, w), dtype=bool)
    best_stats: dict = {}

    for it in range(n_outer_iters):
        c_tp = 1.0 - ALPHA * lam
        c_fp = ALPHA * lam

        cov = np.zeros((h, w), dtype=np.float64)
        blocked = np.zeros((h, w), dtype=bool)
        emitted = np.zeros((h, w), dtype=bool)

        init_gains = c_tp * init_tp[cand_r, cand_c] - c_fp * e_fp[cand_r, cand_c]
        order = np.argsort(-init_gains)
        # Min-heap over negative upper-bound gain: (-gain_ub, tie_idx, r, c)
        heap: list[tuple[float, int, int, int]] = [
            (-float(init_gains[idx]), int(rank), int(cand_r[idx]), int(cand_c[idx]))
            for rank, idx in enumerate(order)
        ]
        heapq.heapify(heap)

        tot_tp = 0.0
        tot_fp = 0.0
        count = 0

        while heap and count < max_dots:
            neg_ub, rank, r, c = heapq.heappop(heap)
            if blocked[r, c]:
                continue
            # Re-evaluate exact marginal expected TP gain given current coverage `cov`
            dtp = 0.0
            for dy, dx, kv in offs:
                nr, nc = r + dy, c + dx
                if 0 <= nr < h and 0 <= nc < w:
                    rem = kv - cov[nr, nc]
                    if rem > 0.0:
                        dtp += p_cov[nr, nc] * rem
            exact_gain = c_tp * dtp - c_fp * e_fp[r, c]
            # FIRST: Minoux (1978) lazy re-insertion! `exact_gain` is ONLY the global maximum
            # over all remaining candidates if `-exact_gain <= heap[0][0]`.
            if heap and (-exact_gain > heap[0][0]):
                heapq.heappush(heap, (-exact_gain, rank, r, c))
                continue

            # SECOND: Now `(r, c)` is provably the maximum-gain remaining pixel across the entire map.
            # Stop at Dinkelbach's exact optimality bar (`exact_gain <= 0`).
            if exact_gain <= 0.0 and count >= min_dots:
                break

            # Accept dot (r, c)!
            emitted[r, c] = True
            count += 1
            tot_tp += dtp
            tot_fp += float(e_fp[r, c])
            for dy, dx, kv in offs:
                nr, nc = r + dy, c + dx
                if 0 <= nr < h and 0 <= nc < w and kv > cov[nr, nc]:
                    cov[nr, nc] = kv
            for dy, dx in excl_offs:
                nr, nc = r + dy, c + dx
                if 0 <= nr < h and 0 <= nc < w:
                    blocked[nr, nc] = True

        exp_dti = tot_tp / (tot_tp + ALPHA * tot_fp + BETA * n_truth_est + EPS)
        lam_trace.append(float(exp_dti))
        best_emitted = emitted
        best_stats = {
            "n_dots": int(count),
            "expected_tp_w": float(tot_tp),
            "expected_fp_w": float(tot_fp),
            "expected_n_truth": float(n_truth_est),
            "expected_dti": float(exp_dti),
            "lambda_used": float(lam),
            "lambda_trace": list(lam_trace),
        }
        if abs(exp_dti - lam) < 1e-3:
            break
        lam = float(exp_dti)

    return best_emitted, best_stats
