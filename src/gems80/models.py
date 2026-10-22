"""Memory-bounded (`< 1.25 GB` RSS) & CPU-vectorized execution of the 3 preregistered
experiments (`E80-1`, `E80-2`, `E80-3`) for Lane H80-DPVC (`GEMSDOE80`)."""
from __future__ import annotations

import gc
import heapq
import time
from typing import Dict, List, Sequence, Tuple

import numpy as np
from scipy import ndimage as ndi
from scipy.stats import rankdata
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.neural_network import MLPClassifier
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from gems80.features import (
    FOLD_MECH_NAMES,
    StaticStack,
    compute_fold_mechanics_2d,
    fill_axial_context_valid,
)

# Exact 25-cell (`5 x 5`) Euclidean < 300 m (`< 3 px`) positive triangular kernel
KERNEL_300M: Tuple[Tuple[int, int, float], ...] = tuple(
    (dy, dx, float(1.0 - np.hypot(dy, dx) / 3.0))
    for dy in range(-2, 3)
    for dx in range(-2, 3)
    if np.hypot(dy, dx) < 3.0 - 1e-9
)
assert len(KERNEL_300M) == 25


def _rank01(values: np.ndarray) -> np.ndarray:
    v = np.asarray(values, dtype=np.float64)
    n = v.size
    if n == 0:
        return np.empty(0, dtype=np.float32)
    if np.ptp(v) == 0.0:
        return np.full(n, 0.5, dtype=np.float32)
    return (rankdata(v, method="average") / float(n)).astype(np.float32)


def _logit(p: np.ndarray) -> np.ndarray:
    c = np.clip(p, 1e-5, 1.0 - 1e-5)
    return np.log(c / (1.0 - c)).astype(np.float32)


def sample_terrain_balanced_negatives(
    neg_pool_flat: np.ndarray,
    pos_flat: np.ndarray,
    slope10_flat: np.ndarray,
    n_neg_total: int,
    hard_frac: float = 0.35,
    rng: np.random.Generator | None = None,
) -> np.ndarray:
    if rng is None:
        rng = np.random.default_rng(80)
    n_neg_total = min(int(n_neg_total), int(neg_pool_flat.size))
    n_hard = int(round(n_neg_total * hard_frac))
    n_unif = n_neg_total - n_hard

    unif = rng.choice(neg_pool_flat, size=n_unif, replace=False)
    slp_pool = slope10_flat[neg_pool_flat]
    p70 = float(np.percentile(slope10_flat[pos_flat], 68.0))
    steep_cand = neg_pool_flat[slp_pool >= p70]
    if steep_cand.size >= n_hard:
        hard = rng.choice(steep_cand, size=n_hard, replace=False)
    else:
        hard = rng.choice(neg_pool_flat, size=n_hard, replace=False)
    return np.unique(np.concatenate([unif, hard]))


def fit_hgb_fast(
    X_tr: np.ndarray,
    y_tr: np.ndarray,
    seed: int,
    max_iter: int = 125,
    max_leaf_nodes: int = 47,
    min_samples_leaf: int = 100,
    l2_regularization: float = 5.0,
) -> HistGradientBoostingClassifier:
    clf = HistGradientBoostingClassifier(
        learning_rate=0.075,
        max_iter=max_iter,
        max_leaf_nodes=max_leaf_nodes,
        min_samples_leaf=min_samples_leaf,
        l2_regularization=l2_regularization,
        early_stopping=True,
        validation_fraction=0.1,
        n_iter_no_change=12,
        random_state=seed,
    )
    clf.fit(X_tr, y_tr)
    return clf


def predict_rows_chunked(
    model: HistGradientBoostingClassifier,
    X_source: np.ndarray,
    row_indices: np.ndarray,
    cols: Sequence[int] | None = None,
    chunk_size: int = 450_000,
) -> np.ndarray:
    n = row_indices.size
    out = np.zeros(n, dtype=np.float32)
    for start in range(0, n, chunk_size):
        end = min(n, start + chunk_size)
        r_slice = row_indices[start:end]
        sub = np.asarray(X_source[r_slice], dtype=np.float32)
        if cols is not None:
            sub = sub[:, cols]
        out[start:end] = model.predict_proba(sub)[:, 1].astype(np.float32)
    return out


def predict_single_chunked(
    model: HistGradientBoostingClassifier,
    X_source: np.ndarray,
    cols: Sequence[int] | None = None,
    chunk_size: int = 750_000,
) -> np.ndarray:
    n = X_source.shape[0]
    out = np.zeros(n, dtype=np.float32)
    for start in range(0, n, chunk_size):
        end = min(n, start + chunk_size)
        sub = (
            np.asarray(X_source[start:end, cols], dtype=np.float32)
            if cols is not None
            else np.asarray(X_source[start:end], dtype=np.float32)
        )
        out[start:end] = model.predict_proba(sub)[:, 1].astype(np.float32)
    return out


def calibrate_terrain_prevalence(
    prob_2d: np.ndarray,
    slope10_2d: np.ndarray,
    train_pos_mask: np.ndarray,
    train_domain_mask: np.ndarray,
    fp: np.ndarray,
    blend: float = 0.65,
) -> np.ndarray:
    out = np.where(fp & np.isfinite(prob_2d), prob_2d, 0.0).astype(np.float32)
    slp_dom = slope10_2d[train_domain_mask]
    edges = np.quantile(slp_dom, np.linspace(0.0, 1.0, 9))
    edges[0] = -1e9
    edges[-1] = 1e9

    global_prev = float(np.mean(train_pos_mask[train_domain_mask])) + 1e-8
    global_p99 = float(np.percentile(out[train_domain_mask], 99.0)) + 1e-8

    calibrated = out.copy()
    for b in range(8):
        in_bin_tr = train_domain_mask & (slope10_2d >= edges[b]) & (slope10_2d < edges[b + 1])
        in_bin_all = fp & (slope10_2d >= edges[b]) & (slope10_2d < edges[b + 1])
        if not np.any(in_bin_tr) or not np.any(in_bin_all):
            continue
        bin_prev = float(np.mean(train_pos_mask[in_bin_tr])) + 1e-8
        bin_p99 = float(np.percentile(out[in_bin_tr], 99.0)) + 1e-8
        target_ratio = np.clip(bin_prev / global_prev, 0.55, 1.45)
        current_ratio = bin_p99 / global_p99
        mult = float((target_ratio / max(current_ratio, 1e-4)) ** blend)
        calibrated[in_bin_all] = out[in_bin_all] * mult

    mx = float(np.max(calibrated[fp]))
    if mx > 0.0:
        calibrated = np.clip(calibrated / max(mx, 1.0), 0.0, 1.0)
    return calibrated.astype(np.float32)


def submodular_pack_300m(
    surface_2d: np.ndarray,
    allowed_2d: np.ndarray,
    target_n: int = 37654,
    support_quantile: float = 0.90,
    gamma: float = 2.0,
    max_candidates: int = 320_000,
    skeleton_prefilter: bool = False,
) -> Tuple[np.ndarray, int]:
    """Exact unblurred 300 m triangular-kernel lazy-greedy submodular dot packer (Minoux 1978),
    with optional 1-px ridge-skeleton prefiltering to eliminate adjacent-pixel FP mass."""
    h, w = surface_2d.shape
    allow = np.asarray(allowed_2d, dtype=bool)
    if skeleton_prefilter:
        mx3 = ndi.maximum_filter(surface_2d, size=3, mode="nearest")
        skel_cand = allow & (surface_2d >= 0.88 * mx3 - 1e-7)
        if int(skel_cand.sum()) >= target_n * 2:
            allow = skel_cand

    vals = surface_2d[allow]
    if vals.size == 0:
        return np.zeros((h, w), dtype=bool), 0
    if vals.size < target_n:
        target_n = int(vals.size)

    ranks = _rank01(vals)
    pos = np.clip((ranks - support_quantile) / max(1.0 - support_quantile, 1e-6), 0.0, 1.0)
    weights_allow = (pos ** gamma).astype(np.float32) + (1e-5 * ranks).astype(np.float32)

    pad = 2
    hp, wp = h + 2 * pad, w + 2 * pad
    w_pad = np.zeros((hp, wp), dtype=np.float32)
    ay, ax = np.nonzero(allow)
    w_pad[ay + pad, ax + pad] = weights_allow

    g0 = np.zeros((h, w), dtype=np.float32)
    for dy, dx, k in KERNEL_300M:
        g0 += k * w_pad[pad + dy : pad + dy + h, pad + dx : pad + dx + w]

    cand_scores = g0[allow]
    if cand_scores.size > max_candidates:
        part = np.argpartition(cand_scores, -max_candidates)[-max_candidates:]
        part = part[np.argsort(cand_scores[part], kind="mergesort")[::-1]]
    else:
        part = np.argsort(cand_scores, kind="mergesort")[::-1]

    cy = (ay[part] + pad).astype(np.int32)
    cx = (ax[part] + pad).astype(np.int32)
    flat_idx = cy * wp + cx
    init_gains = cand_scores[part]

    heap = [(-float(init_gains[i]), int(flat_idx[i])) for i in range(flat_idx.size)]
    heapq.heapify(heap)

    w_list = w_pad.ravel().tolist()
    cov_list = [0.0] * (hp * wp)
    off_k = tuple((dy * wp + dx, k) for dy, dx, k in KERNEL_300M)

    dots = np.zeros((h, w), dtype=bool)
    picked = 0
    pop = heapq.heappop
    push = heapq.heappush

    while picked < target_n and heap:
        _, idx = pop(heap)
        gain = 0.0
        for off, k in off_k:
            j = idx + off
            wj = w_list[j]
            if wj > 0.0:
                rem = k - cov_list[j]
                if rem > 0.0:
                    gain += wj * rem
        if not heap or gain >= -heap[0][0] - 1e-12:
            y = idx // wp - pad
            x = idx % wp - pad
            if not dots[y, x]:
                dots[y, x] = True
                for off, k in off_k:
                    j = idx + off
                    if k > cov_list[j]:
                        cov_list[j] = k
                picked += 1
        elif gain > 0.0:
            push(heap, (-gain, idx))

    return dots, picked


def build_stage1_context_matrix(
    p_phys: np.ndarray,
    p_scarp: np.ndarray,
    p_geotherm: np.ndarray,
    p_subsurf: np.ndarray,
    fp: np.ndarray,
    valid_idx: np.ndarray,
) -> np.ndarray:
    """Compact `(len(valid_idx), 20)` `float16` Stage-1 axial-ridge + corroboration matrix."""
    out = np.empty((valid_idx.size, 20), dtype=np.float16)
    out[:, 0] = _logit(p_phys).ravel()[valid_idx].astype(np.float16)
    out[:, 1] = _logit(p_scarp).ravel()[valid_idx].astype(np.float16)
    out[:, 2] = _logit(p_geotherm).ravel()[valid_idx].astype(np.float16)
    out[:, 3] = _logit(p_subsurf).ravel()[valid_idx].astype(np.float16)

    fill_axial_context_valid(out, 4, p_phys, fp, valid_idx)
    blind_blend = (0.50 * p_geotherm + 0.50 * p_subsurf).astype(np.float32)
    fill_axial_context_valid(out, 10, blind_blend, fp, valid_idx)
    del blind_blend
    gc.collect()

    q_s = _rank01(p_scarp[fp])
    q_g = _rank01(p_geotherm[fp])
    q_u = _rank01(p_subsurf[fp])
    tri = np.sort(np.column_stack([q_s, q_g, q_u]), axis=1)
    out[:, 16] = tri[:, 1].astype(np.float16)
    out[:, 17] = np.cbrt(np.maximum(q_s * q_g * q_u, 0.0)).astype(np.float16)
    out[:, 18] = np.sqrt(np.maximum(q_g * q_u, 0.0)).astype(np.float16)
    out[:, 19] = np.sqrt(np.maximum(q_s * np.maximum(q_g, q_u), 0.0)).astype(np.float16)
    del q_s, q_g, q_u, tri
    gc.collect()
    return out


def build_meta_features_from_parts(s1_part: np.ndarray, mech_part: np.ndarray) -> np.ndarray:
    p_logit = s1_part[:, 0]
    ax_r9 = s1_part[:, 4]
    ax_r21 = s1_part[:, 5]
    blind_r21 = s1_part[:, 11]
    q_sec = s1_part[:, 16]

    ray_s = mech_part[:, 0]
    ray_l = mech_part[:, 1]
    bridge = mech_part[:, 3]
    d_cap = mech_part[:, 5]
    deficit = mech_part[:, 7]
    updip_rel = mech_part[:, 8]

    inter = np.column_stack(
        [
            np.maximum(ax_r9, 0.0) * ray_s,
            np.maximum(ax_r21, 0.0) * (ray_l + bridge),
            np.maximum(blind_r21, 0.0) * deficit,
            q_sec * np.exp(-d_cap / 14.0),
            q_sec * (bridge + 0.5 * ray_s),
            (p_logit + 4.0) * updip_rel,
        ]
    ).astype(np.float32)
    return np.hstack([s1_part, mech_part, inter]).astype(np.float32)


def _build_self_masked_training_dataset(
    stack: StaticStack,
    inv_map: np.ndarray,
    vis_train_domain: np.ndarray,
    full_vis_fold: np.ndarray,
    fp: np.ndarray,
    stage1_context_valid: np.ndarray,
    quad_map: np.ndarray | None = None,
    n_subfolds: int = 3,
    seed: int = 8000,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    vis_in_dom = full_vis_fold & vis_train_domain
    comp_lbl, n_comp = ndi.label(vis_in_dom, structure=np.ones((3, 3), dtype=int))
    comp_ids = np.arange(1, n_comp + 1)
    rng.shuffle(comp_ids)
    buckets = np.array_split(comp_ids, n_subfolds)

    slp_flat = stack.slope10_full.ravel()
    quad_flat = quad_map.ravel() if quad_map is not None else np.zeros(fp.size, dtype=np.int8)
    X_base_list = []
    X_meta_list = []
    y_list = []
    q_list = []

    for cid_part in buckets:
        if cid_part.size == 0:
            continue
        lut = np.zeros(n_comp + 1, dtype=bool)
        lut[cid_part] = True
        hide_mask = lut[comp_lbl]
        if not np.any(hide_mask):
            continue
        hide_buf = ndi.binary_dilation(hide_mask, iterations=3)
        sub_vis = full_vis_fold & (~hide_buf)
        d_sub = ndi.distance_transform_edt(~sub_vis)

        pos_mask = hide_mask & vis_train_domain & (d_sub > 2.15)
        pos_flat = np.flatnonzero(pos_mask)
        if pos_flat.size < 50:
            continue

        neg_mask = vis_train_domain & (~full_vis_fold) & (d_sub > 2.15)
        neg_pool = np.flatnonzero(neg_mask)
        neg_flat = sample_terrain_balanced_negatives(
            neg_pool, pos_flat, slp_flat, n_neg_total=3 * int(pos_flat.size), hard_frac=0.35, rng=rng
        )
        del d_sub

        mech_full = compute_fold_mechanics_2d(
            sub_vis, fp, stack.best_strike_full, stack.tect_bg_full, stack.updip_energy_full
        ).reshape(-1, len(FOLD_MECH_NAMES))
        samp_flat = np.concatenate([pos_flat, neg_flat])
        m_part = mech_full[samp_flat].copy()
        del mech_full
        gc.collect()

        y_sub = np.concatenate(
            [np.ones(pos_flat.size, dtype=np.int8), np.zeros(neg_flat.size, dtype=np.int8)]
        )
        v_rows = inv_map[samp_flat]
        b_part = np.hstack([np.asarray(stack.X_mmap[v_rows], dtype=np.float32), m_part])
        s1_part = np.asarray(stage1_context_valid[v_rows], dtype=np.float32)
        meta_part = build_meta_features_from_parts(s1_part, m_part)

        X_base_list.append(b_part)
        X_meta_list.append(meta_part)
        y_list.append(y_sub)
        q_list.append(quad_flat[samp_flat].astype(np.int8))

    return (
        np.vstack(X_base_list).astype(np.float32),
        np.vstack(X_meta_list).astype(np.float32),
        np.concatenate(y_list).astype(np.int8),
        np.concatenate(q_list).astype(np.int8),
    )


def fit_and_predict_stage1_oof_4quad(
    stack: StaticStack,
    inv_map: np.ndarray,
    quad_map: np.ndarray,
    full_vis_fold: np.ndarray,
    fp: np.ndarray,
    seed: int,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Fit Stage-1 zero-catalogue domain experts (`M_scarp`, `M_geotherm`, `M_updip_subsurf`, `M_phys_all`)
    across the 4 spatial quadrants strictly on `quad_map != q_idx` and predict out-of-fold on
    `quad_map == q_idx` (total predicted rows = 1.0x footprint across all 4 folds)."""
    H, W = fp.shape
    v_idx = stack.valid_idx
    quad_valid = quad_map.ravel()[v_idx]
    d_vis = ndi.distance_transform_edt(~full_vis_fold)

    p_scarp_v = np.zeros(v_idx.size, dtype=np.float32)
    p_geotherm_v = np.zeros(v_idx.size, dtype=np.float32)
    p_subsurf_v = np.zeros(v_idx.size, dtype=np.float32)
    p_all_v = np.zeros(v_idx.size, dtype=np.float32)

    scarp_cols = list(range(0, 18))
    geotherm_cols = list(range(18, 34)) + [11]
    subsurf_cols = list(range(34, 50)) + [11]

    for q_idx in range(4):
        rng = np.random.default_rng(seed + q_idx * 17)
        train_dom = fp & (quad_map != q_idx)
        pos_mask = full_vis_fold & train_dom
        pos_flat = np.flatnonzero(pos_mask)
        neg_pool = np.flatnonzero(train_dom & (~full_vis_fold) & (d_vis > 3.0))

        neg_flat = sample_terrain_balanced_negatives(
            neg_pool,
            pos_flat,
            stack.slope10_full.ravel(),
            n_neg_total=3 * int(pos_flat.size),
            hard_frac=0.40,
            rng=rng,
        )
        samp_flat = np.concatenate([pos_flat, neg_flat])
        y_tr = np.concatenate([np.ones(pos_flat.size, dtype=np.int8), np.zeros(neg_flat.size, dtype=np.int8)])
        v_rows = inv_map[samp_flat]
        X_tr = np.asarray(stack.X_mmap[v_rows], dtype=np.float32)

        m_scarp = fit_hgb_fast(X_tr[:, scarp_cols], y_tr, seed=seed + q_idx * 10 + 1, max_iter=105, max_leaf_nodes=39)
        m_geotherm = fit_hgb_fast(X_tr[:, geotherm_cols], y_tr, seed=seed + q_idx * 10 + 2, max_iter=100, max_leaf_nodes=31)
        m_subsurf = fit_hgb_fast(X_tr[:, subsurf_cols], y_tr, seed=seed + q_idx * 10 + 3, max_iter=105, max_leaf_nodes=39)
        m_all = fit_hgb_fast(X_tr, y_tr, seed=seed + q_idx * 10 + 4, max_iter=120, max_leaf_nodes=55, min_samples_leaf=90)
        del X_tr, y_tr
        gc.collect()

        q_rows = np.flatnonzero(quad_valid == q_idx)
        p_scarp_v[q_rows] = predict_rows_chunked(m_scarp, stack.X_mmap, q_rows, cols=scarp_cols)
        p_geotherm_v[q_rows] = predict_rows_chunked(m_geotherm, stack.X_mmap, q_rows, cols=geotherm_cols)
        p_subsurf_v[q_rows] = predict_rows_chunked(m_subsurf, stack.X_mmap, q_rows, cols=subsurf_cols)
        p_all_v[q_rows] = predict_rows_chunked(m_all, stack.X_mmap, q_rows)
        del m_scarp, m_geotherm, m_subsurf, m_all
        gc.collect()

    del d_vis
    p_scarp = self_reshape(p_scarp_v, fp, v_idx)
    p_geotherm = self_reshape(p_geotherm_v, fp, v_idx)
    p_subsurf = self_reshape(p_subsurf_v, fp, v_idx)
    p_all = self_reshape(p_all_v, fp, v_idx)
    del p_scarp_v, p_geotherm_v, p_subsurf_v, p_all_v
    gc.collect()

    p_phys_raw = (0.42 * p_all + 0.22 * p_scarp + 0.18 * p_geotherm + 0.18 * p_subsurf).astype(np.float32)
    return p_phys_raw, p_scarp, p_geotherm, p_subsurf


def predict_stage2_rows(
    stack: StaticStack,
    s1_ctx_valid: np.ndarray,
    mech_valid: np.ndarray,
    row_indices: np.ndarray,
    m_base1: HistGradientBoostingClassifier,
    m_base2: HistGradientBoostingClassifier,
    m_meta_hgb: HistGradientBoostingClassifier,
    m_meta_mlp,
    m_meta_lr,
) -> np.ndarray:
    """Predict Stage-2 self-masked fracture-transfer + updip multi-physics surface on `row_indices`."""
    n_rows = row_indices.size
    pred_out = np.zeros(n_rows, dtype=np.float32)
    chunk = 450_000
    for start in range(0, n_rows, chunk):
        end = min(n_rows, start + chunk)
        r_slice = row_indices[start:end]
        m_sub = mech_valid[r_slice]
        b_sub = np.hstack([np.asarray(stack.X_mmap[r_slice], dtype=np.float32), m_sub])
        pb = 0.50 * m_base1.predict_proba(b_sub)[:, 1] + 0.50 * m_base2.predict_proba(b_sub)[:, 1]
        del b_sub

        s1_sub = np.asarray(s1_ctx_valid[r_slice], dtype=np.float32)
        meta_sub = build_meta_features_from_parts(s1_sub, m_sub)
        pm_hgb = m_meta_hgb.predict_proba(meta_sub)[:, 1]
        meta_clean = np.nan_to_num(meta_sub, nan=0.0, posinf=0.0, neginf=0.0)
        pm_mlp = m_meta_mlp.predict_proba(meta_clean)[:, 1]
        pm_lr = m_meta_lr.predict_proba(meta_clean)[:, 1]
        del meta_sub, meta_clean

        pred_out[start:end] = (0.48 * pb + 0.32 * pm_hgb + 0.12 * pm_mlp + 0.08 * pm_lr).astype(
            np.float32
        )
    return pred_out


def run_holdout_oof_for_seed(
    stack: StaticStack,
    folds: List[dict],
    quad_map: np.ndarray,
    fp: np.ndarray,
    seed: int = 20,
) -> Dict[str, np.ndarray]:
    """Run out-of-fold predictions across the 4 spatial quadrants for one holdout draw (`seed`),
    returning `{oof_e1_phys, oof_e2_selfmask, oof_e3_hero}`."""
    t_s0 = time.time()
    H, W = fp.shape
    v_idx = stack.valid_idx
    inv_map = np.full(fp.size, -1, dtype=np.int64)
    inv_map[v_idx] = np.arange(v_idx.size, dtype=np.int64)

    # All 4 folds in `buffered_component_draw` share the identical global `visible` mask
    vis_fold = folds[0]["visible"] & fp
    d_vis_fold = ndi.distance_transform_edt(~vis_fold).astype(np.float32)
    w_blind = np.clip((d_vis_fold - 8.0) / 22.0, 0.0, 0.42).astype(np.float32)

    p_phys_raw, p_scarp, p_geotherm, p_subsurf = fit_and_predict_stage1_oof_4quad(
        stack, inv_map, quad_map, vis_fold, fp, seed=seed * 100
    )
    s1_ctx_valid = build_stage1_context_matrix(
        p_phys_raw, p_scarp, p_geotherm, p_subsurf, fp, v_idx
    )
    p_phys_axial_raw = (
        0.55 * p_phys_raw
        + 0.25 * (1.0 / (1.0 + np.exp(-self_reshape(s1_ctx_valid[:, 5], fp, v_idx))))
        + 0.20 * (1.0 / (1.0 + np.exp(-self_reshape(s1_ctx_valid[:, 11], fp, v_idx))))
    ).astype(np.float32)
    del p_phys_raw, p_scarp, p_geotherm, p_subsurf
    gc.collect()
    print(f"  [Holdout Seed {seed}] Stage-1 OOF + axial context complete in {time.time() - t_s0:.1f}s", flush=True)

    t_ds = time.time()
    X_b_all, X_m_all, y_m_all, q_row_all = _build_self_masked_training_dataset(
        stack,
        inv_map,
        fp,
        vis_fold,
        fp,
        s1_ctx_valid,
        quad_map=quad_map,
        n_subfolds=3,
        seed=seed * 1000,
    )
    mech_eval_valid = compute_fold_mechanics_2d(
        vis_fold, fp, stack.best_strike_full, stack.tect_bg_full, stack.updip_energy_full
    ).reshape(-1, len(FOLD_MECH_NAMES))[v_idx].copy()
    gc.collect()
    print(f"  [Holdout Seed {seed}] Self-masked subfold dataset ({X_b_all.shape[0]:,} rows) built in {time.time() - t_ds:.1f}s", flush=True)

    oof_e1_phys = np.zeros((H, W), dtype=np.float32)
    oof_e2_selfmask = np.zeros((H, W), dtype=np.float32)
    oof_e3_hero = np.zeros((H, W), dtype=np.float32)
    quad_valid = quad_map.ravel()[v_idx]

    for q_idx, f in enumerate(folds):
        t_f0 = time.time()
        q_name = f["name"]
        q_eval = f["region"] & fp
        train_dom = fp & (quad_map != q_idx)

        p_phys_cal = calibrate_terrain_prevalence(
            p_phys_axial_raw, stack.slope10_full, vis_fold & train_dom, train_dom, fp, blend=0.65
        )
        oof_e1_phys[q_eval] = p_phys_cal[q_eval]

        tr_mask = q_row_all != q_idx
        X_b_tr = X_b_all[tr_mask]
        X_m_tr = X_m_all[tr_mask]
        y_m_tr = y_m_all[tr_mask]

        m_base1 = fit_hgb_fast(X_b_tr, y_m_tr, seed=seed * 10 + q_idx + 1, max_iter=125, max_leaf_nodes=47)
        m_base2 = fit_hgb_fast(
            X_b_tr, y_m_tr, seed=seed * 10 + q_idx + 2, max_iter=135, max_leaf_nodes=63, min_samples_leaf=80
        )
        del X_b_tr

        m_meta_hgb = fit_hgb_fast(X_m_tr, y_m_tr, seed=seed * 10 + q_idx + 3, max_iter=110, max_leaf_nodes=31)
        X_m_clean = np.nan_to_num(X_m_tr, nan=0.0, posinf=0.0, neginf=0.0)
        del X_m_tr
        m_meta_mlp = make_pipeline(
            StandardScaler(),
            MLPClassifier(
                hidden_layer_sizes=(48, 24),
                max_iter=30,
                alpha=1.5,
                early_stopping=True,
                random_state=seed * 10 + q_idx + 4,
            ),
        )
        m_meta_mlp.fit(X_m_clean, y_m_tr)
        m_meta_lr = make_pipeline(StandardScaler(), LogisticRegression(C=0.2, max_iter=120))
        m_meta_lr.fit(X_m_clean, y_m_tr)
        del X_m_clean, y_m_tr

        # Predict on q_eval plus a 10% stride of train_dom for quantile calibration
        eval_and_cal_rows = np.flatnonzero((quad_valid == q_idx) | ((np.arange(v_idx.size) % 10 == 0) & (quad_valid != q_idx)))
        pred_rows = predict_stage2_rows(
            stack, s1_ctx_valid, mech_eval_valid, eval_and_cal_rows, m_base1, m_base2, m_meta_hgb, m_meta_mlp, m_meta_lr
        )
        del m_base1, m_base2, m_meta_hgb, m_meta_mlp, m_meta_lr

        p_meta_fold_v = np.zeros(v_idx.size, dtype=np.float32)
        p_meta_fold_v[eval_and_cal_rows] = pred_rows
        p_meta_fold = self_reshape(p_meta_fold_v, fp, v_idx)
        cal_dom = np.zeros(fp.size, dtype=bool)
        cal_dom[v_idx[eval_and_cal_rows]] = True
        cal_dom = cal_dom.reshape(H, W) & train_dom

        p_meta_cal = calibrate_terrain_prevalence(
            p_meta_fold, stack.slope10_full, vis_fold & cal_dom, cal_dom, fp, blend=0.55
        )
        oof_e2_selfmask[q_eval] = p_meta_cal[q_eval]

        p_hero = ((1.0 - w_blind) * p_meta_cal + w_blind * p_phys_cal).astype(np.float32)
        p_hero = calibrate_terrain_prevalence(
            p_hero, stack.slope10_full, vis_fold & cal_dom, cal_dom, fp, blend=0.60
        )
        oof_e3_hero[q_eval] = p_hero[q_eval]
        print(f"  [Holdout Seed {seed} Fold {q_name}] Stage-2 OOF complete in {time.time() - t_f0:.1f}s", flush=True)

    del X_b_all, X_m_all, y_m_all, q_row_all, mech_eval_valid, s1_ctx_valid
    gc.collect()
    return {
        "oof_e1_phys": oof_e1_phys,
        "oof_e2_selfmask": oof_e2_selfmask,
        "oof_e3_hero": oof_e3_hero,
    }


def build_full_catalogue_surfaces(
    stack: StaticStack,
    full_cat: np.ndarray,
    fp: np.ndarray,
    seed: int = 8080,
) -> Tuple[np.ndarray, np.ndarray]:
    """Fit full-catalogue Stage-1 & Stage-2 models for the shipped GeoTIFF submission."""
    t0 = time.time()
    v_idx = stack.valid_idx
    inv_map = np.full(fp.size, -1, dtype=np.int64)
    inv_map[v_idx] = np.arange(v_idx.size, dtype=np.int64)

    yy, xx = np.nonzero(fp)
    ym, xm = int(np.median(yy)), int(np.median(xx))
    gy, gx = np.ogrid[: fp.shape[0], : fp.shape[1]]
    quad_full = np.full(fp.shape, -1, dtype=np.int8)
    quad_full[(gy < ym) & (gx < xm) & fp] = 0
    quad_full[(gy < ym) & (gx >= xm) & fp] = 1
    quad_full[(gy >= ym) & (gx < xm) & fp] = 2
    quad_full[(gy >= ym) & (gx >= xm) & fp] = 3

    p_phys_raw, p_scarp, p_geotherm, p_subsurf = fit_and_predict_stage1_oof_4quad(
        stack, inv_map, quad_full, full_cat & fp, fp, seed=seed
    )
    s1_ctx_valid = build_stage1_context_matrix(
        p_phys_raw, p_scarp, p_geotherm, p_subsurf, fp, v_idx
    )
    p_phys_axial = (
        0.55 * p_phys_raw
        + 0.25 * (1.0 / (1.0 + np.exp(-self_reshape(s1_ctx_valid[:, 5], fp, v_idx))))
        + 0.20 * (1.0 / (1.0 + np.exp(-self_reshape(s1_ctx_valid[:, 11], fp, v_idx))))
    ).astype(np.float32)
    p_phys_axial = calibrate_terrain_prevalence(
        p_phys_axial, stack.slope10_full, full_cat & fp, fp, fp, blend=0.65
    )
    del p_phys_raw, p_scarp, p_geotherm, p_subsurf
    gc.collect()

    X_b_tr, X_m_tr, y_m_tr, _ = _build_self_masked_training_dataset(
        stack, inv_map, fp, full_cat & fp, fp, s1_ctx_valid, quad_map=quad_full, n_subfolds=4, seed=seed + 500
    )
    m_base1 = fit_hgb_fast(X_b_tr, y_m_tr, seed=seed + 11, max_iter=135, max_leaf_nodes=47)
    m_base2 = fit_hgb_fast(X_b_tr, y_m_tr, seed=seed + 12, max_iter=145, max_leaf_nodes=63, min_samples_leaf=80)
    del X_b_tr
    gc.collect()

    m_meta_hgb = fit_hgb_fast(X_m_tr, y_m_tr, seed=seed + 13, max_iter=115, max_leaf_nodes=31)
    X_m_clean = np.nan_to_num(X_m_tr, nan=0.0, posinf=0.0, neginf=0.0)
    del X_m_tr
    m_meta_mlp = make_pipeline(
        StandardScaler(),
        MLPClassifier(
            hidden_layer_sizes=(48, 24),
            max_iter=30,
            alpha=1.5,
            early_stopping=True,
            random_state=seed + 14,
        ),
    )
    m_meta_mlp.fit(X_m_clean, y_m_tr)
    m_meta_lr = make_pipeline(StandardScaler(), LogisticRegression(C=0.2, max_iter=120))
    m_meta_lr.fit(X_m_clean, y_m_tr)
    del X_m_clean, y_m_tr
    gc.collect()

    mech_full_valid = compute_fold_mechanics_2d(
        full_cat & fp, fp, stack.best_strike_full, stack.tect_bg_full, stack.updip_energy_full
    ).reshape(-1, len(FOLD_MECH_NAMES))[v_idx].copy()
    all_rows = np.arange(v_idx.size, dtype=np.int64)
    pred_full_v = predict_stage2_rows(
        stack, s1_ctx_valid, mech_full_valid, all_rows, m_base1, m_base2, m_meta_hgb, m_meta_mlp, m_meta_lr
    )
    del m_base1, m_base2, m_meta_hgb, m_meta_mlp, m_meta_lr, s1_ctx_valid, mech_full_valid
    gc.collect()

    p_meta_full = self_reshape(pred_full_v, fp, v_idx)
    p_meta_cal = calibrate_terrain_prevalence(
        p_meta_full, stack.slope10_full, full_cat & fp, fp, fp, blend=0.55
    )
    print(f"  [Full-Catalogue Fit] Complete in {time.time() - t0:.1f}s", flush=True)
    return p_meta_cal, p_phys_axial


def self_reshape(vec_valid: np.ndarray, fp: np.ndarray, valid_idx: np.ndarray) -> np.ndarray:
    out = np.zeros(fp.size, dtype=np.float32)
    out[valid_idx] = np.asarray(vec_valid, dtype=np.float32)
    return out.reshape(fp.shape)
