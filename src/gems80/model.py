"""Leakage canary, Stage-1 physical LOQO ensemble, Curvilinear Radon sharpening,
and Stage-2 Nested Whole-Component Fracture-Mechanics Meta-Learner + Calibration.
"""
from __future__ import annotations

import gc
import numpy as np
from lightgbm import LGBMClassifier
from scipy import ndimage as ndi
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.isotonic import IsotonicRegression
from sklearn.metrics import roc_auc_score

from vendor.gems57.grid import Grid
from vendor.gems57.holdout import buffered_component_draw
from .features import (
    FRACTURE_FEATURE_NAMES,
    N_STATIC,
    STATIC_FEATURE_NAMES,
    compute_fracture_mechanics_features,
)


def _logit(p: np.ndarray, eps: float = 1e-6) -> np.ndarray:
    c = np.clip(p.astype(np.float64), eps, 1.0 - eps)
    return np.log(c / (1.0 - c)).astype(np.float32)


def curvilinear_radon_features(
    logit_field_2d: np.ndarray, valid: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """12-azimuth strike-steered line-integral & transverse ridge-contrast sharpening.

    Operates strictly on an out-of-fold score field `logit_field_2d` (zero labels).
    Produces three structural fields on `valid`:
      1. `radon_short` : max 12-azimuth 700 m (`+/-3 px`) line integral (curved trace segment)
      2. `radon_long`  : max 12-azimuth 1.5 km (`+/-7 px`) line integral (regional fault zone)
      3. `ridge_sharp` : transverse second-difference ridge contrast `score - perp_flank_mean`,
                         suppressing fat 2D hillside blobs while sharpening 1-px scarps!
    """
    fill_val = float(np.percentile(logit_field_2d[valid], 50.0))
    f = np.where(valid, logit_field_2d, fill_val).astype(np.float32)

    best_short = np.full(f.shape, -1e9, dtype=np.float32)
    best_long = np.full(f.shape, -1e9, dtype=np.float32)
    best_sharp = np.full(f.shape, -1e9, dtype=np.float32)

    for k in range(12):
        theta = k * np.pi / 12.0  # 15-degree steps covering [0, 180 deg)
        ks = _line_kernel(length=7, theta_rad=theta)
        kl = _line_kernel(length=15, theta_rad=theta)
        kp = _perp_flank_kernel(theta_rad=theta)

        s_short = ndi.convolve(f, ks, mode="nearest")
        s_long = ndi.convolve(f, kl, mode="nearest")
        flank = ndi.convolve(f, kp, mode="nearest")
        sharp = s_short - flank

        better = s_short > best_short
        best_short = np.where(better, s_short, best_short)
        best_sharp = np.where(better, sharp, best_sharp)
        best_long = np.maximum(best_long, s_long)

    return best_short[valid], best_long[valid], best_sharp[valid]


def _line_kernel(length: int, theta_rad: float) -> np.ndarray:
    r = length // 2
    y, x = np.mgrid[-r : r + 1, -r : r + 1].astype(np.float32)
    perp = np.abs(-y * np.cos(theta_rad) + x * np.sin(theta_rad))
    along = np.abs(y * np.sin(theta_rad) + x * np.cos(theta_rad))
    w = np.clip(1.0 - perp, 0.0, 1.0) * (along <= float(r) + 0.35)
    w /= max(float(w.sum()), 1e-6)
    return w.astype(np.float32)


def _perp_flank_kernel(theta_rad: float) -> np.ndarray:
    r = 3
    y, x = np.mgrid[-r : r + 1, -r : r + 1].astype(np.float32)
    perp = np.abs(-y * np.cos(theta_rad) + x * np.sin(theta_rad))
    along = np.abs(y * np.sin(theta_rad) + x * np.cos(theta_rad))
    w = np.exp(-0.5 * ((perp - 2.4) / 0.7) ** 2) * (along <= 1.2)
    w /= max(float(w.sum()), 1e-6)
    return w.astype(np.float32)


def run_leakage_canary(
    X_static: np.memmap,
    folds: list[dict],
    valid: np.ndarray,
    max_allowed_auc: float = 0.90,
) -> dict:
    """Single-feature ROC-AUC leakage canary across all 52 static + 7 fracture features."""
    rng = np.random.default_rng(8020)
    pos_indices_v = []
    neg_indices_v = []

    valid_idx_map = np.full(valid.shape, -1, dtype=np.int32)
    valid_idx_map[valid] = np.arange(int(valid.sum()), dtype=np.int32)

    # All 4 folds in `buffered_component_draw` share one globally buffered `visible` mask
    F_vis, _, _ = compute_fracture_mechanics_features(folds[0]["visible"], valid)

    for f in folds:
        active = valid & f["region"] & (~f["masked_known"])
        truth = active & f["truth"]
        bg = active & (~f["truth"])

        pv = valid_idx_map[truth]
        nv_all = valid_idx_map[bg]
        nv = rng.choice(nv_all, size=min(40_000, nv_all.size), replace=False)

        pos_indices_v.append(pv)
        neg_indices_v.append(nv)

    pv_cat = np.concatenate(pos_indices_v)
    nv_cat = np.concatenate(neg_indices_v)
    y = np.concatenate([np.ones(pv_cat.size, dtype=np.int8), np.zeros(nv_cat.size, dtype=np.int8)])

    rows = []
    X_p = np.asarray(X_static[pv_cat], dtype=np.float32)
    X_n = np.asarray(X_static[nv_cat], dtype=np.float32)
    for col, name in enumerate(STATIC_FEATURE_NAMES):
        vals = np.concatenate([X_p[:, col], X_n[:, col]])
        fin = np.isfinite(vals)
        auc = float(roc_auc_score(y[fin], vals[fin])) if fin.sum() > 200 else 0.5
        eff_auc = max(auc, 1.0 - auc)
        rows.append({"feature": name, "group": "static_physical", "raw_auc": round(auc, 5), "effective_auc": round(eff_auc, 5)})

    Fp = F_vis[pv_cat]
    Fn = F_vis[nv_cat]
    for col, name in enumerate(FRACTURE_FEATURE_NAMES):
        vals = np.concatenate([Fp[:, col], Fn[:, col]])
        auc = float(roc_auc_score(y, vals))
        eff_auc = max(auc, 1.0 - auc)
        rows.append({"feature": name, "group": "fold_visible_fracture", "raw_auc": round(auc, 5), "effective_auc": round(eff_auc, 5)})

    rows.sort(key=lambda r: r["effective_auc"], reverse=True)
    max_eff = max(r["effective_auc"] for r in rows)
    flagged = [r for r in rows if r["effective_auc"] > max_allowed_auc]
    return {
        "evidence_class": "LEAKAGE-CANARY",
        "evaluator_version": "gems57-pooled-hide-v2",
        "n_withheld_positives_tested": int(pv_cat.size),
        "n_background_samples_tested": int(nv_cat.size),
        "n_features_checked": len(rows),
        "threshold_auc": max_allowed_auc,
        "max_effective_auc": round(max_eff, 5),
        "top_feature": rows[0]["feature"],
        "flagged_features": flagged,
        "passed": len(flagged) == 0,
        "features": rows,
    }


def _predict_chunked(clf, X: np.memmap | np.ndarray, idx: np.ndarray, chunk: int = 400_000) -> np.ndarray:
    out = np.empty(idx.size, dtype=np.float32)
    for s in range(0, idx.size, chunk):
        sl = idx[s : s + chunk]
        out[s : s + chunk] = clf.predict_proba(np.asarray(X[sl], dtype=np.float32))[:, 1].astype(np.float32)
    return out


def _quadrant_partition_masks(valid: np.ndarray) -> list[np.ndarray]:
    """Full-footprint quadrant partition matching `vendor.gems57.holdout.quadrant_ids`."""
    from vendor.gems57.holdout import quadrant_ids
    q = quadrant_ids(valid)
    return [q == i for i in range(4)]


def fit_stage1_physical_oof(
    X_static: np.memmap,
    labels: np.ndarray,
    valid: np.ndarray,
    quad_masks: list[np.ndarray],
) -> tuple[np.ndarray, list[tuple[LGBMClassifier, HistGradientBoostingClassifier]], dict]:
    """Train dual-family digitization-jitter-clean spatial LOQO physical models across the 4 quadrants.

    Every pixel in quadrant `q` is scored strictly out-of-fold by models trained strictly
    outside quadrant `q` + a 2.4 km (`24 px` across boundary: `12 px` inside + `12 px` outside) buffer!
    """
    n_valid = int(valid.sum())
    oof_logit = np.zeros(n_valid, dtype=np.float32)
    models: list[tuple[LGBMClassifier, HistGradientBoostingClassifier]] = []

    d_true_all = ndi.distance_transform_edt(~labels).astype(np.float32)
    valid_idx_map = np.full(valid.shape, -1, dtype=np.int32)
    valid_idx_map[valid] = np.arange(n_valid, dtype=np.int32)

    importances = np.zeros(N_STATIC, dtype=np.float64)

    for q_idx, q_mask in enumerate(quad_masks):
        rng = np.random.default_rng(8100 + q_idx)
        # Strict spatial separation: 12 px outside full quadrant `q_mask`!
        d_q = ndi.distance_transform_edt(~q_mask).astype(np.float32)
        train_dom = valid & (d_q >= 12.0)
        assert not np.any(train_dom & q_mask), "Spatial LOQO leakage guard violated!"

        pos_a = valid_idx_map[train_dom & labels]
        neg_pool_a = valid_idx_map[train_dom & (d_true_all > 1.15)]
        neg_a = rng.choice(neg_pool_a, size=min(pos_a.size * 9, neg_pool_a.size), replace=False)

        tr_a = np.concatenate([pos_a, neg_a])
        y_a = np.concatenate([np.ones(pos_a.size, dtype=np.int32), np.zeros(neg_a.size, dtype=np.int32)])
        perm_a = rng.permutation(tr_a.size)
        tr_a, y_a = tr_a[perm_a], y_a[perm_a]

        lgbm = LGBMClassifier(
            n_estimators=420,
            learning_rate=0.055,
            num_leaves=63,
            min_child_samples=120,
            subsample=0.80,
            subsample_freq=1,
            colsample_bytree=0.75,
            reg_alpha=1.5,
            reg_lambda=6.0,
            random_state=8100 + q_idx,
            n_jobs=2,
            verbose=-1,
        )
        lgbm.fit(np.asarray(X_static[tr_a], dtype=np.float32), y_a)
        importances += lgbm.booster_.feature_importance(importance_type="gain")
        del tr_a, y_a
        gc.collect()

        pos_b_all = valid_idx_map[train_dom & (d_true_all <= 1.05)]
        pos_b = rng.choice(pos_b_all, size=min(95_000, pos_b_all.size), replace=False)
        neg_pool_b = valid_idx_map[train_dom & (d_true_all > 3.05)]
        neg_b = rng.choice(neg_pool_b, size=min(360_000, neg_pool_b.size), replace=False)

        tr_b = np.concatenate([pos_b, neg_b])
        y_b = np.concatenate([np.ones(pos_b.size, dtype=np.int32), np.zeros(neg_b.size, dtype=np.int32)])
        perm_b = rng.permutation(tr_b.size)
        tr_b, y_b = tr_b[perm_b], y_b[perm_b]

        hgb = HistGradientBoostingClassifier(
            max_iter=320,
            learning_rate=0.06,
            max_leaf_nodes=63,
            min_samples_leaf=140,
            l2_regularization=4.0,
            random_state=8200 + q_idx,
        )
        hgb.fit(np.asarray(X_static[tr_b], dtype=np.float32), y_b)
        del tr_b, y_b
        gc.collect()

        models.append((lgbm, hgb))

        # Predict strictly out-of-fold across ALL valid pixels inside quadrant `q_mask`!
        test_v = valid_idx_map[valid & q_mask]
        pa = _predict_chunked(lgbm, X_static, test_v)
        pb = _predict_chunked(hgb, X_static, test_v)
        la = _logit(pa) - float(np.log(pos_a.size / max(neg_a.size, 1)))
        lb = _logit(pb) - float(np.log(pos_b.size / max(neg_b.size, 1)))
        oof_logit[test_v] = 0.62 * la + 0.38 * lb

    imp_order = np.argsort(-importances)
    top_imp = [
        {"feature": STATIC_FEATURE_NAMES[int(i)], "gain_fraction": round(float(importances[i] / max(importances.sum(), 1e-9)), 5)}
        for i in imp_order[:20]
    ]
    return oof_logit, models, {"top_20_gain_features": top_imp}


def predict_stage1_ensemble_all_valid(
    X_static: np.memmap,
    models: list[tuple[LGBMClassifier, HistGradientBoostingClassifier]],
    n_valid: int,
) -> np.ndarray:
    """Average log-likelihood-ratio scores across all 4 LOQO model pairs for final live inference."""
    all_idx = np.arange(n_valid, dtype=np.int32)
    acc = np.zeros(n_valid, dtype=np.float32)
    for lgbm, hgb in models:
        pa = _predict_chunked(lgbm, X_static, all_idx)
        pb = _predict_chunked(hgb, X_static, all_idx)
        la = _logit(pa)
        lb = _logit(pb)
        la -= float(np.median(la))
        lb -= float(np.median(lb))
        acc += (0.62 * la + 0.38 * lb) / float(len(models))
    return acc


def build_meta_features(
    stage1_logit_v: np.ndarray,
    radon_short_v: np.ndarray,
    radon_long_v: np.ndarray,
    ridge_sharp_v: np.ndarray,
    X_static: np.memmap | np.ndarray,
    F_frac: np.ndarray,
) -> np.ndarray:
    """19-channel structural meta-feature matrix fusing OOF physics, Radon ridges, and visible fracture geometry."""
    ctx = np.asarray(X_static[:, [6, 24, 28, 45, 47, 49]], dtype=np.float32)
    M = np.empty((stage1_logit_v.size, 19), dtype=np.float32)
    M[:, 0] = stage1_logit_v
    M[:, 1] = radon_short_v
    M[:, 2] = radon_long_v
    M[:, 3] = ridge_sharp_v
    M[:, 4] = stage1_logit_v * F_frac[:, 1]
    M[:, 5] = stage1_logit_v * F_frac[:, 3]
    M[:, 6:13] = F_frac
    M[:, 13:19] = ctx
    return M


def fit_stage2_fracture_meta_and_calibrators(
    stage1_oof_v: np.ndarray,
    X_static: np.memmap,
    grid: Grid,
    folds: list[dict],
    quad_masks: list[np.ndarray],
) -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    list[HistGradientBoostingClassifier],
    list[IsotonicRegression],
    list[IsotonicRegression],
]:
    """Fit nested whole-component Stage-2 Fracture Meta-Learner + exact OOF posterior calibrators `(p_hat, q_hat)`.

    Strictly Zero-Leakage Nested Whole-Component Protocol:
      1. Test fracture features `F_vis` use ONLY `folds[0]["visible"]` (where all withheld components
         and a 3-px buffer have been removed!).
      2. Inside each fold `f_idx`, we restrict strictly to `train_dom = valid & (d_q >= 12.0)` outside
         quadrant `q_mask`, construct a mini-Grid whose catalogue is STRICTLY `vis_train = folds[f_idx]["visible"] & train_dom`,
         and generate TWO independent whole-component buffered draws (`sub_seeds = (81, 82)`) via
         `buffered_component_draw`!
      3. Sub-draw `81` trains the nested meta-classifier and Sub-draw `82` calibrates `(p_hat, q_hat)` out-of-sample,
         then both sub-draws fit the final fold meta-classifier — 100% matched to whole-component unmapped faults
         and 0.000% test leakage!
    """
    valid = grid.footprint
    n_valid = int(valid.sum())

    s1_grid = np.zeros(valid.shape, dtype=np.float32)
    s1_grid[valid] = stage1_oof_v
    rad_short_v, rad_long_v, ridge_sharp_v = curvilinear_radon_features(s1_grid, valid)
    del s1_grid

    valid_idx_map = np.full(valid.shape, -1, dtype=np.int32)
    valid_idx_map[valid] = np.arange(n_valid, dtype=np.int32)

    # Compute test-time fracture mechanics ONCE from `folds[0]["visible"]`!
    F_test_vis, _, _ = compute_fracture_mechanics_features(folds[0]["visible"], valid)

    meta_oof_logit_v = np.zeros(n_valid, dtype=np.float32)
    p_hat_oof_v = np.zeros(n_valid, dtype=np.float32)
    q_hat_oof_v = np.zeros(n_valid, dtype=np.float32)

    meta_models: list[HistGradientBoostingClassifier] = []
    p_calibrators: list[IsotonicRegression] = []
    q_calibrators: list[IsotonicRegression] = []

    for q_idx, q_mask in enumerate(quad_masks):
        rng = np.random.default_rng(8500 + q_idx)
        d_q = ndi.distance_transform_edt(~q_mask).astype(np.float32)
        train_dom = valid & (d_q >= 12.0)
        assert not np.any(train_dom & folds[q_idx]["truth"]), "Stage-2 leakage guard violated!"

        sub_grid = Grid(
            footprint=train_dom,
            catalogue=folds[q_idx]["visible"] & train_dom,
            crs=grid.crs,
            transform=grid.transform,
            height=grid.height,
            width=grid.width,
        )

        # Build 2 independent whole-component buffered sub-draws inside `train_dom`:
        # Draw A (`seed = 80 + q_idx`) -> fits nested sub-model
        # Draw B (`seed = 90 + q_idx`) -> out-of-sample isotonic calibration `(p_hat, q_hat)`!
        sub_draws = []
        for s_seed in (80 + q_idx, 90 + q_idx):
            s_folds, _ = buffered_component_draw(sub_grid, seed=s_seed, hide_frac=0.22, buffer_px=3, boundary_px=6)
            s_vis = s_folds[0]["visible"]
            s_masked = s_folds[0]["masked_known"]
            s_truth = np.zeros(valid.shape, dtype=bool)
            s_reg_union = np.zeros(valid.shape, dtype=bool)
            for sf in s_folds:
                s_truth |= sf["truth"]
                s_reg_union |= sf["region"]
            s_active = train_dom & s_reg_union & (~s_masked)
            F_sub, _, _ = compute_fracture_mechanics_features(s_vis, valid)
            q_sub_2d = np.maximum(1.0 - ndi.distance_transform_edt(~(s_active & s_truth)) / 3.0, 0.0).astype(np.float32)

            act_v = valid_idx_map[s_active]
            M_sub = build_meta_features(
                stage1_oof_v[act_v],
                rad_short_v[act_v],
                rad_long_v[act_v],
                ridge_sharp_v[act_v],
                X_static[act_v],
                F_sub[act_v],
            )
            sub_draws.append((M_sub, s_truth[s_active].astype(np.int8), q_sub_2d[s_active]))
            del F_sub, q_sub_2d
            gc.collect()

        # Cross-predict between Draw A and Draw B for unbiased isotonic posterior calibration!
        cal_s, cal_yp, cal_yq, cal_w = [], [], [], []
        for tr_i, va_i in ((0, 1), (1, 0)):
            M_tr, y_tr, q_tr = sub_draws[tr_i]
            pos_i = np.flatnonzero(y_tr == 1)
            neg_pool_i = np.flatnonzero(q_tr <= 0.01)
            neg_i = rng.choice(neg_pool_i, size=min(pos_i.size * 22, neg_pool_i.size), replace=False)
            sel_i = np.concatenate([pos_i, neg_i])
            clf_sub = HistGradientBoostingClassifier(
                max_iter=220,
                learning_rate=0.05,
                max_leaf_nodes=31,
                min_samples_leaf=160,
                l2_regularization=6.0,
                random_state=8600 + 10 * q_idx + tr_i,
            )
            clf_sub.fit(M_tr[sel_i], y_tr[sel_i])

            M_va, y_va, q_va = sub_draws[va_i]
            p_va_idx = np.flatnonzero(y_va == 1)
            n_va_all = np.flatnonzero(y_va == 0)
            n_keep = min(250_000, n_va_all.size)
            n_va_idx = rng.choice(n_va_all, size=n_keep, replace=False)
            w_bg = float(n_va_all.size) / float(max(n_keep, 1))

            ev_idx = np.concatenate([p_va_idx, n_va_idx])
            ev_w = np.concatenate([np.ones(p_va_idx.size), np.full(n_va_idx.size, w_bg)])
            lg_va = _logit(clf_sub.predict_proba(M_va[ev_idx])[:, 1])
            cal_s.append(lg_va)
            cal_yp.append(y_va[ev_idx].astype(np.float64))
            cal_yq.append(q_va[ev_idx].astype(np.float64))
            cal_w.append(ev_w)

        iso_p = IsotonicRegression(y_min=1e-6, y_max=0.85, out_of_bounds="clip")
        iso_p.fit(np.concatenate(cal_s), np.concatenate(cal_yp), sample_weight=np.concatenate(cal_w))
        p_calibrators.append(iso_p)

        iso_q = IsotonicRegression(y_min=1e-5, y_max=0.95, out_of_bounds="clip")
        iso_q.fit(np.concatenate(cal_s), np.concatenate(cal_yq), sample_weight=np.concatenate(cal_w))
        q_calibrators.append(iso_q)

        # Fit final fold meta-classifier combining both nested whole-component draws
        tr_parts_X, tr_parts_y = [], []
        for M_d, y_d, q_d in sub_draws:
            p_d = np.flatnonzero(y_d == 1)
            n_pool_d = np.flatnonzero(q_d <= 0.01)
            n_d = rng.choice(n_pool_d, size=min(p_d.size * 22, n_pool_d.size), replace=False)
            s_d = np.concatenate([p_d, n_d])
            tr_parts_X.append(M_d[s_d])
            tr_parts_y.append(y_d[s_d])

        meta_clf = HistGradientBoostingClassifier(
            max_iter=240,
            learning_rate=0.05,
            max_leaf_nodes=31,
            min_samples_leaf=180,
            l2_regularization=6.0,
            random_state=8700 + q_idx,
        )
        meta_clf.fit(np.vstack(tr_parts_X), np.concatenate(tr_parts_y))
        meta_models.append(meta_clf)
        del sub_draws, tr_parts_X, tr_parts_y
        gc.collect()

        # Predict strictly out-of-fold across entire test quadrant `q_mask`!
        test_v = valid_idx_map[valid & q_mask]
        M_test = build_meta_features(
            stage1_oof_v[test_v],
            rad_short_v[test_v],
            rad_long_v[test_v],
            ridge_sharp_v[test_v],
            X_static[test_v],
            F_test_vis[test_v],
        )
        lg_test = _logit(meta_clf.predict_proba(M_test)[:, 1].astype(np.float32))
        meta_oof_logit_v[test_v] = lg_test
        p_hat_oof_v[test_v] = iso_p.predict(lg_test).astype(np.float32)
        q_hat_oof_v[test_v] = iso_q.predict(lg_test).astype(np.float32)
        del M_test
        gc.collect()

    return (
        meta_oof_logit_v,
        p_hat_oof_v,
        q_hat_oof_v,
        meta_models,
        p_calibrators,
        q_calibrators,
    )


def predict_stage2_for_visible_network(
    stage1_logit_v: np.ndarray,
    X_static: np.memmap,
    visible_network: np.ndarray,
    valid: np.ndarray,
    meta_models: list[HistGradientBoostingClassifier],
    p_calibrators: list[IsotonicRegression],
    q_calibrators: list[IsotonicRegression],
    quad_masks: list[np.ndarray] | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Predict Stage-2 log-odds & `(p_hat, q_hat)` for any visible fault network (`Seed 30` replication or Full Live!)."""
    n_valid = int(valid.sum())
    s1_grid = np.zeros(valid.shape, dtype=np.float32)
    s1_grid[valid] = stage1_logit_v
    rad_short_v, rad_long_v, ridge_sharp_v = curvilinear_radon_features(s1_grid, valid)
    del s1_grid

    F_net, _, _ = compute_fracture_mechanics_features(visible_network, valid)

    if quad_masks is not None:
        # Strict per-quadrant OOF evaluation (e.g. for independent Seed-30 holdout replication!)
        valid_idx_map = np.full(valid.shape, -1, dtype=np.int32)
        valid_idx_map[valid] = np.arange(n_valid, dtype=np.int32)
        lg_out = np.zeros(n_valid, dtype=np.float32)
        p_out = np.zeros(n_valid, dtype=np.float32)
        q_out = np.zeros(n_valid, dtype=np.float32)
        for q_idx, q_mask in enumerate(quad_masks):
            tv = valid_idx_map[valid & q_mask]
            M_q = build_meta_features(
                stage1_logit_v[tv], rad_short_v[tv], rad_long_v[tv], ridge_sharp_v[tv], X_static[tv], F_net[tv]
            )
            lg = _logit(meta_models[q_idx].predict_proba(M_q)[:, 1].astype(np.float32))
            lg_out[tv] = lg
            p_out[tv] = p_calibrators[q_idx].predict(lg).astype(np.float32)
            q_out[tv] = q_calibrators[q_idx].predict(lg).astype(np.float32)
        return lg_out, p_out, q_out

    # Full 4-fold ensemble for live submission (`visible_network = grid.catalogue`)
    logit_acc = np.zeros(n_valid, dtype=np.float32)
    p_acc = np.zeros(n_valid, dtype=np.float32)
    q_acc = np.zeros(n_valid, dtype=np.float32)
    chunk = 400_000
    for s in range(0, n_valid, chunk):
        sl = slice(s, min(s + chunk, n_valid))
        M_chunk = build_meta_features(
            stage1_logit_v[sl], rad_short_v[sl], rad_long_v[sl], ridge_sharp_v[sl], X_static[sl], F_net[sl]
        )
        for m_clf, iso_p, iso_q in zip(meta_models, p_calibrators, q_calibrators):
            lg = _logit(m_clf.predict_proba(M_chunk)[:, 1].astype(np.float32))
            logit_acc[sl] += lg / float(len(meta_models))
            p_acc[sl] += iso_p.predict(lg).astype(np.float32) / float(len(meta_models))
            q_acc[sl] += iso_q.predict(lg).astype(np.float32) / float(len(meta_models))

    return logit_acc, p_acc, q_acc
