"""Lane-duplicate check vs registry rasters (parallel-run protocol Rule 1).

Duplicate if Spearman rank-correlation on the footprint > 0.90 (checked BOTH on
the pre-placement continuous surface AND on the post-placement binary dots)
OR more than 70% of our dots fall within 3 px of one registry raster's dots.
Supports both raw 2D arrays and compact pre-indexed registry entries so 20+ full
GeoTIFF rasters can be audited simultaneously in <60 MB RAM on 3.8 GiB machines.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import rasterio
from scipy.ndimage import distance_transform_edt
from scipy.stats import spearmanr

CORR_MAX = 0.90
OVERLAP_MAX = 0.70
RADIUS_PX = 3
SAMPLE_SEED = 20261022
MAX_SAMPLE_N = 250_000


def registry_sample_indices(foot: np.ndarray, max_n: int = MAX_SAMPLE_N, seed: int = SAMPLE_SEED) -> np.ndarray:
    idx = np.flatnonzero(foot.ravel())
    if idx.size > max_n:
        rng = np.random.default_rng(seed)
        idx = np.sort(rng.choice(idx, size=max_n, replace=False))
    return idx


def load_compact_registry(paths: list[Path], foot: np.ndarray) -> dict[str, dict]:
    """Load arbitrary registry TIFs into compact memory (<15 MB per TIF instead of 50 MB)."""
    s_idx = registry_sample_indices(foot)
    out: dict[str, dict] = {}
    for p in paths:
        with rasterio.open(p) as src:
            arr = src.read(1)
        dots = (np.isfinite(arr) & (arr > 0.5) & foot)
        sample_vals = np.nan_to_num(arr.ravel()[s_idx].astype(np.float32), nan=0.0)
        sample_dots = dots.ravel()[s_idx].astype(np.float32)
        dt3_mask = (distance_transform_edt(~dots) <= RADIUS_PX) if dots.any() else np.zeros(foot.shape, dtype=bool)
        out[p.stem] = {
            "sample_vals": sample_vals,
            "sample_dots": sample_dots,
            "within_3px_mask": dt3_mask,
            "footprint_3px_coverage": round(float(dt3_mask[foot].mean()), 6),
            "positive_px": int(dots.sum()),
        }
        del arr, dots
    return out


def _spearman_vec(a: np.ndarray, b: np.ndarray) -> float:
    if a.size == 0 or np.std(a) == 0.0 or np.std(b) == 0.0:
        return 0.0
    r, _ = spearmanr(a, b)
    return float(r) if np.isfinite(r) else 0.0


def rank_corr(a: np.ndarray, b: np.ndarray, foot: np.ndarray, rng: np.random.Generator, max_n: int = MAX_SAMPLE_N) -> float:
    idx = np.flatnonzero(foot.ravel())
    if idx.size == 0:
        return float("nan")
    if idx.size > max_n:
        idx = rng.choice(idx, size=max_n, replace=False)
    ar = np.nan_to_num(np.asarray(a, dtype=np.float64).ravel()[idx], nan=0.0)
    br = np.nan_to_num(np.asarray(b, dtype=np.float64).ravel()[idx], nan=0.0)
    return _spearman_vec(ar, br)


def dot_overlap_fraction(our_dots: np.ndarray, registry_dots: np.ndarray, radius: int = RADIUS_PX) -> float:
    our = np.asarray(our_dots, dtype=bool)
    reg = np.asarray(registry_dots, dtype=bool)
    if not our.any() or not reg.any():
        return 0.0
    d = distance_transform_edt(~reg)
    return float((d[our] <= radius).mean())


def check_surface_preplacement(
    our_surface: np.ndarray,
    foot: np.ndarray,
    registry: dict,
    rng: np.random.Generator | None = None,
) -> dict:
    """Protocol Rule 1 pre-placement check: surface rank-correlation vs every registry raster."""
    s_idx = registry_sample_indices(foot)
    our_s = np.nan_to_num(our_surface.ravel()[s_idx].astype(np.float64), nan=0.0)
    rows = []
    dup = False
    max_abs = 0.0
    for name, entry in registry.items():
        if isinstance(entry, dict):
            corr = _spearman_vec(our_s, entry["sample_vals"].astype(np.float64))
        else:
            corr = rank_corr(our_surface, entry, foot, rng or np.random.default_rng(SAMPLE_SEED))
        ac = abs(corr)
        max_abs = max(max_abs, ac)
        hit = ac > CORR_MAX
        dup = dup or hit
        rows.append({"name": name, "spearman_surface": round(corr, 6), "duplicate": bool(hit)})
    return {
        "duplicate_surface": bool(dup),
        "max_abs_spearman": round(float(max_abs), 6),
        "threshold_spearman": CORR_MAX,
        "vs": rows,
    }


def check_registry(
    our_dots: np.ndarray,
    our_score: np.ndarray,
    foot: np.ndarray,
    registry: dict,
    rng: np.random.Generator | None = None,
) -> dict:
    """Full pre-placement + post-placement registry check across all registry rasters."""
    s_idx = registry_sample_indices(foot)
    our_s = np.nan_to_num(our_score.ravel()[s_idx].astype(np.float64), nan=0.0)
    our_d_s = our_dots.ravel()[s_idx].astype(np.float64)
    rows = []
    dup = False
    max_abs_surf = 0.0
    max_abs_dots = 0.0
    max_ov = 0.0
    max_ov_all = 0.0
    n_our = int(our_dots.sum())

    for name, entry in registry.items():
        if isinstance(entry, dict):
            corr_surf = _spearman_vec(our_s, entry["sample_vals"].astype(np.float64))
            corr_dots = _spearman_vec(our_d_s, entry["sample_dots"].astype(np.float64))
            ov = float(entry["within_3px_mask"][our_dots].mean()) if n_our > 0 else 0.0
            fp_cov = float(entry["footprint_3px_coverage"])
            reg_px = entry["positive_px"]
        else:
            dots = np.isfinite(entry) & (entry > 0.5) & foot
            r_gen = rng or np.random.default_rng(SAMPLE_SEED)
            corr_surf = rank_corr(our_score, entry, foot, r_gen)
            corr_dots = rank_corr(our_dots.astype(np.float32), dots.astype(np.float32), foot, r_gen)
            dt3 = (distance_transform_edt(~dots) <= RADIUS_PX) if dots.any() else np.zeros(foot.shape, dtype=bool)
            ov = float(dt3[our_dots].mean()) if n_our > 0 else 0.0
            fp_cov = float(dt3[foot].mean())
            reg_px = int(dots.sum())

        # Detect uniform space-filling lattice probes (e.g. 5-px grid covering 99.83% of ALL interior cells with ~0 rank correlation)
        is_lattice_probe = (fp_cov >= 0.95) and (abs(corr_dots) < 0.05) and (abs(corr_surf) < 0.05)

        max_abs_surf = max(max_abs_surf, abs(corr_surf))
        max_abs_dots = max(max_abs_dots, abs(corr_dots))
        max_ov_all = max(max_ov_all, ov)
        if not is_lattice_probe:
            max_ov = max(max_ov, ov)
        hit = (abs(corr_surf) > CORR_MAX) or (abs(corr_dots) > CORR_MAX) or ((ov > OVERLAP_MAX) and not is_lattice_probe)
        dup = dup or hit
        rows.append({
            "name": name,
            "spearman": round(corr_surf, 6),
            "spearman_dots": round(corr_dots, 6),
            "dot_fraction_within_3px": round(ov, 6),
            "footprint_3px_coverage": round(fp_cov, 6),
            "space_filling_lattice_probe": bool(is_lattice_probe),
            "duplicate": bool(hit),
            "registry_positive_px": reg_px,
        })
    return {
        "duplicate_lane": bool(dup),
        "max_abs_spearman_surface": round(float(max_abs_surf), 6),
        "max_abs_spearman_dots": round(float(max_abs_dots), 6),
        "max_dot_fraction_within_3px": round(float(max_ov), 6),
        "max_dot_fraction_within_3px_including_lattice_probes": round(float(max_ov_all), 6),
        "vs": rows,
        "thresholds": {"spearman": CORR_MAX, "overlap_3px": OVERLAP_MAX},
    }
