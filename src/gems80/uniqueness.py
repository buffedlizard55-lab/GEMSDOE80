"""Lane-duplicate check vs registry rasters (parallel-run protocol).

Duplicate if Spearman rank-correlation on the footprint > 0.90
OR more than 70% of our dots fall within 3 px of one registry raster's dots.
"""

from __future__ import annotations

import numpy as np
from scipy.ndimage import distance_transform_edt
from scipy.stats import spearmanr

CORR_MAX = 0.90
OVERLAP_MAX = 0.70
RADIUS_PX = 3


def rank_corr(a: np.ndarray, b: np.ndarray, foot: np.ndarray, rng: np.random.Generator, max_n: int = 250_000) -> float:
    idx = np.flatnonzero(foot.ravel())
    if idx.size == 0:
        return float("nan")
    if idx.size > max_n:
        idx = rng.choice(idx, size=max_n, replace=False)
    ar = np.asarray(a, dtype=np.float64).ravel()[idx]
    br = np.asarray(b, dtype=np.float64).ravel()[idx]
    if np.std(ar) == 0 or np.std(br) == 0:
        return 0.0
    r, _ = spearmanr(ar, br)
    return float(r)


def dot_overlap_fraction(our_dots: np.ndarray, registry_dots: np.ndarray, radius: int = RADIUS_PX) -> float:
    our = np.asarray(our_dots, dtype=bool)
    reg = np.asarray(registry_dots, dtype=bool)
    n = int(our.sum())
    if n == 0:
        return 0.0
    if not reg.any():
        return 0.0
    d = distance_transform_edt(~reg)
    return float((d[our] <= radius).mean())


def check_registry(our_dots: np.ndarray, our_score: np.ndarray, foot: np.ndarray, registry: dict[str, np.ndarray],
                   rng: np.random.Generator) -> dict:
    rows = []
    dup = False
    for name, arr in registry.items():
        dots = np.isfinite(arr) & (arr > 0.5)
        corr = rank_corr(our_score, arr, foot, rng)
        ov = dot_overlap_fraction(our_dots, dots)
        hit = (corr is not None and abs(corr) > CORR_MAX) or (ov > OVERLAP_MAX)
        dup = dup or hit
        rows.append({
            "name": name,
            "spearman": corr,
            "dot_fraction_within_3px": ov,
            "duplicate": bool(hit),
            "registry_positive_px": int(dots.sum()),
        })
    return {"duplicate_lane": bool(dup), "vs": rows, "thresholds": {"spearman": CORR_MAX, "overlap_3px": OVERLAP_MAX}}
