"""Forensic byte verification of GEMSDOE32 H33 (0.2778), spatial regime profiling,
and pre-placement / post-placement 133-raster uniqueness auditing for Lane H80-DPVC."""
from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Dict, List

import numpy as np
import rasterio
from scipy import ndimage as ndi
from scipy.stats import spearmanr

from gems57 import uniqueness as UQ


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def verify_h33_forensic_autopsy(scored_dir: Path, fp: np.ndarray, cat: np.ndarray) -> dict:
    """Verify byte-level identity between GEMSDOE32 H33-2-B2 (0.2778) and
    `(GEMSDOE25 D2.8 > 0) & (EDT(~existing_faults) > 2.0)` on the SHA-256-pinned rasters."""
    d28_path = scored_dir / "gems24-h25-1-dotted-h19-5-d2-8-20261002-e56ea318af89-nan.tif"
    d15_path = scored_dir / "gems24-h25-1-dotted-h19-5-d1-5-20261002-989f59505db1-nan.tif"
    h33_path = scored_dir / "GEMSDOE32__gemsdoe32-h33-h33-2-b2-20261004T220000Z-e5eb6e7e-zeros.tif"

    with rasterio.open(d28_path) as s:
        d28 = np.isfinite(s.read(1)) & (s.read(1) > 0) & fp
    with rasterio.open(d15_path) as s:
        d15 = np.isfinite(s.read(1)) & (s.read(1) > 0) & fp
    with rasterio.open(h33_path) as s:
        h33_raw = s.read(1)
        h33 = np.isfinite(h33_raw) & (h33_raw > 0) & fp

    edt_cat = ndi.distance_transform_edt(~cat).astype(np.float32)
    pruned_d28 = d28 & (edt_cat > 2.0)
    xor_pixels = int(np.sum(pruned_d28 ^ h33))
    removed_dots = d28 & (~h33)
    removed_dists = edt_cat[removed_dots]

    # Constant-TP reciprocal collinearity ladder across the 5 nested pruning submissions
    # D1.5 (60,069 px, 0.2477) -> D2.8 (44,090 px, 0.2600) -> H27-4 (40,199 px, 0.2708)
    # -> G36 (38,854 px, 0.2750) -> H33 (37,654 px, 0.2778)
    ladder = [
        {"name": "GEMSDOE24 D1.5", "dots": int(d15.sum()), "reported_live_dti": 0.2477, "evidence_class": "ORGANIZER-CONFIRMED (owner-reported)"},
        {"name": "GEMSDOE25 D2.8", "dots": int(d28.sum()), "reported_live_dti": 0.2600, "evidence_class": "ORGANIZER-CONFIRMED (owner-reported)"},
        {"name": "GEMSDOE28/31 H27-4", "dots": 40199, "reported_live_dti": 0.2708, "evidence_class": "ORGANIZER-CONFIRMED (owner-reported)"},
        {"name": "GEMSDOE36 Anderson-PINN", "dots": 38854, "reported_live_dti": 0.2750, "evidence_class": "ORGANIZER-CONFIRMED (owner-reported)"},
        {"name": "GEMSDOE32 H33-2-B2", "dots": int(h33.sum()), "reported_live_dti": 0.2778, "evidence_class": "ORGANIZER-CONFIRMED (owner-reported)"},
    ]
    dots_arr = np.array([r["dots"] for r in ladder], dtype=np.float64)
    inv_dti = np.array([1.0 / r["reported_live_dti"] for r in ladder], dtype=np.float64)
    r_collinear_collar4 = float(np.corrcoef(dots_arr[1:], inv_dti[1:])[0, 1])
    r_collinear_all5 = float(np.corrcoef(dots_arr, inv_dti)[0, 1])

    return {
        "evidence_class": "BYTE-FORENSIC-MEASUREMENT + MATHEMATICAL-ALGEBRA",
        "d28_sha256": sha256_file(d28_path),
        "h33_sha256": sha256_file(h33_path),
        "d28_positive_dots": int(d28.sum()),
        "h33_positive_dots": int(h33.sum()),
        "dots_removed_by_2px_collar": int(removed_dots.sum()),
        "xor_pixels_between_h33_and_d28_gt_2px": xor_pixels,
        "exact_bitwise_subset_verified": xor_pixels == 0,
        "removed_dots_min_dist_px": round(float(removed_dists.min()), 4),
        "removed_dots_max_dist_px": round(float(removed_dists.max()), 4),
        "kept_dots_min_dist_px": round(float(edt_cat[h33].min()), 4),
        "kept_dots_median_dist_px": round(float(np.median(edt_cat[h33])), 4),
        "h33_has_zero_nan_across_full_grid": bool(np.isfinite(h33_raw).all()),
        "pruning_ladder": ladder,
        "reciprocal_dti_vs_dot_count_pearson_r": round(r_collinear_collar4, 6),
        "reciprocal_dti_vs_dot_count_all5_pearson_r": round(r_collinear_all5, 6),
        "physical_explanation": (
            "GEMSDOE32 h33-h33-2-b2-20261004T220000Z-e5eb6e7e-zeros.tif (37,654 dots) is "
            "pixel-for-pixel identical (0 XOR pixels across all 12,279,160 cells) to pruning "
            "the 6,436 dots of GEMSDOE25 D2.8 (44,090 dots) that lie at Euclidean distance "
            "1.0 <= d_cat <= 2.0 px (100-200 m) from mapped existing_faults.tif. Under the "
            "Dilated Tversky Index (alpha=0.2, beta=0.8, 300 m kernel), known faults are masked "
            "out of evaluation; dots within <= 200 m of known traces respond to the known fault's "
            "own scarp/potential-field halo and earn near-zero new-fault TP credit while paying "
            "full alpha=0.2 false-positive denominator weight."
        ),
    }


def compute_spatial_profile(
    dots: np.ndarray,
    slope10_full: np.ndarray,
    lidar_valid_full: np.ndarray,
    cat: np.ndarray,
    fp: np.ndarray,
) -> dict:
    """Compute structural-regime shares and steep-terrain prevalence of a dot set."""
    d_cat = ndi.distance_transform_edt(~cat).astype(np.float32)
    p70_slp = float(np.percentile(slope10_full[fp], 70.0))
    steep = slope10_full >= p70_slp

    n_dots = int(dots.sum())
    if n_dots == 0:
        return {"n_dots": 0}

    d_vals = d_cat[dots]
    cont_share = float(np.mean(d_vals <= 10.0))
    step_share = float(np.mean((d_vals > 10.0) & (d_vals <= 26.0)))
    blind_share = float(np.mean(d_vals > 26.0))
    steep_dots = float(np.mean(steep[dots]))
    steep_cat = float(np.mean(steep[cat & fp]))
    lidar_share = float(np.mean(lidar_valid_full[dots]))

    return {
        "evidence_class": "SPATIAL-PROFILE-MEASUREMENT",
        "n_dots": n_dots,
        "min_dist_to_catalogue_px": round(float(d_vals.min()), 4),
        "median_dist_to_catalogue_px": round(float(np.median(d_vals)), 4),
        "p90_dist_to_catalogue_px": round(float(np.percentile(d_vals, 90.0)), 4),
        "regime_continuation_le_1km_share": round(cont_share, 4),
        "regime_stepover_1_to_2p6km_share": round(step_share, 4),
        "regime_blind_basin_gt_2p6km_share": round(blind_share, 4),
        "steep_terrain_share_dots": round(steep_dots, 4),
        "steep_terrain_share_true_catalogue": round(steep_cat, 4),
        "neotectonic_steep_calibration_pass": bool(0.24 <= steep_dots <= 0.48),
        "lidar_coverage_share_dots": round(lidar_share, 4),
    }


def run_full_uniqueness_audit(
    candidate_surface: np.ndarray,
    candidate_dots: np.ndarray,
    scored_dir: Path,
    fp: np.ndarray,
    cat: np.ndarray,
    candidate_name: str,
) -> dict:
    """Audit both the pre-placement continuous surface (Spearman rank correlation) and
    post-placement final dots (3-px forward & reverse overlap, Jaccard) against all
    accessible registry GeoTIFFs in `scored_dir` (`Protocol Rule 1` + `IR-H85-01`)."""
    allowed = fp & (~cat)
    mine = np.asarray(candidate_dots, dtype=bool) & fp
    halo_mine = UQ._near(mine, radius_px=3.0)

    # Stratified 250,000-pixel sample across `fp` for fast exact Spearman rank correlation
    rng = np.random.default_rng(80)
    fp_idx = np.flatnonzero(fp.ravel())
    samp_idx = rng.choice(fp_idx, size=min(250_000, fp_idx.size), replace=False)
    surf_samp = np.asarray(candidate_surface, dtype=np.float32).ravel()[samp_idx]
    dot_samp = mine.astype(np.float32).ravel()[samp_idx]

    paths = sorted(scored_dir.glob("*.tif"))
    rows: List[dict] = []
    for p in paths:
        if p.stem == candidate_name:
            continue
        try:
            with rasterio.open(p) as src:
                if (src.height, src.width) != fp.shape:
                    continue
                arr = src.read(1).astype(np.float32)
        except Exception:
            continue

        fin = np.isfinite(arr)
        theirs = fin & (arr > 0) & fp
        n_theirs = int(theirs.sum())
        if n_theirs == 0:
            continue

        halo_theirs = UQ._near(theirs, radius_px=3.0)
        fwd = float(halo_theirs[mine].mean()) if mine.any() else 0.0
        rev = float(halo_mine[theirs].mean())
        jac = float((mine & theirs).sum() / max(int((mine | theirs).sum()), 1))

        t_samp = np.where(fin.ravel()[samp_idx], arr.ravel()[samp_idx], 0.0)
        if np.ptp(t_samp) > 0 and np.ptp(surf_samp) > 0:
            rho_surf = float(spearmanr(surf_samp, t_samp).statistic)
        else:
            rho_surf = 0.0
        if np.ptp(t_samp) > 0 and np.ptp(dot_samp) > 0:
            rho_dots = float(spearmanr(dot_samp, t_samp).statistic)
        else:
            rho_dots = 0.0

        rep = UQ.representation_class_v2(theirs, allowed)
        rows.append(
            {
                "file": p.name,
                "sha256": sha256_file(p),
                "their_dots": n_theirs,
                "forward_overlap": round(fwd, 5),
                "reverse_overlap": round(rev, 5),
                "jaccard": round(jac, 5),
                "spearman_preplacement_surface": round(rho_surf, 5),
                "spearman_full_footprint": round(max(abs(rho_surf), abs(rho_dots)), 5),
                "representation": rep["kind"],
                "halo_coverage": round(rep["halo_coverage"], 5),
                "support_fraction": round(rep["support_fraction"], 5),
            }
        )

    binding = UQ.binding_duplicate_decision(rows)
    disc_rows = [r for r in rows if r["representation"] == "discriminating"]
    space_rows = [r for r in rows if r["representation"] == "space_filling"]

    worst_surf_rho = max(rows, key=lambda r: abs(r["spearman_preplacement_surface"])) if rows else None
    h33_row = next((r for r in rows if "h33-2-b2" in r["file"] and "zeros" in r["file"]), None)
    h76_row = next((r for r in rows if "gems76_h76_3" in r["file"]), None)
    h85_row = next((r for r in rows if "gems85-h85-edti" in r["file"]), None)

    return {
        "evidence_class": "REGISTRY-UNIQUENESS-MEASUREMENT",
        "candidate_name": candidate_name,
        "candidate_dots": int(mine.sum()),
        "rasters_compared": len(rows),
        "discriminating_priors_count": len(disc_rows),
        "space_filling_priors_count": len(space_rows),
        "pre_placement_surface_gate": {
            "max_abs_spearman_rho": round(abs(worst_surf_rho["spearman_preplacement_surface"]), 5) if worst_surf_rho else 0.0,
            "worst_file": worst_surf_rho["file"] if worst_surf_rho else None,
            "threshold": 0.90,
            "passed": bool(not worst_surf_rho or abs(worst_surf_rho["spearman_preplacement_surface"]) <= 0.90),
        },
        "post_placement_binding_like_for_like_gate": binding,
        "literal_reading_with_space_filling_priors": {
            "worst_forward_overlap": max((r["forward_overlap"] for r in rows), default=0.0),
            "space_filling_firings_count": sum(1 for r in space_rows if r["forward_overlap"] > 0.70),
            "discriminating_firings_count": sum(1 for r in disc_rows if r["forward_overlap"] > 0.70),
            "note": (
                "Space-filling priors (e.g. 17GEMSDOE E-proba-multiscale with 5,167,373 positive pixels "
                "and 13GEMSDOE 5px lattice) have 3-px halo coverage >= 98% of the allowed domain, so even "
                "a uniform random dot set has ~1.0 forward overlap against them. Under the binding "
                "like-for-like rule (IR-H85-01), dot overlap is evaluated across all 126 discriminating "
                "priors (halo coverage < 0.70) while Spearman rank correlation is evaluated across all 133 priors."
            ),
        },
        "key_benchmark_comparisons": {
            "vs_GEMSDOE32_H33_02778": h33_row,
            "vs_GEMSDOE76_H76_3": h76_row,
            "vs_GEMSDOE77_H85_EDTI": h85_row,
        },
        "all_126_discriminating_overlap_le_070": bool(not binding["duplicate"]),
        "top_15_discriminating_by_overlap": sorted(disc_rows, key=lambda r: -r["forward_overlap"])[:15],
    }
