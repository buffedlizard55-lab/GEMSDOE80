"""Memory-Mapped 34-Channel Label-Free Static Multi-Physics Stack (`data/cache_h80v2/`).

Stores footprint-indexed (`N_fp = 5,167,373`) standardized `float16` features on disk so
peak RAM in a 3.8 GiB sandbox never exceeds ~650 MB while retaining all 5 physical domains:
  - Domain 1: 1m USGS 3DEP LiDAR scarplet micro-geomorphology (native NaN outside tiles!)
  - Domain 2: Seamless 100m DEM breaklines, Piedmont bajada index & Pluvial Lake shoreline rejection
  - Domain 3: Up-dip-pulled Gravity/Magnetic/MT/Basement potential-field worms & ridges
  - Domain 4: Exact 3D Crustal-Stress Cauchy Slip/Dilation Tendency (Morris-Ferrill-Simpson)
  - Domain 5: Geothermal upflow plumes (GDR 1391 geothermometers + ScienceBase heat flow +
              3D Euler contacts + GeoDAWN K/Th boiling-zone & demagnetization conduits)
"""

from __future__ import annotations

import gc
import json
from pathlib import Path

import numpy as np
import rasterio
from scipy.ndimage import gaussian_filter, maximum_filter, minimum_filter

from .geomechanics import (
    build_regional_stress_and_tendency,
    compute_strike_worm,
    fill_median,
    joint_strike_field,
    multiscale_hessian_ridge,
    robust_zscore,
    updip_pullback,
)
from .geothermal import build_geothermal_channels

STATIC_FEATURE_NAMES: tuple[str, ...] = (
    # Domain 1: 1m LiDAR Scarp Micro-Geomorphology (8 channels; NaN outside tile coverage)
    "lid_ex_max",
    "lid_step_max",
    "lid_tect_ratio",
    "lid_crest_curv",
    "lid_coh100",
    "lid_chain_vote",
    "lid_stress_align",
    "lid_valid",
    # Domain 2: Seamless 100m Topography, Piedmont Habitat & Shoreline Decontamination (6 channels)
    "elev_det",
    "elev_slope",
    "elev_lrm_breakline",
    "elev_hess_ridge",
    "piedmont_index",
    "shoreline_mimic_penalty",
    # Domain 3: Up-Dip-Pulled Potential-Field Worms & Subsurface Basement/MT Conduits (8 channels)
    "grav_hg_local",
    "grav_updip_pull",
    "grav_strike_worm",
    "mag_tmi_hg_local",
    "mag_updip_pull",
    "mag_strike_worm",
    "base_depth_step",
    "mt_cond_edge",
    # Domain 4: 3D Crustal Geomechanics & Strain (6 channels)
    "stress_ts",
    "stress_td",
    "stress_crit",
    "stress_aphi",
    "orient_coherence",
    "geod_strain_combo",
    # Domain 5: Geothermal Upflow, Heat Flow, Euler Contacts & Alteration + 4-Line Gate (6 channels)
    "gdr_geotherm_plume",
    "sb_heatflow_anom",
    "euler_contact_density",
    "rad_kth_alter",
    "rad_uth_travertine",
    "hydro_demag_conduit",
)


def build_or_load_static_stack(
    data_dir: Path,
    foot: np.ndarray,
    transform,
) -> tuple[np.memmap, dict[str, np.ndarray], list[str]]:
    """Return `(mmap_f16, aux_maps, feature_names)` where `mmap_f16` has shape `(N_fp, 34)`."""
    cache_dir = data_dir / "cache_h80v2"
    cache_dir.mkdir(parents=True, exist_ok=True)
    mmap_path = cache_dir / "static_f16.mmap"
    aux_path = cache_dir / "aux_fields.npz"
    meta_path = cache_dir / "meta.json"

    fp_idx = np.flatnonzero(foot.ravel())
    n_fp = int(fp_idx.size)
    n_ch = len(STATIC_FEATURE_NAMES)

    if mmap_path.is_file() and aux_path.is_file() and meta_path.is_file():
        meta = json.loads(meta_path.read_text())
        if meta.get("n_fp") == n_fp and meta.get("names") == list(STATIC_FEATURE_NAMES):
            mm = np.memmap(mmap_path, dtype=np.float16, mode="r", shape=(n_fp, n_ch))
            aux = dict(np.load(aux_path))
            return mm, aux, list(STATIC_FEATURE_NAMES)

    mm = np.memmap(mmap_path, dtype=np.float16, mode="w+", shape=(n_fp, n_ch))

    def put(col_name: str, arr2d: np.ndarray, keep_nan_mask: np.ndarray | None = None, standardize: bool = True) -> None:
        j = STATIC_FEATURE_NAMES.index(col_name)
        if standardize:
            v_mask = foot if keep_nan_mask is None else (foot & keep_nan_mask)
            z = robust_zscore(arr2d, v_mask, clip_lo=-6.0, clip_hi=9.0)
            if keep_nan_mask is not None:
                vec = z.ravel()[fp_idx].astype(np.float16)
                vec[~keep_nan_mask.ravel()[fp_idx]] = np.float16(np.nan)
                mm[:, j] = vec
                return
            mm[:, j] = z.ravel()[fp_idx].astype(np.float16)
        else:
            mm[:, j] = arr2d.ravel()[fp_idx].astype(np.float16)

    feat_path = data_dir / "training_features.tif"
    nodata_val = np.float32(-3.4028234663852886e38)

    def read_fband(b: int) -> np.ndarray:
        with rasterio.open(feat_path) as src:
            a = src.read(b).astype(np.float32)
        return np.where((a == nodata_val) | ~np.isfinite(a), np.nan, a)

    # Load key 100m bands carefully with bounded RAM
    rtp = read_fband(2)
    tmi_hg = read_fband(3)
    geod_inv = read_fband(4)
    geod_shear = read_fband(7)
    geod_dil = read_fband(8)
    tmi_vg = read_fband(9)
    det_elev = read_fband(12)
    iso_grav = read_fband(13)
    depth_base = read_fband(15)
    cond_surf = read_fband(17)
    grav_hg = read_fband(18)
    elev_slope = read_fband(19)

    # Upward-continued 150m TMI for deep magnetic basement structure
    ext_tif = data_dir / "external" / "geodawn_extensions_u8.tif"
    if ext_tif.is_file():
        with rasterio.open(ext_tif) as se:
            tmi_up150 = se.read(4).astype(np.float32)
        gy_u, gx_u = np.gradient(gaussian_filter(fill_median(tmi_up150, foot), sigma=1.5, mode="nearest"))
        deep_mag_hg = np.hypot(gy_u, gx_u).astype(np.float32)
        del tmi_up150, gy_u, gx_u
    else:
        deep_mag_hg = fill_median(tmi_hg, foot)

    # 1m LiDAR scarp features (native NaN outside the 706 tiles!)
    lid_tif = data_dir / "external" / "lidar_scarp_features_u8.tif"
    if lid_tif.is_file():
        with rasterio.open(lid_tif) as sl:
            l_ex = sl.read(1).astype(np.float32)
            l_step = sl.read(3).astype(np.float32)
            l_lapn = sl.read(4).astype(np.float32)
            l_lapp = sl.read(5).astype(np.float32)
            l_down = sl.read(6).astype(np.float32)
            l_up = sl.read(7).astype(np.float32)
            l_coh = sl.read(10).astype(np.float32)
            l_strike = sl.read(11).astype(np.float32)
            l_val_mask = (sl.read(12) > 0) & foot
    else:
        l_ex = l_step = l_lapn = l_lapp = l_down = l_up = l_coh = l_strike = np.zeros(foot.shape, dtype=np.float32)
        l_val_mask = np.zeros(foot.shape, dtype=bool)

    lid_tect = (l_down - 0.65 * l_up) / (l_down + l_up + 4.0)
    lid_curv = np.sqrt(np.maximum(l_lapn * l_lapp, 0.0))
    lid_chain = compute_strike_worm(l_step * (0.3 + l_coh / 255.0), l_val_mask, half_len=5, flank=3)

    # Multi-physics orientation & 3D crustal geomechanics
    gy_b, gx_b = np.gradient(gaussian_filter(fill_median(depth_base, foot), sigma=1.5, mode="nearest"))
    base_step = np.hypot(gy_b, gx_b).astype(np.float32)
    del gy_b, gx_b

    grav_pull = updip_pullback(grav_hg, iso_grav, foot, shift_px=3.0)
    mag_pull = updip_pullback(deep_mag_hg, iso_grav, foot, shift_px=2.5)

    cand_strike_az, r_orient = joint_strike_field(
        [
            (det_elev, 1.3),
            (np.where(l_val_mask, l_step, np.nan), 1.2),
            (iso_grav, 1.0),
            (rtp, 0.9),
            (depth_base, 1.0),
        ],
        foot,
    )

    stress = build_regional_stress_and_tendency(
        data_dir / "external" / "sb_slip_tendency_in_footprint.csv",
        transform,
        foot.shape,
        cand_strike_az,
        r_orient,
        foot,
    )

    # 1m LiDAR strike alignment with optimal regional SHmax strike
    lid_az_rad = (np.pi * l_strike / 255.0) % np.pi
    lid_align = np.abs(np.cos(lid_az_rad - stress["shmax_az"])).astype(np.float32)

    put("lid_ex_max", l_ex, keep_nan_mask=l_val_mask)
    put("lid_step_max", l_step, keep_nan_mask=l_val_mask)
    put("lid_tect_ratio", lid_tect, keep_nan_mask=l_val_mask)
    put("lid_crest_curv", lid_curv, keep_nan_mask=l_val_mask)
    put("lid_coh100", l_coh, keep_nan_mask=l_val_mask)
    put("lid_chain_vote", lid_chain, keep_nan_mask=l_val_mask)
    put("lid_stress_align", lid_align, keep_nan_mask=l_val_mask)
    put("lid_valid", l_val_mask.astype(np.float32), standardize=False)

    del l_ex, l_lapn, l_lapp, l_down, l_up, l_coh, l_strike, lid_tect, lid_curv, lid_chain, lid_az_rad, lid_align
    gc.collect()

    # Domain 2: Topography, Piedmont Habitat & Pluvial Paleoshoreline Rejection
    elev_f = fill_median(det_elev, foot)
    lrm_resid = elev_f - gaussian_filter(elev_f, sigma=5.0, mode="nearest")
    gy_l, gx_l = np.gradient(lrm_resid)
    lrm_break = np.hypot(gy_l, gx_l).astype(np.float32)
    del lrm_resid, gy_l, gx_l

    elev_hess = multiscale_hessian_ridge(det_elev, foot, sigmas=(1.2, 2.5), polarity="both")

    # Piedmont bajada habitat: close (0.5 - 6 km) to regional mountain front relief without sitting high inside mountain crest or flat central playa
    local_relief_3km = maximum_filter(elev_f, size=25) - minimum_filter(elev_f, size=25)
    near_front = gaussian_filter((local_relief_3km > np.percentile(local_relief_3km[foot], 72)).astype(np.float32), sigma=28.0, mode="nearest")
    slope_f = fill_median(elev_slope, foot)
    not_rugged_core = np.exp(-np.maximum(slope_f - np.percentile(slope_f[foot], 78), 0.0) / (np.std(slope_f[foot]) + 1e-6))
    not_dead_playa = 1.0 - np.exp(-slope_f / (0.18 * np.median(slope_f[foot]) + 1e-6))
    piedmont_idx = (near_front * not_rugged_core * (0.35 + 0.65 * not_dead_playa)).astype(np.float32)

    # Pluvial shoreline mimic score: sharp narrow bench along smooth contour where subsurface gravity/basement step is absent
    subsurf_support = 0.5 * robust_zscore(grav_pull, foot, -3.0, 6.0) + 0.5 * robust_zscore(base_step, foot, -3.0, 6.0)
    shoreline_penalty = (
        np.maximum(robust_zscore(lrm_break, foot, -3.0, 6.0), 0.0)
        * np.exp(-np.maximum(subsurf_support, 0.0))
        * (np.abs(elev_f) < np.percentile(np.abs(elev_f[foot]), 55))
    ).astype(np.float32)

    put("elev_det", det_elev)
    put("elev_slope", elev_slope)
    put("elev_lrm_breakline", lrm_break)
    put("elev_hess_ridge", elev_hess)
    put("piedmont_index", piedmont_idx)
    put("shoreline_mimic_penalty", shoreline_penalty)

    del lrm_break, elev_hess, local_relief_3km, near_front, not_rugged_core, not_dead_playa, piedmont_idx, subsurf_support, shoreline_penalty, l_step
    gc.collect()

    # Domain 3: Potential-Field Worms & Subsurface Conduits
    grav_worm = compute_strike_worm(grav_hg, foot, half_len=7, flank=3)
    mag_worm = compute_strike_worm(tmi_hg, foot, half_len=7, flank=3)
    cond_f = fill_median(cond_surf, foot)
    gy_c, gx_c = np.gradient(gaussian_filter(cond_f, sigma=1.5, mode="nearest"))
    mt_edge = np.hypot(gy_c, gx_c).astype(np.float32)
    del gy_c, gx_c

    grav_f = fill_median(grav_hg, foot)
    grav_loc = grav_f - gaussian_filter(grav_f, sigma=25.0, mode="nearest")
    mag_f = fill_median(tmi_hg, foot)
    mag_loc = mag_f - gaussian_filter(mag_f, sigma=25.0, mode="nearest")

    put("grav_hg_local", grav_loc)
    put("grav_updip_pull", grav_pull)
    put("grav_strike_worm", grav_worm)
    put("mag_tmi_hg_local", mag_loc)
    put("mag_updip_pull", mag_pull)
    put("mag_strike_worm", mag_worm)
    put("base_depth_step", base_step)
    put("mt_cond_edge", mt_edge)

    del grav_worm, mag_worm, mt_edge, grav_loc, mag_loc, grav_pull, mag_pull, base_step, deep_mag_hg
    gc.collect()

    # Domain 4: 3D Crustal Geomechanics & Regional Strain
    strain_combo = (
        0.45 * robust_zscore(geod_inv, foot)
        + 0.35 * robust_zscore(geod_shear, foot)
        + 0.20 * robust_zscore(geod_dil, foot)
    )
    put("stress_ts", stress["stress_ts"])
    put("stress_td", stress["stress_td"])
    put("stress_crit", stress["stress_crit"])
    put("stress_aphi", stress["stress_aphi"])
    put("orient_coherence", r_orient)
    put("geod_strain_combo", strain_combo)
    del strain_combo, geod_inv, geod_shear, geod_dil
    gc.collect()

    # Domain 5: Geothermal Upflow & Hydrothermal Alteration Inversion
    geoth = build_geothermal_channels(data_dir, transform, foot.shape, rtp, tmi_vg, cond_surf, foot)
    for k_name, arr_g in geoth.items():
        put(k_name, arr_g)

    mm.flush()
    np.savez_compressed(
        aux_path,
        cand_strike_az=cand_strike_az.astype(np.float32),
        stress_td=stress["stress_td"].astype(np.float32),
    )
    meta_path.write_text(json.dumps({"n_fp": n_fp, "names": list(STATIC_FEATURE_NAMES)}, indent=2) + "\n")

    mm_read = np.memmap(mmap_path, dtype=np.float16, mode="r", shape=(n_fp, n_ch))
    return mm_read, {"cand_strike_az": cand_strike_az, "stress_td": stress["stress_td"]}, list(STATIC_FEATURE_NAMES)
