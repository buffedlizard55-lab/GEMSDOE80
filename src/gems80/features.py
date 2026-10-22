"""Multi-physics, updip 60-degree normal-fault basement-to-surface advection,
hydrothermal magnetite-destruction / clay-cap / geothermometer vent conduits,
TIGER/BLM/Lahontan-deconfounded 1m/10m tectonic geomorphology, and fold-proper
fracture-transfer features for Lane H80-DPVC (`GEMSDOE80`).

Memory-engineered for `< 1.1 GB` peak container RSS on a 3.85 GB sandbox:
- Each physical domain computes its channels inside an isolated helper scope and
  releases all `(3730, 3292)` intermediates via `gc.collect()`;
- Static features are persisted on disk as clipped `float16` (`5,167,373 x 50`),
  lossless for 255-bin gradient boosted trees.
"""
from __future__ import annotations

import gc
import io
import json
import struct
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import rasterio
from rasterio.warp import transform as warp_transform
from scipy import ndimage as ndi
from scipy.spatial import cKDTree
from skimage.morphology import skeletonize
from sklearn.metrics import roc_auc_score

GRID_SHAPE = (3730, 3292)
ORIGIN_X = 243350.0
ORIGIN_Y = 4508550.0
PIXEL_M = 100.0
COT_60_DEG = float(1.0 / np.sqrt(3.0))  # 0.57735027

LIDAR_SCALES: Dict[int, Tuple[float, float]] = {
    1: (0.0, 2.5),
    2: (0.0, 0.8),
    3: (0.0, 6.0),
    4: (0.0, 2.5),
    5: (0.0, 2.5),
    6: (0.0, 1.5),
    7: (0.0, 1.5),
    8: (0.0, 2.5),
    9: (0.0, 60.0),     # True 1m local relief (m)
    10: (0.0, 1.0),     # 100 m scarp strike coherence
    11: (0.0, np.pi),   # Local scarp strike (radians)
}

SCARP_STATIC_NAMES: List[str] = [
    "lid_deconf_piedmont_step",   # 0
    "lid_deconf_scarp_dipole",    # 1
    "lid_tect_vs_fluvial",        # 2
    "lid_antislope_scarp",        # 3
    "lid_step_anom11",            # 4
    "lid_coh100",                 # 5
    "lid_worm_1100m",             # 6
    "t10_deconf_piedmont_break",  # 7
    "t10_aspect_planar",          # 8
    "t10_slope_local_anom",       # 9
    "t10_worm_1500m",             # 10
    "t10_relief_norm",            # 11
    "chi_knick_raised",           # 12
    "chi_ksn_anom",               # 13
    "chi_interfluve_scarp",       # 14
    "lahontan_shoreline_idx",     # 15
    "det_elev_hp9",               # 16
    "det_slope_worm_1500m",       # 17
]

GEOTHERM_STATIC_NAMES: List[str] = [
    "mt_clay_cap_feeder",         # 18
    "mt_cond_grad_hp",            # 19
    "mt_cond_worm_1500m",         # 20
    "mt_cond_raw",                # 21
    "hydro_demag_trough",         # 22
    "hydro_shallow_demag_u150",   # 23
    "rad_tc_grad_f32",            # 24
    "rad_sinter_leach_edge",      # 25
    "rad_worm_tc_1500m",          # 26
    "rad_k_th_adularia",          # 27
    "rad_u_th_redox_grad",        # 28
    "gdr_vent_aniso_near",        # 29
    "gdr_vent_aniso_far",         # 30
    "gdr_vent_x_subsurf_edge",    # 31
    "volc_vent_conduit",          # 32
    "geotherm_multiline_idx",     # 33
]

SUBSURF_UPDIP_STATIC_NAMES: List[str] = [
    "verduzco_grav_tilt_edge",    # 34
    "verduzco_grav_worm_2000m",   # 35
    "verduzco_mag_tilt_edge",     # 36
    "verduzco_mag_worm_2000m",    # 37
    "updip60_grav_step",          # 38  <-- H80-A Novel Updip 60-deg Advected Gravity Step
    "updip60_grav_tilt",          # 39  <-- H80-A Novel Updip 60-deg Advected Tilt Zero-Crossing
    "updip60_dbase_grad",         # 40  <-- H80-A Novel Updip 60-deg Advected Basement Fault Step
    "updip60_dipole_asymmetry",   # 41  <-- H80-A Novel Updip vs Downdip Basement Asymmetry
    "updip60_x_piedmont_scarp",   # 42  <-- H80-A Novel Updip Step x Deconfounded Piedmont Scarp
    "crossgrad_parallel_gm",      # 43
    "crossgrad_ortho_gm",         # 44
    "dbase_grad_norm",            # 45
    "dbase_curv_horst",           # 46
    "grav_hp8",                   # 47
    "strain_transtension",        # 48
    "seismic_micro_line",         # 49
]

STATIC_FEATURE_NAMES: List[str] = (
    SCARP_STATIC_NAMES + GEOTHERM_STATIC_NAMES + SUBSURF_UPDIP_STATIC_NAMES
)
assert len(STATIC_FEATURE_NAMES) == 50

FOLD_MECH_NAMES: List[str] = [
    "vis_ray_collinear_short",    # 50
    "vis_ray_collinear_long",     # 51
    "vis_wingcrack_splay",        # 52
    "vis_relay_bridge",           # 53
    "vis_strike_corridor",        # 54
    "vis_dist_capped",            # 55
    "vis_host_len_transfer",      # 56
    "vis_completeness_deficit",   # 57
    "vis_updip_relay_coupler",    # 58  <-- H80-A Fold-Proper Relay x Updip Basement Coupler
]

ALL_BASE_FEATURE_NAMES: List[str] = STATIC_FEATURE_NAMES + FOLD_MECH_NAMES
assert len(ALL_BASE_FEATURE_NAMES) == 59


def _clean_band(arr: np.ndarray, fp: np.ndarray) -> np.ndarray:
    ok = fp & np.isfinite(arr) & (arr > -1e20) & (arr < 1e20)
    med = float(np.median(arr[ok]))
    return np.where(ok, arr, med).astype(np.float32)


def _robust_pos_z(arr: np.ndarray, fp: np.ndarray) -> np.ndarray:
    v = arr[fp & np.isfinite(arr)]
    p50 = float(np.percentile(v, 50.0))
    iqr = float(np.percentile(v, 84.0) - p50) + 1e-6
    z = np.clip((np.nan_to_num(arr, nan=p50) - p50) / iqr, 0.0, 6.0)
    return z.astype(np.float32)


def _grad_mag_f32(arr: np.ndarray) -> np.ndarray:
    gy = np.empty_like(arr, dtype=np.float32)
    gx = np.empty_like(arr, dtype=np.float32)
    gy[1:-1, :] = 0.5 * (arr[2:, :] - arr[:-2, :])
    gy[0, :] = arr[1, :] - arr[0, :]
    gy[-1, :] = arr[-1, :] - arr[-2, :]
    gx[:, 1:-1] = 0.5 * (arr[:, 2:] - arr[:, :-2])
    gx[:, 0] = arr[:, 1] - arr[:, 0]
    gx[:, -1] = arr[:, -1] - arr[:, -2]
    return np.hypot(gy, gx).astype(np.float32)


def _grad_yx_f32(arr: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    gy = np.empty_like(arr, dtype=np.float32)
    gx = np.empty_like(arr, dtype=np.float32)
    gy[1:-1, :] = 0.5 * (arr[2:, :] - arr[:-2, :])
    gy[0, :] = arr[1, :] - arr[0, :]
    gy[-1, :] = arr[-1, :] - arr[-2, :]
    gx[:, 1:-1] = 0.5 * (arr[:, 2:] - arr[:, :-2])
    gx[:, 0] = arr[:, 1] - arr[:, 0]
    gx[:, -1] = arr[:, -1] - arr[:, -2]
    return gy, gx


def _strike_worm(
    arr: np.ndarray,
    length_px: int = 15,
    flank_px: int = 3,
    n_angles: int = 6,
) -> Tuple[np.ndarray, np.ndarray]:
    best = np.full(arr.shape, -1e9, dtype=np.float32)
    best_theta = np.zeros(arr.shape, dtype=np.float32)
    base = np.nan_to_num(arr, nan=0.0).astype(np.float32)
    for k in range(n_angles):
        deg = k * (180.0 / n_angles)
        rot = ndi.rotate(base, deg, reshape=False, order=1, mode="nearest")
        line = ndi.uniform_filter1d(rot, size=length_px, axis=0, mode="nearest")
        flank = 0.5 * (np.roll(line, flank_px, axis=1) + np.roll(line, -flank_px, axis=1))
        w = ndi.rotate((line - flank).astype(np.float32), -deg, reshape=False, order=1, mode="nearest")
        upd = w > best
        best = np.where(upd, w, best)
        best_theta = np.where(upd, np.deg2rad(deg), best_theta)
    return best, best_theta


def updip_advect_field(
    field: np.ndarray,
    uy_downdip: np.ndarray,
    ux_downdip: np.ndarray,
    shift_px: np.ndarray,
) -> np.ndarray:
    """Advect a deep subsurface gradient/tilt field updip to the predicted surface
    rupture trace by sampling `field` at `(y + shift_px * uy_downdip, x + shift_px * ux_downdip)`
    (nearest/bilinear grid lookup, zero float64 overhead)."""
    H, W = field.shape
    yy, xx = np.indices((H, W), dtype=np.float32)
    sy = np.clip(np.rint(yy + shift_px * uy_downdip).astype(np.int32), 0, H - 1)
    sx = np.clip(np.rint(xx + shift_px * ux_downdip).astype(np.int32), 0, W - 1)
    return field[sy, sx].astype(np.float32)


def fill_axial_context_valid(
    out_mat: np.ndarray,
    col_offset: int,
    prob_map: np.ndarray,
    fp: np.ndarray,
    valid_idx: np.ndarray,
) -> None:
    """Write 6 exact 8-azimuth 2D multi-scale axial-ridge channels directly into
    `out_mat[:, col_offset : col_offset + 6]` (`float16`)."""
    p = np.where(fp & np.isfinite(prob_map), prob_map, 0.0).astype(np.float32)
    dirs = (
        (1, 0, 0, 3, 0, 2),
        (2, 1, -1, 2, -1, 2),
        (1, 1, -2, 2, -1, 1),
        (1, 2, -2, 1, -2, 1),
        (0, 1, -3, 0, -2, 0),
        (1, -2, 2, 1, 2, 1),
        (1, -1, 2, 2, 1, 1),
        (2, -1, 1, 2, 1, 2),
    )
    pad = 22
    pp = np.pad(p, pad, mode="edge")
    H, W = p.shape

    best_r9 = np.full((H, W), -1e9, dtype=np.float32)
    best_r21 = np.full((H, W), -1e9, dtype=np.float32)
    best_m21 = np.full((H, W), -1e9, dtype=np.float32)
    best_nms = np.full((H, W), -1e9, dtype=np.float32)
    acc = np.zeros((H + 2 * pad, W + 2 * pad), dtype=np.float32)
    c0, c1 = pad, pad + H
    d0, d1 = pad, pad + W

    for dy, dx, ny3, nx3, ny2, nx2 in dirs:
        step_len = float(np.hypot(dy, dx))
        k9 = max(1, int(round(4.0 / step_len)))
        k21 = max(2, int(round(10.0 / step_len)))

        line9 = pp[c0:c1, d0:d1].copy()
        for s in range(1, k9 + 1):
            line9 += pp[c0 + s * dy : c1 + s * dy, d0 + s * dx : d1 + s * dx]
            line9 += pp[c0 - s * dy : c1 - s * dy, d0 - s * dx : d1 - s * dx]
        line9 *= 1.0 / float(2 * k9 + 1)

        line21 = line9 * float(2 * k9 + 1)
        for s in range(k9 + 1, k21 + 1):
            line21 += pp[c0 + s * dy : c1 + s * dy, d0 + s * dx : d1 + s * dx]
            line21 += pp[c0 - s * dy : c1 - s * dy, d0 - s * dx : d1 - s * dx]
        line21 *= 1.0 / float(2 * k21 + 1)

        acc[c0:c1, d0:d1] = line9
        flank3 = 0.5 * (
            acc[c0 + ny3 : c1 + ny3, d0 + nx3 : d1 + nx3]
            + acc[c0 - ny3 : c1 - ny3, d0 - nx3 : d1 - nx3]
        )
        flank2 = 0.5 * (
            acc[c0 + ny2 : c1 + ny2, d0 + nx2 : d1 + nx2]
            + acc[c0 - ny2 : c1 - ny2, d0 - nx2 : d1 - nx2]
        )
        np.maximum(best_r9, line9 - flank3, out=best_r9)
        np.maximum(best_nms, line9 - flank2, out=best_nms)
        del line9, flank3, flank2

        acc[c0:c1, d0:d1] = line21
        ny4, nx4 = int(round(ny3 * 1.34)), int(round(nx3 * 1.34))
        flank4 = 0.5 * (
            acc[c0 + ny4 : c1 + ny4, d0 + nx4 : d1 + nx4]
            + acc[c0 - ny4 : c1 - ny4, d0 - nx4 : d1 - nx4]
        )
        np.maximum(best_r21, line21 - flank4, out=best_r21)
        np.maximum(best_m21, line21, out=best_m21)
        del line21, flank4

    del pp, acc
    out_mat[:, col_offset + 0] = best_r9.ravel()[valid_idx].astype(np.float16)
    out_mat[:, col_offset + 1] = best_r21.ravel()[valid_idx].astype(np.float16)
    out_mat[:, col_offset + 2] = best_m21.ravel()[valid_idx].astype(np.float16)
    out_mat[:, col_offset + 3] = best_nms.ravel()[valid_idx].astype(np.float16)
    del best_r9, best_r21, best_m21, best_nms
    gc.collect()

    g15 = ndi.gaussian_filter(p, sigma=1.5)
    g8 = ndi.gaussian_filter(p, sigma=8.0)
    out_mat[:, col_offset + 4] = (g15 - g8).ravel()[valid_idx].astype(np.float16)
    del g15, g8
    out_mat[:, col_offset + 5] = ndi.gaussian_filter(p, sigma=25.0).ravel()[valid_idx].astype(np.float16)
    gc.collect()


def _anisotropic_strike_plume(
    impulse: np.ndarray,
    sigma_along: float,
    sigma_across: float,
    strike_deg_cw_from_north: float,
) -> np.ndarray:
    rot = ndi.rotate(impulse, -strike_deg_cw_from_north, reshape=False, order=1, mode="constant")
    filt = ndi.gaussian_filter(rot, sigma=(sigma_along, sigma_across))
    return ndi.rotate(filt, strike_deg_cw_from_north, reshape=False, order=1, mode="constant").astype(
        np.float32
    )


def _rasterize_geothermal_evidence(ext_dir: Path, fp: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Rasterize GDR-1391 well/spring temperatures, quartz/chalcedony/Na-K-Ca geothermometers,
    2-meter shallow temperature probe anomalies (`F2mDAB`), paleo-geothermal sinter/travertine
    deposits, and Quaternary volcanic vents — strictly excluding `dist_known_fault_px`."""
    H, W = fp.shape
    vent_imp = np.zeros((H, W), dtype=np.float32)
    volc_imp = np.zeros((H, W), dtype=np.float32)

    ws = pd.read_csv(ext_dir / "gdr_wellspring_in_footprint.csv")
    temp = pd.to_numeric(ws["temp_c"], errors="coerce").fillna(0.0).to_numpy()
    qtz = pd.to_numeric(ws["geothermquartz_c"], errors="coerce").fillna(0.0).to_numpy()
    chal = pd.to_numeric(ws["geothermchalc_c"], errors="coerce").fillna(0.0).to_numpy()
    cat = pd.to_numeric(ws["geothermcat_c"], errors="coerce").fillna(0.0).to_numpy()
    w_therm = (
        np.clip((temp - 24.0) / 45.0, 0.0, 3.0)
        + 0.45 * np.clip((np.maximum(qtz, chal) - 70.0) / 80.0, 0.0, 2.0)
        + 0.35 * np.clip((cat - 85.0) / 90.0, 0.0, 1.8)
    )
    sel = w_therm > 0.05
    rr = np.clip(ws["row"].to_numpy()[sel].astype(int), 0, H - 1)
    cc = np.clip(ws["col"].to_numpy()[sel].astype(int), 0, W - 1)
    np.add.at(vent_imp, (rr, cc), w_therm[sel].astype(np.float32))

    z2m = ext_dir / "2m_temperature_probe_INGENIOUS_regional_data.zip"
    if z2m.exists():
        with zipfile.ZipFile(z2m) as zf:
            dbf = zf.read(
                "2m_temperature_probe_INGENIOUS_regional_data/"
                "2m_temperature_probe_n83geo/2m_temperature_probe_n83geo.dbf"
            )
        nrec, hdrlen, reclen = struct.unpack("<IHH", dbf[4:12])
        fields = []
        off = 32
        while dbf[off] != 0x0D:
            fn = dbf[off : off + 11].split(b"\x00")[0].decode("latin1")
            fl = dbf[off + 16]
            fields.append((fn, fl))
            off += 32
        lons, lats, anoms = [], [], []
        pos = hdrlen
        for _ in range(nrec):
            if dbf[pos] == 0x2A:
                pos += reclen
                continue
            p = pos + 1
            row_d = {}
            for fn, fl in fields:
                if fn in ("LatDegreeN", "LongDegree", "F2mDAB"):
                    row_d[fn] = dbf[p : p + fl].decode("latin1").strip()
                p += fl
            pos += reclen
            try:
                lat = float(row_d["LatDegreeN"])
                lon = float(row_d["LongDegree"])
                anom = float(row_d["F2mDAB"])
            except (KeyError, ValueError):
                continue
            if np.isfinite(lat) and np.isfinite(lon) and np.isfinite(anom) and anom >= 1.5:
                lons.append(lon)
                lats.append(lat)
                anoms.append(min(anom / 3.0, 4.0))
        if lons:
            xs, ys = warp_transform("EPSG:4269", "EPSG:32611", lons, lats)
            r2 = ((ORIGIN_Y - np.asarray(ys)) // PIXEL_M).astype(int)
            c2 = ((np.asarray(xs) - ORIGIN_X) // PIXEL_M).astype(int)
            inb = (r2 >= 0) & (r2 < H) & (c2 >= 0) & (c2 < W)
            np.add.at(vent_imp, (r2[inb], c2[inb]), np.asarray(anoms, dtype=np.float32)[inb])

    zpal = ext_dir / "paleo_geothermal_regional.zip"
    if zpal.exists():
        with zipfile.ZipFile(zpal) as zf:
            txt = zf.read("Paleo_geothermal_final.csv").decode("cp1252")
        dfp = pd.read_csv(io.StringIO(txt))
        lons = pd.to_numeric(dfp["Long_deg_NAD83"], errors="coerce").to_numpy()
        lats = pd.to_numeric(dfp["Lat_deg_NAD83"], errors="coerce").to_numpy()
        okp = np.isfinite(lons) & np.isfinite(lats)
        if np.any(okp):
            xs, ys = warp_transform("EPSG:4269", "EPSG:32611", lons[okp].tolist(), lats[okp].tolist())
            rp = ((ORIGIN_Y - np.asarray(ys)) // PIXEL_M).astype(int)
            cp = ((np.asarray(xs) - ORIGIN_X) // PIXEL_M).astype(int)
            inb = (rp >= 0) & (rp < H) & (cp >= 0) & (cp < W)
            np.add.at(vent_imp, (rp[inb], cp[inb]), 2.5)

    vv_path = ext_dir / "gdr_volcanic_vents_in_footprint.csv"
    if vv_path.exists():
        vv = pd.read_csv(vv_path)
        rv = np.clip(vv["row"].to_numpy().astype(int), 0, H - 1)
        cv = np.clip(vv["col"].to_numpy().astype(int), 0, W - 1)
        np.maximum.at(volc_imp, (rv, cv), 1.5)

    return vent_imp, volc_imp


@dataclass
class StaticStack:
    valid_idx: np.ndarray          # flat indices of `fp` (`5,167,373` pixels)
    X_mmap: np.memmap              # shape `(len(valid_idx), 50)`, dtype `float16`
    feature_names: List[str]
    slope10_full: np.ndarray       # shape `(3730, 3292)` float32
    lidar_valid_full: np.ndarray   # shape `(3730, 3292)` bool
    best_strike_full: np.ndarray   # shape `(3730, 3292)` float32
    tect_bg_full: np.ndarray       # shape `(3730, 3292)` float32
    updip_energy_full: np.ndarray  # shape `(3730, 3292)` float32 (updip 60-deg basement-step energy)


def _put_f16(X_mm: np.memmap, col: int, arr2d: np.ndarray, valid_idx: np.ndarray) -> None:
    vec = arr2d.ravel()[valid_idx].astype(np.float32)
    vec = np.where(np.isfinite(vec), np.clip(vec, -60000.0, 60000.0), np.nan)
    X_mm[:, col] = vec.astype(np.float16)


def _compute_scarp_domain(
    X_mm: np.memmap, data_dir: Path, ext_dir: Path, fp: np.ndarray, valid_idx: np.ndarray
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Compute columns 0..17 inside isolated scope with TIGER road & BLM mine deconfounding."""
    # Anthropogenic deconfounding gate from US Census TIGER roads and BLM mining claims
    anthro_gate = np.ones(fp.shape, dtype=np.float32)
    road_p = ext_dir / "tiger_road_distance_m.tif"
    if road_p.exists():
        with rasterio.open(road_p) as sr:
            d_road = np.maximum(sr.read(1).astype(np.float32), 0.0)
        anthro_gate *= (1.0 - 0.55 * np.exp(-((d_road / 160.0) ** 2))).astype(np.float32)
        del d_road
    mine_p = ext_dir / "blm_closed_claim_distance_m.tif"
    if mine_p.exists():
        with rasterio.open(mine_p) as sm:
            d_mine = np.maximum(sm.read(1).astype(np.float32), 0.0)
        anthro_gate *= (1.0 - 0.30 * np.exp(-((d_mine / 140.0) ** 2))).astype(np.float32)
        del d_mine

    with rasterio.open(ext_dir / "lidar_scarp_features_u8.tif") as s:
        lid_valid = (s.read(12) > 0) & fp

        def read_lid(band: int) -> np.ndarray:
            lo, hi = LIDAR_SCALES[band]
            return lo + (s.read(band).astype(np.float32) / 255.0) * (hi - lo)

        l_step = read_lid(3)
        l_lapn = read_lid(4)
        l_lapp = read_lid(5)
        l_down = read_lid(6)
        l_up = read_lid(7)
        l_cross = read_lid(8)
        l_relief = read_lid(9)
        l_coh = read_lid(10)
        l_strike = read_lid(11)

    nan_mask = ~lid_valid
    lid_piedmont = (l_step / np.sqrt(l_relief + 9.0)) * anthro_gate
    lid_w, _ = _strike_worm(np.where(lid_valid, lid_piedmont * (0.4 + l_coh), 0.0), length_px=11, flank_px=3)

    _put_f16(X_mm, 0, np.where(nan_mask, np.nan, lid_piedmont), valid_idx)
    _put_f16(X_mm, 1, np.where(nan_mask, np.nan, np.sqrt(np.maximum(l_lapn * l_lapp, 0.0)) * anthro_gate), valid_idx)
    _put_f16(
        X_mm,
        2,
        np.where(nan_mask, np.nan, np.maximum(l_down, 1.35 * l_up) * (0.30 + l_coh) * anthro_gate / (l_cross + 0.08)),
        valid_idx,
    )
    _put_f16(X_mm, 3, np.where(nan_mask, np.nan, l_up * (0.25 + l_coh) * anthro_gate / (l_cross + 0.10)), valid_idx)
    _put_f16(
        X_mm,
        4,
        np.where(nan_mask, np.nan, (l_step - ndi.uniform_filter(np.where(lid_valid, l_step, 0.0), size=11)) * anthro_gate),
        valid_idx,
    )
    _put_f16(X_mm, 5, np.where(nan_mask, np.nan, l_coh), valid_idx)
    _put_f16(X_mm, 6, np.where(nan_mask, np.nan, lid_w), valid_idx)
    del l_step, l_lapn, l_lapp, l_down, l_up, l_cross, l_relief, l_coh, lid_piedmont, lid_w
    gc.collect()

    with rasterio.open(ext_dir / "topo_10m_derived_100m_u8.tif") as s:
        t_slp90 = _clean_band(s.read(1).astype(np.float32), fp)
        t_slp_std = _clean_band(s.read(4).astype(np.float32), fp)
        t_rel = _clean_band(s.read(5).astype(np.float32), fp)
        t_curv = _clean_band(s.read(6).astype(np.float32), fp)
        t_asp_coh = _clean_band(s.read(7).astype(np.float32), fp) / 255.0
        t_line = _clean_band(s.read(8).astype(np.float32), fp)

    slp_bg13 = ndi.uniform_filter(t_slp90, size=13)
    t10_slp_anom = ((t_slp90 - slp_bg13) / (slp_bg13 + 12.0)) * anthro_gate
    t10_w, strike_t10 = _strike_worm(t10_slp_anom, length_px=15, flank_px=3)

    piedmont_break_10m = ((0.45 * t_slp_std + 0.35 * t_curv + 0.20 * t_line) / np.sqrt(t_rel + 12.0)) * anthro_gate
    _put_f16(X_mm, 7, piedmont_break_10m, valid_idx)
    _put_f16(X_mm, 8, t_asp_coh * t_curv / np.sqrt(t_rel + 12.0), valid_idx)
    _put_f16(X_mm, 9, t10_slp_anom, valid_idx)
    _put_f16(X_mm, 10, t10_w, valid_idx)
    _put_f16(X_mm, 11, np.log1p(np.maximum(t_rel, 0.0)), valid_idx)
    del t_slp_std, t_rel, t_curv, t_line, slp_bg13
    gc.collect()

    with rasterio.open(ext_dir / "dem10_chi_100m_u8.tif") as s:
        knick_n = s.read(2).astype(np.float32)
        raised_n = s.read(3).astype(np.float32)
    with rasterio.open(ext_dir / "dem10_chi_100m_i16.tif") as s:
        ksn = _clean_band(s.read(3).astype(np.float32), fp)
        log_area = _clean_band(s.read(4).astype(np.float32) / 1000.0, fp)
        zmean_m = _clean_band(s.read(5).astype(np.float32), fp)

    _put_f16(X_mm, 12, np.log1p(np.maximum(knick_n + 2.0 * raised_n, 0.0)), valid_idx)
    _put_f16(X_mm, 13, ksn - ndi.uniform_filter(ksn, size=17), valid_idx)
    _put_f16(X_mm, 14, t10_slp_anom * np.clip(4.2 - log_area, 0.0, 3.0), valid_idx)
    lahontan = (
        np.exp(-(((zmean_m - 1334.0) / 16.0) ** 2))
        + 0.75 * np.exp(-(((zmean_m - 1298.0) / 14.0) ** 2))
        + 0.65 * np.exp(-(((zmean_m - 1212.0) / 15.0) ** 2))
    ) * t_asp_coh
    _put_f16(X_mm, 15, lahontan, valid_idx)
    del knick_n, raised_n, ksn, log_area, zmean_m, lahontan, t10_slp_anom, t_asp_coh, anthro_gate
    gc.collect()

    with rasterio.open(data_dir / "training_features.tif") as src:
        b_det_elev = _clean_band(src.read(12), fp)
        b_det_slp = _clean_band(src.read(19), fp)
    det_slp_bg = ndi.uniform_filter(b_det_slp, size=15)
    det_slp_ratio = (b_det_slp - det_slp_bg) / (det_slp_bg + 2.5)
    det_slp_w, _ = _strike_worm(det_slp_ratio, length_px=15, flank_px=3)
    _put_f16(X_mm, 16, b_det_elev - ndi.gaussian_filter(b_det_elev, sigma=9.0), valid_idx)
    _put_f16(X_mm, 17, det_slp_w, valid_idx)
    X_mm.flush()

    t10_z = _robust_pos_z(t10_w, fp)
    pied_z = _robust_pos_z(piedmont_break_10m, fp)
    del b_det_elev, b_det_slp, det_slp_bg, det_slp_ratio, det_slp_w, piedmont_break_10m
    gc.collect()
    return t_slp90, lid_valid, l_strike, strike_t10, t10_z, pied_z


def _compute_geotherm_and_subsurf_domains(
    X_mm: np.memmap,
    data_dir: Path,
    ext_dir: Path,
    fp: np.ndarray,
    valid_idx: np.ndarray,
    lid_valid: np.ndarray,
    l_strike: np.ndarray,
    strike_t10: np.ndarray,
    t10_z: np.ndarray,
    pied_z: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute columns 18..49 including H80-A's Updip 60-Deg Basement-to-Surface Advection."""
    with rasterio.open(data_dir / "training_features.tif") as src:
        b_rtp = _clean_band(src.read(2), fp)
        b_tmi_hg = _clean_band(src.read(3), fp)
        b_grav_slp = _clean_band(src.read(5), fp)
        b_rad_tc = _clean_band(src.read(6), fp)
        b_tmi_vg = _clean_band(src.read(9), fp)
        b_grav_vg = _clean_band(src.read(11), fp)
        b_dbase = _clean_band(src.read(15), fp)
        b_cond = _clean_band(src.read(17), fp)

    # Verduzco–Miller–Singh depth-invariant tilt zero-crossing edges
    tdr_grav = np.arctan2(b_grav_vg, np.maximum(b_grav_slp, 1e-4)).astype(np.float32)
    egrav_tilt = (_grad_mag_f32(tdr_grav) * np.exp(-((tdr_grav / 0.5) ** 2))).astype(np.float32)
    del tdr_grav, b_grav_vg
    w_egrav, strike_grav = _strike_worm(egrav_tilt, length_px=20, flank_px=4)

    tdr_mag = np.arctan2(b_tmi_vg, np.maximum(b_tmi_hg, 1e-4)).astype(np.float32)
    emag_tilt = (_grad_mag_f32(tdr_mag) * np.exp(-((tdr_mag / 0.5) ** 2))).astype(np.float32)
    del tdr_mag, b_tmi_vg
    w_emag, _ = _strike_worm(emag_tilt, length_px=20, flank_px=4)

    _put_f16(X_mm, 34, egrav_tilt, valid_idx)
    _put_f16(X_mm, 35, w_egrav, valid_idx)
    _put_f16(X_mm, 36, emag_tilt, valid_idx)
    _put_f16(X_mm, 37, w_emag, valid_idx)
    del emag_tilt, w_emag
    gc.collect()

    # H80-A Core Innovation: Updip 60-Deg Normal-Fault Basement-to-Surface Advection
    db_smooth = ndi.gaussian_filter(np.maximum(b_dbase, 0.0), sigma=3.0)
    dby, dbx = _grad_yx_f32(db_smooth)
    db_grad_mag = np.hypot(dby, dbx) + 1e-5
    # Unit vector pointing down-dip (into deeper basin fill)
    uy_down = (dby / db_grad_mag).astype(np.float32)
    ux_down = (dbx / db_grad_mag).astype(np.float32)
    # Horizontal offset from mid-basement scarp (0.5 * z_b * cot(60 deg)) in 100m pixels
    shift_px = np.clip((0.5 * db_smooth * COT_60_DEG) / PIXEL_M, 0.5, 12.0).astype(np.float32)
    db_step_norm = (db_grad_mag / np.sqrt(db_smooth + 100.0)).astype(np.float32)

    updip_grav_step = updip_advect_field(b_grav_slp, uy_down, ux_down, shift_px)
    updip_grav_tilt = updip_advect_field(egrav_tilt, uy_down, ux_down, shift_px)
    updip_dbase_grad = updip_advect_field(db_step_norm, uy_down, ux_down, shift_px)
    downdip_dbase_grad = updip_advect_field(db_step_norm, -uy_down, -ux_down, shift_px)
    updip_dipole_asym = (updip_dbase_grad - downdip_dbase_grad).astype(np.float32)

    updip_energy = (
        0.40 * _robust_pos_z(updip_grav_step, fp)
        + 0.35 * _robust_pos_z(updip_grav_tilt, fp)
        + 0.25 * _robust_pos_z(updip_dbase_grad, fp)
    ).astype(np.float32)

    _put_f16(X_mm, 38, updip_grav_step, valid_idx)
    _put_f16(X_mm, 39, updip_grav_tilt, valid_idx)
    _put_f16(X_mm, 40, updip_dbase_grad, valid_idx)
    _put_f16(X_mm, 41, updip_dipole_asym, valid_idx)
    _put_f16(X_mm, 42, updip_energy * pied_z, valid_idx)
    _put_f16(X_mm, 45, db_step_norm, valid_idx)
    _put_f16(X_mm, 46, -ndi.laplace(db_smooth), valid_idx)
    del db_smooth, dby, dbx, db_grad_mag, uy_down, ux_down, shift_px
    del updip_grav_step, updip_grav_tilt, updip_dbase_grad, downdip_dbase_grad, updip_dipole_asym, db_step_norm, b_dbase
    gc.collect()

    # Geothermal Domain (`18..33`): MT clay-cap feeder, shallow demagnetization, radiometrics, GDR vents
    cond_grad = _grad_mag_f32(ndi.gaussian_filter(b_cond, sigma=2.0))
    mt_feeder = np.maximum(-(b_cond - ndi.uniform_filter(b_cond, size=31)), 0.0) * cond_grad
    w_mt, _ = _strike_worm(cond_grad, length_px=15, flank_px=3)
    _put_f16(X_mm, 18, mt_feeder, valid_idx)
    _put_f16(X_mm, 19, cond_grad - ndi.uniform_filter(cond_grad, size=21), valid_idx)
    _put_f16(X_mm, 20, w_mt, valid_idx)
    _put_f16(X_mm, 21, b_cond, valid_idx)
    del w_mt, b_cond
    gc.collect()

    with rasterio.open(ext_dir / "geodawn_extensions_u8.tif") as s:
        ext_k_th = _clean_band(s.read(1).astype(np.float32), fp)
        ext_u_th = _clean_band(s.read(3).astype(np.float32), fp)
        ext_up150 = _clean_band(s.read(4).astype(np.float32), fp)

    demag_trough = np.maximum(-(b_rtp - ndi.uniform_filter(b_rtp, size=17)), 0.0) * np.log1p(
        np.maximum(b_tmi_hg, 0.0)
    )
    z_u150 = (ext_up150 - ndi.uniform_filter(ext_up150, 25)) / (float(np.std(ext_up150[fp])) + 1e-5)
    z_rtp = (b_rtp - ndi.uniform_filter(b_rtp, 25)) / (float(np.std(b_rtp[fp])) + 1e-5)
    _put_f16(X_mm, 22, demag_trough, valid_idx)
    _put_f16(X_mm, 23, z_u150 - z_rtp, valid_idx)
    del ext_up150, z_u150, z_rtp

    rad_grad_f32 = _grad_mag_f32(ndi.gaussian_filter(b_rad_tc, sigma=1.2))
    w_rad, _ = _strike_worm(rad_grad_f32, length_px=15, flank_px=3)
    _put_f16(X_mm, 24, rad_grad_f32, valid_idx)
    _put_f16(X_mm, 25, np.maximum(-(b_rad_tc - ndi.uniform_filter(b_rad_tc, size=21)), 0.0) * rad_grad_f32, valid_idx)
    _put_f16(X_mm, 26, w_rad, valid_idx)
    _put_f16(X_mm, 27, ext_k_th - ndi.uniform_filter(ext_k_th, size=21), valid_idx)
    _put_f16(X_mm, 28, _grad_mag_f32(ndi.gaussian_filter(ext_u_th, sigma=1.5)), valid_idx)
    del b_rad_tc, w_rad, ext_k_th, ext_u_th
    gc.collect()

    vent_imp, volc_imp = _rasterize_geothermal_evidence(ext_dir, fp)
    vent_near = 0.68 * _anisotropic_strike_plume(vent_imp, 14.0, 3.5, 22.0) + 0.32 * _anisotropic_strike_plume(
        vent_imp, 14.0, 3.5, -35.0
    )
    vent_far = 0.68 * _anisotropic_strike_plume(vent_imp, 40.0, 12.0, 22.0) + 0.32 * _anisotropic_strike_plume(
        vent_imp, 40.0, 12.0, -35.0
    )
    volc_conduit = 0.70 * _anisotropic_strike_plume(volc_imp, 18.0, 4.0, 20.0) + 0.30 * _anisotropic_strike_plume(
        volc_imp, 18.0, 4.0, -35.0
    )
    del vent_imp, volc_imp

    subsurf_edge_z = (
        0.35 * _robust_pos_z(cond_grad, fp)
        + 0.40 * updip_energy
        + 0.25 * _robust_pos_z(rad_grad_f32, fp)
    )
    geotherm_multi = (
        (0.65 * _robust_pos_z(vent_near, fp) + 0.35 * _robust_pos_z(vent_far, fp))
        * (0.45 * _robust_pos_z(mt_feeder, fp) + 0.35 * _robust_pos_z(demag_trough, fp) + 0.20 * updip_energy + 0.25)
    ).astype(np.float32)

    _put_f16(X_mm, 29, vent_near, valid_idx)
    _put_f16(X_mm, 30, vent_far, valid_idx)
    _put_f16(X_mm, 31, vent_near * subsurf_edge_z, valid_idx)
    _put_f16(X_mm, 32, volc_conduit, valid_idx)
    _put_f16(X_mm, 33, geotherm_multi, valid_idx)
    vent_far_z = _robust_pos_z(vent_far, fp)
    del cond_grad, egrav_tilt, rad_grad_f32, mt_feeder, demag_trough, vent_near, vent_far, volc_conduit, subsurf_edge_z, geotherm_multi
    gc.collect()

    with rasterio.open(data_dir / "training_features.tif") as src:
        b_geod_2inv = _clean_band(src.read(4), fp)
        b_shear = _clean_band(src.read(7), fp)
        b_dilate = _clean_band(src.read(8), fp)
        b_deq = _clean_band(src.read(10), fp)
        b_grav = _clean_band(src.read(13), fp)
        b_ieq = _clean_band(src.read(16), fp)

    gs = ndi.gaussian_filter(b_grav, sigma=2.0)
    ms = ndi.gaussian_filter(b_rtp, sigma=2.0)
    gg_y, gg_x = _grad_yx_f32(gs)
    gm_y, gm_x = _grad_yx_f32(ms)
    ng = np.hypot(gg_y, gg_x) + 1e-6
    nm = np.hypot(gm_y, gm_x) + 1e-6
    amp_gm = np.sqrt(ng * nm)
    _put_f16(X_mm, 43, (np.abs(gg_y * gm_y + gg_x * gm_x) / (ng * nm)) * amp_gm, valid_idx)
    _put_f16(X_mm, 44, (np.abs(gg_x * gm_y - gg_y * gm_x) / (ng * nm)) * amp_gm, valid_idx)
    _put_f16(X_mm, 47, b_grav - ndi.gaussian_filter(b_grav, sigma=8.0), valid_idx)
    del gs, ms, gg_y, gg_x, gm_y, gm_x, ng, nm, amp_gm, b_rtp, b_grav, b_grav_slp, b_tmi_hg
    gc.collect()

    strain_tt = (
        0.45 * _robust_pos_z(b_shear - ndi.gaussian_filter(b_shear, 35.0), fp)
        + 0.35 * _robust_pos_z(b_dilate - ndi.gaussian_filter(b_dilate, 35.0), fp)
        + 0.20 * _robust_pos_z(b_geod_2inv, fp)
    ).astype(np.float32)
    leq = np.log1p(np.maximum(b_ieq, 0.0))
    seismic_line = (leq - ndi.uniform_filter(leq, size=25)) - 0.25 * np.log1p(np.maximum(b_deq, 0.0))
    _put_f16(X_mm, 48, strain_tt, valid_idx)
    _put_f16(X_mm, 49, seismic_line, valid_idx)
    X_mm.flush()

    best_strike = np.where(lid_valid, l_strike, 0.65 * strike_t10 + 0.35 * strike_grav).astype(np.float32)
    tect_bg = ndi.gaussian_filter(
        0.30 * t10_z + 0.25 * _robust_pos_z(w_egrav, fp) + 0.25 * updip_energy + 0.20 * (strain_tt + 0.5 * vent_far_z),
        sigma=25.0,
    ).astype(np.float32)
    return best_strike, tect_bg, updip_energy


def build_or_load_static_stack(
    data_dir: Path,
    ext_dir: Path,
    cache_dir: Path,
    fp: np.ndarray,
    force_rebuild: bool = False,
) -> StaticStack:
    cache_dir.mkdir(parents=True, exist_ok=True)
    mmap_path = cache_dir / "gems80_static_v1.f16.mmap"
    aux_path = cache_dir / "gems80_static_v1_aux.npz"
    meta_path = cache_dir / "gems80_static_v1_meta.json"

    valid_idx = np.flatnonzero(fp).astype(np.int64)
    n_pix = int(valid_idx.size)
    n_feat = len(STATIC_FEATURE_NAMES)

    if (
        not force_rebuild
        and mmap_path.exists()
        and aux_path.exists()
        and meta_path.exists()
        and json.loads(meta_path.read_text()).get("feature_names") == STATIC_FEATURE_NAMES
    ):
        X_mm = np.memmap(mmap_path, dtype=np.float16, mode="r", shape=(n_pix, n_feat))
        aux = np.load(aux_path)
        return StaticStack(
            valid_idx=valid_idx,
            X_mmap=X_mm,
            feature_names=list(STATIC_FEATURE_NAMES),
            slope10_full=aux["slope10_full"],
            lidar_valid_full=aux["lidar_valid_full"],
            best_strike_full=aux["best_strike_full"],
            tect_bg_full=aux["tect_bg_full"],
            updip_energy_full=aux["updip_energy_full"],
        )

    X_mm = np.memmap(mmap_path, dtype=np.float16, mode="w+", shape=(n_pix, n_feat))
    t_slp90, lid_valid, l_strike, strike_t10, t10_z, pied_z = _compute_scarp_domain(
        X_mm, data_dir, ext_dir, fp, valid_idx
    )
    gc.collect()
    best_strike, tect_bg, updip_energy = _compute_geotherm_and_subsurf_domains(
        X_mm, data_dir, ext_dir, fp, valid_idx, lid_valid, l_strike, strike_t10, t10_z, pied_z
    )
    del l_strike, strike_t10, t10_z, pied_z
    del X_mm
    gc.collect()

    np.savez_compressed(
        aux_path,
        slope10_full=t_slp90.astype(np.float32),
        lidar_valid_full=lid_valid,
        best_strike_full=best_strike,
        tect_bg_full=tect_bg,
        updip_energy_full=updip_energy,
    )
    meta_path.write_text(
        json.dumps(
            {
                "version": "gems80_static_v1_f16",
                "n_pixels": n_pix,
                "n_features": n_feat,
                "feature_names": STATIC_FEATURE_NAMES,
            },
            indent=2,
        )
    )
    X_mm_ro = np.memmap(mmap_path, dtype=np.float16, mode="r", shape=(n_pix, n_feat))
    return StaticStack(
        valid_idx=valid_idx,
        X_mmap=X_mm_ro,
        feature_names=list(STATIC_FEATURE_NAMES),
        slope10_full=t_slp90.astype(np.float32),
        lidar_valid_full=lid_valid,
        best_strike_full=best_strike,
        tect_bg_full=tect_bg,
        updip_energy_full=updip_energy,
    )


def compute_fold_mechanics_2d(
    vis_mask: np.ndarray,
    fp: np.ndarray,
    best_strike_full: np.ndarray,
    tect_bg_full: np.ndarray,
    updip_energy_full: np.ndarray,
) -> np.ndarray:
    """Compute all 9 fold-proper fracture-transfer features strictly from `vis_mask`."""
    H, W = vis_mask.shape
    vis = vis_mask & fp

    comp_lbl, _ = ndi.label(vis, structure=np.ones((3, 3), dtype=int))
    comp_sizes = np.bincount(comp_lbl.ravel()).astype(np.float32)
    comp_sizes[0] = 0.0

    skel = skeletonize(vis)
    sr, sc = np.nonzero(skel)
    if len(sr) == 0:
        return np.zeros((H, W, len(FOLD_MECH_NAMES)), dtype=np.float32)

    dist_px, (nr_idx, nc_idx) = ndi.distance_transform_edt(~vis, return_indices=True)
    dist_f32 = dist_px.astype(np.float32)
    host_len = comp_sizes[comp_lbl[nr_idx, nc_idx]]
    del dist_px, nr_idx, nc_idx

    skel_tree = cKDTree(np.column_stack([sr, sc]).astype(np.float32))
    skel_comp = comp_lbl[sr, sc]

    deg = ndi.convolve(skel.astype(np.uint8), np.ones((3, 3), dtype=np.uint8), mode="constant")
    tip_mask = skel & (deg <= 2)
    tr, tc = np.nonzero(tip_mask)
    tip_comp = comp_lbl[tr, tc]
    keep = comp_sizes[tip_comp] >= 4.0
    tr, tc, tip_comp = tr[keep], tc[keep], tip_comp[keep]

    ray_short = np.zeros((H, W), dtype=np.float32)
    ray_long = np.zeros((H, W), dtype=np.float32)
    splay = np.zeros((H, W), dtype=np.float32)
    bridge = np.zeros((H, W), dtype=np.float32)
    corridor = np.zeros((H, W), dtype=np.float32)

    R = 44
    yy, xx = np.mgrid[-R : R + 1, -R : R + 1].astype(np.float32)
    rr_d = np.hypot(yy, xx) + 1e-6
    uy = yy / rr_d
    ux = xx / rr_d

    tip_ty = np.zeros(len(tr), dtype=np.float32)
    tip_tx = np.zeros(len(tr), dtype=np.float32)
    tip_w = np.log1p(comp_sizes[tip_comp]).astype(np.float32)

    nbr_lists = skel_tree.query_ball_point(np.column_stack([tr, tc]), r=11.0)
    valid_tips = []
    for i, nbrs in enumerate(nbr_lists):
        nbrs_arr = np.asarray(nbrs, dtype=np.int64)
        same = nbrs_arr[skel_comp[nbrs_arr] == tip_comp[i]]
        dy = tr[i] - sr[same]
        dx = tc[i] - sc[same]
        dd = np.hypot(dy, dx)
        ring = same[(dd >= 3.0) & (dd <= 11.0)]
        if len(ring) == 0:
            ring = same[dd >= 1.0]
        if len(ring) == 0:
            continue
        vy = float(tr[i] - np.mean(sr[ring]))
        vx = float(tc[i] - np.mean(sc[ring]))
        vn = float(np.hypot(vy, vx))
        if vn < 1e-4:
            continue
        tip_ty[i] = vy / vn
        tip_tx[i] = vx / vn
        valid_tips.append(i)

    for i in valid_tips:
        r0, c0 = int(tr[i]), int(tc[i])
        y0, y1 = max(0, r0 - R), min(H, r0 + R + 1)
        x0, x1 = max(0, c0 - R), min(W, c0 + R + 1)
        py0, py1 = y0 - (r0 - R), y1 - (r0 - R)
        px0, px1 = x0 - (c0 - R), x1 - (c0 - R)

        cos_a = uy[py0:py1, px0:px1] * tip_ty[i] + ux[py0:py1, px0:px1] * tip_tx[i]
        d = rr_d[py0:py1, px0:px1]
        w = tip_w[i]

        np.maximum(
            ray_short[y0:y1, x0:x1],
            w * np.exp(-d / 11.0) * np.exp(-(((1.0 - cos_a) / 0.055) ** 2)),
            out=ray_short[y0:y1, x0:x1],
        )
        np.maximum(
            ray_long[y0:y1, x0:x1],
            w * np.exp(-d / 26.0) * np.exp(-(((1.0 - cos_a) / 0.075) ** 2)),
            out=ray_long[y0:y1, x0:x1],
        )
        np.maximum(
            splay[y0:y1, x0:x1],
            w * np.exp(-d / 15.0) * np.exp(-(((0.85 - cos_a) / 0.07) ** 2)),
            out=splay[y0:y1, x0:x1],
        )
        np.maximum(
            corridor[y0:y1, x0:x1],
            w * np.exp(-d / 12.0) * (np.abs(cos_a) ** 3),
            out=corridor[y0:y1, x0:x1],
        )

    if len(valid_tips) >= 2:
        vt = np.asarray(valid_tips, dtype=np.int64)
        pts = np.column_stack([tr[vt], tc[vt]]).astype(np.float32)
        pairs = cKDTree(pts).query_pairs(r=48.0, output_type="ndarray")
        for ia, ib in pairs:
            a, b = int(vt[ia]), int(vt[ib])
            if tip_comp[a] == tip_comp[b]:
                continue
            dy = float(tr[b] - tr[a])
            dx = float(tc[b] - tc[a])
            gap = float(np.hypot(dy, dx))
            if gap < 3.0:
                continue
            gy, gx = dy / gap, dx / gap
            dot_a = tip_ty[a] * gy + tip_tx[a] * gx
            dot_b = -(tip_ty[b] * gy + tip_tx[b] * gx)
            anti = -(tip_ty[a] * tip_ty[b] + tip_tx[a] * tip_tx[b])
            align = max(0.0, min(dot_a, dot_b))
            relay_score = 0.0
            if align > 0.30:
                relay_score = align * (tip_w[a] + tip_w[b]) * 0.5
            elif anti > 0.70 and dot_a > -0.25 and dot_b > -0.25 and gap <= 26.0:
                relay_score = 0.75 * anti * (tip_w[a] + tip_w[b]) * 0.5
            if relay_score <= 0.05:
                continue

            pad = 7
            y0 = max(0, min(int(tr[a]), int(tr[b])) - pad)
            y1 = min(H, max(int(tr[a]), int(tr[b])) + pad + 1)
            x0 = max(0, min(int(tc[a]), int(tc[b])) - pad)
            x1 = min(W, max(int(tc[a]), int(tc[b])) + pad + 1)
            gyy, gxx = np.mgrid[y0:y1, x0:x1].astype(np.float32)
            proj = np.clip(((gyy - tr[a]) * gy + (gxx - tc[a]) * gx) / gap, 0.0, 1.0)
            perp = np.hypot(gyy - (tr[a] + proj * dy), gxx - (tc[a] + proj * dx))
            tube = relay_score * np.exp(-((perp / 4.2) ** 2)) * np.exp(-gap / 36.0)
            np.maximum(bridge[y0:y1, x0:x1], tube, out=bridge[y0:y1, x0:x1])

    vis_dens = ndi.gaussian_filter(vis.astype(np.float32), sigma=28.0)
    scale_bg = float(np.percentile(vis_dens[fp], 90.0)) / (
        float(np.percentile(tect_bg_full[fp], 90.0)) + 1e-6
    )
    deficit = np.maximum(tect_bg_full * scale_bg - vis_dens, 0.0).astype(np.float32)
    updip_relay = ((ray_short + 0.7 * bridge + 0.5 * splay) * (0.25 + updip_energy_full)).astype(np.float32)

    return np.stack(
        [
            ray_short,
            ray_long,
            splay,
            bridge,
            corridor,
            np.minimum(dist_f32, 45.0).astype(np.float32),
            (np.log1p(host_len) * np.exp(-dist_f32 / 18.0)).astype(np.float32),
            deficit,
            updip_relay,
        ],
        axis=-1,
    )


def run_leakage_canary_audit(
    stack: StaticStack,
    folds: List[dict],
    fp: np.ndarray,
    max_allowed_auc: float = 0.90,
    seed: int = 80,
) -> dict:
    """Test every single feature (all 50 static + 9 fold-visible mechanics = 59 features)
    alone on the spatially-blocked holdout before trusting any model result (`Protocol Rule 4`)."""
    rng = np.random.default_rng(seed)
    inv_map = np.full(fp.size, -1, dtype=np.int64)
    inv_map[stack.valid_idx] = np.arange(stack.valid_idx.size, dtype=np.int64)

    pos_flat_parts = []
    neg_flat_parts = []
    mech_pos_parts = []
    mech_neg_parts = []

    for f in folds:
        reg = f["region"]
        truth = f["truth"] & reg
        neg_pool = reg & (~f["visible"]) & (~f["hidden_all"])
        p_idx = np.flatnonzero(truth)
        n_pool = np.flatnonzero(neg_pool)
        n_idx = rng.choice(n_pool, size=min(len(n_pool), 4 * len(p_idx)), replace=False)

        mech2d = compute_fold_mechanics_2d(
            f["visible"], fp, stack.best_strike_full, stack.tect_bg_full, stack.updip_energy_full
        ).reshape(-1, len(FOLD_MECH_NAMES))
        pos_flat_parts.append(p_idx)
        neg_flat_parts.append(n_idx)
        mech_pos_parts.append(mech2d[p_idx])
        mech_neg_parts.append(mech2d[n_idx])
        del mech2d
        gc.collect()

    pos_flat = np.concatenate(pos_flat_parts)
    neg_flat = np.concatenate(neg_flat_parts)
    y = np.concatenate([np.ones(len(pos_flat)), np.zeros(len(neg_flat))])

    rows_pos = inv_map[pos_flat]
    rows_neg = inv_map[neg_flat]
    X_pos = np.hstack([np.asarray(stack.X_mmap[rows_pos], dtype=np.float32), np.vstack(mech_pos_parts)])
    X_neg = np.hstack([np.asarray(stack.X_mmap[rows_neg], dtype=np.float32), np.vstack(mech_neg_parts)])
    X_all = np.vstack([X_pos, X_neg])

    results = []
    violators = []
    for j, name in enumerate(ALL_BASE_FEATURE_NAMES):
        col = X_all[:, j]
        fin = np.isfinite(col)
        if np.sum(fin) < 200 or np.ptp(col[fin]) == 0.0:
            auc = 0.5
        else:
            auc = float(roc_auc_score(y[fin], col[fin]))
        sym_auc = float(max(auc, 1.0 - auc))
        rec = {
            "feature_index": j,
            "feature": name,
            "group": (
                "deconfounded_scarp_geomorph"
                if j < 18
                else "geothermal_vent_mt_rad"
                if j < 34
                else "updip60_subsurface_tilt_strain"
                if j < 50
                else "fold_visible_fracture_mechanics"
            ),
            "roc_auc": round(auc, 5),
            "symmetric_auc": round(sym_auc, 5),
            "finite_fraction": round(float(np.mean(fin)), 5),
            "status": "PASS" if sym_auc <= max_allowed_auc else "LEAKAGE_FAIL",
        }
        results.append(rec)
        if sym_auc > max_allowed_auc:
            violators.append(rec)

    results_sorted = sorted(results, key=lambda r: r["symmetric_auc"], reverse=True)
    summary = {
        "evidence_class": "LEAKAGE-CANARY-AUDIT",
        "evaluator_version": "gems57-pooled-hide-v2",
        "n_withheld_positives_tested": int(len(pos_flat)),
        "n_off_fault_negatives_tested": int(len(neg_flat)),
        "threshold_symmetric_auc": max_allowed_auc,
        "max_symmetric_auc": results_sorted[0]["symmetric_auc"],
        "top_feature": results_sorted[0]["feature"],
        "n_features_tested": len(results),
        "n_violations": len(violators),
        "all_passed": len(violators) == 0,
        "features_ranked": results_sorted,
    }
    if violators:
        raise RuntimeError(f"LEAKAGE CANARY VIOLATION (> {max_allowed_auc}): {violators}")
    return summary
