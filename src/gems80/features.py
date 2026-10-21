"""Zero-leakage multi-physics feature engineering for GEMSDOE80 (Lane H80: ADPT-80).

Core scientific innovations in this module:
1. **360° Directed Normal-Fault Half-Space Tomography & Step Monotonicity (`Group B`)**:
   Evaluates antisymmetric half-space steps across 8 directed azimuths covering [0, 2*pi)
   at short scarp (200 m, 400 m) and range-front piedmont (800 m, 1600 m) baselines,
   paired with the **Step Monotonicity Ratio** `M(x, phi; h) in [0, 1]` (equal to 1 on a
   monotone normal-fault scarp and ~0 on incised stream valleys and symmetric ridges)
   and **Cross-Physics Hanging-Wall / Footwall Polarity Concordance** (co-polarized
   elevation drop + hanging-wall basement cover thickening + isostatic gravity low +
   magnetic tilt zero-crossing + smectite/brine conductivity increase).
2. **Fluvial Interfluve vs Knickpoint Separation & Pleistocene Lake Lahontan Guard**:
   Uses upstream drainage area `log10(A)` and normalized channel steepness `k_sn` to
   separate un-incised alluvial-fan interfluve scarps and fault-crossing wash knickpoints
   from ordinary fluvial banks, plus orthometric elevation distance from Pleistocene
   Lake Lahontan (~1332 m) and Lake Bonneville (~1552 m) highstand benches.
3. **Label-Free 10 m -> 1 m LiDAR Cross-Sensor Super-Resolution Bridge**:
   Trains a purely physical regressor inside the 75.3% interior-eroded 1 m LiDAR domain
   (never touching `existing_faults.tif` or any labels) to impute 1 m scarp dipole
   equivalent across the 24.7% non-LiDAR gap without tile-edge seam artifacts.
4. **Strictly Fold-Visible Elastic Fracture-Mechanics Continuation Features**:
   Computed fresh for any visible fault network `vis` (never seeing hidden truth).
"""
from __future__ import annotations

import csv
import gc
import io
import struct
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from pyproj import Transformer
from scipy import ndimage as ndi
from sklearn.ensemble import HistGradientBoostingRegressor

STATIC_FEATURE_NAMES: tuple[str, ...] = (
    # Group A: Seamless 1m LiDAR + 10m DEM + Fluvial Chi / k_sn Geomorphology (14)
    "lid_ex_max_clean",
    "lid_step_max_clean",
    "lid_scarp_dipole_clean",
    "lid_down_up_max_clean",
    "lid_cross_contour_clean",
    "lid_coh100_clean",
    "superres_scarp_seamless",
    "topo10_slope_p90",
    "topo10_slope_std",
    "topo10_curv_prof",
    "topo10_aspect_coh",
    "topo10_hs_lineament",
    "chi_step_max",
    "chi_ksn_knickpoint",
    # Group B: 360-deg Signed Dip-Polarity Tomography, Step Monotonicity & Lake Guard (12)
    "bipolar_step_200m",
    "bipolar_mono_200m",
    "bipolar_step_400m",
    "bipolar_mono_400m",
    "bipolar_step_800m",
    "bipolar_mono_800m",
    "hw_basement_thickening",
    "hw_gravity_drop",
    "hw_mag_tilt_drop",
    "hw_cond_increase",
    "bipolar_normal_fault_index",
    "along_strike_step_1km",
    # Group C: Fluvial / Pluvial Geomorphic Guards & Official Elevation/Potential Bands (12)
    "fluvial_log10_area",
    "interfluve_scarp_index",
    "lahontan_highstand_guard",
    "det_elev_hp3",
    "det_elev_hp10",
    "det_elev_slope",
    "rtp_hp8",
    "tmi_hg_abs",
    "tc_zero_prox",
    "iso_grav_hp8",
    "iso_grav_hg_abs",
    "crossgrad_gm_parallel",
    # Group D: Deep Crustal Mechanics, Alteration & GDR #1391 Hydrothermal Conduits (8)
    "depth_to_base",
    "depth_base_grad",
    "cond_surf_val",
    "cond_surf_grad",
    "crustal_strain_comb",
    "seismic_density_grad",
    "hydro_alteration_maglow",
    "gdr_thermal_conduit_plume",
    # Group E: Multi-Scale 2D Piedmont Range-Front Structural Context (6)
    "ctx_slope_600m",
    "ctx_slope_2km",
    "ctx_piedmont_transition_2km",
    "ctx_depth_2km",
    "ctx_relief_std_2km",
    "ctx_scarp_chain_1km",
)

FRACTURE_FEATURE_NAMES: tuple[str, ...] = (
    "vis_dist_capped",
    "vis_axial_tip_lobe",
    "vis_wing_crack_lobe",
    "vis_relay_stepover",
    "vis_host_len_prox",
    "vis_accomm_density_2km",
    "vis_lateral_flank_ring",
)

ALL_FEATURE_NAMES: tuple[str, ...] = STATIC_FEATURE_NAMES + FRACTURE_FEATURE_NAMES
N_STATIC = len(STATIC_FEATURE_NAMES)     # 52
N_FRACTURE = len(FRACTURE_FEATURE_NAMES) # 7

# 8 unit normal vectors covering full 360 deg: (dy, dx)
DIRS_8: tuple[tuple[float, float], ...] = tuple(
    (float(np.sin(k * np.pi / 4.0)), float(np.cos(k * np.pi / 4.0))) for k in range(8)
)


def _clean_band(arr: np.ndarray, valid: np.ndarray) -> np.ndarray:
    out = arr.astype(np.float32, copy=True)
    good = valid & np.isfinite(out) & (out > -1e25) & (out < 1e25)
    med = float(np.median(out[good])) if good.any() else 0.0
    out[~good] = med
    return out


def _robust_z_pos(arr: np.ndarray, valid: np.ndarray, clip: float = 8.0) -> np.ndarray:
    vals = arr[valid]
    med = float(np.median(vals))
    q25, q75 = np.percentile(vals, [25.0, 75.0])
    iqr_sd = max(float(q75 - q25) / 1.349, 1e-6)
    return np.clip((arr - med) / iqr_sd, -clip, clip).astype(np.float32)


def _shift2d(arr: np.ndarray, dy: int, dx: int) -> np.ndarray:
    """Fast edge-clamped integer 2D shift without wrap-around seam artifacts."""
    h, w = arr.shape
    out = np.empty_like(arr)
    sy0 = max(0, -dy)
    sy1 = min(h, h - dy)
    dy0 = max(0, dy)
    dy1 = min(h, h + dy)

    sx0 = max(0, -dx)
    sx1 = min(w, w - dx)
    dx0 = max(0, dx)
    dx1 = min(w, w + dx)

    out[dy0:dy1, dx0:dx1] = arr[sy0:sy1, sx0:sx1]
    if dy0 > 0:
        out[:dy0, dx0:dx1] = arr[0:1, sx0:sx1]
    if dy1 < h:
        out[dy1:, dx0:dx1] = arr[h - 1 : h, sx0:sx1]
    if dx0 > 0:
        out[:, :dx0] = out[:, dx0 : dx0 + 1]
    if dx1 < w:
        out[:, dx1:] = out[:, dx1 - 1 : dx1]
    return out


def _directed_step_and_monotonicity(
    field_smooth: np.ndarray, h_px: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute max signed drop S(x) >= 0, Step Monotonicity Ratio M(x) in [0, 1], and argmax dir index.

    For normal vector n_a = (dy_a, dx_a):
      up_a  = field(x - h * n_a)   (uplifted footwall side)
      dn_a  = field(x + h * n_a)   (down-thrown hanging-wall side)
      step  = up_a - dn_a          (> 0 when hanging wall is lower than footwall)
      mono  = |up_a - dn_a| / (|up_a - c| + |c - dn_a| + eps)
              -> 1.0 on a monotone normal-fault scarp; ~0.0 on stream valleys and symmetric ridges!
    """
    best_step = np.full(field_smooth.shape, -1e9, dtype=np.float32)
    best_mono = np.zeros(field_smooth.shape, dtype=np.float32)
    best_dir = np.zeros(field_smooth.shape, dtype=np.uint8)
    c = field_smooth
    eps = np.float32(1e-4)

    for a_idx, (uy, ux) in enumerate(DIRS_8):
        dy = int(round(h_px * uy))
        dx = int(round(h_px * ux))
        up = _shift2d(c, -dy, -dx)
        dn = _shift2d(c, dy, dx)
        step = up - dn
        denom = np.abs(up - c) + np.abs(c - dn) + eps
        mono = np.abs(step) / denom
        better = step > best_step
        best_step = np.where(better, step, best_step)
        best_mono = np.where(better, mono, best_mono)
        best_dir = np.where(better, np.uint8(a_idx), best_dir)
    return np.maximum(best_step, 0.0), best_mono, best_dir


def _directed_contrast_at_phi(field_smooth: np.ndarray, phi_idx: np.ndarray, h_px: int, sign: float) -> np.ndarray:
    """Evaluate sign * (field(x + h * n_{phi*}) - field(x - h * n_{phi*})) at pixel's down-thrown azimuth phi*."""
    out = np.zeros_like(field_smooth, dtype=np.float32)
    for a_idx, (uy, ux) in enumerate(DIRS_8):
        mask = phi_idx == a_idx
        if not mask.any():
            continue
        dy = int(round(h_px * uy))
        dx = int(round(h_px * ux))
        hw = _shift2d(field_smooth, dy, dx)
        fw = _shift2d(field_smooth, -dy, -dx)
        out[mask] = sign * (hw[mask] - fw[mask])
    return out


def _along_strike_continuity(step_field: np.ndarray, phi_idx: np.ndarray, half_len: int = 4) -> np.ndarray:
    """Integrate scarp step strength along orthogonal strike tangent t = (-dx, dy) over 2*half_len + 1 px."""
    # Tangent index is (phi_idx + 2) % 4 among the 4 undirected axes
    axis_idx = phi_idx % 4
    out = np.zeros_like(step_field, dtype=np.float32)
    for ax in range(4):
        mask = axis_idx == ax
        if not mask.any():
            continue
        # Normal is DIRS_8[ax] = (uy, ux) -> strike tangent is (-ux, uy)
        uy, ux = DIRS_8[ax]
        ty, tx = int(round(-ux)), int(round(uy))
        acc = step_field.copy()
        for k in range(1, half_len + 1):
            acc += _shift2d(step_field, k * ty, k * tx)
            acc += _shift2d(step_field, -k * ty, -k * tx)
        out[mask] = acc[mask] / float(2 * half_len + 1)
    return out


def _load_gdr_hydrothermal_grid(
    ext_dir: Path, shape: tuple[int, int], transform: rasterio.Affine
) -> np.ndarray:
    """Rasterize verified DOE GDR #1391 thermal springs/geothermometers, 2m probes, and paleo-sinter."""
    grid = np.zeros(shape, dtype=np.float32)

    # 1. Well & spring geothermometry + measured thermal springs
    ws_path = ext_dir / "gdr_wellspring_in_footprint.csv"
    if ws_path.is_file():
        df = pd.read_csv(ws_path)
        for col in ("temp_c", "geothermquartz_c", "geothermchalc_c", "geothermcat_c", "row", "col"):
            df[col] = pd.to_numeric(df[col], errors="coerce")
        df = df.dropna(subset=["row", "col"])
        r = df["row"].to_numpy(dtype=np.int32)
        c = df["col"].to_numpy(dtype=np.int32)
        inb = (r >= 0) & (r < shape[0]) & (c >= 0) & (c < shape[1])
        t_surf = np.nan_to_num(df["temp_c"].to_numpy(dtype=np.float32), nan=0.0)
        t_geo = np.nanmax(
            np.column_stack([
                np.nan_to_num(df["geothermquartz_c"].to_numpy(dtype=np.float32), nan=0.0),
                np.nan_to_num(df["geothermchalc_c"].to_numpy(dtype=np.float32), nan=0.0),
                np.nan_to_num(df["geothermcat_c"].to_numpy(dtype=np.float32), nan=0.0),
            ]),
            axis=1,
        )
        weight = np.clip((t_surf - 22.0) / 35.0, 0.0, 2.5) + np.clip((t_geo - 85.0) / 80.0, 0.0, 2.5)
        keep = inb & (weight > 0.05)
        np.maximum.at(grid, (r[keep], c[keep]), weight[keep])

    # 2. UNR 2-meter temperature probes (NAD83 lon/lat shapefile inside zip)
    probe_zip = ext_dir / "2m_temperature_probe_INGENIOUS_regional_data.zip"
    if probe_zip.is_file():
        with zipfile.ZipFile(probe_zip) as zf:
            shp_names = [n for n in zf.namelist() if n.endswith(".shp")]
            dbf_names = [n for n in zf.namelist() if n.endswith(".dbf")]
            if shp_names and dbf_names:
                shp_bytes = zf.read(shp_names[0])
                dbf_bytes = zf.read(dbf_names[0])
                xs, ys = _parse_point_shp(shp_bytes)
                records = _parse_simple_dbf(dbf_bytes)
                if len(xs) == len(records) and len(xs) > 0:
                    tr = Transformer.from_crs("EPSG:4269", "EPSG:32611", always_xy=True)
                    ux, uy = tr.transform(np.asarray(xs, dtype=np.float64), np.asarray(ys, dtype=np.float64))
                    cols = np.floor((ux - transform.c) / transform.a).astype(np.int32)
                    rows = np.floor((uy - transform.f) / transform.e).astype(np.int32)
                    for i, rec in enumerate(records):
                        rr, cc = int(rows[i]), int(cols[i])
                        if 0 <= rr < shape[0] and 0 <= cc < shape[1]:
                            anom = 0.0
                            for k, val in rec.items():
                                if "2m" in k.lower() or "anom" in k.lower() or "temp" in k.lower():
                                    try:
                                        v = float(val)
                                        if np.isfinite(v) and v > 1.2 and v < 80.0:
                                            anom = max(anom, min(2.5, v / 3.0))
                                    except (ValueError, TypeError):
                                        pass
                            if anom > 0.1:
                                grid[rr, cc] = max(float(grid[rr, cc]), anom + 0.5)

    # 3. Paleo-geothermal sinter / travertine / tufa deposits (CSV inside zip)
    sinter_zip = ext_dir / "paleo_geothermal_regional.zip"
    if sinter_zip.is_file():
        with zipfile.ZipFile(sinter_zip) as zf:
            csv_names = [n for n in zf.namelist() if n.endswith(".csv")]
            if csv_names:
                text = zf.read(csv_names[0]).decode("utf-8", errors="ignore")
                reader = csv.DictReader(io.StringIO(text))
                lons, lats = [], []
                for row in reader:
                    lat = row.get("lat") or row.get("LAT") or row.get("Latitude") or row.get("Y")
                    lon = row.get("lon") or row.get("LON") or row.get("Longitude") or row.get("X")
                    try:
                        la, lo = float(lat), float(lon)
                        if -125.0 < lo < -110.0 and 34.0 < la < 44.0:
                            lons.append(lo)
                            lats.append(la)
                    except (TypeError, ValueError):
                        pass
                if lons:
                    tr = Transformer.from_crs("EPSG:4326", "EPSG:32611", always_xy=True)
                    ux, uy = tr.transform(np.asarray(lons), np.asarray(lats))
                    cols = np.floor((ux - transform.c) / transform.a).astype(np.int32)
                    rows = np.floor((uy - transform.f) / transform.e).astype(np.int32)
                    inb = (rows >= 0) & (rows < shape[0]) & (cols >= 0) & (cols < shape[1])
                    np.maximum.at(grid, (rows[inb], cols[inb]), 2.5)

    return grid


def _parse_point_shp(buf: bytes) -> tuple[list[float], list[float]]:
    xs, ys = [], []
    off = 100
    n = len(buf)
    while off + 8 <= n:
        _, content_words = struct.unpack(">2i", buf[off : off + 8])
        off += 8
        rec_len = content_words * 2
        if off + rec_len > n or rec_len < 20:
            break
        shp_type = struct.unpack("<i", buf[off : off + 4])[0]
        if shp_type in (1, 11, 21):
            x, y = struct.unpack("<2d", buf[off + 4 : off + 20])
            xs.append(x)
            ys.append(y)
        off += rec_len
    return xs, ys


def _parse_simple_dbf(buf: bytes) -> list[dict[str, str]]:
    n_rec, header_len, rec_len = struct.unpack("<IHH", buf[4:12])
    fields = []
    pos = 32
    while pos + 32 <= header_len and buf[pos] != 0x0D:
        name = buf[pos : pos + 11].split(b"\x00", 1)[0].decode("ascii", errors="ignore")
        flen = buf[pos + 16]
        fields.append((name, flen))
        pos += 32
    rows = []
    off = header_len
    for _ in range(n_rec):
        if off + rec_len > len(buf):
            break
        p = off + 1
        d = {}
        for name, flen in fields:
            d[name] = buf[p : p + flen].decode("latin1", errors="ignore").strip()
            p += flen
        rows.append(d)
        off += rec_len
    return rows


def build_or_load_static_features(
    valid: np.ndarray, data_dir: Path, ext_dir: Path, cache_path: Path
) -> tuple[np.memmap, np.ndarray]:
    """Compute or load the 52-channel zero-catalogue static feature matrix on `valid`."""
    meta_path = cache_path.with_suffix(".strike.npy")
    n_valid = int(valid.sum())
    if cache_path.is_file() and meta_path.is_file():
        mm = np.memmap(cache_path, dtype=np.float32, mode="r", shape=(n_valid, N_STATIC))
        phi_idx = np.load(meta_path)
        if phi_idx.shape == valid.shape:
            return mm, phi_idx

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    mm = np.memmap(cache_path, dtype=np.float32, mode="w+", shape=(n_valid, N_STATIC))

    # --- 1. Read 10m DEM + Fluvial Chi / k_sn / Orthometric Elevation (100% AOI) ---
    with rasterio.open(ext_dir / "topo_10m_derived_100m_u8.tif") as ds:
        t10_slp90 = ds.read(1).astype(np.float32) / 255.0
        t10_slpstd = ds.read(4).astype(np.float32) / 255.0
        t10_relief = ds.read(5).astype(np.float32) / 255.0
        t10_curv = ds.read(6).astype(np.float32) / 255.0
        t10_coh = ds.read(7).astype(np.float32) / 255.0
        t10_line = ds.read(8).astype(np.float32) / 255.0

    with rasterio.open(ext_dir / "dem10_chi_100m_i16.tif") as ds:
        chi_step = np.where(ds.read(1) == -32768, 0.0, ds.read(1).astype(np.float32) / 1000.0)
        chi_ksn = np.where(ds.read(3) == -32768, 0.0, ds.read(3).astype(np.float32))
        log10_area = np.where(ds.read(4) == -32768, 3.0, ds.read(4).astype(np.float32) / 1000.0)
        zmean_m = _clean_band(np.where(ds.read(5) == -32768, np.nan, ds.read(5).astype(np.float32)), valid)
        transform = ds.transform

    # --- 2. Read 1m LiDAR Scarp stack (interior eroded by 2 px to eliminate tile-edge seams!) ---
    with rasterio.open(ext_dir / "lidar_scarp_features_u8.tif") as ds:
        lid_valid_raw = (ds.read(12) > 0) & valid
        lid_interior = ndi.binary_erosion(lid_valid_raw, iterations=2)
        b_ex = ds.read(1).astype(np.float32) / 255.0
        b_step = ds.read(3).astype(np.float32) / 255.0
        b_dipole = 0.5 * (ds.read(4).astype(np.float32) + ds.read(5).astype(np.float32)) / 255.0
        b_du = np.maximum(ds.read(6), ds.read(7)).astype(np.float32) / 255.0
        b_cross = ds.read(8).astype(np.float32) / 255.0
        b_coh = ds.read(10).astype(np.float32) / 255.0

    for idx_col, raw_arr in enumerate((b_ex, b_step, b_dipole, b_du, b_cross, b_coh)):
        tmp = np.where(lid_interior, raw_arr, np.nan).astype(np.float32)
        mm[:, idx_col] = tmp[valid]
    lid_composite_target = 0.35 * b_ex + 0.30 * b_dipole + 0.20 * b_cross + 0.15 * b_coh
    del b_ex, b_step, b_dipole, b_du, b_cross, b_coh
    gc.collect()

    mm[:, 7] = t10_slp90[valid]
    mm[:, 8] = t10_slpstd[valid]
    mm[:, 9] = t10_curv[valid]
    mm[:, 10] = t10_coh[valid]
    mm[:, 11] = t10_line[valid]
    mm[:, 12] = np.log1p(np.maximum(chi_step, 0.0)).astype(np.float32)[valid]
    mm[:, 13] = np.log1p(np.maximum(chi_ksn, 0.0)).astype(np.float32)[valid]

    # --- 3. Read official 19-band geophysical + detrended elevation channels ---
    with rasterio.open(data_dir / "training_features.tif") as ds:
        rtp = _clean_band(ds.read(2), valid)
        tmi_hg = _clean_band(ds.read(3), valid)
        geod_inv = _clean_band(ds.read(4), valid)
        tc = _clean_band(ds.read(6), valid)
        geod_shr = _clean_band(ds.read(7), valid)
        geod_dil = _clean_band(ds.read(8), valid)
        deq = _clean_band(ds.read(10), valid)
        det_elev = _clean_band(ds.read(12), valid)
        iso_grav = _clean_band(ds.read(13), valid)
        depth_base = _clean_band(ds.read(15), valid)
        ieq = _clean_band(ds.read(16), valid)
        cond_surf = _clean_band(ds.read(17), valid)
        grav_hg = _clean_band(ds.read(18), valid)
        det_slope = _clean_band(ds.read(19), valid)

    # --- 4. Group B: 360-deg Signed Dip-Polarity Normal-Fault Half-Space Tomography ---
    z_sm1 = ndi.gaussian_filter(det_elev, sigma=1.0)
    z_sm2 = ndi.gaussian_filter(det_elev, sigma=2.2)
    z_sm4 = ndi.gaussian_filter(det_elev, sigma=4.5)

    step2, mono2, phi2 = _directed_step_and_monotonicity(z_sm1, h_px=2)
    step4, mono4, phi4 = _directed_step_and_monotonicity(z_sm2, h_px=4)
    step8, mono8, phi8 = _directed_step_and_monotonicity(z_sm4, h_px=8)

    mm[:, 14] = np.log1p(step2).astype(np.float32)[valid]
    mm[:, 15] = mono2[valid]
    mm[:, 16] = np.log1p(step4).astype(np.float32)[valid]
    mm[:, 17] = mono4[valid]
    mm[:, 18] = np.log1p(step8).astype(np.float32)[valid]
    mm[:, 19] = mono8[valid]

    # Signed hanging-wall vs footwall subsurface contrasts at down-thrown azimuth phi4!
    dep_sm = ndi.gaussian_filter(depth_base, sigma=3.0)
    grv_sm = ndi.gaussian_filter(iso_grav, sigma=3.0)
    tc_sm = ndi.gaussian_filter(tc, sigma=2.0)
    cnd_sm = ndi.gaussian_filter(cond_surf, sigma=3.0)

    hw_dep = _directed_contrast_at_phi(dep_sm, phi4, h_px=8, sign=+1.0)  # + if hanging wall deeper!
    hw_grv = _directed_contrast_at_phi(grv_sm, phi4, h_px=8, sign=-1.0)  # + if hanging wall lower density!
    hw_tc = _directed_contrast_at_phi(tc_sm, phi4, h_px=5, sign=-1.0)    # + if magnetic tilt drops into basin!
    hw_cnd = _directed_contrast_at_phi(cnd_sm, phi4, h_px=6, sign=+1.0)  # + if hanging wall more conductive!

    mm[:, 20] = hw_dep[valid]
    mm[:, 21] = hw_grv[valid]
    mm[:, 22] = hw_tc[valid]
    mm[:, 23] = hw_cnd[valid]

    # Combined 360-deg Co-Polarized Normal-Fault Tomography Index
    z_scarp = np.maximum(0.0, _robust_z_pos(np.log1p(step4) * mono4, valid))
    z_dep = np.maximum(0.0, _robust_z_pos(hw_dep, valid))
    z_grv = np.maximum(0.0, _robust_z_pos(hw_grv, valid))
    z_tc = np.maximum(0.0, _robust_z_pos(hw_tc, valid))
    bipolar_nf = z_scarp * (0.45 + 0.25 * z_dep + 0.20 * z_grv + 0.10 * z_tc)
    mm[:, 24] = bipolar_nf.astype(np.float32)[valid]

    strike_cont = _along_strike_continuity(np.log1p(step2) * mono2, phi2, half_len=4)
    mm[:, 25] = strike_cont.astype(np.float32)[valid]
    np.save(meta_path, phi2)

    del z_sm1, z_sm2, z_sm4, dep_sm, grv_sm, tc_sm, cnd_sm, z_scarp, z_dep, z_grv, z_tc
    gc.collect()

    # --- 5. Group C: Fluvial Interfluve Index, Lake Lahontan Guard & Potential Fields ---
    mm[:, 26] = log10_area.astype(np.float32)[valid]
    interfluve_w = np.clip((4.6 - log10_area) / 1.4, 0.0, 1.0)
    interfluve_scarp = (np.log1p(step2) * mono2) * interfluve_w + 0.25 * np.log1p(np.maximum(chi_ksn, 0.0)) * (1.0 - interfluve_w)
    mm[:, 27] = interfluve_scarp.astype(np.float32)[valid]

    # Distance in elevation from Pleistocene Lake Lahontan (~1332 m) & Bonneville (~1552 m) highstands
    d_lahontan = np.minimum(np.abs(zmean_m - 1332.0), np.abs(zmean_m - 1270.0))
    d_bonneville = np.abs(zmean_m - 1552.0)
    pluvial_prox = np.exp(-0.5 * (np.minimum(d_lahontan, d_bonneville) / 14.0) ** 2).astype(np.float32)
    mm[:, 28] = pluvial_prox[valid]

    det_hp3 = det_elev - ndi.gaussian_filter(det_elev, sigma=3.0)
    det_hp10 = det_elev - ndi.gaussian_filter(det_elev, sigma=10.0)
    rtp_hp8 = rtp - ndi.gaussian_filter(rtp, sigma=8.0)
    grav_hp8 = iso_grav - ndi.gaussian_filter(iso_grav, sigma=8.0)

    mm[:, 29] = det_hp3.astype(np.float32)[valid]
    mm[:, 30] = det_hp10.astype(np.float32)[valid]
    mm[:, 31] = det_slope.astype(np.float32)[valid]
    mm[:, 32] = rtp_hp8.astype(np.float32)[valid]
    mm[:, 33] = np.log1p(np.abs(tmi_hg)).astype(np.float32)[valid]
    tc_p90 = float(np.percentile(np.abs(tc[valid]), 90)) + 1e-6
    mm[:, 34] = np.clip(1.0 - np.abs(tc) / tc_p90, 0.0, 1.0).astype(np.float32)[valid]
    mm[:, 35] = grav_hp8.astype(np.float32)[valid]
    mm[:, 36] = np.log1p(np.abs(grav_hg)).astype(np.float32)[valid]

    gy_g, gx_g = np.gradient(ndi.gaussian_filter(iso_grav, sigma=2.5))
    gy_m, gx_m = np.gradient(ndi.gaussian_filter(rtp, sigma=2.5))
    ng = np.hypot(gy_g, gx_g) + 1e-8
    nm = np.hypot(gy_m, gx_m) + 1e-8
    cos_gm = ((gy_g * gy_m + gx_g * gx_m) / (ng * nm)) ** 2
    mm[:, 37] = (cos_gm * np.log1p(ng * 10.0) * np.log1p(nm * 0.1)).astype(np.float32)[valid]
    del gy_g, gx_g, gy_m, gx_m, ng, nm, cos_gm
    gc.collect()

    # --- 6. Group D: Deep Crustal Mechanics, Hydrothermal Alteration & GDR Conduits ---
    gy_d, gx_d = np.gradient(ndi.gaussian_filter(depth_base, sigma=2.0))
    gy_c, gx_c = np.gradient(ndi.gaussian_filter(cond_surf, sigma=2.0))
    mm[:, 38] = depth_base.astype(np.float32)[valid]
    mm[:, 39] = np.hypot(gy_d, gx_d).astype(np.float32)[valid]
    mm[:, 40] = cond_surf.astype(np.float32)[valid]
    mm[:, 41] = np.hypot(gy_c, gx_c).astype(np.float32)[valid]
    del gy_d, gx_d, gy_c, gx_c

    strain_comb = (
        np.abs(_robust_z_pos(geod_inv, valid))
        + np.abs(_robust_z_pos(geod_shr, valid))
        + 0.5 * np.abs(_robust_z_pos(geod_dil, valid))
    )
    gy_eq, gx_eq = np.gradient(ndi.gaussian_filter(ieq, sigma=3.0))
    mm[:, 42] = strain_comb.astype(np.float32)[valid]
    mm[:, 43] = (np.hypot(gy_eq, gx_eq) + 0.1 * np.exp(-np.maximum(deq, 0.0) / 25.0)).astype(np.float32)[valid]
    del geod_inv, geod_shr, geod_dil, deq, ieq, gy_eq, gx_eq

    with rasterio.open(ext_dir / "geodawn_rad_u8.tif") as ds:
        rad_k = ds.read(1).astype(np.float32)
        rad_th = ds.read(2).astype(np.float32)
    with rasterio.open(ext_dir / "geodawn_extensions_u8.tif") as ds:
        tmi_up150 = ds.read(4).astype(np.float32)

    k_th_z = np.maximum(0.0, _robust_z_pos(rad_k / (rad_th + 8.0), valid))
    mag_low_z = np.maximum(0.0, _robust_z_pos(-rtp_hp8, valid))
    up_diff = np.abs(_robust_z_pos(rtp_hp8, valid) - _robust_z_pos(tmi_up150 - ndi.gaussian_filter(tmi_up150, 8.0), valid))
    mm[:, 44] = (mag_low_z * (0.4 + 0.6 * k_th_z) + 0.25 * up_diff).astype(np.float32)[valid]
    del rad_k, rad_th, tmi_up150, k_th_z, mag_low_z, up_diff, rtp, iso_grav, rtp_hp8, grav_hp8
    gc.collect()

    gdr_pts = _load_gdr_hydrothermal_grid(ext_dir, valid.shape, transform)
    gdr_plume = 0.65 * ndi.gaussian_filter(gdr_pts, sigma=3.5) + 0.35 * ndi.gaussian_filter(gdr_pts, sigma=11.0)
    mm[:, 45] = np.log1p(gdr_plume * 15.0).astype(np.float32)[valid]
    del gdr_pts, gdr_plume

    # --- 7. Group E: Multi-Scale 2D Piedmont Range-Front Structural Context ---
    slp_600 = ndi.gaussian_filter(det_slope, sigma=6.0)
    slp_2km = ndi.gaussian_filter(det_slope, sigma=18.0)
    # Signed upslope minus downslope slope difference across 1.6 km: peaks right on the piedmont range-front break!
    pied_trans = _directed_contrast_at_phi(slp_600, phi4, h_px=12, sign=-1.0)
    dep_2km = ndi.gaussian_filter(depth_base, sigma=18.0)
    rel_sq = ndi.gaussian_filter(det_elev * det_elev, sigma=15.0) - ndi.gaussian_filter(det_elev, sigma=15.0) ** 2
    rel_std_2km = np.sqrt(np.maximum(rel_sq, 0.0))
    scarp_chain_1km = ndi.gaussian_filter(strike_cont, sigma=4.0)

    mm[:, 46] = slp_600.astype(np.float32)[valid]
    mm[:, 47] = slp_2km.astype(np.float32)[valid]
    mm[:, 48] = pied_trans.astype(np.float32)[valid]
    mm[:, 49] = dep_2km.astype(np.float32)[valid]
    mm[:, 50] = rel_std_2km.astype(np.float32)[valid]
    mm[:, 51] = scarp_chain_1km.astype(np.float32)[valid]

    del slp_600, slp_2km, pied_trans, dep_2km, rel_sq, rel_std_2km, scarp_chain_1km, det_elev, det_slope, depth_base
    gc.collect()

    # --- 8. Fill `superres_scarp_seamless` (col 6) using ONLY zero-label 10m + physical features! ---
    # Trained strictly on interior-LiDAR pixels (`lid_interior[valid]`) with ZERO access to any fault catalogue.
    predictor_cols = [7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 25, 26, 27, 29, 31, 46, 48]
    lid_mask_v = lid_interior[valid]
    y_super = lid_composite_target[valid]
    rng = np.random.default_rng(8001)
    tr_idx = rng.choice(np.flatnonzero(lid_mask_v), size=min(220_000, int(lid_mask_v.sum())), replace=False)
    sr_reg = HistGradientBoostingRegressor(
        max_iter=120, learning_rate=0.07, max_leaf_nodes=31, min_samples_leaf=150, l2_regularization=3.0, random_state=8001
    )
    sr_reg.fit(np.asarray(mm[tr_idx][:, predictor_cols], dtype=np.float32), y_super[tr_idx])
    super_pred = np.empty(n_valid, dtype=np.float32)
    step = 600_000
    for s in range(0, n_valid, step):
        super_pred[s : s + step] = sr_reg.predict(np.asarray(mm[s : s + step, predictor_cols], dtype=np.float32))
    # Inside interior-valid 1m LiDAR, blend 75% actual 1m LiDAR composite + 25% super-resolved 10m;
    # outside 1m LiDAR (the 24.7% non-LiDAR gap + 2-px tile seams), use quantile-matched super-resolved 10m!
    gap_v = ~lid_mask_v
    q_grid = np.linspace(0.0, 1.0, 1001)
    super_cal = super_pred.copy()
    super_cal[gap_v] = np.interp(
        super_pred[gap_v],
        np.quantile(super_pred[gap_v], q_grid),
        np.quantile(y_super[lid_mask_v], q_grid),
    ).astype(np.float32)
    super_cal[lid_mask_v] = 0.75 * y_super[lid_mask_v] + 0.25 * super_pred[lid_mask_v]
    mm[:, 6] = super_cal
    mm.flush()

    del lid_composite_target, super_pred, super_cal, t10_slp90, t10_slpstd, t10_relief, t10_curv, t10_coh, t10_line
    gc.collect()
    return np.memmap(cache_path, dtype=np.float32, mode="r", shape=(n_valid, N_STATIC)), phi2


def compute_component_tips(catalogue: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Extract principal endpoints and outward unit tangents for each 8-connected trace >= 4 px."""
    cat = np.asarray(catalogue, bool)
    comp, n = ndi.label(cat, structure=np.ones((3, 3), dtype=int))
    if n == 0:
        z = np.zeros(0, dtype=np.float32)
        return z, z, z, z, z
    ys, xs = np.nonzero(cat)
    ids = comp[ys, xs]
    counts = np.bincount(ids, minlength=n + 1).astype(np.float64)
    keep_ids = np.flatnonzero(counts >= 4)
    if keep_ids.size == 0:
        z = np.zeros(0, dtype=np.float32)
        return z, z, z, z, z

    cy = np.bincount(ids, weights=ys, minlength=n + 1) / np.maximum(counts, 1.0)
    cx = np.bincount(ids, weights=xs, minlength=n + 1) / np.maximum(counts, 1.0)
    dy = ys - cy[ids]
    dx = xs - cx[ids]
    syy = np.bincount(ids, weights=dy * dy, minlength=n + 1)
    sxx = np.bincount(ids, weights=dx * dx, minlength=n + 1)
    sxy = np.bincount(ids, weights=dy * dx, minlength=n + 1)
    theta = 0.5 * np.arctan2(2.0 * sxy, sxx - syy + 1e-12)
    vy = np.cos(theta)  # note: sxx - syy uses x first! Be exact with eigenvector below:
    # Exact principal eigenvector of [[syy, sxy], [sxy, sxx]]:
    ang = 0.5 * np.arctan2(2.0 * sxy, syy - sxx + 1e-12)
    ey = np.cos(ang)
    ex = np.sin(ang)

    proj = dy * ey[ids] + dx * ex[ids]
    order = np.argsort(ids, kind="mergesort")
    ids_s = ids[order]
    proj_s = proj[order]
    ys_s = ys[order]
    xs_s = xs[order]
    splits = np.flatnonzero(np.diff(ids_s)) + 1
    groups_idx = np.split(np.arange(ids_s.size), splits)

    ty_list, tx_list, oy_list, ox_list, len_list = [], [], [], [], []
    for g in groups_idx:
        cid = int(ids_s[g[0]])
        cnt = float(counts[cid])
        if cnt < 4.0:
            continue
        p_g = proj_s[g]
        i_min = g[int(np.argmin(p_g))]
        i_max = g[int(np.argmax(p_g))]
        # Local tangent from pixels within 8 px of each tip for curved faults!
        for tip_idx, sign in ((i_min, -1.0), (i_max, +1.0)):
            ty0, tx0 = float(ys_s[tip_idx]), float(xs_s[tip_idx])
            d2 = (ys_s[g] - ty0) ** 2 + (xs_s[g] - tx0) ** 2
            near = g[(d2 > 0.5) & (d2 <= 64.0)]
            if near.size >= 2:
                vy_loc = float(ty0 - np.mean(ys_s[near]))
                vx_loc = float(tx0 - np.mean(xs_s[near]))
                nrm = float(np.hypot(vy_loc, vx_loc))
                if nrm > 1e-4:
                    uy, ux = vy_loc / nrm, vx_loc / nrm
                else:
                    uy, ux = float(sign * ey[cid]), float(sign * ex[cid])
            else:
                uy, ux = float(sign * ey[cid]), float(sign * ex[cid])
            ty_list.append(ty0)
            tx_list.append(tx0)
            oy_list.append(uy)
            ox_list.append(ux)
            len_list.append(cnt)

    return (
        np.asarray(ty_list, dtype=np.float32),
        np.asarray(tx_list, dtype=np.float32),
        np.asarray(oy_list, dtype=np.float32),
        np.asarray(ox_list, dtype=np.float32),
        np.asarray(len_list, dtype=np.float32),
    )


def compute_fracture_mechanics_features(
    visible: np.ndarray, valid: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute strictly fold-proper 7 fracture-mechanics channels on `valid`, plus 2D `d_vis` & `tip_cone`.

    Zero-leakage guarantee: depends ONLY on `visible` (the fold's un-hidden catalogue).
    """
    from scipy.spatial import cKDTree

    vis = np.asarray(visible, dtype=bool)
    d_vis, inds = ndi.distance_transform_edt(~vis, return_indices=True)
    d_vis = d_vis.astype(np.float32)

    comp, _ = ndi.label(vis, structure=np.ones((3, 3), dtype=int))
    sizes = np.bincount(comp.ravel())
    sizes[0] = 0
    nearest_len = sizes[comp[inds[0], inds[1]]].astype(np.float32)
    del comp, sizes, inds

    ty, tx, oy, ox, tlen = compute_component_tips(vis)
    n_valid = int(valid.sum())
    F = np.zeros((n_valid, N_FRACTURE), dtype=np.float32)
    F[:, 0] = np.minimum(d_vis[valid], 50.0)

    tip_cone_2d = np.zeros(vis.shape, dtype=bool)
    if ty.size > 0:
        tree = cKDTree(np.column_stack([ty, tx]))
        vy_all, vx_all = np.nonzero(valid)
        k_q = min(2, ty.size)
        # Chunked KD-tree query to keep peak RAM tiny (<120 MB)
        ax_lobe = np.zeros(n_valid, dtype=np.float32)
        wing_lobe = np.zeros(n_valid, dtype=np.float32)
        relay_bridge = np.zeros(n_valid, dtype=np.float32)
        tip_cone_v = np.zeros(n_valid, dtype=bool)
        w_scale = np.log1p(tlen) / np.log1p(25.0)

        step = 500_000
        sin25, cos25 = float(np.sin(np.deg2rad(26.0))), float(np.cos(np.deg2rad(26.0)))
        for s in range(0, n_valid, step):
            pts = np.column_stack([vy_all[s : s + step], vx_all[s : s + step]])
            dists, idxs = tree.query(pts, k=k_q, distance_upper_bound=42.0)
            if k_q == 1:
                dists = dists[:, None]
                idxs = idxs[:, None]
            for j_col in range(k_q):
                idx_j = idxs[:, j_col]
                ok = idx_j < ty.size
                if not ok.any():
                    continue
                idx_ok = idx_j[ok]
                dy = pts[ok, 0] - ty[idx_ok]
                dx = pts[ok, 1] - tx[idx_ok]
                u = dy * oy[idx_ok] + dx * ox[idx_ok]       # along outward tip tangent
                v = -dy * ox[idx_ok] + dx * oy[idx_ok]      # cross-strike signed offset
                sc = w_scale[idx_ok]

                # 1. Axial continuation lobe ahead of fault tip (u > 0, tight |v|)
                lobe_a = np.where(
                    u > 0.4,
                    sc * np.exp(-u / 18.0) * np.exp(-0.5 * (v / 3.2) ** 2),
                    0.0,
                )
                ax_lobe[s : s + step][ok] = np.maximum(ax_lobe[s : s + step][ok], lobe_a)

                # Tip-cone exemption where unmapped trace continues right off a visible tip (`1.05 < d_vis <= 2.1`)
                cone_hit = (u > 0.8) & (u < 6.0) & (np.abs(v) <= 1.8)
                tip_cone_v[s : s + step][ok] |= cone_hit

                # 2. +/- 26-deg Andersonian tensile wing-crack splays at fault tips
                u_w1 = u * cos25 + v * sin25
                v_w1 = -u * sin25 + v * cos25
                u_w2 = u * cos25 - v * sin25
                v_w2 = u * sin25 + v * cos25
                w1 = np.where(u_w1 > 0.5, sc * np.exp(-u_w1 / 14.0) * np.exp(-0.5 * (v_w1 / 3.5) ** 2), 0.0)
                w2 = np.where(u_w2 > 0.5, sc * np.exp(-u_w2 / 14.0) * np.exp(-0.5 * (v_w2 / 3.5) ** 2), 0.0)
                wing_lobe[s : s + step][ok] = np.maximum(wing_lobe[s : s + step][ok], np.maximum(w1, w2))

            # 3. En-echelon relay ramp / stepover bridge between TWO distinct fault tips!
            if k_q >= 2:
                both = (idxs[:, 0] < ty.size) & (idxs[:, 1] < ty.size) & (idxs[:, 0] != idxs[:, 1])
                if both.any():
                    i0, i1 = idxs[both, 0], idxs[both, 1]
                    d0, d1 = dists[both, 0], dists[both, 1]
                    # Facing / sub-parallel stepover geometry: dot(o0, o1) <= 0.25 or overlapping relay
                    anti_align = np.clip(0.65 - 0.35 * (oy[i0] * oy[i1] + ox[i0] * ox[i1]), 0.15, 1.0)
                    bridge = anti_align * np.exp(-(d0 + d1) / 22.0) * np.minimum(w_scale[i0], w_scale[i1])
                    relay_bridge[s : s + step][both] = bridge

        F[:, 1] = ax_lobe
        F[:, 2] = wing_lobe
        F[:, 3] = relay_bridge
        tip_cone_2d[valid] = tip_cone_v

    F[:, 4] = (np.log1p(nearest_len[valid]) * np.exp(-d_vis[valid] / 14.0)).astype(np.float32)
    accomm = ndi.gaussian_filter(vis.astype(np.float32), sigma=16.0)
    F[:, 5] = accomm[valid].astype(np.float32)
    F[:, 6] = ((d_vis[valid] > 0.0) & (d_vis[valid] <= 2.2) & (~tip_cone_2d[valid])).astype(np.float32)
    return F, d_vis, tip_cone_2d
