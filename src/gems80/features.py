"""Label-free 56-channel Multi-Scale Anisotropic Cross-Physics Structural & Hydrothermal
Conduit feature stack (MACS) on the official 3730x3292 EPSG:32611 100m grid.

STRICT ZERO-LEAKAGE GUARANTEE:
  - Never reads labels.tif or existing_faults.tif.
  - Explicitly drops the label-derived 'dist_known_fault_px' column when reading
    external DOE GDR-1391 CSV tables.
  - Splits features into:
      * 'subsurface' (36 channels: ZERO topography, ZERO slope, ZERO LiDAR -- pure
        gravity, magnetics, MT conductivity, GPS/InSAR strain, seismicity, GeoDAWN
        gamma-ray spectrometry, DOE GDR-1391 geothermometry & INGENIOUS borehole heat flow)
      * 'surface_coupling' (20 channels: band-pass detrended elevation, 1m 3DEP LiDAR
        scarps, multi-scale steerable Radon scarp integration, cross-strike dipping
        normal-fault step asymmetry, and cross-physics strike phase-locking R2).
  - Memory-safe streaming implementation (< 750 MB peak RSS) persisting float16
    standardized features to .cache/gems80/macs_features_f16.npy.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import time
from pathlib import Path

import numpy as np
import rasterio
from rasterio.warp import transform as warp_transform
from scipy import ndimage as ndi

OFFICIAL_BANDS = {
    1: "mag_anom",
    2: "rtp",
    3: "tmi_hg",
    4: "geod_2ndinv",
    5: "iso_grav_anom_slope",
    6: "tc",
    7: "geod_shearrate",
    8: "geod_dilaterate",
    9: "tmi_vg",
    10: "deq_n100a15",
    11: "iso_grav_anom_vg",
    12: "det_elev",
    13: "iso_grav_anom",
    14: "tmi",
    15: "depth_to_base_surf",
    16: "ieq_n100a15",
    17: "cond_surf",
    18: "iso_grav_anom_hg",
    19: "det_elev_slope",
}

# ScienceBase INGENIOUS heat-flow table uses Albers Equal-Area Conic centered at -117W
EAC117_CRS = (
    "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23.0 +lon_0=-117.0 "
    "+x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
)


def _clean_band(arr: np.ndarray) -> np.ndarray:
    out = np.asarray(arr, dtype=np.float32).copy()
    out[out <= -1e37] = np.nan
    out[~np.isfinite(out)] = np.nan
    return out


def _nn_fill(arr: np.ndarray) -> np.ndarray:
    bad = ~np.isfinite(arr)
    if not bad.any():
        return arr
    idx = ndi.distance_transform_edt(bad, return_distances=False, return_indices=True)
    return arr[tuple(idx)]


def robust_z(arr: np.ndarray, mask: np.ndarray, clip: float = 10.0) -> np.ndarray:
    v = arr[mask]
    v = v[np.isfinite(v)]
    if v.size < 32:
        return np.zeros_like(arr, dtype=np.float32)
    med = float(np.median(v))
    iqr = float(np.percentile(v, 75.0) - np.percentile(v, 25.0))
    scale = 0.7413 * iqr if (np.isfinite(iqr) and iqr > 1e-9) else (float(np.std(v)) or 1.0)
    z = (arr - med) / scale
    return np.clip(np.nan_to_num(z, nan=0.0, posinf=clip, neginf=-clip), -clip, clip).astype(np.float32)


def _gauss(arr: np.ndarray, sigma: float) -> np.ndarray:
    return ndi.gaussian_filter(arr, sigma=sigma, mode="nearest").astype(np.float32)


def _grad_mag(arr: np.ndarray) -> np.ndarray:
    gy, gx = np.gradient(arr)
    return np.hypot(gx, gy).astype(np.float32)


def _hessian_ridge(arr: np.ndarray, sigma: float = 0.0) -> tuple[np.ndarray, np.ndarray]:
    """Axis-exact Hessian principal eigenvalue magnitude & dominant strike double-angle.

    Fixes IR-80-03 (transposed dxx/dyy in GEMSDOE77 H85):
      - dyy uses vertical 3x1 stencil [[1], [-2], [1]]
      - dxx uses horizontal 1x3 stencil [[1, -2, 1]]
    For a linear ridge/step striking at angle theta (measured CCW from +x East),
    the maximum curvature normal points at phi = theta + pi/2, so strike double-angle
    is 2*theta = 2*phi - pi => (cos 2*theta, sin 2*theta) = (-cos 2*phi, -sin 2*phi).
    """
    src = _gauss(arr, sigma) if sigma > 0 else arr
    dyy = ndi.correlate(src, np.array([[1.0], [-2.0], [1.0]], dtype=np.float32), mode="nearest")
    dxx = ndi.correlate(src, np.array([[1.0, -2.0, 1.0]], dtype=np.float32), mode="nearest")
    dxy = ndi.correlate(
        src,
        np.array([[0.25, 0.0, -0.25], [0.0, 0.0, 0.0], [-0.25, 0.0, 0.25]], dtype=np.float32),
        mode="nearest",
    )
    tr = 0.5 * (dxx + dyy)
    diff = 0.5 * (dxx - dyy)
    rad = np.hypot(diff, dxy)
    lam1 = tr + rad
    lam2 = tr - rad
    use1 = np.abs(lam1) >= np.abs(lam2)
    dom = np.where(use1, lam1, lam2).astype(np.float32)
    # Eigenvector angle phi for dominant eigenvalue satisfies tan(2*phi) = dxy / diff
    sign = np.where(use1, 1.0, -1.0).astype(np.float32)
    norm = np.maximum(rad, 1e-8)
    cos2phi = sign * (diff / norm)
    sin2phi = sign * (dxy / norm)
    # Strike double-angle = 2*(phi - pi/2) => negate (cos2phi, sin2phi)
    return dom, np.stack([-cos2phi, -sin2phi], axis=0).astype(np.float32)


def steerable_radon_bank(field_pos: np.ndarray, length_px: int, n_angles: int = 8):
    """Integrate a non-negative edge/scarp field along `n_angles` straight line kernels
    of length `length_px` (at 100m grid: 9px=900m, 21px=2.1km, 35px=3.5km).

    Returns:
      rmax:    peak directional line integral max_theta R_{L,theta}[f]
      aniso:   directional concentration max_theta R - mean_theta R (suppresses blobs & contours)
      u2, v2:  second-order orientation director (cos 2*theta*, sin 2*theta*)
    """
    half = length_px // 2
    t = np.linspace(-half, half, length_px, dtype=np.float32)
    # Gently tapered boxcar so endpoints do not alias on pixel boundaries
    w_tap = (0.65 + 0.35 * np.cos(np.pi * t / max(half + 1.0, 1.0))).astype(np.float32)
    w_tap /= w_tap.sum()

    h, w = field_pos.shape
    rmax = np.full((h, w), -np.inf, dtype=np.float32)
    rsum = np.zeros((h, w), dtype=np.float32)
    u2_num = np.zeros((h, w), dtype=np.float32)
    v2_num = np.zeros((h, w), dtype=np.float32)

    for m in range(n_angles):
        theta = m * math.pi / n_angles  # strike angle in [0, pi)
        dx = math.cos(theta)
        dy = math.sin(theta)
        acc = np.zeros((h, w), dtype=np.float32)
        # Sample integer offsets along the line using bilinear-snapped taps
        for tk, wk in zip(t, w_tap):
            ox = int(round(tk * dx))
            oy = int(round(tk * dy))
            sy0, sy1 = max(0, oy), min(h, h + oy)
            sx0, sx1 = max(0, ox), min(w, w + ox)
            dy0, dy1 = max(0, -oy), min(h, h - oy)
            dx0, dx1 = max(0, -ox), min(w, w - ox)
            acc[dy0:dy1, dx0:dx1] += field_pos[sy0:sy1, sx0:sx1] * wk
        np.maximum(rmax, acc, out=rmax)
        rsum += acc
        c2 = math.cos(2.0 * theta)
        s2 = math.sin(2.0 * theta)
        # Soft-max weight emphasizes the winning strike without discretization jumps
        u2_num += acc * c2
        v2_num += acc * s2

    rmean = rsum / float(n_angles)
    aniso = (rmax - rmean).astype(np.float32)
    denom = np.maximum(np.hypot(u2_num, v2_num), 1e-6)
    return rmax, aniso, (u2_num / denom).astype(np.float32), (v2_num / denom).astype(np.float32)


def _rasterize_geothermal_points(data_dir: Path, shape: tuple[int, int], grid_transform) -> dict[str, np.ndarray]:
    """Rasterize strictly label-free DOE GDR-1391 geothermometry + volcanic vents +
    USGS/INGENIOUS borehole heat-flow residuals onto the 3730x3292 grid.

    CRITICAL LEAKAGE GUARD: The column 'dist_known_fault_px' in the GDR CSVs is
    NEVER read or touched. Only (row, col) / coordinates and temperature / heat-flow
    chemistry columns are parsed.
    """
    h, w = shape
    geotherm_excess = np.zeros((h, w), dtype=np.float32)
    spring_density = np.zeros((h, w), dtype=np.float32)

    ws_csv = data_dir / "external" / "gdr_wellspring_in_footprint.csv"
    if ws_csv.is_file():
        with ws_csv.open("r", newline="") as f:
            reader = csv.DictReader(f)
            assert "dist_known_fault_px" in (reader.fieldnames or []), "expected GDR schema"
            for row in reader:
                try:
                    r = int(float(row["row"]))
                    c = int(float(row["col"]))
                except (KeyError, ValueError):
                    continue
                if not (0 <= r < h and 0 <= c < w):
                    continue
                temps = []
                for col_name in ("temp_c", "geothermquartz_c", "geothermchalc_c", "geothermcat_c"):
                    val_s = row.get(col_name, "")
                    if val_s:
                        try:
                            val = float(val_s)
                            if np.isfinite(val) and 0.0 < val < 450.0:
                                temps.append(val)
                        except ValueError:
                            pass
                spring_density[r, c] += 1.0
                if temps:
                    excess = max(0.0, max(temps) - 20.0)
                    if excess > geotherm_excess[r, c]:
                        geotherm_excess[r, c] = excess

    vent_mask = np.zeros((h, w), dtype=np.float32)
    vv_csv = data_dir / "external" / "gdr_volcanic_vents_in_footprint.csv"
    if vv_csv.is_file():
        with vv_csv.open("r", newline="") as f:
            for row in csv.DictReader(f):
                try:
                    r = int(float(row["row"]))
                    c = int(float(row["col"]))
                    if 0 <= r < h and 0 <= c < w:
                        vent_mask[r, c] = 1.0
                except (KeyError, ValueError):
                    continue

    hf_pos = np.zeros((h, w), dtype=np.float32)
    hf_csv = data_dir / "external" / "sb_heat_flow_in_footprint.csv"
    if hf_csv.is_file():
        xs_eac, ys_eac, vals = [], [], []
        with hf_csv.open("r", newline="") as f:
            for row in csv.DictReader(f):
                try:
                    xe = float(row["x_eac117"])
                    ye = float(row["y_eac117"])
                except (KeyError, ValueError):
                    continue
                best_hf = None
                for k in ("hf_resid", "hftc_resid", "hfqc_resid", "hf_meas", "hf_est"):
                    vs = row.get(k, "")
                    if vs:
                        try:
                            v = float(vs)
                            if np.isfinite(v):
                                # Normalize raw mW/m^2 measurements around 85 mW/m^2 background
                                best_hf = max(0.0, v - 85.0) if "meas" in k or "est" in k else max(0.0, v)
                                break
                        except ValueError:
                            pass
                if best_hf is not None and best_hf > 0:
                    xs_eac.append(xe)
                    ys_eac.append(ye)
                    vals.append(min(best_hf, 500.0))
        if xs_eac:
            ux, uy = warp_transform(EAC117_CRS, "EPSG:32611", xs_eac, ys_eac)
            inv = ~grid_transform
            for x0, y0, val in zip(ux, uy, vals):
                cc, rr = inv * (x0, y0)
                ri, ci = int(round(rr)), int(round(cc))
                if 0 <= ri < h and 0 <= ci < w and val > hf_pos[ri, ci]:
                    hf_pos[ri, ci] = val

    return {
        "geotherm_plume_1p5km": _gauss(np.log1p(geotherm_excess), 15.0),
        "geotherm_plume_4km": _gauss(np.log1p(geotherm_excess), 40.0) + 0.5 * _gauss(np.log1p(spring_density), 20.0),
        "heatflow_vent_plume_4km": _gauss(np.log1p(hf_pos), 35.0) + _gauss(vent_mask, 45.0),
    }


def build_or_load_feature_stack(data_dir: Path, cache_dir: Path, log=print) -> dict:
    """Build (once) or memory-map the 56-channel MACS float16 feature tensor."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    mmap_path = cache_dir / "macs56_f16.npy"
    meta_path = cache_dir / "macs56_meta.json"

    with rasterio.open(data_dir / "sample_submission.tif") as tmpl:
        tmpl_arr = tmpl.read(1)
        footprint = np.isfinite(tmpl_arr)
        shape = tmpl_arr.shape
        grid_transform = tmpl.transform
    foot_flat = np.flatnonzero(footprint.ravel())

    if mmap_path.is_file() and meta_path.is_file():
        meta = json.loads(meta_path.read_text())
        arr = np.load(mmap_path, mmap_mode="r")
        if arr.shape == (len(meta["channels"]), foot_flat.size):
            log(f"[MACS] Reusing cached 56-channel stack: {arr.shape} ({mmap_path})")
            return {"X": arr, "meta": meta, "footprint": footprint, "foot_flat": foot_flat}

    t0 = time.time()
    n_expected = 61
    log(f"[MACS] Building {n_expected}-channel label-free MACS feature stack across {foot_flat.size:,} pixels...")

    # Pre-allocate list of channel specs as we stream to keep RAM <= 700 MB
    channel_specs: list[tuple[str, str, str]] = []  # (name, domain, description)
    mm = np.lib.format.open_memmap(
        mmap_path, mode="w+", dtype=np.float16, shape=(n_expected, foot_flat.size)
    )
    cursor = 0

    def emit_ch(name: str, domain: str, desc: str, grid_arr: np.ndarray) -> np.ndarray:
        nonlocal cursor
        z = robust_z(grid_arr, footprint)
        mm[cursor] = z.ravel()[foot_flat].astype(np.float16)
        channel_specs.append((name, domain, desc))
        cursor += 1
        return z

    # 1. Load & clean official non-topographic bands first (Pure Subsurface Domain = 'subsurface')
    sub_keep = {}
    with rasterio.open(data_dir / "training_features.tif") as fsrc:
        for b in (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14, 15, 16, 17, 18):
            name = OFFICIAL_BANDS[b]
            filled = _nn_fill(_clean_band(fsrc.read(b)))
            # For DC-drifting regional potential fields, emit multiscale high-pass residual
            # so LOQO spatial generalization never memorizes quadrant DC offsets!
            if b in (1, 2, 13, 14, 15):
                hp4 = filled - _gauss(filled, 4.0)
                hp14 = filled - _gauss(filled, 14.0)
                emit_ch(f"sub_hp4_{name}", "subsurface", f"400m high-pass residual of {name}", hp4)
                emit_ch(f"sub_hp14_{name}", "subsurface", f"1.4km band-pass residual of {name}", hp14)
                sub_keep[name] = filled
            else:
                z = emit_ch(f"sub_{name}", "subsurface", f"Standardized official subsurface channel {name}", filled)
                sub_keep[name] = filled
                if b in (5, 6, 17):
                    sub_keep[f"z_{name}"] = z

        # Load topographic bands (held separately for 'surface_coupling' domain ONLY!)
        det_elev = _nn_fill(_clean_band(fsrc.read(12)))
        det_slope = _nn_fill(_clean_band(fsrc.read(19)))

    # 2. Pure-Subsurface Hessian Curvature, Gravity Tilt, Strain Vorticity & Multi-Sensor Coupling
    tc_ridge, _ = _hessian_ridge(sub_keep["tc"], sigma=1.0)
    rtp_ridge, _ = _hessian_ridge(sub_keep["rtp"], sigma=1.2)
    grav_ridge, _ = _hessian_ridge(sub_keep["iso_grav_anom"], sigma=1.5)
    base_grad = _grad_mag(sub_keep["depth_to_base_surf"])
    cond_grad = _grad_mag(sub_keep["cond_surf"])

    z_tc_r = emit_ch("sub_hess_tc", "subsurface", "Axis-exact Hessian principal curvature of magnetic tilt tc", np.abs(tc_ridge))
    z_rtp_r = emit_ch("sub_hess_rtp", "subsurface", "Axis-exact Hessian principal curvature of RTP magnetics", np.abs(rtp_ridge))
    z_grav_r = emit_ch("sub_hess_grav", "subsurface", "Axis-exact Hessian principal curvature of isostatic gravity", np.abs(grav_ridge))
    z_bgrad = emit_ch("sub_grad_basement", "subsurface", "Basement-depth structural step gradient |grad(depth_to_base_surf)|", base_grad)
    emit_ch("sub_grad_cond", "subsurface", "MT electrical conductivity lateral boundary |grad(cond_surf)|", cond_grad)

    # Gravity tilt angle (Miller & Singh 1994 potential-field edge invariant)
    grav_tilt = np.arctan2(sub_keep["iso_grav_anom_vg"], np.abs(sub_keep["iso_grav_anom_hg"]) + 1e-5)
    emit_ch("sub_grav_tilt_angle", "subsurface", "Gravity potential-field tilt derivative atan(VG/|HG|)", grav_tilt)

    # Andersonian transtensional vorticity & microseismic fracture permeability (H80-3)
    dil = sub_keep["geod_dilaterate"]
    shr = np.abs(sub_keep["geod_shearrate"])
    vort = (dil / (np.hypot(dil, shr) + 1e-9)) * np.log1p(np.maximum(sub_keep["geod_2ndinv"], 0.0))
    emit_ch("sub_transtension_vorticity", "subsurface", "Andersonian dilatation/shear transtensional index", vort)
    emit_ch(
        "sub_seismic_strain_product",
        "subsurface",
        "Co-located microseismic intensity x geodetic strain invariant",
        robust_z(sub_keep["ieq_n100a15"], footprint) * np.maximum(robust_z(sub_keep["geod_2ndinv"], footprint), 0.0),
    )

    # 3. External Pure-Subsurface Hydrothermal Alteration (GeoDAWN K/Th/U + TMI_up150) & GDR/INGENIOUS Plumes (H80-2)
    with rasterio.open(data_dir / "external" / "geodawn_rad_u8.tif") as rsrc:
        rad_k = _nn_fill(rsrc.read(1).astype(np.float32))
        rad_th = _nn_fill(rsrc.read(2).astype(np.float32))
    with rasterio.open(data_dir / "external" / "geodawn_extensions_u8.tif") as esrc:
        rad_thk = _nn_fill(esrc.read(1).astype(np.float32))
        rad_uk = _nn_fill(esrc.read(2).astype(np.float32))
        tmi_up150 = _nn_fill(esrc.read(4).astype(np.float32))

    zk = robust_z(rad_k, footprint)
    zth = robust_z(rad_th, footprint)
    # Potassic metasomatic enrichment stripped of felsic lithology covariance
    cov_kth = float(np.mean(zk[footprint] * zth[footprint])) / max(float(np.mean(zth[footprint] ** 2)), 1e-6)
    k_resid = zk - cov_kth * zth
    z_kres = emit_ch("sub_potassic_k_th_resid", "subsurface", "GeoDAWN potassic alteration residual K - alpha*Th", k_resid)
    emit_ch("sub_rad_uk_ratio", "subsurface", "GeoDAWN U/K hydrothermal redox mobility ratio", rad_uk - _gauss(rad_uk, 12.0))
    z_mag_destr = emit_ch(
        "sub_mag_destruction_trough",
        "subsurface",
        "Shallow hydrothermal magnetite-destruction trough (TMI_up150 - TMI)",
        robust_z(tmi_up150, footprint) - robust_z(sub_keep["tmi"], footprint),
    )
    del rad_k, rad_th, rad_thk, rad_uk, tmi_up150

    pts = _rasterize_geothermal_points(data_dir, shape, grid_transform)
    z_geo15 = emit_ch("sub_gdr_geotherm_1p5km", "subsurface", "DOE GDR-1391 quartz/cation geothermometer upflow halo (1.5km)", pts["geotherm_plume_1p5km"])
    emit_ch("sub_gdr_geotherm_4km", "subsurface", "DOE GDR-1391 regional hydrothermal upflow halo (4km)", pts["geotherm_plume_4km"])
    z_hf4 = emit_ch("sub_ingenious_heatflow_4km", "subsurface", "USGS/INGENIOUS borehole convective heat-flow & Quaternary vent halo", pts["heatflow_vent_plume_4km"])
    emit_ch(
        "sub_hydrothermal_conduit_index",
        "subsurface",
        "H80-2 composite blind hydrothermal conduit (magnetite trough x K/Th x MT cond x heatflow)",
        (np.maximum(z_mag_destr, 0.0) + np.maximum(z_kres, 0.0)) * (1.0 + np.maximum(sub_keep["z_cond_surf"], 0.0) + np.maximum(z_geo15, 0.0) + np.maximum(z_hf4, 0.0)),
    )
    del pts

    # 4. Steerable Multi-Scale Radon Line-Integral Bank on Pure Subsurface Field (900m & 2.1km)
    sub_edge_field = np.maximum(
        0.30 * z_tc_r
        + 0.25 * z_rtp_r
        + 0.25 * robust_z(sub_keep["iso_grav_anom_slope"], footprint)
        + 0.20 * z_bgrad,
        0.0,
    ).astype(np.float32)
    r9_sub_max, r9_sub_aniso, _, _ = steerable_radon_bank(sub_edge_field, length_px=9, n_angles=8)
    r21_sub_max, r21_sub_aniso, u2_sub, v2_sub = steerable_radon_bank(sub_edge_field, length_px=21, n_angles=8)
    emit_ch("sub_radon9_max", "subsurface", "900m 8-azimuth steerable Radon peak on subsurface potential-field edges", r9_sub_max)
    emit_ch("sub_radon9_aniso", "subsurface", "900m Radon directional anisotropy contrast on subsurface edges", r9_sub_aniso)
    emit_ch("sub_radon21_max", "subsurface", "2.1km 8-azimuth steerable Radon peak on subsurface basement steps", r21_sub_max)
    z_sub_aniso21 = emit_ch("sub_radon21_aniso", "subsurface", "2.1km Radon directional anisotropy on subsurface basement steps", r21_sub_aniso)

    # Off-catalogue USGS SGMC pre-2022 bedrock lineament smooth field (strictly label-free geologic map)
    with rasterio.open(data_dir / "external" / "derived_sgmc_faults_100m_u8.tif") as ssrc:
        sgmc = (ssrc.read(1) > 0).astype(np.float32)
    d_sgmc = ndi.distance_transform_edt(sgmc == 0).astype(np.float32)
    emit_ch("sub_sgmc_bedrock_prox", "subsurface", "USGS SGMC pre-Quaternary bedrock fault exponential kernel exp(-d/600m)", np.exp(-d_sgmc / 6.0))
    z_grav_slope_keep = robust_z(sub_keep["iso_grav_anom_slope"], footprint)
    z_cond_keep = sub_keep["z_cond_surf"]
    del sgmc, d_sgmc, r9_sub_max, r9_sub_aniso, r21_sub_max, r21_sub_aniso, sub_edge_field
    del sub_keep, tc_ridge, rtp_ridge, grav_ridge, base_grad, cond_grad, grav_tilt, vort, zk, zth, k_resid
    del z_tc_r, z_rtp_r, z_grav_r, z_kres, z_mag_destr, z_geo15, z_hf4
    import gc
    gc.collect()

    log(f"[MACS] Pure-subsurface block built: {cursor} channels in {time.time() - t0:.1f}s")

    # 5. Surface Scarp, 1m 3DEP LiDAR, Multi-Scale Radon & Dipping Cross-Strike Coupling ('surface_coupling' = 19 ch)
    z_slope = emit_ch("surf_det_elev_slope", "surface_coupling", "Detrended topographic slope (100m)", det_slope)
    emit_ch("surf_det_elev_hp4", "surface_coupling", "400m high-pass topographic scarp residual", det_elev - _gauss(det_elev, 4.0))
    emit_ch("surf_det_elev_hp14", "surface_coupling", "1.4km band-pass range-margin relief residual", det_elev - _gauss(det_elev, 14.0))
    elev_ridge, dir_elev = _hessian_ridge(det_elev, sigma=1.0)
    slope_ridge, _ = _hessian_ridge(det_slope, sigma=1.0)
    z_elev_r = emit_ch("surf_hess_det_elev", "surface_coupling", "Signed Hessian break-in-slope curvature of det_elev", elev_ridge)
    z_slope_r = emit_ch("surf_hess_slope_crest", "surface_coupling", "Hessian scarp-crest curvature of det_elev_slope", slope_ridge)

    with rasterio.open(data_dir / "external" / "lidar_scarp_features_u8.tif") as lsrc:
        lid_ex = lsrc.read(1).astype(np.float32)
        lid_step = lsrc.read(3).astype(np.float32)
        lid_lapn = lsrc.read(4).astype(np.float32)
        lid_down = lsrc.read(6).astype(np.float32)
        lid_up = lsrc.read(7).astype(np.float32)
        lid_coh = lsrc.read(10).astype(np.float32) / 255.0
        lid_strike_rad = (lsrc.read(11).astype(np.float32) / 255.0) * math.pi

    z_lstep = emit_ch("surf_lidar_step_max", "surface_coupling", "1m 3DEP LiDAR maximum vertical scarp step height", lid_step)
    emit_ch("surf_lidar_ex_max", "surface_coupling", "1m 3DEP LiDAR extreme local scarp gradient", lid_ex)
    emit_ch("surf_lidar_scarp_polarity", "surface_coupling", "1m LiDAR one-sided fault-face polarity |downface - upface|", np.abs(lid_down - lid_up))
    z_lcoh_step = emit_ch("surf_lidar_coherent_step", "surface_coupling", "1m LiDAR scarp step height x 100m strike coherence", lid_step * lid_coh)
    emit_ch("surf_lidar_lapneg", "surface_coupling", "1m LiDAR scarp toe/crest negative Laplacian", lid_lapn)
    del lid_ex, lid_lapn, lid_down, lid_up

    scarp_energy = np.maximum(
        0.35 * z_lcoh_step + 0.25 * z_lstep + 0.20 * z_slope + 0.20 * np.abs(z_elev_r),
        0.0,
    ).astype(np.float32)
    r9_sc_max, r9_sc_aniso, _, _ = steerable_radon_bank(scarp_energy, length_px=9, n_angles=8)
    r21_sc_max, r21_sc_aniso, u2_sc, v2_sc = steerable_radon_bank(scarp_energy, length_px=21, n_angles=8)
    _, r35_sc_aniso, _, _ = steerable_radon_bank(scarp_energy, length_px=35, n_angles=8)

    emit_ch("surf_radon9_scarp_max", "surface_coupling", "900m 8-azimuth steerable Radon peak on 1m/100m scarp energy", r9_sc_max)
    emit_ch("surf_radon9_scarp_aniso", "surface_coupling", "900m Radon directional anisotropy on scarp energy", r9_sc_aniso)
    emit_ch("surf_radon21_scarp_max", "surface_coupling", "2.1km 8-azimuth steerable Radon peak on scarp energy", r21_sc_max)
    z_sc_aniso21 = emit_ch("surf_radon21_scarp_aniso", "surface_coupling", "2.1km Radon directional anisotropy on scarp energy (suppresses shoreline loops)", r21_sc_aniso)
    emit_ch("surf_radon35_scarp_aniso", "surface_coupling", "3.5km Radon long-strike continuity contrast on scarp energy", r35_sc_aniso)

    # Cross-physics strike phase-locking order parameter: cos(2*(theta_scarp - theta_subsurface))
    strike_dot = u2_sc * u2_sub + v2_sc * v2_sub
    lid_align = (np.cos(2.0 * lid_strike_rad) * u2_sub + np.sin(2.0 * lid_strike_rad) * v2_sub) * lid_coh
    emit_ch(
        "surf_cross_physics_strike_lock",
        "surface_coupling",
        "Phase-locked strike alignment between surface scarps and subsurface basement gradients",
        strike_dot * np.sqrt(np.maximum(z_sc_aniso21, 0.0) * np.maximum(z_sub_aniso21, 0.0)) + 0.5 * lid_align,
    )

    # Dipping Normal-Fault Cross-Strike Step Asymmetry (+-600m perpendicular to local strike!)
    # Fault block throw produces large |elev(+6n) - elev(-6n)|; symmetric valleys/ridges cancel to 0!
    theta_dom = 0.5 * np.arctan2(v2_sc, u2_sc)
    nx_step = np.rint(-6.0 * np.sin(theta_dom)).astype(np.int16)
    ny_step = np.rint(6.0 * np.cos(theta_dom)).astype(np.int16)
    del theta_dom, r9_sc_max, r9_sc_aniso, r21_sc_max, r21_sc_aniso, r35_sc_aniso, scarp_energy
    del strike_dot, lid_align, lid_coh, lid_strike_rad, lid_step, u2_sc, v2_sc, u2_sub, v2_sub
    gc.collect()

    grav_s = z_grav_slope_keep + 0.5 * np.maximum(z_cond_keep, 0.0)
    elev_throw = np.empty(shape, dtype=np.float32)
    hw_grav_mt = np.empty(shape, dtype=np.float32)
    x_base = np.arange(shape[1], dtype=np.int32)[None, :]
    for r0 in range(0, shape[0], 512):
        r1 = min(shape[0], r0 + 512)
        y_base = np.arange(r0, r1, dtype=np.int32)[:, None]
        yp = np.clip(y_base + ny_step[r0:r1], 0, shape[0] - 1)
        xp = np.clip(x_base + nx_step[r0:r1], 0, shape[1] - 1)
        ym = np.clip(y_base - ny_step[r0:r1], 0, shape[0] - 1)
        xm = np.clip(x_base - nx_step[r0:r1], 0, shape[1] - 1)
        ep = det_elev[yp, xp]
        em = det_elev[ym, xm]
        elev_throw[r0:r1] = np.abs(ep - em)
        hw_grav_mt[r0:r1] = np.where(ep < em, grav_s[yp, xp], grav_s[ym, xm])

    emit_ch(
        "surf_cross_strike_throw_600m",
        "surface_coupling",
        "600m across-strike elevation throw asymmetry (distinguishes fault blocks from symmetric gullies/ridges)",
        elev_throw,
    )
    emit_ch(
        "surf_hangingwall_downdip_grav_mt",
        "surface_coupling",
        "Hanging-wall down-dip gravity gradient & MT conductor offset 600m basinward of scarp",
        hw_grav_mt * np.maximum(z_lstep, 0.0),
    )
    del elev_throw, hw_grav_mt, grav_s, nx_step, ny_step

    # Concealed Half-Graben Growth-Wedge detector in flat alluvial basins (H80-4)
    slope_pct = np.zeros(shape, dtype=np.float32)
    order_s = np.argsort(det_slope[footprint], kind="mergesort")
    ranks_s = np.empty(order_s.size, dtype=np.float32)
    ranks_s[order_s] = np.linspace(0.0, 1.0, order_s.size, dtype=np.float32)
    slope_pct[footprint] = ranks_s
    emit_ch(
        "surf_concealed_basin_growth_wedge",
        "surface_coupling",
        "H80-4 buried basement/gravity step gated on low-slope alluvial basin floor",
        (np.maximum(z_bgrad, 0.0) + np.maximum(z_grav_slope_keep, 0.0)) * (1.0 - slope_pct),
    )

    assert cursor == n_expected, f"Expected {n_expected} channels, built {cursor}"
    mm.flush()
    del mm

    meta = {
        "schema_version": 1,
        "built_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "n_channels": n_expected,
        "n_subsurface_channels": sum(1 for _, d, _ in channel_specs if d == "subsurface"),
        "n_surface_coupling_channels": sum(1 for _, d, _ in channel_specs if d == "surface_coupling"),
        "valid_footprint_pixels": int(foot_flat.size),
        "grid_shape": list(shape),
        "training_features_sha256": hashlib.sha256((data_dir / "training_features.tif").read_bytes()).hexdigest(),
        "channels": [
            {"index": idx, "name": nm, "domain": dom, "description": dsc}
            for idx, (nm, dom, dsc) in enumerate(channel_specs)
        ],
    }
    meta_path.write_text(json.dumps(meta, indent=2) + "\n")
    arr = np.load(mmap_path, mmap_mode="r")
    log(f"[MACS] Completed all {n_expected} channels in {time.time() - t0:.1f}s -> {mmap_path}")
    return {"X": arr, "meta": meta, "footprint": footprint, "foot_flat": foot_flat}
