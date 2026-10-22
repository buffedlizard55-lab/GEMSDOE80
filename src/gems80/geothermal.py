"""Hydrothermal Upflow Plume, Geothermometer & Potassic/Demag Alteration Inversion (H80-G3).

Strictly label-free geothermal & alteration channel synthesis from official public sources:
  1. DOE GDR 1391 (DOI: 10.15121/1881483): well/spring measured temperatures & silica/cation
     chemical geothermometers (`dist_known_fault_px` column explicitly stripped & verified absent).
  2. DOE GDR 1391 Quaternary volcanic vents (`dist_known_fault_px` stripped).
  3. USGS ScienceBase Great Basin Heat Flow release (DOI: 10.5066/P9BZPVUC, DeAngelo et al. 2022):
     in-footprint heat-flow wells & positive convective heat-flow residuals (`hf_resid`, `hf_used`).
  4. Reid et al. (1990) 3D Magnetic Euler Deconvolution shallow contact solutions (`depth <= 450 m`).
  5. USGS GeoDAWN airborne radiometrics (DOI: 10.5066/P93LGLVQ): local high-pass K/Th potassic
     boiling-zone metasomatism, U/Th carbonate sinter/travertine anomaly, and magnetite-destruction
     hydrothermal troughs (`-rtp` / `-tmi_vg` x `cond_surf`).
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from scipy.ndimage import gaussian_filter

from .geomechanics import fill_median, multiscale_hessian_ridge


def build_geothermal_channels(
    data_dir: Path,
    transform,
    shape: tuple[int, int],
    rtp: np.ndarray,
    tmi_vg: np.ndarray,
    cond_surf: np.ndarray,
    valid: np.ndarray,
) -> dict[str, np.ndarray]:
    h, w = shape
    ext = data_dir / "external"

    # 1. GDR 1391 Well/Spring temperatures + Quartz/Chalcedony/Na-K-Ca geothermometers
    core_impulse = np.zeros(shape, dtype=np.float32)
    ws_csv = ext / "gdr_wellspring_in_footprint.csv"
    if ws_csv.is_file():
        df_ws = pd.read_csv(ws_csv)
        if "dist_known_fault_px" in df_ws.columns:
            df_ws = df_ws.drop(columns=["dist_known_fault_px"])
        assert "dist_known_fault_px" not in df_ws.columns

        rows = np.clip(np.round(df_ws["row"].to_numpy()).astype(np.int32), 0, h - 1)
        cols = np.clip(np.round(df_ws["col"].to_numpy()).astype(np.int32), 0, w - 1)
        t_meas = np.nan_to_num(df_ws["temp_c"].to_numpy(dtype=np.float32), nan=0.0)
        gt_max = np.nanmax(
            np.column_stack([
                np.nan_to_num(df_ws["geothermquartz_c"].to_numpy(dtype=np.float32), nan=0.0),
                np.nan_to_num(df_ws["geothermchalc_c"].to_numpy(dtype=np.float32), nan=0.0),
                np.nan_to_num(df_ws["geothermcat_c"].to_numpy(dtype=np.float32), nan=0.0),
            ]),
            axis=1,
        )
        # Thermal anomaly weight: direct surface discharge excess (>22 C) + deep reservoir geothermometer (>85 C)
        excess = np.maximum(t_meas - 22.0, 0.0) / 40.0 + 0.65 * np.maximum(gt_max - 85.0, 0.0) / 80.0
        pos = excess > 0.0
        np.maximum.at(core_impulse, (rows[pos], cols[pos]), np.clip(excess[pos], 0.0, 5.0))

    vv_csv = ext / "gdr_volcanic_vents_in_footprint.csv"
    if vv_csv.is_file():
        df_vv = pd.read_csv(vv_csv)
        if "dist_known_fault_px" in df_vv.columns:
            df_vv = df_vv.drop(columns=["dist_known_fault_px"])
        assert "dist_known_fault_px" not in df_vv.columns
        vr = np.clip(np.round(df_vv["row"].to_numpy()).astype(np.int32), 0, h - 1)
        vc = np.clip(np.round(df_vv["col"].to_numpy()).astype(np.int32), 0, w - 1)
        np.maximum.at(core_impulse, (vr, vc), 1.5)

    # Multi-scale upflow conduit core (sigma=4.5 px = 450 m) + shallow outflow plume (sigma=18 px = 1.8 km)
    gdr_plume = (
        2.5 * gaussian_filter(core_impulse, sigma=4.5, mode="constant")
        + 1.0 * gaussian_filter(core_impulse, sigma=18.0, mode="constant")
    )

    def _as_float(v) -> float:
        if v is None or isinstance(v, bool):
            return float("nan")
        try:
            return float(v)
        except (TypeError, ValueError):
            return float("nan")

    # 2. USGS ScienceBase Conductive & Convective Residual Heat Flow (DOI: 10.5066/P9BZPVUC)
    hf_impulse = np.zeros(shape, dtype=np.float32)
    hf_json = ext / "sb_heat_flow_in_footprint.json"
    if hf_json.is_file():
        payload = json.loads(hf_json.read_text())
        for rec in payload.get("records", []):
            parts = rec.get("parts") or []
            if not parts or not parts[0]:
                continue
            # Use single-point wells from USGS_gbHeatFlowWells_wEstimates.shp / Wells.shp
            if len(parts) == 1 and len(parts[0]) == 1:
                x_m, y_m = parts[0][0]
                c_idx = int(round((x_m - transform.c) / transform.a))
                r_idx = int(round((y_m - transform.f) / transform.e))
                if 0 <= r_idx < h and 0 <= c_idx < w:
                    hf_res = _as_float(rec.get("hf_resid"))
                    hf_meas = _as_float(rec.get("hf_meas"))
                    val = 0.0
                    if np.isfinite(hf_res) and hf_res > 12.0:
                        val += min(np.log1p((hf_res - 12.0) / 25.0), 4.0)
                    if np.isfinite(hf_meas) and hf_meas > 85.0:
                        val += min(np.log1p((hf_meas - 85.0) / 40.0), 3.0)
                    if val > hf_impulse[r_idx, c_idx]:
                        hf_impulse[r_idx, c_idx] = float(val)
    sb_heatflow = (
        2.2 * gaussian_filter(hf_impulse, sigma=5.0, mode="constant")
        + 1.0 * gaussian_filter(hf_impulse, sigma=20.0, mode="constant")
    )

    # 3. Shallow low-dispersion 3D Magnetic Euler Deconvolution Contact Clusters (Reid et al. 1990)
    euler_field = np.zeros(shape, dtype=np.float32)
    euler_csv = ext / "h31_1_euler_clusters.csv"
    if euler_csv.is_file():
        dfe = pd.read_csv(euler_csv)
        good = (
            (dfe["median_depth_m"] <= 500.0)
            & (dfe["depth_mad_m"] <= 80.0)
            & (dfe["n_solutions"] >= 6)
        )
        dfe = dfe.loc[good]
        er = np.clip(np.round(dfe["row"].to_numpy()).astype(np.int32), 0, h - 1)
        ec = np.clip(np.round(dfe["col"].to_numpy()).astype(np.int32), 0, w - 1)
        ew = np.clip(np.log1p(dfe["n_solutions"].to_numpy(dtype=np.float32)) * (550.0 - dfe["median_depth_m"].to_numpy(dtype=np.float32)) / 400.0, 0.1, 4.0)
        np.maximum.at(euler_field, (er, ec), ew)
        euler_field = gaussian_filter(euler_field, sigma=2.2, mode="constant")

    # 4. GeoDAWN Airborne Radiometric Potassic (K/Th) & Travertine/Redox (U/Th) Residuals
    rad_kth = np.zeros(shape, dtype=np.float32)
    rad_uth = np.zeros(shape, dtype=np.float32)
    rad_tif = ext / "geodawn_rad_u8.tif"
    ext_tif = ext / "geodawn_extensions_u8.tif"
    if rad_tif.is_file() and ext_tif.is_file():
        with rasterio.open(rad_tif) as sr:
            k_band = sr.read(1).astype(np.float32)
            th_band = sr.read(2).astype(np.float32)
            rad_ok = valid & (k_band > 0) & (th_band > 0)
        with rasterio.open(ext_tif) as se:
            thk_u8 = se.read(1).astype(np.float32)
            uth_u8 = se.read(3).astype(np.float32)
        # High K relative to Th marks boiling-zone adularia/sericite alteration; subtract 3 km regional lithology
        k_over_th = np.where(rad_ok, (k_band + 8.0) / (th_band + 8.0) - 0.004 * thk_u8, np.nan)
        k_f = fill_median(k_over_th, rad_ok)
        rad_kth = np.where(rad_ok, k_f - gaussian_filter(k_f, sigma=25.0, mode="nearest"), 0.0).astype(np.float32)

        u_f = fill_median(np.where(rad_ok, uth_u8, np.nan), rad_ok)
        gy_u, gx_u = np.gradient(gaussian_filter(u_f, sigma=1.2, mode="nearest"))
        rad_uth = np.where(
            rad_ok,
            0.6 * np.maximum(u_f - gaussian_filter(u_f, sigma=20.0, mode="nearest"), 0.0) + 0.4 * np.hypot(gy_u, gx_u),
            0.0,
        ).astype(np.float32)

    # 5. Hydrothermal Magnetite-Destruction Trough x Clay-Cap Electrical Conductivity
    demag_trough = multiscale_hessian_ridge(rtp, valid, sigmas=(1.5, 3.0), polarity="trough")
    cond_f = fill_median(cond_surf, valid)
    cond_anom = np.maximum(cond_f - gaussian_filter(cond_f, sigma=30.0, mode="nearest"), 0.0)
    hydro_demag = demag_trough * (0.4 + 0.6 * np.tanh(cond_anom / (np.std(cond_anom[valid]) + 1e-6)))

    return {
        "gdr_geotherm_plume": np.where(valid, gdr_plume, 0.0).astype(np.float32),
        "sb_heatflow_anom": np.where(valid, sb_heatflow, 0.0).astype(np.float32),
        "euler_contact_density": np.where(valid, euler_field, 0.0).astype(np.float32),
        "rad_kth_alter": np.where(valid, rad_kth, 0.0).astype(np.float32),
        "rad_uth_travertine": np.where(valid, rad_uth, 0.0).astype(np.float32),
        "hydro_demag_conduit": np.where(valid, hydro_demag, 0.0).astype(np.float32),
    }
