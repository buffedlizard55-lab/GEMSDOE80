#!/usr/bin/env python3
"""
Generate a GEMS Prize Challenge submission GeoTIFF.

Pipeline stage: inference -> submission raster.

Design (unsupervised, leakage-safe):
  * No training labels are used to build the surface, so there is no label
    leakage into the submission.
  * It REQUIRES data/sample_submission.tif (official DrivenData template) and
    data/training_features.tif. Without them it exits 2 and writes nothing: the
    grid is never guessed and no synthetic structure is ever generated. The
    feature grid must match the template's shape, CRS and geotransform exactly,
    otherwise it exits 3.

The probability surface combines label-free geophysical/geological bands of
training_features.tif (found by their band-description metadata, robust to band
order):
  - surface conductivity / depth to conductive base          (geothermal upflow)
  - dilatation rate, shear strain rate, 2nd strain invariant (crustal deformation)
  - isostatic gravity anomaly + slope                        (density contrasts)
  - magnetics: RTP/TMI anomaly, vertical & horizontal slope  (structural fabric)
  - earthquake density                                       (seismicity on structures)
Each band is percentile-normalized to [0,1] on valid pixels and combined into a
single [0,1] fault-probability layer; pixels where features are nodata are
written as NaN (null outside the valid footprint).

Outputs:
  submissions/<name>.tif        single-band float32 GeoTIFF, EPSG:32611, values in [0,1]
  submissions/<name>.tif.json   sidecar: sha256, grid report, validity flags
  submissions/<name>.zip        zip of the GeoTIFF (DrivenData accepts zips)

Links (official, manual review):
  https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/
  https://www.drivendata.org/competitions/306/competition-doe-gems/data/      (login required)
  https://github.com/drivendataorg/gems-prize-reference-solution
  https://epsg.io/32611
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import rasterio
from rasterio.transform import from_bounds

REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "data"
CANDIDATES = REPO / "artifacts" / "candidates"  # git-ignored; never the Pages folder
NOTE_MAX = 140

# Competition format (official DrivenData specification).
CRS = "EPSG:32611"          # UTM zone 11N
RESOLUTION_M = 100.0        # 100 m resolution
DTYPE = "float32"
NODATA = np.nan             # null outside the valid footprint
VAL_MIN, VAL_MAX = 0.0, 1.0

def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def band_index_by_keywords(descriptions, keywords):
    """Return 0-based band index whose description contains all keywords."""
    low = [d.lower() for d in descriptions]
    for i, d in enumerate(low):
        if all(k in d for k in keywords):
            return i
    return None


def robust_normalize(arr: np.ndarray, valid: np.ndarray, lo: float = 1.0, hi: float = 99.0):
    """Min-max normalize to [0,1] using percentiles on valid pixels."""
    out = np.zeros(arr.shape, dtype=np.float32)
    if not valid.any():
        return out, False
    vals = arr[valid]
    if vals.size == 0:
        return out, False
    p_lo, p_hi = np.percentile(vals, [lo, hi])
    if not np.isfinite(p_lo) or not np.isfinite(p_hi) or p_hi <= p_lo:
        return out, False
    n = np.clip((arr - p_lo) / (p_hi - p_lo), 0.0, 1.0)
    out[valid] = n[valid]
    return out, True


def build_surface_from_features(feature_path: Path):
    """Build [0,1] structural-propensity surface from feature bands."""
    with rasterio.open(feature_path) as src:
        meta = src.meta.copy()
        descriptions = [src.tags(i + 1).get("description", "") for i in range(src.count)]
        data = src.read()  # (bands, rows, cols)
        nodata = src.nodata
        transform = src.transform
        crs = src.crs
        height, width = src.height, src.width

    valid = np.ones((height, width), dtype=bool)
    for b in range(data.shape[0]):
        band = data[b]
        if nodata is not None:
            valid &= ~(band == nodata)
        valid &= np.isfinite(band)

    terms = []  # (name, normalized [0,1] array, ok)
    # conductivity (and depth-to-conductive-base)
    idx = band_index_by_keywords(descriptions, ["conductivity"])
    if idx is not None:
        n, ok = robust_normalize(data[idx], valid)
        terms.append(("surface_conductivity", n, ok))
    idx = band_index_by_keywords(descriptions, ["depth", "conductive"])
    if idx is not None:
        n, ok = robust_normalize(data[idx], valid)
        terms.append(("depth_to_conductive_base", n, ok))
    # strain rate magnitude
    strain_bands = []
    for kw in [["dilatation"], ["shear", "strain"], ["second", "invariant", "strain"]]:
        idx = band_index_by_keywords(descriptions, kw)
        if idx is not None:
            strain_bands.append(idx)
    if strain_bands:
        mag = np.zeros((height, width), dtype=np.float32)
        for idx in strain_bands:
            mag += np.abs(data[idx])
        n, ok = robust_normalize(mag, valid)
        terms.append(("strain_rate_magnitude", n, ok))
    # gravity magnitude (anomaly + slope)
    grav = []
    for kw in [["isostatic", "gravity", "anomaly"], ["slope", "isostatic"]]:
        idx = band_index_by_keywords(descriptions, kw)
        if idx is not None:
            grav.append(idx)
    if grav:
        mag = np.zeros((height, width), dtype=np.float32)
        for idx in grav:
            mag += np.abs(data[idx])
        n, ok = robust_normalize(mag, valid)
        terms.append(("gravity_magnitude", n, ok))
    # magnetics gradient magnitude (RTP/TMI anomaly + vertical + horizontal slope)
    mag_bands = []
    for kw in [["reduced-to-pole"], ["total", "magnetic", "intensity"], ["vertical", "slope"], ["horizontal", "slope"]]:
        idx = band_index_by_keywords(descriptions, kw)
        if idx is not None:
            mag_bands.append(idx)
    if mag_bands:
        mag = np.zeros((height, width), dtype=np.float32)
        for idx in mag_bands:
            mag += np.abs(data[idx])
        n, ok = robust_normalize(mag, valid)
        terms.append(("magnetics_gradient_magnitude", n, ok))
    # earthquake density
    idx = band_index_by_keywords(descriptions, ["density", "earthquake"])
    if idx is not None:
        n, ok = robust_normalize(data[idx], valid)
        terms.append(("earthquake_density", n, ok))

    surface = np.zeros((height, width), dtype=np.float32)
    used = 0
    for name, n, ok in terms:
        if ok:
            surface += n
            used += 1
    if used:
        surface = surface / used
    surface = np.clip(surface, VAL_MIN, VAL_MAX)
    # gentle emphasis on higher-confidence areas (still in [0,1])
    surface = surface ** 1.25
    surface = np.where(valid, surface, np.nan)
    return surface, transform, crs, height, width, terms


def write_geotiff(path: Path, surface: np.ndarray, transform, crs, height: int, width: int):
    profile = dict(
        driver="GTiff",
        dtype=DTYPE,
        count=1,
        height=height,
        width=width,
        transform=transform,
        crs=crs,
        nodata=NODATA,
        compress="lzw",
        tiled=False,  # GeoTIFF tiles must be multiples of 16; strips avoid that constraint
    )
    try:
        with rasterio.open(path, "w", **profile) as dst:
            dst.write(surface.astype(np.float32), 1)
            dst.set_band_description(1, "fault_probability")
    except Exception:
        # never leave a partial or unreadable raster where a submission could pick it up
        if path.exists():
            path.unlink()
        raise


def make_sidecar(name, path, surface, transform, crs, height, width, illustrative, terms):
    valid = ~np.isnan(surface)
    vals = surface[valid] if valid.any() else np.array([], dtype=np.float32)
    report = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "submission_file": path.name,
        "submission_name": name,
        "format": {
            "driver": "GTiff",
            "dtype": DTYPE,
            "count": 1,
            "crs": str(crs),
            "resolution_m": RESOLUTION_M,
            "height": height,
            "width": width,
            "transform": [transform.a, transform.b, transform.c, transform.d, transform.e, transform.f],
            "nodata": "NaN",
            "value_range_in_[0,1]": True,
        },
        "valid_pixel_count": int(valid.sum()),
        "valid_pixel_fraction": float(valid.mean()),
        "value_min": float(vals.min()) if vals.size else None,
        "value_max": float(vals.max()) if vals.size else None,
        "illustrative_fallback": illustrative,
        "note": (
            "Unsupervised geothermal structural-propensity prior; no training labels used; "
            "not copied from any previous GEMSDOE submission."
        ),
        "template_match_pending": illustrative,
    }
    return report


def write_zip(tif_path: Path):
    import zipfile
    zip_path = tif_path.with_suffix(".zip")
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as z:
        z.write(tif_path, arcname=tif_path.name)
    return zip_path


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--template", default=str(DATA / "sample_submission.tif"))
    ap.add_argument("--features", default=str(DATA / "training_features.tif"))
    ap.add_argument("--out", default=None)
    ap.add_argument("--name", default=None, help="unique submission name (default: timestamped)")
    ap.add_argument("--note", default="", help=f"short note, at most {NOTE_MAX} characters")
    args = ap.parse_args()
    if len(args.note) > NOTE_MAX:
        print(f"[BLOCKED] --note is {len(args.note)} characters; the limit is {NOTE_MAX}.")
        return 4

    template = Path(args.template)
    features = Path(args.features)

    # Refuse rather than fabricate: no template -> no grid; no features -> no evidence.
    missing = [str(p) for p in (template, features) if not p.is_file()]
    if missing:
        print("[BLOCKED] No submission written. Missing input(s):", ", ".join(missing))
        print("[BLOCKED] DrivenData data needs login; place competition files in data/ manually (see scripts/download_competition_data.sh).")
        print("[BLOCKED] This script never invents a grid or synthetic structure.")
        return 2

    CANDIDATES.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    name = args.name or f"gemsdoe80-structural-prior-v1-{stamp}"
    out = Path(args.out) if args.out else CANDIDATES / f"{name}.tif"

    with rasterio.open(template) as src:
        t_h, t_w, t_crs, t_tr = src.height, src.width, src.crs, src.transform
    surface, transform, crs, height, width, terms = build_surface_from_features(features)
    if (height, width) != (t_h, t_w) or transform != t_tr or str(crs) != str(t_crs):
        print("[FAIL] features grid does not match the official template grid; refusing to write.")
        print(f"       features: {height}x{width} {crs} {transform}")
        print(f"       template: {t_h}x{t_w} {t_crs} {t_tr}")
        return 3

    write_geotiff(out, surface, transform, crs, height, width)
    report = make_sidecar(name, out, surface, transform, crs, height, width, False, terms)
    report["note_for_submission_form"] = args.note
    report["status"] = "CANDIDATE-UNGATED: holdout gate not run. Do not download or submit."
    report["download_ready"] = False
    report["sha256"] = sha256_file(out)
    sidecar_path = out.with_suffix(".tif.json")
    sidecar_path.write_text(json.dumps(report, indent=2))
    zip_path = write_zip(out)
    report["zip_sha256"] = sha256_file(zip_path)
    sidecar_path.write_text(json.dumps(report, indent=2))

    print(f"[OK] Wrote {out}  sha256={report['sha256']}")
    print(f"[OK] Wrote {sidecar_path}")
    print(f"[OK] Wrote {zip_path}")
    print("[STATUS] CANDIDATE-UNGATED: this file is NOT approved for download or submission.")
    print("[NEXT] validate_submission.py (format), then the holdout gate; only then copy to docs/submissions/.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
