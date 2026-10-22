"""Official competition grid I/O. Geometry is read from sample_submission.tif, never invented."""

from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
import rasterio

NODATA_F32 = np.float32(-3.4028234663852886e38)
EXPECTED_HEIGHT = 3730
EXPECTED_WIDTH = 3292
EXPECTED_TRANSFORM = (100.0, 0.0, 243350.0, 0.0, -100.0, 4508550.0)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_template(path: Path) -> dict:
    with rasterio.open(path) as src:
        arr = src.read(1)
        foot = np.isfinite(arr) & (src.read_masks(1) > 0)
        meta = {
            "driver": "GTiff",
            "dtype": "float32",
            "count": 1,
            "height": src.height,
            "width": src.width,
            "crs": src.crs,
            "transform": src.transform,
            "compress": "deflate",
            "tiled": True,
            "blockxsize": 256,
            "blockysize": 256,
        }
        report = {
            "height": src.height,
            "width": src.width,
            "crs": str(src.crs),
            "transform": [float(src.transform.a), float(src.transform.b), float(src.transform.c),
                          float(src.transform.d), float(src.transform.e), float(src.transform.f)],
            "nodata": src.nodata,
            "footprint_px": int(foot.sum()),
            "dtype": src.dtypes[0],
        }
    return {"arr": arr, "foot": foot, "meta": meta, "report": report}


def load_labels(path: Path, foot: np.ndarray) -> np.ndarray:
    with rasterio.open(path) as src:
        a = src.read(1)
    return np.isfinite(a) & (a >= 1) & foot


def load_band(path: Path, band: int):
    with rasterio.open(path) as src:
        a = src.read(band).astype(np.float32)
        desc = src.descriptions[band - 1] if src.descriptions else f"band{band}"
        tags = src.tags(band)
    valid = np.isfinite(a) & (a > NODATA_F32 * np.float32(0.5))
    return a, valid, desc, tags


def feature_descriptions(path: Path) -> list[str]:
    with rasterio.open(path) as src:
        return [d or f"band{i+1}" for i, d in enumerate(src.descriptions)]
