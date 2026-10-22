"""Automated GeoTIFF validator tests for GEMSDOE80 (Lane H80-DPVC)."""
import hashlib
import json
import zipfile
from pathlib import Path

import numpy as np
import rasterio

ROOT = Path(__file__).resolve().parents[1]
TIF_STEM = "gems80-h80-dpvc-updip-vent-submod-41600-20261021-zeros"


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def test_shipped_tif_and_zip_automated_validator():
    sub_tif = ROOT / "submissions" / f"{TIF_STEM}.tif"
    sub_zip = ROOT / "submissions" / f"{TIF_STEM}.zip"
    sub_json = ROOT / "submissions" / f"{TIF_STEM}.json"

    dl_tif = ROOT / "docs" / "downloads" / f"{TIF_STEM}.tif"
    dl_zip = ROOT / "docs" / "downloads" / f"{TIF_STEM}.zip"
    dl_json = ROOT / "docs" / "downloads" / f"{TIF_STEM}.json"

    for p in (sub_tif, sub_zip, sub_json, dl_tif, dl_zip, dl_json):
        assert p.exists(), f"Missing artifact: {p}"

    # Byte-exact identity between submissions/ and docs/downloads/
    assert _sha256(sub_tif) == _sha256(dl_tif)
    assert _sha256(sub_zip) == _sha256(dl_zip)

    receipt = json.loads(sub_json.read_text())
    val = receipt["validator"]
    assert val["all_checks_passed"] is True
    assert all(bool(v) for v in val["checks"].values())
    assert len(val["checks"]) == 15
    assert receipt["sha256"] == _sha256(dl_tif)
    assert receipt["zip_sha256"] == _sha256(dl_zip)
    assert receipt["template_report"]["ok"] is True

    # Verify GeoTIFF directly against sample_submission.tif and labels.tif
    with rasterio.open(ROOT / "data" / "bridge" / "sample_submission.tif") as tpl:
        tpl_arr = tpl.read(1)
        foot = np.isfinite(tpl_arr)
        tpl_crs = tpl.crs
        tpl_transform = tpl.transform
        tpl_shape = (tpl.height, tpl.width)

    with rasterio.open(ROOT / "data" / "bridge" / "labels.tif") as lab:
        cat = lab.read(1) == 1

    with rasterio.open(dl_tif) as src:
        assert src.count == 1
        assert src.dtypes[0] == "float32"
        assert (src.height, src.width) == (3730, 3292) == tpl_shape
        assert src.crs == tpl_crs
        assert "32611" in str(src.crs)
        assert src.transform == tpl_transform
        assert src.nodata is None

        arr = src.read(1)
        assert arr.dtype == np.float32
        assert int(np.isnan(arr).sum()) == 0, "NaN detected in submission GeoTIFF"
        assert int(np.isinf(arr).sum()) == 0, "Inf detected in submission GeoTIFF"
        assert float(arr.min()) >= 0.0 and float(arr.max()) <= 1.0, "Values outside [0, 1]"
        assert float(np.abs(arr[~foot]).max()) == 0.0, "Non-zero value outside footprint"
        assert int(((arr > 0.0) & cat).sum()) == 0, "Dots placed on mapped catalogue faults"
        assert int((arr > 0.5).sum()) == 41600

    # Verify ZIP archive contains the exact TIF
    with zipfile.ZipFile(dl_zip, "r") as zf:
        names = zf.namelist()
        assert f"{TIF_STEM}.tif" in names
        data = zf.read(f"{TIF_STEM}.tif")
        assert hashlib.sha256(data).hexdigest() == _sha256(dl_tif)
