"""Regression: generator round trip on a non-multiple-of-16 grid (tiling bug fixed 2026-10-21)."""
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

rasterio = pytest.importorskip("rasterio")
from rasterio.transform import from_origin  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
GEN = REPO / "scripts" / "generate_submission.py"
VAL = REPO / "scripts" / "validate_submission.py"
TR = from_origin(243350, 4508550, 100, 100)


def _write(path, arr, desc=None):
    arr = np.asarray(arr, dtype=np.float32)
    count = 1 if arr.ndim == 2 else arr.shape[0]
    h, w = arr.shape[-2:]
    with rasterio.open(path, "w", driver="GTiff", height=h, width=w, count=count, dtype="float32",
                       crs="EPSG:32611", transform=TR) as ds:
        for b in range(count):
            ds.write(arr if arr.ndim == 2 else arr[b], b + 1)
            if desc:
                ds.update_tags(b + 1, description=desc[b])


def test_generator_roundtrip_non_multiple_of_16(tmp_path):
    rng = np.random.default_rng(1)
    tpl = tmp_path / "sample.tif"
    feat = tmp_path / "features.tif"
    out = tmp_path / "out.tif"
    _write(tpl, np.zeros((30, 40)))
    _write(feat, rng.normal(5, 1, size=(3, 30, 40)),
           desc=["surface conductivity", "dilatation rate", "earthquake density"])
    g = subprocess.run([sys.executable, str(GEN), "--template", str(tpl), "--features", str(feat),
                        "--out", str(out)], capture_output=True, text=True)
    assert g.returncode == 0, g.stdout + g.stderr
    v = subprocess.run([sys.executable, str(VAL), str(out), "--template", str(tpl)],
                       capture_output=True, text=True)
    assert v.returncode == 0, v.stdout + v.stderr
    with rasterio.open(out) as d:
        data = d.read(1)
        assert d.dtypes == ("float32",) and d.count == 1
        assert d.crs.to_epsg() == 32611
        assert data.shape == (30, 40)
        assert np.nanmin(data) >= 0 and np.nanmax(data) <= 1


def test_generator_refuses_grid_mismatch_and_leaves_no_file(tmp_path):
    tpl = tmp_path / "sample.tif"
    feat = tmp_path / "features.tif"
    out = tmp_path / "out.tif"
    _write(tpl, np.zeros((30, 40)))
    _write(feat, np.zeros((2, 20, 40)))
    g = subprocess.run([sys.executable, str(GEN), "--template", str(tpl), "--features", str(feat),
                        "--out", str(out)], capture_output=True, text=True)
    assert g.returncode == 3
    assert not out.exists()
