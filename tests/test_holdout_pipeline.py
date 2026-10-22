"""End-to-end smoke test of the internal holdout + leakage canary on SYNTHETIC rasters.

The rasters are random toy data, not GEMS competition data. The test checks plumbing and
the canary's behaviour: a band copied from the labels must be flagged (AUC near 1), a random
band must not (AUC near 0.5).
"""
from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]

try:
    import numpy as np
    import rasterio
    from rasterio.transform import from_origin
except ImportError:  # pragma: no cover
    np = None


@unittest.skipIf(np is None, "numpy/rasterio required")
class HoldoutSmokeTests(unittest.TestCase):
    def test_canary_flags_label_copy_and_not_noise(self) -> None:
        rng = np.random.default_rng(1)
        H, W = 120, 120
        gt = np.zeros((H, W), dtype=np.float32)
        for i in range(6):  # six straight synthetic fault traces
            gt[10 + 15 * i, 5:115] = 1
        feat_leak = gt + 0.01 * rng.standard_normal((H, W)).astype(np.float32)
        feat_noise = rng.random((H, W)).astype(np.float32)
        pred = np.clip(gt * 0.6 + 0.1 * rng.random((H, W)), 0, 1).astype(np.float32)
        tr = from_origin(0, 12000, 100, 100)
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)

            def w(path, arr, desc=None):
                with rasterio.open(path, "w", driver="GTiff", height=H, width=W, count=arr.shape[0] if arr.ndim == 3 else 1,
                                   dtype="float32", crs="EPSG:32611", transform=tr) as ds:
                    if arr.ndim == 3:
                        for b in range(arr.shape[0]):
                            ds.write(arr[b], b + 1)
                            ds.update_tags(b + 1, description=desc[b])
                    else:
                        ds.write(arr, 1)
            w(root / "template.tif", np.zeros((H, W), np.float32))
            w(root / "pred.tif", pred)
            w(root / "labels.tif", gt)
            w(root / "features.tif", np.stack([feat_leak, feat_noise]), desc=["leak_copy", "noise"])
            res = subprocess.run(
                [sys.executable, str(REPO / "scripts" / "validate_submission.py"), str(root / "pred.tif"),
                 "--template", str(root / "template.tif"), "--labels", str(root / "labels.tif"),
                 "--features", str(root / "features.tif")],
                capture_output=True, text=True, check=False,
            )
            report = json.loads(res.stdout)
        self.assertTrue(report["format_valid"], report["errors"])
        self.assertTrue(report["holdout"]["available"])
        self.assertIn("ci95", report["holdout"])
        aucs = {a["description"]: a["auc"] for a in report["leakage_canary"]["per_band_auc"]}
        self.assertGreater(aucs["leak_copy"], 0.99)
        self.assertLess(abs(aucs["noise"] - 0.5), 0.1)
        self.assertTrue(report["leakage_canary"]["leakage_flag"])
        self.assertFalse(report["organizer_confirmed"])


if __name__ == "__main__":
    unittest.main()
