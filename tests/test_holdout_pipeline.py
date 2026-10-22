"""Smoke-test hide-and-recover leakage canaries on SYNTHETIC rasters.

These toy arrays are not GEMS data. The test verifies that per-feature AUC includes both
withheld positives and negative pixels in the scoring domain. It does not score a static
full-catalogue prediction as if it were a rebuilt holdout result.
"""
from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import rasterio
from rasterio.transform import from_origin

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from gems80.holdout import spatial_block_dti_ci, withhold_segments  # noqa: E402
from run_h80_experiments import leakage_canary  # noqa: E402
from evaluate_holdout import evaluate  # noqa: E402


class HoldoutCanarySmokeTests(unittest.TestCase):
    def test_canary_flags_label_copy_and_not_noise_on_full_domain(self) -> None:
        rng = np.random.default_rng(1)
        h, w = 120, 120
        foot = np.ones((h, w), dtype=bool)
        labels = np.zeros((h, w), dtype=bool)
        for i in range(6):
            labels[10 + 15 * i, 5:115] = True
        leak = labels.astype(np.float32) + 0.01 * rng.standard_normal((h, w)).astype(np.float32)
        noise = rng.random((h, w), dtype=np.float32)
        valid = np.ones((h, w), dtype=bool)

        split = withhold_segments(labels, foot, hold_frac=0.3, buffer_px=3, seed=7)
        canary = leakage_canary(
            {"label_copy": (leak, valid), "noise": (noise, valid)},
            split["hidden"], foot, split["visible"],
        )
        aucs = {row["feature"]: row["auc"] for row in canary["per_feature"]}
        self.assertGreater(aucs["label_copy"], 0.99)
        self.assertLess(abs(aucs["noise"] - 0.5), 0.1)
        self.assertTrue(canary["leakage_flag"])

    def test_evaluator_rebuilds_from_visible_faults_and_writes_no_tif(self) -> None:
        h = w = 120
        labels = np.zeros((h, w), dtype=np.float32)
        for i in range(6):
            labels[10 + 15 * i, 5:115] = 1
        rng = np.random.default_rng(2)
        features = rng.normal(size=(17, h, w)).astype(np.float32)
        # The H80-C model's geological fields are finite and nonconstant.
        features[7] = rng.normal(size=(h, w)).astype(np.float32)
        features[1] = rng.normal(size=(h, w)).astype(np.float32)
        features[5] = rng.normal(size=(h, w)).astype(np.float32)
        transform = from_origin(0, 12000, 100, 100)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            template_path = root / "template.tif"
            labels_path = root / "labels.tif"
            features_path = root / "features.tif"
            profile = dict(driver="GTiff", height=h, width=w, count=1, dtype="float32",
                           crs="EPSG:32611", transform=transform)
            for path, array in ((template_path, np.zeros((h, w), dtype=np.float32)),
                                (labels_path, labels)):
                with rasterio.open(path, "w", **profile) as dst:
                    dst.write(array, 1)
            with rasterio.open(features_path, "w", **{**profile, "count": 17}) as dst:
                dst.write(features)
                for band in range(1, 18):
                    dst.set_band_description(band, f"synthetic_band_{band}")
            result = evaluate(template_path, labels_path, features_path)
        self.assertEqual(result["label"], "HOLDOUT-DTI")
        self.assertEqual(result["submission_tif_written"], False)
        self.assertGreater(result["withheld_positive_count"], 0)
        self.assertTrue(result["leakage_canary"]["available"])

    def test_spatial_block_ci_is_deterministic_and_requires_enough_blocks(self) -> None:
        small = np.ones((120, 120), dtype=bool)
        dots = np.zeros_like(small)
        truth = np.zeros_like(small)
        truth[10:110:10, 10:110] = True
        lo, hi, n_blocks = spatial_block_dti_ci(dots, truth, small, block_px=50, n_boot=100, seed=4)
        self.assertIsNone(lo)
        self.assertIsNone(hi)
        self.assertLess(n_blocks, 20)

        valid = np.ones((1000, 1000), dtype=bool)
        truth = np.zeros_like(valid)
        truth[100:900:10, 100:900] = True
        pred = np.zeros_like(valid)
        pred[100:900:10, 101:901] = truth[100:900:10, 100:900]
        first = spatial_block_dti_ci(pred, truth, valid, block_px=50, n_boot=100, seed=4)
        second = spatial_block_dti_ci(pred, truth, valid, block_px=50, n_boot=100, seed=4)
        self.assertGreaterEqual(first[2], 20)
        self.assertEqual(first, second)
        self.assertLessEqual(0.0, first[0])
        self.assertLessEqual(first[0], first[1])
        self.assertLessEqual(first[1], 1.0)


if __name__ == "__main__":
    unittest.main()
