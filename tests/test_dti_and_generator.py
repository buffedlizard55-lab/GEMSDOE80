"""Unit tests for the shared binary DTI implementation and generator guard.

All arrays are synthetic; these are correctness tests, not competition scores.
"""
from __future__ import annotations

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from gems80.metric import dti_binary  # noqa: E402


class DTITests(unittest.TestCase):
    def test_single_gt_pixel_one_neighbour_hand_value(self) -> None:
        gt = np.zeros((7, 7), dtype=bool); gt[3, 3] = True
        pred = np.zeros_like(gt); pred[3, 4] = True
        r = dti_binary(pred, gt)
        tp = 2 / 3
        fn = 1 - tp
        fp = 1 - 2 / 3
        self.assertAlmostEqual(r["TP_w"], tp, places=9)
        self.assertAlmostEqual(r["FN_w"], fn, places=9)
        self.assertAlmostEqual(r["FP_w"], fp, places=9)
        self.assertAlmostEqual(r["DTI"], tp / (tp + 0.2 * fp + 0.8 * fn), places=9)

    def test_perfect_prediction_scores_one_and_empty_scores_zero(self) -> None:
        gt = np.zeros((9, 9), dtype=bool); gt[2:7, 4] = True
        self.assertAlmostEqual(dti_binary(gt.copy(), gt)["DTI"], 1.0, places=9)
        self.assertEqual(dti_binary(np.zeros_like(gt), gt)["DTI"], 0.0)

    def test_fp_outside_kernel_has_full_penalty(self) -> None:
        gt = np.zeros((11, 11), dtype=bool); gt[5, 5] = True
        pred = np.zeros_like(gt); pred[0, 0] = True
        r = dti_binary(pred, gt)
        self.assertAlmostEqual(r["FP_w"], 1.0, places=9)
        self.assertAlmostEqual(r["FN_w"], 1.0, places=9)
        self.assertEqual(r["TP_w"], 0.0)

    def test_fn_counts_only_gt_pixels(self) -> None:
        gt = np.zeros((20, 20), dtype=bool); gt[10, 10] = True
        r = dti_binary(np.zeros_like(gt), gt)
        self.assertAlmostEqual(r["FN_w"], 1.0, places=9)


class GeneratorGuardTests(unittest.TestCase):
    def test_generator_refuses_without_official_template(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            result = subprocess.run(
                [sys.executable, str(REPO / "scripts" / "generate_submission.py"),
                 "--template", str(Path(tmp) / "missing_sample.tif"),
                 "--features", str(Path(tmp) / "missing_features.tif"),
                 "--out", str(Path(tmp) / "should_not_exist.tif")],
                capture_output=True, text=True, check=False,
            )
            self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
            self.assertFalse((Path(tmp) / "should_not_exist.tif").exists())
            self.assertIn("BLOCKED", result.stdout)


if __name__ == "__main__":
    unittest.main()
