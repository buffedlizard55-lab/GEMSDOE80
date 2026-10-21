"""Unit tests for the distance-weighted Tversky implementation and the generator guard.

All expected values are computed by hand from the official definition
(problem description, "Performance metric": TP_w, FP_w, FN_w, DTI, R = 300 m = 3 px).
None of these tests is a competition score; the data are small synthetic rasters.
"""
from __future__ import annotations

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

try:
    import numpy as np
    import validate_submission as vs
except ImportError:  # pragma: no cover
    np = None


@unittest.skipIf(np is None, "numpy/scipy required")
class DTITests(unittest.TestCase):
    def test_single_gt_pixel_one_neighbour_hand_value(self) -> None:
        # GT at centre, prediction p=0.5 one pixel east. k(1) = 1 - 1/3 = 2/3.
        gt = np.zeros((7, 7)); gt[3, 3] = 1
        p = np.zeros((7, 7)); p[3, 4] = 0.5
        r = vs.dti_weighted(p, gt)
        tp = 0.5 * (2 / 3)
        fn = 1 - tp
        fp = 0.5 * (1 - 2 / 3)
        self.assertAlmostEqual(r["TP_w"], tp, places=9)
        self.assertAlmostEqual(r["FN_w"], fn, places=9)
        self.assertAlmostEqual(r["FP_w"], fp, places=9)
        self.assertAlmostEqual(r["DTI"], tp / (tp + 0.2 * fp + 0.8 * fn), places=9)  # = 0.370370...

    def test_perfect_prediction_scores_one_and_empty_scores_zero(self) -> None:
        gt = np.zeros((9, 9)); gt[2:7, 4] = 1
        self.assertAlmostEqual(vs.dti_weighted(gt.copy(), gt)["DTI"], 1.0, places=9)
        self.assertEqual(vs.dti_weighted(np.zeros((9, 9)), gt)["DTI"], 0.0)

    def test_fp_outside_kernel_has_full_penalty(self) -> None:
        gt = np.zeros((11, 11)); gt[5, 5] = 1
        p = np.zeros((11, 11)); p[0, 0] = 1.0  # far away: k=0 -> FP_w contribution = 1
        r = vs.dti_weighted(p, gt)
        self.assertAlmostEqual(r["FP_w"], 1.0, places=9)
        self.assertAlmostEqual(r["FN_w"], 1.0, places=9)
        self.assertEqual(r["TP_w"], 0.0)

    def test_fn_counts_only_gt_pixels(self) -> None:
        # Regression test for the bug where FN_w summed (1 - 0) over every background pixel.
        gt = np.zeros((20, 20)); gt[10, 10] = 1
        r = vs.dti_weighted(np.zeros((20, 20)), gt)
        self.assertAlmostEqual(r["FN_w"], 1.0, places=9)


@unittest.skipIf(np is None, "numpy/rasterio required")
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
