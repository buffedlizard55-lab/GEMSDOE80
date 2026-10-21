"""Official-metric self-checks (DrivenData page 967 worked example)."""
from __future__ import annotations

import math
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from gems80.metric import dti_from_components, worked_example_dti


class MetricTests(unittest.TestCase):
    def test_worked_example(self):
        # TP_w=3.00, FP_w=1.89, FN_w=2.00 -> 3/(3+0.2*1.89+0.8*2)=0.602651...
        got = worked_example_dti()
        expected = 3.0 / (3.0 + 0.2 * 1.89 + 0.8 * 2.0)
        self.assertTrue(math.isclose(got, expected, rel_tol=1e-12))
        self.assertTrue(math.isclose(got, 0.60, rel_tol=0, abs_tol=0.003))

    def test_perfect_match(self):
        self.assertTrue(math.isclose(dti_from_components(10.0, 0.0, 0.0), 1.0, rel_tol=0, abs_tol=1e-12))

    def test_no_prediction(self):
        # all FN: DTI = 0 / (0 + 0 + 0.8*|G|) = 0
        self.assertEqual(dti_from_components(0.0, 0.0, 10.0), 0.0)


if __name__ == "__main__":
    unittest.main()
