"""Verify Dilated Tversky Index (DTI) kernel and metric identities in vendor/gems57."""
import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "vendor"))
sys.path.insert(0, str(ROOT / "src"))

from gems57 import metric, evaluate_holdout


def test_dti_parameters_and_version():
    assert metric.ALPHA == 0.2
    assert metric.BETA == 0.8
    assert abs(metric.ALPHA + metric.BETA - 1.0) < 1e-12
    assert metric.PIXEL_M == 100.0
    assert evaluate_holdout.VERSION == "gems57-pooled-hide-v2"


def test_official_cross_example_tp_fn():
    # Competition problem page worked example: 5-px vertical truth crossed by 5-px horizontal prediction.
    gt = np.zeros((9, 9), dtype=bool)
    gt[2:7, 4] = True
    pr = np.zeros((9, 9), dtype=np.float64)
    pr[4, 2:7] = 1.0
    valid = np.ones((9, 9), dtype=bool)

    res = metric.dti_exact(pr, gt, valid=valid)
    assert abs(res["tp"] - 3.0) < 1e-4
    assert abs(res["fn"] - 2.0) < 1e-4
    assert round(res["dti"], 2) == 0.60


def test_tp_plus_fn_equals_truth_count():
    rng = np.random.default_rng(42)
    pr = (rng.random((32, 32)) < 0.04).astype(np.float64)
    gt = rng.random((32, 32)) < 0.05
    if not gt.any():
        gt[0, 0] = True
    valid = np.ones((32, 32), dtype=bool)

    res = metric.dti_exact(pr, gt, valid=valid)
    assert abs((res["tp"] + res["fn"]) - float(gt.sum())) < 1e-6
