"""Regression tests for audited uncertainty, persistent budget and honest provenance."""
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from gems80.holdout import bootstrap_dti_ci, withhold_segments, score_holdout
from gems80.io_grid import load_template

ROOT = Path(__file__).resolve().parents[1]


def test_ci_perfect_and_empty_predictions():
    valid = np.ones((40, 40), bool)
    truth = np.zeros_like(valid)
    truth[::4, ::4] = True
    comp = np.zeros_like(valid, dtype=int)
    ci = bootstrap_dti_ci(truth, truth, comp, np.array([1]), valid)
    assert np.allclose(ci, (1, 1))
    assert bootstrap_dti_ci(np.zeros_like(truth), truth, comp, np.array([1]), valid) == (0, 0)


def test_ci_deterministic_with_multiplicity_and_no_component_deduplication():
    rng = np.random.default_rng(5)
    foot = np.ones((40, 40), bool)
    pred = rng.random(foot.shape) > .8
    truth = rng.random(foot.shape) > .9
    # Component ids are deliberately irrelevant: spatial contributions are resampled.
    a = bootstrap_dti_ci(pred, truth, np.ones_like(pred), np.array([1]), foot)
    b = bootstrap_dti_ci(pred, truth, np.zeros_like(pred), np.array([99]), foot)
    assert a == b
    assert a[0] < a[1]
    assert 0 <= a[0] <= a[1] <= 1


def test_invalid_split_and_uncalibrated_promotion():
    foot = np.ones((20, 20), bool)
    cat = np.eye(20, dtype=bool)
    with pytest.raises(ValueError):
        withhold_segments(cat, foot, buffer_px=-1)
    split = withhold_segments(cat, foot)
    result = score_holdout(cat, split, foot)
    assert not result['promotion_eligible']
    assert 'NOT spatially-blocked' in result['split_method']


def test_exhausted_budget_stops_before_data_access():
    proc = subprocess.run([sys.executable, 'scripts/run_h80_experiments.py'],
                          cwd=ROOT, capture_output=True, text=True)
    assert proc.returncode == 2
    assert 'budget exhausted' in proc.stdout


def test_published_card_fail_closed():
    for path in [ROOT/'research/run-card.json', ROOT/'docs/run-card.json']:
        card = json.loads(path.read_text())
        assert card['verdict'] == 'negative'
        assert not card['submission']['promotion_eligible']
        assert not card['submission']['organizer_confirmed']


def test_uniqueness_empty_registry_and_identical_dots():
    from gems80.uniqueness import check_registry
    foot = np.ones((20, 20), bool)
    dots = np.eye(20, dtype=bool)
    rng = np.random.default_rng(7)
    with pytest.raises(ValueError, match="Registry empty"):
        check_registry(dots, dots, foot, {}, rng)
    result = check_registry(dots, dots, foot, {'identical': dots}, rng)
    assert result['duplicate_lane']
    assert result['vs'][0]['final_dot_spearman'] == pytest.approx(1)


def test_template_numeric_nodata_is_not_footprint(tmp_path):
    import rasterio
    from rasterio.transform import from_origin
    path = tmp_path/'template.tif'
    with rasterio.open(path, 'w', driver='GTiff', height=2, width=2,
                       count=1, dtype='float32', nodata=-9999, crs='EPSG:32611',
                       transform=from_origin(0, 200, 100, 100)) as ds:
        ds.write(np.array([[0, -9999], [0, 0]], dtype='float32'), 1)
    assert load_template(path)['foot'].sum() == 3
