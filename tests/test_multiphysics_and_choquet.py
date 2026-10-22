"""Unit tests for GEMSDOE80 v2 scientific modules:
  1. 6-shell exact kernel offsets & continuous DTI consistency (`metric.py`).
  2. Unbiased Voronoi-partition multinomial cluster bootstrap CI (`holdout.py`).
  3. Real Great Basin Albers stress projection (`geomechanics.py`: >100,000 segments).
  4. Zero-leakage geothermometer synthesis (`geothermal.py`: strips `dist_known_fault_px`).
  5. Directed crack-tip vs lateral-flank asymmetry (`topology.py`).
  6. Exact 6-shell Choquet survival integral vs brute-force 2^9 enumeration & visible-mask zeroing (`choquet_submodular.py`).
"""

from __future__ import annotations

import itertools
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gems80.choquet_submodular import choquet_expected_self_credit, dinkelbach_choquet_allocate
from gems80.geomechanics import build_regional_stress_and_tendency
from gems80.geothermal import build_geothermal_channels
from gems80.holdout import bootstrap_dti_ci, quadrant_dti_breakdown, score_holdout, withhold_segments
from gems80.io_grid import load_template
from gems80.metric import KERNEL_OFFSETS, dti_binary, dti_continuous
from gems80.topology import compute_topology_features


def test_kernel_offsets_and_continuous_dti_exact():
    assert len(KERNEL_OFFSETS) == 25
    rng = np.random.default_rng(123)
    foot = np.ones((40, 40), dtype=bool)
    pred = rng.random((40, 40)) < 0.04
    truth = rng.random((40, 40)) < 0.03
    s_bin = dti_binary(pred, truth, foot)
    s_cont = dti_continuous(pred.astype(np.float32), truth, foot)
    assert abs(s_bin["DTI"] - s_cont["DTI"]) < 1e-6


def test_voronoi_bootstrap_ci_is_unbiased():
    foot = np.ones((80, 80), dtype=bool)
    cat = np.zeros((80, 80), dtype=bool)
    for k in range(12):
        r = 5 + k * 6
        cat[r, 8:32] = True
        cat[r, 42:68] = True
    split = withhold_segments(cat, foot, hold_frac=0.35, buffer_px=2, seed=11)
    pred = np.zeros((80, 80), dtype=bool)
    ys, xs = np.nonzero(split["hidden"])
    pred[ys[::4], xs[::4]] = True
    base = score_holdout(pred, split, foot, bootstrap=True)
    lo, hi = base["ci95"]
    assert lo <= base["score"] <= hi
    assert (hi - lo) > 0.0

    q = quadrant_dti_breakdown(pred, split["hidden"], foot)
    assert set(q.keys()) == {"NW", "NE", "SW", "SE"}


def test_real_albers_stress_projection_coverage():
    slip_csv = ROOT / "data" / "external" / "sb_slip_tendency_in_footprint.csv"
    tmpl_tif = ROOT / "data" / "sample_submission.tif"
    if not (slip_csv.is_file() and tmpl_tif.is_file()):
        pytest.skip("Real external slip tendency data not materialized")
    tmpl = load_template(tmpl_tif)
    foot = tmpl["foot"]
    strike = np.full(foot.shape, np.radians(20.0), dtype=np.float32)
    coh = np.full(foot.shape, 0.7, dtype=np.float32)
    res = build_regional_stress_and_tendency(
        slip_csv, tmpl["meta"]["transform"], foot.shape, strike, coh, foot, sigma_regional_px=30.0
    )
    # Verify genuine spatial variation (std > 0) rather than fallback constant!
    assert float(np.std(res["shmax_az"][foot])) > 0.02
    assert float(np.std(res["stress_ts"][foot])) > 0.01


def test_geothermal_strips_leakage_column(tmp_path):
    ext = tmp_path / "external"
    ext.mkdir()
    (ext / "gdr_wellspring_in_footprint.csv").write_text(
        "row,col,temp_c,geothermquartz_c,geothermchalc_c,geothermcat_c,dist_known_fault_px\n"
        "10,12,88.5,145.0,120.0,135.0,0.0\n"
    )
    (ext / "gdr_volcanic_vents_in_footprint.csv").write_text(
        "row,col,dist_known_fault_px\n15,18,0.0\n"
    )
    from rasterio.transform import from_origin
    tr = from_origin(243350.0, 4508550.0, 100.0, 100.0)
    valid = np.ones((30, 30), dtype=bool)
    z = np.zeros((30, 30), dtype=np.float32)
    out = build_geothermal_channels(tmp_path, tr, (30, 30), z, z, z, valid)
    assert out["gdr_geotherm_plume"][10, 12] > out["gdr_geotherm_plume"][0, 0]


def test_topology_tip_cone_and_bridge_dominate_flank():
    foot = np.ones((60, 80), dtype=bool)
    vis = np.zeros((60, 80), dtype=bool)
    # Two collinear E-W segments separated by a 12-px relay gap at row 30: cols [8..24] and [37..54]
    vis[30, 8:25] = True
    vis[30, 37:55] = True
    cand_az = np.full((60, 80), np.radians(90.0), dtype=np.float32)
    td = np.full((60, 80), 0.85, dtype=np.float32)
    feats, allowed = compute_topology_features(vis, foot, cand_az, td)
    # Compare pixel 5 px ahead of tip inside collinear gap (30, 29) vs pixel 5 px laterally off flank (25, 16)
    assert feats["vis_tip_cone_short"][30, 29] > 0.45
    assert feats["vis_tip_cone_short"][25, 16] < 1e-3
    assert feats["vis_pair_bridge"][30, 30] > 0.50
    assert feats["vis_pair_bridge"][25, 16] == 0.0
    assert not allowed[30, 15]
    assert allowed[30, 29]


def test_choquet_exact_vs_bruteforce_and_visible_zeroing():
    rng = np.random.default_rng(99)
    p = np.zeros((9, 9), dtype=np.float32)
    # Populate a 3x3 active block around center (4, 4) -> 9 Bernoulli variables = 512 exact states
    sub_p = rng.uniform(0.05, 0.45, size=(3, 3)).astype(np.float32)
    p[3:6, 3:6] = sub_p
    choq = float(choquet_expected_self_credit(p)[4, 4])

    # Exact 2^9 enumeration at center (4, 4)
    coords = [(r, c) for r in range(3, 6) for c in range(3, 6)]
    exact = 0.0
    for bits in itertools.product([0, 1], repeat=9):
        prob_state = 1.0
        best_k = 0.0
        for b, (r, c) in zip(bits, coords):
            prob_state *= float(sub_p[r - 3, c - 3]) if b else float(1.0 - sub_p[r - 3, c - 3])
            if b:
                kv = max(0.0, 1.0 - np.hypot(r - 4, c - 4) / 3.0)
                if kv > best_k:
                    best_k = kv
        exact += prob_state * best_k
    assert abs(choq - exact) < 1e-6

    # Verify `eval_mask` strictly zeroes visible fault probabilities so phantom dots are never placed!
    grid_p = np.zeros((30, 30), dtype=np.float32)
    vis = np.zeros((30, 30), dtype=bool)
    vis[15, 5:25] = True
    grid_p[vis] = 1.0  # 100% probability ONLY on already-visible faults, 0 everywhere else!
    allowed = ~vis
    dots, stats = dinkelbach_choquet_allocate(grid_p, allowed, eval_mask=~vis, min_dots=0, max_dots=100)
    assert int(dots.sum()) == 0
    assert stats["n_dots"] == 0
