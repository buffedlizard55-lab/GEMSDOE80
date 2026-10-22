"""Complete automated unit & integration test suite for Lane MACS-61 (G80)."""
from __future__ import annotations

import hashlib
import json
import os
import re
import sys
import zipfile
from pathlib import Path

import numpy as np
import rasterio
from scipy.ndimage import distance_transform_edt

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
if str(REPO / "src") not in sys.path:
    sys.path.insert(0, str(REPO / "src"))
if str(REPO / "vendor") not in sys.path:
    sys.path.insert(0, str(REPO / "vendor"))

from gems57 import submission_writer as sw
from gems57 import uniqueness as uq
from gems57.grid import Grid
from gems80.macs61_emitters import celf_submodular_ranked
from gems80.macs61_features import _hessian_ridge, steerable_radon_bank
from scripts.validate_submission import validate_geotiff

MACS61_STEM = "gems80-macs61-triexpert-submod300-31492dots-20261021T235000Z-379576aa-zeros"
MACS61_SHA = "7bc62f5837739efb6a98aebc34cb908d153840d7f9a75f306c5636ec2e0d8853"
MACS61_ZIP_SHA = "b25361699d13eb90ab907ab3ad21053bf7c12c85384fffc6f5d45fce791fcbbc"


def _data_dir() -> Path:
    for p in (REPO / "data", REPO / "data" / "bridge"):
        if (p / "sample_submission.tif").is_file():
            return p
    raise FileNotFoundError("Missing sample_submission.tif in data/ and data/bridge/")


def test_macs61_data_placement_and_official_grid():
    rep = json.loads((REPO / "reports" / "prepare_data_report.json").read_text())
    assert rep["status"] == "READY"
    assert [rep["grid"]["height"], rep["grid"]["width"]] == [3730, 3292]
    assert rep["grid"]["crs"] == "EPSG:32611"
    assert rep["grid"]["valid_footprint_pixels"] == 5_167_373
    assert rep["labels"]["known_fault_pixels"] == 60_988
    assert len(rep["official_feature_bands"]) == 19


def test_macs61_hessian_ridge_axis_exact_fix_ir80_03():
    """Verify exact N-S ridge axis recovery after IR-80-03 fix (cos2theta = -1 for N-S, +1 for E-W)."""
    y, x = np.mgrid[0:64, 0:64].astype(np.float32)
    ns_ridge = np.exp(-0.5 * ((x - 32.0) / 2.0) ** 2).astype(np.float32)
    dom_ns, dir_ns = _hessian_ridge(ns_ridge, sigma=1.5)
    assert abs(float(dom_ns[32, 32])) > 0.05
    assert float(dir_ns[0, 32, 32]) < -0.98
    assert abs(float(dir_ns[1, 32, 32])) < 0.05

    ew_ridge = np.exp(-0.5 * ((y - 32.0) / 2.0) ** 2).astype(np.float32)
    dom_ew, dir_ew = _hessian_ridge(ew_ridge, sigma=1.5)
    assert abs(float(dom_ew[32, 32])) > 0.05
    assert float(dir_ew[0, 32, 32]) > 0.98
    assert abs(float(dir_ew[1, 32, 32])) < 0.05


def test_macs61_steerable_radon_straight_fault_beats_isotropic_blob():
    """Verify multi-scale Radon straightness boost suppresses isotropic mounds relative to straight scarps."""
    field = np.zeros((160, 160), dtype=np.float32)
    field[20:140, 40] = 1.0
    yy, xx = np.mgrid[0:160, 0:160]
    blob = np.exp(-0.5 * (((yy - 80) ** 2 + (xx - 120) ** 2) / 9.0)).astype(np.float32)
    field += blob

    r9, _, _, _ = steerable_radon_bank(field, length_px=9)
    r21, _, _, _ = steerable_radon_bank(field, length_px=21)
    r35, _, _, _ = steerable_radon_bank(field, length_px=35)
    straight_score = float(r9[80, 40] * (1.0 + 0.35 * r21[80, 40] + 0.35 * r35[80, 40]))
    blob_score = float(r9[80, 120] * (1.0 + 0.35 * r21[80, 120] + 0.35 * r35[80, 120]))
    assert straight_score > 1.6 * blob_score


def test_macs61_celf_submodular_enforces_exact_min_separation():
    """Verify CELF submodular 300 m emitter returns exact binary [0, 1] dots with zero collisions."""
    H, W = 90, 90
    allowed = np.ones((H, W), dtype=bool)
    allowed[:6, :] = False
    yy, xx = np.mgrid[0:H, 0:W]
    prob = (0.05 + 0.55 * np.exp(-0.5 * ((xx - 45.0) / 2.0) ** 2)).astype(np.float32)
    rank_map, _ = celf_submodular_ranked(prob, allowed, max_picks=18, min_sep_px=2.85)
    dots = ((rank_map > 0) & (rank_map <= 18)).astype(np.float32)
    assert dots.dtype == np.float32
    assert set(np.unique(dots)).issubset({0.0, 1.0})
    assert int((dots == 1.0).sum()) == 18
    assert int(dots[~allowed].sum()) == 0
    ry, rx = np.nonzero(dots == 1.0)
    pts = np.column_stack([ry, rx]).astype(np.float64)
    d2 = ((pts[:, None, :] - pts[None, :, :]) ** 2).sum(axis=-1)
    np.fill_diagonal(d2, 999.0)
    assert float(np.sqrt(d2.min())) >= 2.85 - 1e-6


def test_macs61_shared_template_repairs_ir80_05_and_ir80_06(tmp_path):
    """Verify upstream repairs IR-80-05 (symmetric blanket uniqueness) and IR-80-06 (0644 permissions)."""
    allowed = np.ones((40, 40), dtype=bool)
    dense = np.ones((40, 40), dtype=np.float32)
    sparse = np.zeros((40, 40), dtype=np.float32)
    sparse[5, 5] = 1.0
    sparse[25, 25] = 1.0

    # Symmetric dot_overlap distinguishes unidirectional blanket coverage from mutual near-duplication
    fwd = uq.dot_overlap(sparse > 0.5, dense > 0.5)
    rev = uq.dot_overlap(dense > 0.5, sparse > 0.5)
    assert fwd == 1.0 and rev < 0.05

    tpl = _data_dir() / "sample_submission.tif"
    with rasterio.open(tpl) as src:
        foot = np.isfinite(src.read(1))
    full_sparse = np.zeros(foot.shape, dtype=np.float32)
    fy, fx = np.nonzero(foot)
    full_sparse[fy[0], fx[0]] = 1.0
    out_tif = tmp_path / "out.tif"
    sw.write_submission(
        out_tif,
        full_sparse,
        tpl,
        foot,
        note="GEMS80 permcheck",
        name="gems80-permcheck",
    )
    mode = os.stat(out_tif).st_mode & 0o777
    assert (mode & 0o444) == 0o444


def test_macs61_leakage_canary_and_monotonic_loqo_ablation():
    canary = json.loads((REPO / "reports" / "leakage_canary.json").read_text())
    assert canary["n_channels_tested"] == 61
    assert canary["leakage_detected"] is False
    assert canary["all_channels_clean"] is True
    assert canary["global_max_symmetric_auc"] < 0.90

    ab = json.loads((REPO / "reports" / "ablation_experiments.json").read_text())
    s = ab["pooled_dual_seed"]["scores"]
    assert s["ctrl1_h33_analog_loqo"]["dti"] < s["ctrl2_official19_loqo_d28"]["dti"]
    assert s["ctrl3_h85_repaired_loqo"]["dti"] < s["exp1_macs61_single_d28"]["dti"]
    assert s["exp1_macs61_single_d28"]["dti"] < s["exp2_macs61_dual_d28"]["dti"]
    assert s["exp2_macs61_dual_d28"]["dti"] < s["exp3a_macs61_dual_c35_skel"]["dti"]
    assert s["exp3a_macs61_dual_c35_skel"]["dti"] > 0.115


def test_macs61_promoted_submission_portal_ready_and_orthogonal_to_all_peers():
    rc = json.loads((REPO / "research" / "macs61-run-card.json").read_text())
    rc_docs = json.loads((REPO / "docs" / "research" / "run-card.json").read_text())
    assert rc == rc_docs
    assert rc["verdict"] == "promote"
    assert rc["budget"]["weekly_submission_slots_spent"] == 0
    assert 1 <= len(rc["submission"]["note"]) <= 140

    dl_tif = REPO / "docs" / "downloads" / f"{MACS61_STEM}.tif"
    dl_zip = REPO / "docs" / "downloads" / f"{MACS61_STEM}.zip"
    sub_tif = REPO / "submissions" / f"{MACS61_STEM}.tif"
    dsub_tif = REPO / "docs" / "submissions" / f"{MACS61_STEM}.tif"
    for p in (dl_tif, dl_zip, sub_tif, dsub_tif):
        assert p.is_file(), f"Missing {p}"

    sha_tif = hashlib.sha256(dl_tif.read_bytes()).hexdigest()
    sha_zip = hashlib.sha256(dl_zip.read_bytes()).hexdigest()
    assert sha_tif == MACS61_SHA == rc["raster_sha256"]
    assert sha_tif == hashlib.sha256(sub_tif.read_bytes()).hexdigest()
    assert sha_tif == hashlib.sha256(dsub_tif.read_bytes()).hexdigest()
    assert sha_zip == MACS61_ZIP_SHA == rc["raster_zip_sha256"]

    v = validate_geotiff(dl_tif)
    assert v["passed"] is True
    assert all(v["checks"].values())
    assert v["positive_pixels"] == 31_492
    assert v["catalogue_overlap_pixels"] == 0
    assert v["nan_count_total"] == 0

    ddir = _data_dir()
    with rasterio.open(ddir / "sample_submission.tif") as ds:
        foot = np.isfinite(ds.read(1))
    with rasterio.open(ddir / "existing_faults.tif") as ds:
        cat = (ds.read(1) > 0) & foot
    allowed = foot & (~cat)
    with rasterio.open(dl_tif) as ds:
        ours = ds.read(1)
    od = (ours > 0.5) & allowed
    dt_ours = distance_transform_edt(~od)

    for peer in sorted((REPO / "docs" / "downloads").glob("*-zeros.tif")):
        if peer == dl_tif:
            continue
        with rasterio.open(peer) as ds:
            parr = ds.read(1)
        pd = (parr > 0.5) & allowed
        fwd = float((distance_transform_edt(~pd)[od] <= 3.0).mean())
        rev = float((dt_ours[pd] <= 3.0).mean())
        assert max(fwd, rev) <= 0.70, f"Overlap violation vs {peer.name}: fwd={fwd}, rev={rev}"

    for page_name in ("index.html", "executive-summary.html", "macs61-lane.html"):
        txt = (REPO / "docs" / page_name).read_text()
        assert MACS61_SHA in txt
        for ref in re.findall(r'(?:href|src)="([^"#?]+)"', txt):
            if ref.startswith(("http://", "https://", "mailto:")):
                continue
            assert (REPO / "docs" / ref).resolve().is_file(), f"Broken ref {ref} in {page_name}"
