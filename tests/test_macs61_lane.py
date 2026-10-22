"""Comprehensive end-to-end automated verification suite for GEMSDOE80."""

from __future__ import annotations

import hashlib
import json
import re
import stat
import sys
import zipfile
from pathlib import Path

import numpy as np
import rasterio

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "vendor"))
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))

from gems57 import submission_writer as SW  # noqa: E402
from gems57 import uniqueness as UQ  # noqa: E402
from gems80.emitters import celf_submodular_ranked  # noqa: E402
from gems80.features import _hessian_ridge, steerable_radon_bank  # noqa: E402
from validate_submission import validate_submission  # noqa: E402


def _sha256(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def test_data_placement_and_official_grid():
    rep = json.loads((REPO / "reports" / "prepare_data_report.json").read_text())
    assert rep["status"] == "READY"
    assert [rep["grid"]["height"], rep["grid"]["width"]] == [3730, 3292]
    assert rep["grid"]["crs"] == "EPSG:32611"
    assert rep["grid"]["valid_footprint_pixels"] == 5_167_373
    assert rep["labels"]["known_fault_pixels"] == 60_988
    assert len(rep["official_feature_bands"]) == 19


def test_hessian_ridge_axis_exact_fix_ir80_03():
    yy, xx = np.mgrid[0:81, 0:81].astype(np.float32)
    # East-West horizontal ridge at y = 40 => strike theta = 0 => cos(2*theta) = +1
    horiz = np.exp(-0.5 * ((yy - 40.0) / 2.0) ** 2).astype(np.float32)
    dom_h, dir_h = _hessian_ridge(horiz, sigma=0.5)
    assert abs(float(dom_h[40, 40])) > abs(float(dom_h[10, 40])) * 5.0
    assert float(dir_h[0, 40, 40]) > 0.95  # cos(2*0) == +1

    # North-South vertical ridge at x = 40 => strike theta = pi/2 => cos(2*theta) = -1
    vert = np.exp(-0.5 * ((xx - 40.0) / 2.0) ** 2).astype(np.float32)
    _, dir_v = _hessian_ridge(vert, sigma=0.5)
    assert float(dir_v[0, 40, 40]) < -0.95  # cos(pi) == -1


def test_steerable_radon_straight_fault_beats_isotropic_blob():
    yy, xx = np.mgrid[0:65, 0:65].astype(np.float32)
    line = np.exp(-0.5 * ((xx - 32.0) / 1.2) ** 2).astype(np.float32)
    blob = np.exp(-0.5 * (((xx - 32.0) ** 2 + (yy - 32.0) ** 2) / 9.0)).astype(np.float32)
    _, aniso_line, _, _ = steerable_radon_bank(line, length_px=21, n_angles=8)
    _, aniso_blob, _, _ = steerable_radon_bank(blob, length_px=21, n_angles=8)
    assert float(aniso_line[32, 32]) > float(aniso_blob[32, 32]) * 2.5


def test_celf_submodular_enforces_exact_min_separation():
    yy, xx = np.mgrid[0:50, 0:50].astype(np.float32)
    post = np.exp(-0.5 * ((xx - 25.0) / 3.0) ** 2).astype(np.float32)
    allowed = np.ones((50, 50), dtype=bool)
    rk, gains = celf_submodular_ranked(post, allowed, max_picks=12, candidate_pool=2500, min_sep_px=2.85)
    py, px = np.nonzero(rk > 0)
    assert py.size == 12
    assert np.all(gains > 0)
    # Pairwise Euclidean distance between all placed dots must be >= 2.85 px!
    for i in range(py.size):
        for j in range(i + 1, py.size):
            dist = float(np.hypot(py[i] - py[j], px[i] - px[j]))
            assert dist >= 2.85 - 1e-6


def test_shared_template_repairs_ir80_05_and_ir80_06(tmp_path):
    # IR-80-05: 516,738-px blanket with asymmetric containment (fwd=0.99, rev=0.18 => sym=0.305 <= 0.70)
    # must NOT falsely block a distinct 31,492-dot submission, while true copy (fwd=0.85, rev=0.80) MUST block!
    rows_clean = [
        {"file": "normal_dot.tif", "their_dots": 37654, "forward_overlap": 0.46, "reverse_overlap": 0.38, "spearman_full_footprint": 0.01, "halo_coverage": 0.15},
        {"file": "blanket_10pct.tif", "their_dots": 516738, "forward_overlap": 0.92, "reverse_overlap": 0.18, "spearman_full_footprint": 0.21, "halo_coverage": 0.53},
    ]
    dec_clean = UQ.binding_duplicate_decision(rows_clean)
    assert dec_clean["duplicate"] is False
    assert dec_clean["unrepaired_v2_worst_blanket_overlap"] == 0.92

    rows_true_copy = rows_clean + [
        {"file": "copied_blanket.tif", "their_dots": 150000, "forward_overlap": 0.85, "reverse_overlap": 0.80, "spearman_full_footprint": 0.82, "halo_coverage": 0.40}
    ]
    assert UQ.binding_duplicate_decision(rows_true_copy)["duplicate"] is True

    # IR-80-06: submission_writer writes world-readable 0644 files
    sample = REPO / "data" / "sample_submission.tif"
    with rasterio.open(sample) as s:
        fp = np.isfinite(s.read(1))
    arr = np.zeros(fp.shape, dtype=np.float32)
    fy, fx = np.nonzero(fp)
    arr[fy[:100], fx[:100]] = 1.0
    out_t = tmp_path / "test_perm.tif"
    SW.write_submission(out_t, arr, sample, fp, note="unit test perm", name="test-perm")
    assert stat.S_IMODE(out_t.stat().st_mode) == 0o644
    assert stat.S_IMODE(out_t.with_suffix(".zip").stat().st_mode) == 0o644


def test_leakage_canary_and_monotonic_loqo_ablation():
    canary = json.loads((REPO / "reports" / "leakage_canary.json").read_text())
    assert canary["n_channels_tested"] == 61
    assert canary["leakage_detected"] is False
    assert canary["all_channels_clean"] is True
    assert canary["global_max_symmetric_auc"] <= 0.90
    for ch in canary["channels"]:
        assert "dist_known_fault" not in ch["name"]
        assert ch["max_fold_symmetric_auc"] <= 0.90

    abl = json.loads((REPO / "reports" / "ablation_experiments.json").read_text())
    assert abl["experiments_run"] == 3
    scores = abl["pooled_dual_seed"]["scores"]
    assert scores["exp1_macs61_single_d28"]["dti"] > scores["ctrl3_h85_repaired_loqo"]["dti"]
    assert scores["exp2_macs61_dual_d28"]["dti"] > scores["exp1_macs61_single_d28"]["dti"]
    assert scores["exp3a_macs61_dual_c35_skel"]["dti"] > scores["exp2_macs61_dual_d28"]["dti"]


def test_promoted_submission_portal_ready_unique_and_links_valid():
    rc = json.loads((REPO / "research" / "run-card.json").read_text())
    rc_docs = json.loads((REPO / "docs" / "research" / "run-card.json").read_text())
    assert rc == rc_docs
    assert rc["verdict"] == "promote"
    assert 1 <= rc["submission"]["note_length_chars"] <= 140

    tif_path = REPO / "docs" / "submissions" / rc["submission"]["file"]
    zip_path = REPO / "docs" / "submissions" / rc["submission"]["zip_file"]
    assert tif_path.is_file() and zip_path.is_file()
    assert _sha256(tif_path) == rc["raster_sha256"]
    assert _sha256(zip_path) == rc["raster_zip_sha256"]

    with zipfile.ZipFile(zip_path) as zf:
        assert zf.namelist() == [tif_path.name]
        assert zf.testzip() is None

    val = validate_submission(
        tif_path,
        REPO / "data" / "sample_submission.tif",
        REPO / "data" / "existing_faults.tif",
    )
    assert val["all_passed"] is True
    assert val["download_ready"] is True
    assert val["portal_ready_without_manual_inspection"] is True
    assert val["n_checks_passed"] == val["n_checks_total"] >= 18

    # Verify missing file graceful failure path in validator too
    missing_val = validate_submission(REPO / "docs" / "submissions" / "does_not_exist.tif")
    assert missing_val["all_passed"] is False

    # Verify zero broken local links in GitHub Pages HTML files & SHA-256 consistency in README
    for html_name in ("index.html", "executive-summary.html"):
        html = (REPO / "docs" / html_name).read_text()
        assert rc["raster_sha256"][:12] in html
        assert rc["submission"]["name"] in html
        for target in re.findall(r'(?:href|src)="(\./[^"#?]+)"', html):
            resolved = (REPO / "docs" / target[2:]).resolve()
            assert resolved.is_file() and resolved.stat().st_size > 0, f"Broken relative asset in {html_name}: {target}"

    readme = (REPO / "README.md").read_text()
    assert rc["raster_sha256"] in readme
    assert rc["submission"]["name"] in readme
    assert "<prompt>" in readme and "</prompt>" in readme
