"""Comprehensive automated verification suite for GEMSDOE80 (Passes 1, 2, and 3)."""
from __future__ import annotations

import hashlib
import json
import re
import sys
from pathlib import Path

import numpy as np
import rasterio

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))
if str(ROOT / "vendor") not in sys.path:
    sys.path.insert(0, str(ROOT / "vendor"))

from gems57 import evaluate_holdout as eh
from gems57.grid import load_grid
from scripts.restore_data import DATA_DIR
from scripts.validate_submission import validate_geotiff

EXPECTED_EVAL_SHA = "ad86640b46f3c5fc6b1aceec4c4bfc82d03e35008b37ab2a8ad6076e6810ff3a"
EXPECTED_TIF_NAME = "gems80-adpt-kkt-lagrangian-36500d-20261022T012049Z-5bd11c85-zeros.tif"


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# ==============================================================================
# PASS 1: Shared Evaluator Integrity & Complete 12-Check GeoTIFF Validation
# ==============================================================================
def test_pass1_shared_evaluator_unmodified() -> None:
    hashes = eh.implementation_hashes()
    assert eh.VERSION == "gems57-pooled-hide-v2"
    assert hashes["evaluate_holdout.py"] == EXPECTED_EVAL_SHA


def test_pass1_geotiff_strict_range_no_nan_and_grid_match() -> None:
    sub_tif = ROOT / "submissions" / EXPECTED_TIF_NAME
    dl_tif = ROOT / "docs" / "downloads" / EXPECTED_TIF_NAME
    assert sub_tif.is_file(), f"Missing submission GeoTIFF: {sub_tif}"
    assert dl_tif.is_file(), f"Missing downloadable GitHub Pages GeoTIFF: {dl_tif}"
    assert _sha256(sub_tif) == _sha256(dl_tif), "Submissions and docs/downloads GeoTIFFs must be byte-identical"

    report = validate_geotiff(dl_tif)
    assert report["passed"] is True, f"Validator failed: {report}"
    assert all(report["checks"].values()), f"Failed checks: {report['checks']}"
    assert report["band_count"] == 1
    assert report["dtype"] == "float32"
    assert report["shape"] == [3730, 3292]
    assert report["crs"] == "EPSG:32611"
    assert report["min_value"] == 0.0
    assert report["max_value"] == 1.0
    assert report["nan_count_inside_footprint"] == 0
    assert report["nan_count_total"] == 0
    assert report["positive_pixels"] == 36500
    assert report["catalogue_overlap_pixels"] == 0


# ==============================================================================
# PASS 2: Leakage Canary, Capped 3-Experiment Budget & 92-Repo Uniqueness Audit
# ==============================================================================
def test_pass2_leakage_canary_all_176_channels() -> None:
    canary = json.loads((ROOT / "evidence" / "leakage_canary.json").read_text())
    assert canary["passed"] is True
    assert canary["n_features_checked"] == 176
    assert canary["sgmc_channels_included"] == 0
    assert canary["prior_rasters_in_features"] == 0
    assert canary["max_effective_auc"] <= 0.90
    assert len(canary["flagged_features"]) == 0
    assert len(canary["all_features"]) == 176
    for feat in canary["all_features"]:
        assert 0.50 <= feat["symmetric_auc"] <= 0.90, f"Canary ceiling exceeded on {feat}"


def test_pass2_experiments_cap_and_statistical_victory() -> None:
    exp = json.loads((ROOT / "evidence" / "experiments_summary.json").read_text())
    assert exp["experiments_executed"] == 3
    assert exp["experiment_budget_cap"] == 3
    assert exp["winner_arm"] == "H80-E3-LagrangianUnique"

    s20 = exp["arms_seed20"]["H80-E3-LagrangianUnique"]["pooled_holdout_dti"]
    s30 = exp["arms_seed30"]["H80-E3-LagrangianUnique"]["pooled_holdout_dti"]
    b1_20 = exp["arms_seed20"]["B1-H33-FoldFair"]["pooled_holdout_dti"]
    b1_30 = exp["arms_seed30"]["B1-H33-FoldFair"]["pooled_holdout_dti"]

    assert s20 > 0.178 and s20 > b1_20 + 0.080
    assert s30 > 0.200 and s30 > b1_30 + 0.100
    assert exp["paired_deltas_seed20"]["B1-H33-FoldFair"]["ci95"][0] > 0.070
    assert exp["paired_deltas_seed30"]["B1-H33-FoldFair"]["ci95"][0] > 0.080


def test_pass2_live_recheck_against_all_1081_registry_rasters() -> None:
    from gems80.fast_uniqueness import RegistryUniquenessScanner

    exp = json.loads((ROOT / "evidence" / "experiments_summary.json").read_text())
    e3 = exp["arms_seed20"]["H80-E3-LagrangianUnique"]
    for kind_key in ("surface_uniqueness", "dots_uniqueness"):
        u = e3[kind_key]
        assert u["repos_scanned"] == 92
        assert u["n_registry_rasters_scanned"] == 1081
        assert u["n_discriminating_priors"] == 1001
        assert u["binding_duplicate"] is False
        assert u["discriminating_overlap_firings_count"] == 0
        assert u["max_discriminating_overlap_3px"] <= 0.70
        assert u["max_discriminating_spearman"] <= 0.90

    # Re-scan the actual physical GeoTIFF bytes live against the 1,081-raster index when cache is mounted,
    # or against all shipped GEMSDOE80 rasters in docs/downloads/ on fresh CI checkouts
    from gems80.fast_uniqueness import PORT_DIR, REG_INDEX_PATH

    data_root = DATA_DIR if (DATA_DIR / "existing_faults.tif").is_file() else (
        ROOT / "data" if (ROOT / "data" / "existing_faults.tif").is_file() else ROOT / "data" / "bridge"
    )
    grid = load_grid(data_root)
    dl_tif = ROOT / "docs" / "downloads" / EXPECTED_TIF_NAME
    with rasterio.open(dl_tif) as src:
        live_arr = src.read(1)
    if REG_INDEX_PATH.is_file() and (PORT_DIR / "unit_ranks.npy").is_file():
        scanner = RegistryUniquenessScanner(grid.footprint & (~grid.catalogue))
        live_check = scanner.check(live_arr, "dots")
        assert live_check["binding_duplicate"] is False
        assert live_check["discriminating_overlap_firings_count"] == 0
        assert live_check["max_discriminating_overlap_3px"] <= 0.70
        assert live_check["max_discriminating_spearman"] <= 0.90
    else:
        from scipy.ndimage import distance_transform_edt

        allowed = grid.footprint & (~grid.catalogue)
        cand_dots = (live_arr > 0.5) & allowed
        assert int(cand_dots.sum()) == 36500
        for other in sorted((ROOT / "docs" / "downloads").glob("*-zeros.tif")):
            if other == dl_tif:
                continue
            with rasterio.open(other) as osrc:
                oarr = osrc.read(1)
            odots = (oarr > 0.5) & allowed
            ov = float((distance_transform_edt(~odots)[cand_dots] <= 3.0).mean())
            assert ov <= 0.70, f"Overlap {ov:.4f} > 0.70 vs {other.name}"


# ==============================================================================
# PASS 3: Website Links, Run Card Sync, SHA-256 & <= 140-Char UI Note Audit
# ==============================================================================
def test_pass3_docs_site_links_assets_and_metadata_sync() -> None:
    rc = json.loads((ROOT / "evidence" / "run_card.json").read_text())
    dl_tif = ROOT / "docs" / "downloads" / EXPECTED_TIF_NAME
    actual_sha = _sha256(dl_tif)

    assert rc["verdict"] == "promote"
    assert rc["submission"]["sha256"] == actual_sha
    assert 1 <= len(rc["submission"]["submission_note"]) <= 140

    html_pages = [
        "index.html",
        "executive-summary.html",
        "gemsdoe32-diagnosis.html",
        "hypotheses-and-methods.html",
        "experiments-and-uniqueness.html",
    ]
    for page_name in html_pages:
        page_path = ROOT / "docs" / page_name
        assert page_path.is_file(), f"Missing page {page_path}"
        html = page_path.read_text()

        # Verify every relative href and src resolves to an actual file in `docs/`
        refs = re.findall(r'(?:href|src)="([^"#?]+)"', html)
        assert len(refs) > 0
        for ref in refs:
            if ref.startswith(("http://", "https://", "mailto:")):
                continue
            target = (ROOT / "docs" / ref).resolve()
            assert target.is_file(), f"Broken local reference '{ref}' inside {page_name} -> {target}"

    # Verify SHA-256 and direct download link appear at the top of both index.html and executive-summary.html
    for hero_page in ("index.html", "executive-summary.html"):
        text = (ROOT / "docs" / hero_page).read_text()
        assert actual_sha in text
        assert f"downloads/{EXPECTED_TIF_NAME}" in text

    readme = (ROOT / "README.md").read_text()
    assert actual_sha in readme
    assert EXPECTED_TIF_NAME in readme
    assert "HOLDOUT-DTI" in readme and "ORGANIZER-CONFIRMED" in readme
    assert "https://www.drivendata.org/competitions/306/competition-doe-gems/" in readme
    assert "https://github.com/drivendataorg/gems-prize-reference-solution" in readme
    assert "https://gdr.openei.org/submissions/1391" in readme
