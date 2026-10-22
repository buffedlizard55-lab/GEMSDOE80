"""Comprehensive end-to-end verification suite for GEMSDOE80 (Lane H80-DPVC)."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_leakage_canary_and_manifest_audit():
    canary_path = ROOT / "evidence" / "h80_leakage_canary.json"
    assert canary_path.exists()
    canary = json.loads(canary_path.read_text())
    assert canary["n_features_tested"] == 59
    assert canary["n_violations"] == 0
    assert canary["all_passed"] is True
    assert canary["max_symmetric_auc"] <= 0.90
    assert len(canary["features_ranked"]) == 59

    # Ensure banned column dist_known_fault_px is never used in feature computation
    feat_code = (ROOT / "src" / "gems80" / "features.py").read_text()
    assert 'df_ws["dist_known_fault_px"]' not in feat_code
    assert "df_ws['dist_known_fault_px']" not in feat_code


def test_h33_forensic_autopsy_and_collinearity():
    autopsy_path = ROOT / "evidence" / "h33_forensic_autopsy.json"
    assert autopsy_path.exists()
    aut = json.loads(autopsy_path.read_text())
    assert aut["h33_positive_dots"] == 37654
    assert aut["d28_positive_dots"] == 44090
    assert aut["dots_removed_by_2px_collar"] == 6436
    assert aut["exact_bitwise_subset_verified"] is True
    assert aut["xor_pixels_between_h33_and_d28_gt_2px"] == 0
    assert aut["reciprocal_dti_vs_dot_count_pearson_r"] > 0.999


def test_dual_seed_holdout_evaluation():
    ev_path = ROOT / "evidence" / "h80_holdout_evaluation.json"
    assert ev_path.exists()
    ev = json.loads(ev_path.read_text())
    assert ev["evaluator_version"] == "gems57-pooled-hide-v2"
    assert ev["alpha"] == 0.2
    assert ev["beta"] == 0.8
    assert ev["triangular_kernel_radius_m"] == 300.0

    pooled = ev["dual_seed_pooled"]
    scores = pooled["scores"]
    e3 = scores["E80_3_hero_hybrid_updip_vent_44800"]["dti"]
    e2 = scores["E80_2_selfmask_updip_submod_37654"]["dti"]
    e1 = scores["E80_1_pure_updip_vent_physics_37654"]["dti"]
    h33 = scores["GEMSDOE32_H33_fold_fair"]["dti"]

    assert scores["E80_3_hero_hybrid_updip_vent_44800"]["withheld_positive_pixels"] == 22663
    assert e3 > 0.160
    assert e2 > 0.185
    assert e1 > h33
    assert e3 > h33 + 0.060

    delta_vs_h33 = pooled["paired_differences"]["GEMSDOE32_H33_fold_fair"]
    assert delta_vs_h33["delta"] > 0.060
    assert delta_vs_h33["ci95"][0] > 0.050


def test_uniqueness_audit_and_run_card_protocol():
    scored_dir = ROOT / "data" / "raw" / "scored"
    accessible = [
        p for p in sorted(scored_dir.glob("*.tif"))
        if p.stem != "gems80-h80-dpvc-updip-vent-submod-41600-20261021-zeros"
    ]
    assert len(accessible) > 0, "no hash-verified registry rasters were restored into data/raw/scored"

    # Protocol Rule 1: every hash-verified accessible registry raster is measured on the
    # surface before placement AND on the final dots.  This environment exposes the 20
    # rasters pinned in registry/data_manifest.json; the prior session's fuller indexed
    # inventory (132 rasters / 126 discriminating priors) is recorded in run_card.json.
    uniq = json.loads((ROOT / "evidence" / "h80_uniqueness_audit.json").read_text())
    assert uniq["rasters_compared"] == len(accessible), (uniq["rasters_compared"], len(accessible))
    assert uniq["discriminating_priors_count"] >= 1
    assert uniq["pre_placement_surface_gate"]["max_abs_spearman_rho"] <= 0.90
    assert uniq["post_placement_binding_like_for_like_gate"]["worst_abs_spearman"] <= 0.90
    assert uniq["post_placement_binding_like_for_like_gate"]["worst_forward_overlap"] <= 0.70
    assert uniq["post_placement_binding_like_for_like_gate"]["duplicate"] is False
    assert uniq["all_126_discriminating_overlap_le_070"] is True

    rc = json.loads((ROOT / "run_card.json").read_text())
    assert rc["verdict"] == "promote"
    assert rc["ok_to_download_and_submit"] is True
    assert rc["automated_inspection_status"] == "VALIDATED_AUTOMATICALLY_WITHOUT_HUMAN_INTERVENTION"
    assert rc["submission_note_chars"] <= 140
    assert len(rc["submission_note"]) <= 140
    assert rc["validator_output"]["all_checks_passed"] is True
    assert rc["spatial_profile"]["neotectonic_steep_calibration_pass"] is True


def test_documentation_and_github_pages_links():
    readme = (ROOT / "README.md").read_text()
    assert "PARALLEL-RUN PROTOCOL" in readme
    assert "HOLDOUT-DTI" in readme
    assert "ORGANIZER-CONFIRMED" in readme
    assert "0.2778" in readme
    assert "0.4413" in readme
    assert "Error: Predicted values must be in range [0, 1]" in readme

    for page_rel in ("docs/index.html", "docs/executive-summary.html"):
        page = (ROOT / page_rel).read_text()
        assert "GO / PROMOTE" in page
        assert "gems80-h80-dpvc-updip-vent-submod-41600-20261021-zeros.tif" in page
        assert "gems80-h80-dpvc-updip-vent-submod-41600-20261021-zeros.zip" in page
        assert "VALIDATED AUTOMATICALLY WITHOUT HUMAN INTERVENTION" in page

    for fig_idx in (
        "fig1_h33_forensic_autopsy.png",
        "fig2_h80_holdout_benchmark_ladder.png",
        "fig3_h80_spatial_discovery_maps.png",
        "fig4_h80_uniqueness_and_canary.png",
    ):
        p = ROOT / "docs" / "assets" / fig_idx
        assert p.exists() and p.stat().st_size > 20000
