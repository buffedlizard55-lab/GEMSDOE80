#!/usr/bin/env python3
"""Generate publication-grade empirical diagnostic charts in `docs/assets/`."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gems80.holdout import withhold_segments
from gems80.io_grid import load_band, load_labels, load_template


def plot_ablation_and_landscape(exp_summary: dict, out_path: Path) -> None:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14.5, 5.8), dpi=160)
    fig.patch.set_facecolor("#ffffff")

    # Panel A: Pre-registered HOLDOUT-DTI progression with 95% CIs
    labels = [
        "GEMSDOE80 v1\n(H80-C prior)",
        "Random Control\n(37,654 NMS)",
        "Exp 1 / 3\n34-ch MPhys OOF\n+ NMS",
        "Exp 2 / 3\nTwo-Regime OOF\n(Tip+Bridge) + NMS",
        "Exp 3 / 3 (Selected)\nTwo-Regime OOF +\n6-Shell Choquet",
    ]
    exps = exp_summary["experiments"]
    rand_h = exp_summary["reference_random_control"]
    scores = [
        0.0507,
        rand_h["score"],
        exps[0]["holdout"]["score"],
        exps[1]["holdout"]["score"],
        exps[2]["holdout"]["score"],
    ]
    lows = [
        0.0507,
        rand_h["ci95"][0],
        exps[0]["holdout"]["ci95"][0],
        exps[1]["holdout"]["ci95"][0],
        exps[2]["holdout"]["ci95"][0],
    ]
    highs = [
        0.0507,
        rand_h["ci95"][1],
        exps[0]["holdout"]["ci95"][1],
        exps[1]["holdout"]["ci95"][1],
        exps[2]["holdout"]["ci95"][1],
    ]
    colors = ["#94a3b8", "#cbd5e1", "#0284c7", "#0f766e", "#15803d"]
    x = np.arange(len(labels))
    bars = ax1.bar(x, scores, width=0.58, color=colors, edgecolor="#0f172a", linewidth=0.9, zorder=3)
    yerr = [np.array(scores) - np.array(lows), np.array(highs) - np.array(scores)]
    ax1.errorbar(x[1:], scores[1:], yerr=[yerr[0][1:], yerr[1][1:]], fmt="none", ecolor="#0f172a", elinewidth=1.6, capsize=5, zorder=4)

    for i, (b, sc) in enumerate(zip(bars, scores)):
        top_y = highs[i] + 0.006
        tag = f"{sc:.4f}\n(+{(sc/0.0507 - 1)*100:.0f}%)" if i >= 2 else f"{sc:.4f}"
        ax1.text(
            b.get_x() + b.get_width() / 2.0,
            top_y,
            tag,
            ha="center",
            va="bottom",
            fontsize=9.2,
            fontweight="bold",
            color="#0f172a",
            bbox=dict(boxstyle="round,pad=0.16", facecolor="#ffffff", alpha=0.90, edgecolor="none"),
            zorder=6,
        )

    ax1.axhline(0.1185, color="#b45309", linestyle="--", linewidth=1.3, zorder=2, label="Best prior off-cat registry raster on same split (0.1185)")
    ax1.set_ylim(0, 0.285)
    ax1.set_xticks(x)
    ax1.set_xticklabels(labels, fontsize=8.7)
    ax1.set_ylabel("HOLDOUT-DTI (alpha=0.2, beta=0.8, R=300 m)", fontsize=10, fontweight="bold")
    ax1.set_title("A. Strict 3-Experiment Pre-Registered Holdout Progression\n(18,313 withheld positives, 957 components, 1,000-draw Voronoi CI)", fontsize=10.8, fontweight="bold")
    ax1.grid(axis="y", linestyle=":", alpha=0.55, zorder=0)
    ax1.legend(loc="upper left", fontsize=8.5, frameon=True)

    # Panel B: Verified Public ORGANIZER-CONFIRMED Leaderboard & GEMSDOE Architecture Ladder
    lb_names = [
        "GEMSDOE43\nNon-topo HGB",
        "GEMSDOE66\nLiDAR scarps",
        "19GEMSDOE\nDense multiline",
        "GEMSDOE70\nSubmod OOF",
        "GEMSDOE24\nPoisson d=2.8",
        "GEMSDOE32\nH33-2-B2 (user)",
        "LB #5\nNo Fault...",
        "LB #2\ncfld__",
        "LB #1\nxiaofanhu",
    ]
    lb_vals = [0.0424, 0.1152, 0.1922, 0.2478, 0.2600, 0.2778, 0.3539, 0.3904, 0.4413]
    lb_cols = ["#94a3b8", "#64748b", "#38bdf8", "#0284c7", "#0369a1", "#1d4ed8", "#d97706", "#ea580c", "#b91c1c"]
    y_pos = np.arange(len(lb_names))
    ax2.barh(y_pos, lb_vals, height=0.62, color=lb_cols, edgecolor="#0f172a", linewidth=0.85, zorder=3)
    for yi, v in zip(y_pos, lb_vals):
        ax2.text(v + 0.007, yi, f"{v:.4f}", va="center", fontsize=8.8, fontweight="bold", color="#0f172a")
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels(lb_names, fontsize=8.5)
    ax2.set_xlim(0, 0.505)
    ax2.set_xlabel("ORGANIZER-CONFIRMED Public DW-Tversky Score (DrivenData #306)", fontsize=10, fontweight="bold")
    ax2.set_title("B. Verified Public Leaderboard & Sibling Architecture Audit\n(Audited 2026-10-22; live GEMSDOE80 v2 score UNSCORED until upload)", fontsize=10.8, fontweight="bold")
    ax2.grid(axis="x", linestyle=":", alpha=0.55, zorder=0)

    fig.tight_layout(pad=1.6)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=160)
    plt.close(fig)


def plot_nevada_zoom(out_path: Path) -> None:
    data = ROOT / "data"
    tmpl = load_template(data / "sample_submission.tif")
    foot = tmpl["foot"]
    cat = load_labels(data / "existing_faults.tif", foot)
    split = withhold_segments(cat, foot)

    maps_h = np.load(data / "cache_h80v2" / "holdout_experiment_maps.npz")
    p_phys = maps_h["p_phys"].astype(np.float32)
    p_fused = maps_h["p_fused"].astype(np.float32)
    dots_e3 = np.unpackbits(maps_h["dots_e3"])[: foot.size].reshape(foot.shape).astype(bool)

    elev, _, _, _ = load_band(data / "training_features.tif", 12)

    # Automatically locate the richest 26 km x 26 km (260x260 px) geothermal fault relay window
    from scipy.ndimage import uniform_filter
    hit_density = uniform_filter((split["hidden"] & foot).astype(np.float32), size=220)
    vis_density = uniform_filter((split["visible"] & foot).astype(np.float32), size=220)
    combo = hit_density * vis_density
    cr, cc = np.unravel_index(int(np.argmax(combo)), foot.shape)
    half = 130
    r0, r1 = max(0, cr - half), min(foot.shape[0], cr + half)
    c0, c1 = max(0, cc - half), min(foot.shape[1], cc + half)

    tr = tmpl["meta"]["transform"]
    extent_km = [
        (tr.c + c0 * tr.a) / 1000.0,
        (tr.c + c1 * tr.a) / 1000.0,
        (tr.f + r1 * tr.e) / 1000.0,
        (tr.f + r0 * tr.e) / 1000.0,
    ]

    fig, axes = plt.subplots(2, 2, figsize=(13.2, 12.0), dpi=160)
    fig.patch.set_facecolor("#ffffff")
    elev_sub = np.nan_to_num(elev[r0:r1, c0:c1], nan=0.0)
    vis_sub = split["visible"][r0:r1, c0:c1]
    hid_sub = split["hidden"][r0:r1, c0:c1]
    phys_sub = p_phys[r0:r1, c0:c1]
    fused_sub = p_fused[r0:r1, c0:c1]
    dots_sub = dots_e3[r0:r1, c0:c1]

    # (a) Ground-truth hide-and-recover split over detrended relief
    ax = axes[0, 0]
    ax.imshow(elev_sub, cmap="terrain", extent=extent_km, alpha=0.82)
    vy, vx = np.nonzero(vis_sub)
    hy, hx = np.nonzero(hid_sub)
    ax.scatter(extent_km[0] + vx * 0.1, extent_km[3] - vy * 0.1, s=3.0, c="#0f172a", label="Visible catalogue faults")
    ax.scatter(extent_km[0] + hx * 0.1, extent_km[3] - hy * 0.1, s=4.5, c="#dc2626", label="Withheld hidden faults (target)")
    ax.set_title("A. NW Great Basin Relay Sector (26×26 km, EPSG:32611)\nBlack = Visible Catalogue, Red = Withheld Truth Segments", fontsize=10, fontweight="bold")
    ax.legend(loc="lower right", fontsize=8.2, framealpha=0.92)

    # (b) Pure 34-ch Label-Free Multi-Physics OOF Posterior (`M_phys`)
    ax = axes[0, 1]
    im1 = ax.imshow(np.log10(np.maximum(phys_sub, 1e-4)), cmap="magma", vmin=-3.5, vmax=-0.8, extent=extent_km)
    ax.scatter(extent_km[0] + hx * 0.1, extent_km[3] - hy * 0.1, s=2.2, c="#38bdf8", alpha=0.55, label="Withheld truth overlay")
    ax.set_title("B. Regime B Label-Free Multi-Physics OOF log10 P_phys\n(1m LiDAR + 3D Stress + Up-Dip Worms + Heat Flow)", fontsize=10, fontweight="bold")
    ax.legend(loc="lower right", fontsize=8.2, framealpha=0.9)
    fig.colorbar(im1, ax=ax, fraction=0.044, pad=0.03, label="log10 P_phys(x)")

    # (c) Two-Regime Calibrated Posterior (`p_fused`)
    ax = axes[1, 0]
    im2 = ax.imshow(np.log10(np.maximum(fused_sub, 1e-4)), cmap="inferno", vmin=-3.5, vmax=-0.5, extent=extent_km)
    ax.scatter(extent_km[0] + vx * 0.1, extent_km[3] - vy * 0.1, s=1.8, c="#22d3ee", alpha=0.65, label="Visible Walker Lane & range-front faults")
    ax.set_title("C. Two-Regime Fused Posterior log10 P_fused (Exp 2/3)\nMode-II Crack-Tip Cones + Fault-Arrest + Collinear Bridges", fontsize=10, fontweight="bold")
    ax.legend(loc="lower right", fontsize=8.2, framealpha=0.9)
    fig.colorbar(im2, ax=ax, fraction=0.044, pad=0.03, label="log10 P_fused(x)")

    # (d) Exact 6-Shell Choquet Dinkelbach Emitted Binary Dots vs Withheld Truth
    ax = axes[1, 1]
    ax.imshow(elev_sub, cmap="Greys", alpha=0.35, extent=extent_km)
    dy, dx = np.nonzero(dots_sub)
    ax.scatter(extent_km[0] + vx * 0.1, extent_km[3] - vy * 0.1, s=2.2, c="#475569", alpha=0.65, label="Visible faults (0 emitted)")
    ax.scatter(extent_km[0] + hx * 0.1, extent_km[3] - hy * 0.1, s=6.5, c="#ef4444", alpha=0.85, label="Withheld truth faults")
    ax.scatter(extent_km[0] + dx * 0.1, extent_km[3] - dy * 0.1, s=9.0, facecolors="none", edgecolors="#16a34a", linewidths=1.1, label="Choquet binary dots {1.0}")
    ax.set_title("D. 6-Shell Choquet Dinkelbach Binary Allocation (Exp 3)\nAdaptive Spacing Along High-Posterior Continuations", fontsize=10, fontweight="bold")
    ax.legend(loc="lower right", fontsize=8.0, framealpha=0.93)

    for a in axes.ravel():
        a.set_xlabel("UTM Zone 11N Easting (km)", fontsize=9)
        a.set_ylabel("UTM Zone 11N Northing (km)", fontsize=9)

    fig.tight_layout(pad=1.5)
    fig.savefig(out_path, dpi=160)
    plt.close(fig)


def plot_features_and_uniqueness(exp_summary: dict, live_audit: dict, out_path: Path) -> None:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14.5, 5.8), dpi=160)
    fig.patch.set_facecolor("#ffffff")

    # Left: Top 14 features from M_struct + holdout AUC canary bars
    auc_map = {r["feature"]: r["effective_auc"] for r in exp_summary["leakage_canary_all_47_features"]["all_features"]}
    top_feats = exp_summary["experiments"][1]["model_diagnostics"]["top_features_struct"][:14][::-1]
    fnames = [f[0] for f in top_feats]
    fgains = np.array([f[1] for f in top_feats], dtype=np.float64) / 1000.0
    faucs = [auc_map.get(n, 0.5) for n in fnames]

    y = np.arange(len(fnames))
    bar_cols = ["#0f766e" if n.startswith("vis_") else "#0369a1" for n in fnames]
    ax1.barh(y, fgains, height=0.64, color=bar_cols, edgecolor="#0f172a", linewidth=0.75)
    for yi, g, a in zip(y, fgains, faucs):
        ax1.text(g + 0.5, yi, f"{g:.1f}k  (AUC={a:.3f})", va="center", fontsize=8.1, color="#0f172a")
    ax1.set_yticks(y)
    ax1.set_yticklabels(fnames, fontsize=8.5)
    ax1.set_xlim(0, max(fgains) * 1.36)
    ax1.set_xlabel("LightGBM Split Gain (×10³) + Holdout Canary AUC (All <= 0.789 < 0.90)", fontsize=9.5, fontweight="bold")
    ax1.set_title("A. Top 14 Two-Regime Feature Gains & Zero-Leakage Audit\n(Teal = Fold-Visible Topology, Blue = Label-Free Multi-Physics)", fontsize=10.4, fontweight="bold")
    ax1.grid(axis="x", linestyle=":", alpha=0.55)

    # Right: 20-Raster Uniqueness Scatter on Primary Live Submission
    vs_rows = live_audit["post_placement_uniqueness"]["vs"]
    xs_corr = [abs(r["spearman"]) for r in vs_rows]
    ys_ov = [r["dot_fraction_within_3px"] for r in vs_rows]
    is_lat = [r["space_filling_lattice_probe"] for r in vs_rows]

    ax2.axvspan(0.0, 0.90, ymin=0.0, ymax=0.70 / 1.05, color="#dcfce7", alpha=0.45, label="Allowed Unique Lane Zone (|rho|<=0.90 & 3px<=0.70)")
    ax2.axvline(0.90, color="#dc2626", linestyle="--", linewidth=1.3, label="Rank-correlation ceiling (0.90)")
    ax2.axhline(0.70, color="#ea580c", linestyle="--", linewidth=1.3, label="3-px overlap ceiling (0.70)")

    for xc, yo, lat, row in zip(xs_corr, ys_ov, is_lat, vs_rows):
        if lat:
            ax2.scatter([xc], [yo], s=95, marker="X", c="#9333ea", edgecolors="#0f172a", zorder=5,
                        label="5-px covering lattice probe (covers 99.83% of ALL Nevada)")
            ax2.annotate("13gems 5-px lattice probe\n(rho=+0.0005, 99.8% map coverage)", (xc, yo), xytext=(0.08, 0.93),
                         arrowprops=dict(arrowstyle="->", color="#9333ea"), fontsize=8.2, fontweight="bold", color="#581c87")
        elif "h33-2-b2" in row["name"]:
            ax2.scatter([xc], [yo], s=95, marker="D", c="#dc2626", edgecolors="#0f172a", zorder=5, label="GEMSDOE32 H33-2-B2 (0.2778)")
            ax2.annotate(f"H33-2-B2 (0.2778)\n(|rho|={xc:.3f}, ov={yo:.3f})", (xc, yo), xytext=(0.18, yo - 0.08),
                         arrowprops=dict(arrowstyle="->", color="#dc2626"), fontsize=8.3, fontweight="bold", color="#991b1b")
        else:
            ax2.scatter([xc], [yo], s=52, c="#0284c7", alpha=0.85, edgecolors="#0f172a", zorder=4)

    ax2.set_xlim(-0.02, 1.00)
    ax2.set_ylim(0.0, 1.05)
    ax2.set_xlabel("|Spearman Rank Correlation| vs Registry Raster", fontsize=9.8, fontweight="bold")
    ax2.set_ylabel("Fraction of Emitted Dots Within 3 px (300 m) of Registry Dots", fontsize=9.5, fontweight="bold")
    ax2.set_title("B. Parallel-Run Uniqueness Gate Across All 20 Registry TIFs\n(Max Non-Lattice Overlap = 0.6516 <= 0.70, Max |Spearman| = 0.0952 <= 0.90)", fontsize=10.3, fontweight="bold")
    ax2.grid(True, linestyle=":", alpha=0.55)
    ax2.legend(loc="center right", fontsize=7.9, framealpha=0.94)

    fig.tight_layout(pad=1.6)
    fig.savefig(out_path, dpi=160)
    plt.close(fig)


def main() -> int:
    exp_summary = json.loads((ROOT / "research" / "h80_experiments_summary.json").read_text())
    live_audit = json.loads(
        (ROOT / "docs" / "downloads" / "gemsdoe80-h80g-choquet-tworegime-20261022T000000Z.audit.json").read_text()
    )
    assets = ROOT / "docs" / "assets"
    plot_ablation_and_landscape(exp_summary, assets / "ablation_ladder.png")
    plot_nevada_zoom(assets / "nevada_zoom_anatomy.png")
    plot_features_and_uniqueness(exp_summary, live_audit, assets / "feature_and_uniqueness_audit.png")
    print("[OK] Generated all 3 publication diagnostic figures in docs/assets/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
