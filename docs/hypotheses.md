# Candidate geological hypotheses (GEMSDOE80)

Status: H80-C was run on the spatially-blocked holdout (experiment 3 of 3).
H80-A and H80-B were run and not selected. H80-D and H80-E were **not** run
(budget = 3 experiments). None copies a previous GEMSDOE submission.

Competition (official):
- Task: pixel-wise fault probability, GeoDAWN region, NW Great Basin.
- Metric: distance-weighted Tversky, α=0.2, β=0.8, triangular kernel R=300 m.
- Test labels are faults **excluded** from the public USGS / INGENIOUS catalogue.
- https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/

## Why GEMSDOE32 H33-2-B2 reported 0.2778 — and whether we can beat it

H33-2-B2 is a **metric-geometry** result, not a new detector. The GEMSDOE32
audit describes it as the 0.2708 dotted base with every dot at
d(catalogue) ≤ 2 px deleted (37,654 unit dots, 0 on the catalogue, all-finite
zeros encoding, official 3730×3292 grid).

Under DTI, TP_w is a covering: each truth pixel takes the **max** kernel-weighted
prediction within 3 px. Extra mass that does not raise that max is FP at cost
α=0.2. Sparse, unit, off-catalogue dots raise credit density. A solid surface
of ~121k pixels scored worse than ~44k dots on the group's own Monte Carlo.

The user lists 0.2778 as the GEMSDOE32 download score. GEMSDOE32's own page
labelled that file **UNSCORED** and a **projected** 0.2747. This checkout treats
0.2778 as **not ORGANIZER-CONFIRMED**. The official public leaderboard read
2026-10-21 lists #1 xiaofanhu **0.4413**, which proves higher scores exist.

**Can we beat 0.2778?** In principle yes (leader 0.4413). Doing so requires
higher credit density on the *hidden* expert faults, which are off-catalogue.
Copying H33's dots is forbidden (lane protocol) and would not be unique.
H80 targets a *different* physical locus (strain-transfer residual, not the
H19-5 ridge). Whether that locus matches the hidden set is unknown until an
organizer receipt exists. Holdout against withheld *catalogue* segments is a
proxy that **cannot** reward a genuinely new fault (same limitation GEMSDOE32
flagged as IR-32-PROXY-01).

## Ranked candidates

| Rank | ID | Layers | Signature | Why off-catalogue | Non-fault mimic | vs this repo / GEMSDOE family | Expected DTI | Cost | Data |
|------|----|--------|-----------|-------------------|-----------------|-------------------------------|--------------|------|------|
| 1 | **H80-C (selected)** | dilatation (b8), RTP (b2), tilt/tc (b6) | along-fabric \|t·∇dil\| residual vs distance-to-visible-catalogue × tc; dil>0; buf=4; NMS sep=3 | Transfer/accommodation zones along magnetic fabric, not the mapped trace | Lithologic contacts; geodetic interpolation seams | Not H19-5 ridge, not H33 B=2 prune. Max 3 px overlap vs those rasters = 0.157 | unknown on live; HOLDOUT-DTI 0.0507 | low | training_features.tif |
| 2 | H80-A (run) | dilatation, tmi_hg | residual \|∇dil\| × tmi_hg | Strain-gradient halo after removing catalogue-distance expectation | Same | Unique (overlap 0.164 vs H33) | HOLDOUT-DTI 0.0470 | low | training_features.tif |
| 3 | H80-B (run) | 2nd invariant, tmi_hg, conductivity | residual \|∇2ndinv\| × (HG, cond) | Strain-magnitude gradient + clay/conductive corroboration | Alteration without slip | Unique (overlap 0.166 vs H33) | HOLDOUT-DTI 0.0473 | low | training_features.tif |
| 4 | H80-D (not run) | earthquake density (b16) residual after catalogue-distance kernel | seismicity not explained by mapped traces | Blind seismogenic strands | Induced / volcanic clusters | Not run (budget) | — | low | training_features.tif |
| 5 | H80-E (not run) | depth-to-basement (b15) horizontal gradient residual | basement steps offset from the catalogue | Buried basin-bounding faults | Stratigraphic density contacts | GEMSDOE61 already tried basement-step (low live claims); skipped | — | low | training_features.tif |

## Validation (protocol)

1. Whole-segment holdout (~30% of catalogue components, 3 px buffer); catalogue-based distance computed from **visible** faults only; visible pixels masked in the prediction; pooled DTI α=0.2 β=0.8 R=300 m.
2. Leakage canary: per-feature AUC vs withheld labels. Max AUC observed **0.607** (strain 2nd invariant) — under the 0.90 leakage bar.
3. Lane check vs GEMSDOE32 H33-2-B2, H19-5 d2.8, and the 0.2600 anchor: Spearman |r| ≤ 0.01, 3 px overlap ≤ 0.166. Thresholds 0.90 / 0.70. **Not a duplicate.**
4. No weekly submission slot spent.

Official links:
- https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/
- https://github.com/drivendataorg/gems-prize-reference-solution
- https://www.usgs.gov/data/geodawn-airborne-magnetic-and-radiometric-surveys-northwestern-great-basin-nevada-and
- https://gbcge.org/current-projects/ingenious/
- https://epsg.io/32611
- https://en.wikipedia.org/wiki/Tversky_index
