# Candidate geological hypotheses (GEMSDOE80 — H81 update, 2026-10-21)

Status: H80-C was run (HOLDOUT-DTI 0.0507) and packaged as the current download.
No weekly submission slot spent. H81 is the next lane — 5 newly ranked
hypotheses that were **not** tried before, grounded in the geothermal structural
literature. They differ from every GEMSDOE dot-family and from H80's
strain-residual along magnetic fabric.

Competition (official):
- Task: pixel-wise fault probability, GeoDAWN region, NW Great Basin.
- Metric: distance-weighted Tversky, α=0.2, β=0.8, triangular kernel R=300 m.
- Test labels are faults **excluded** from the public USGS / INGENIOUS catalogue.
- https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/

---

## Why GEMSDOE32 H33-2-B2 reported 0.2778 — and whether we can beat it

H33-2-B2 is a **metric-geometry** result, not a new detector. The GEMSDOE32 audit
describes it as the 0.2708 dotted base with every dot at d(catalogue) ≤ 2 px
deleted (37,654 unit dots, 0 on the catalogue, all-finite zeros encoding,
official 3730×3292 grid).

Under DTI, TP_w is a covering: each truth pixel takes the **max** kernel-weighted
prediction within 3 px. Extra mass that does not raise that max is FP at cost
α=0.2. Sparse, unit, off-catalogue dots raise credit density. A solid surface
of ~121k pixels scored worse than ~44k dots on the group's own Monte Carlo.
GEMSDOE32 measured the marginal break-even credit per dot as **0.0548**;
theory gives `0.2·DTI ≈ 0.052` at 0.26 — two independent routes to the same bar.
Dots below that bar hurt the score. H33-2-B2 keeps 37,654 dots whose mean credit
on the hidden truth model was just above the bar — a flank prune that thins the
field where catalogue distance makes credit certain to be FN.

The user lists 0.2778 as the GEMSDOE32 download score. GEMSDOE32 labels that
file **UNSCORED** with `receipt: null` and a **projected 0.2747** from a
4-fold live-mirror (official DTI, not holdout). The official public
leaderboard read 2026-10-21 lists #1 **xiaofanhu 0.4413** (26 submissions)
and #48 **extradr19 0.2778** — the latter is a coincidence; the leaderboard
exposes no file hash. 0.4413 proves higher scores exist.

**Can we beat 0.2778 on live?** In principle yes (leader 0.4413). Doing so
requires higher credit density on the *hidden* expert faults, which are
off-catalogue by construction (problem description § "Labels": test faults are
not in the USGS database). Copying H33's dots is forbidden (lane protocol) and
would not be unique. H80-C (strain-transfer residual) is already unique vs
that family (Spearman ≈ 0, 3 px overlap 0.15) but its HOLDOUT-DTI is only
0.0507 on the catalogue-hide-and-recover proxy. The proxy **cannot reward a
genuinely new fault** — it withholds catalogue segments (IR-32-PROXY-01). We need
a detector whose physics predicts where blind faults *should* be, not where the
catalogue already is. The H81 candidates below do that.

**What controls live DTI?**  
`DTI = TP_w / (TP_w + 0.2 FP_w + 0.8 FN_w)` with `TP_w = Σ_g max_x p(x)k(d)` for
`k(d)=max(1-d/300,0)`. If truth mass is `|G|` and we emit `N` unit dots with
mean credit `c = TP_w/N`, then `DTI ≈ c·N / (0.2N + 0.8|G|)` when every truth pixel
is within range (the dilute limit). Marginal benefit of one dot is
`d DTI/dN ≈ (c - 0.2 DTI) / (0.2N + 0.8|G|)`. This is the **metric-geometry
budget**: the field's geometry sets `c`, the budget sets the prune.

---

## Ranked candidates (not tried before in this repo)

| Rank | ID | Layers (training_features band #) | Physical signature (transform) | Why it should catch faults **missing** from USGS/INGENIOUS, not those already in it | Named non-fault mimic | How it differs from this repo / GEMSDOE family | Expected DTI gain | Cost | External data |
|------|----|----------------------------------|--------------------------------|----------------------------------------------------------|-----------------------|--------------------------------------------------|-------------------|------|--------------|
| **1** | **H81-A tip-linkage / relay-ramp extrapolation** | labels only (visible fault endpoints) + geod_dilaterate (08) as dilational gate | Endpoint detection on skeletonized visible traces (morphological hit-or-miss); link nearby endpoints (gap 1.0–3.5 km, relay-ramp width) with a Brechtel-style straight linkage axis; score = `exp(-d_link/σ)·[dil>0]`; distance-residual vs visible-catalogue median is *not* applied — geometry itself is the residual | 32% of Great Basin geothermal systems sit in step-overs / relay ramps between overlapping normal faults, 25% at fault terminations where horse-tailing generates closely spaced splays (Faulds et al. 2021, GRC 1034804; Faulds & Hinz 2015 OSTI 1724082). By definition the breaching faults that host permeability sit **between** mapped traces, not on them. A linkage line between two nearby mapped tips predicts exactly that zone. The test set is hidden faults, so it should be enriched in relay-linkage segments. | Fluvial channel braid terminations; linear dune alignments; cadastral grids | Pure catalogue geometry — never used in H80 (which was gradient magnitude along magnetic fabric × tilt) nor in H33/H19-5 (ridge surfaces pruned by catalogue distance). No magnetic/RTP/tilt, no Sobel magnitude. Overlap vs H80-C expected < 0.4 | **+0.01–0.03 on HOLDOUT** (holdout tip recovery is directly measurable) | **low** — labels + dil only | none — all bands in training_features.tif |
| 2 | H81-B earthquake-density residual (seismogenic blind) | ieq_n100a15 (16), deq_n100a15 (10), geod_shearrate (07), geod_dilaterate (08) | Kernel density of earthquake intensity (16) residualized against distance-to-visible-catalogue (median-in-bin subtraction, max_bin=40); gated by `dil>0` and `shear > p55`; multiply by inverse distance-to-earthquake (10) normalized to [0,1] | Microseismic swarms track active blind faults that have not yet produced a mappable surface scarp (Kreemer & Young 2022; NGL MAGNET). INGENIOUS labels include subsurface evidence, but the USGS surface catalogue does not. Off-catalogue seismicity halos are the signature of unmapped active strands. | Mining-blast clusters; induced seismicity near Clay- ton Valley brine operations | H80 did not use earthquake density as primary; previous high scorers used only strain/magnetic, not seismic corroboration | +0.005–0.015 | low | none |
| 3 | H81-C conductive clay-cap / alteration channel | cond_surf (17), det_elev_slope (19), iso_grav_anom_slope (05), iso_grav_anom_hg (18) | Laplacian of conductivity (17) (altered clay cap = high conductivity halo edge); corroborate where topographic slope (19) and gravity HG (18) are high (basin architecture); gate on conductivity within 2–98 pctile and distance-residual | Hydrothermal alteration creates conductive smectite clay caps and silicification halos that extend laterally beyond the fault plane (Glen et al. 2022 GRC 1034630; USGS GeoDAWN rad overview). The surface conductivity anomaly marks the *fluid pathway*, not the mapped trace, so it is offset from the catalogue. | Pedogenic caliche / salt pan conductivity; agricultural irrigation | H80 used conductivity only as secondary multiply in H80-B; never as Laplacian edge | +0.005–0.015 | low | none (conductivity is band 17) |
| 4 | H81-D gravity-magnetic concordance (basin-bounding basement step) | iso_grav_anom_hg (18), tmi_hg (03), depth_to_base_surf (15), tmi_vg (09) | Joint gradient product `|∇Hgrav|·|∇Hmag|` (edges that are both gravity and magnetic); divide by depth-to-basement horizontal gradient magnitude; keep where joint edge > dist-to-visible-catalogue expectation | Buried basin-bounding normal faults offset the basement-cover contact; the surface trace is hidden under Quaternary alluvium. Concordant gravity + magnetic edges at a basement step are offset 0.5–2 km from the published trace (Giddens & Faulds 2025 SGW Figure 7). The test set should be enriched in such pediment faults. | Dike / sill contacts; stratigraphic pinchouts | H80 never used gravity HG jointly; GEMSDOE61 basement-step was depth alone (failed). This adds the joint-gradient filter | +0.003–0.01 | low–mid | none |
| 5 | H81-E lidar micro-scarp curvature orthogonal to extension | det_elev_slope (19) + det_elev (12) at 100 m as proxy; 1 m DEM curvature if obtainable | Profile curvature perpendicular to geodetic extension azimuth (regional strain tensor from Kreemer); `κ = ∂²z/∂n²` where n ⟂ extension; NMS on curvature maxima | USGS Quaternary mapping misses sub-meter scarps (<1 ka) resolved in 3DEP 1 m lidar (Giddens et al. 2025: lidar hillshade used to map Quaternary faults). Such scarps are exactly the hidden faults the competition experts added. | Fluvial terrace risers, beach ridges of Lake Lahontan, anthropogenic cut/fill | H80 used no topographic curvature; previous lidar-scarp entries (GEMSDOE64,67) used lidar as direct feature, not curvature-orthogonal-to-strain; those scored <0.15, but none combined with extension direction | unknown (proxy) to +0.01 | **high** if 1 m needed | USGS 3DEP 1 m DEM (free, official) — `https://apps.nationalmap.gov/lidar-explorer/` or S3 `usgs-lidar`; ~706 of 716 tiles listed in competition PDF (OCR-recovered). Check obtainable before proposing as viable — **it is obtainable** (GeoDAWN § says 3DEP lidar coordinated, similar extent) |

### Why the ranking

1. **Expected DTI improvement** dominates cost. H81-A directly targets the most
   common geothermal structural setting (step-over, Faulds 32% + terminations
   25%) with a mechanism that by construction cannot reproduce a mapped trace.
   It should recover withheld tips in the hide-and-recover split — a falsifiable
   prediction.
2. **Implementation cost** is second. H81-A needs only morphological ops on the
   label raster plus one strain gate, no external data, no training, CPU-only.
   H81-E would need 1 m DEM tiling and strain azimuth — viable but expensive
   (budget 2 h).
3. **Lane uniqueness** — all five are distinct from H80's along-fabric
   dilatation gradient and from the GEMSDOE dotted-ridge family (whose signature
   is ridge detection on RTp/tmi_hg/tilt). None repeats H33's metric prune.

### Validation (protocol, must pass before spending a slot)

1. Whole-segment holdout (~30% of catalogue components, 3 px buffer); catalogue-
   based linkage hypotheses derived **only** from visible faults; visible pixels
   masked in the prediction; pooled DTI α=0.2 β=0.8 R=300 m triangular kernel.
2. Leakage canary: per-feature AUC vs withheld labels. Max AUC > 0.90 = leakage
   until proven otherwise. H81-A's features are geometric, so leakage risk is low.
3. Lane check vs GEMSDOE32 H33-2-B2, H19-5 d2.8, and the 0.2600 anchor:
   Spearman |r| ≤ 0.90 and 3 px overlap ≤ 0.70. Failure = log duplicate and stop
   that arm (parallel-run protocol).
4. No weekly submission slot spent until the top candidate beats the current
   holdout best (H80-C 0.0507) on the **same** split.

Official links:
- https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/  (metric + format)
- https://github.com/drivendataorg/gems-prize-reference-solution
- https://www.usgs.gov/data/geodawn-airborne-magnetic-and-radiometric-surveys-northwestern-great-basin-nevada-and  (GeoDAWN, DOI 10.5066/P93LGLVQ)
- https://gbcge.org/current-projects/ingenious/  (INGENIOUS, DOE DE-EE0009254)
- https://epsg.io/32611
- https://en.wikipedia.org/wiki/Tversky_index
- Faulds et al. 2021 — Assessment of Favorable Structural Settings (OSTi 1724082 / GRC) — step-over 32%, terminations 25%, intersections 22%
- Giddens & Faulds 2025 — Analysis of Step-Overs in Normal Fault Zones (SGW) — relay ramp width 2.8 km mean, 3.3 km for high-T
- Glen et al. 2022 — Geophysical Characterization of Structure and Geology (GRC 1034630) — derivative/filtering delineate buried faults
