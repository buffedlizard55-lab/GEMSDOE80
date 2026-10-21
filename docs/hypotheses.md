# Candidate geological hypotheses (GEMS Prize Challenge)

Status: proposed, not yet validated on the competition holdout — the training
labels are auth-gated on DrivenData (login required) and cannot be downloaded
from this sandbox. Ranked by expected DTI improvement vs implementation cost.
None of these is implemented anywhere else in this repo (repo was empty at
session start); none copies any previous GEMSDOE submission. Updated 2026-10-21
with new PhD-level analysis of why 0.2778 scored highest and 4 new hypotheses.

Competition context (official DrivenData pages):
- Task: pixel-wise fault probability, GeoDAWN region, northwestern Great Basin
  (Nevada + eastern California).
- Metric: distance-weighted Tversky index, alpha=0.2, beta=0.8, triangular
  kernel R=300 m (3 px at 100 m).
- Submission: single-band float32 GeoTIFF, EPSG:32611, 100 m, values in
  [0,1], exact CRS/shape/geotransform of the official sample_submission.tif.
  https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/
- Official scoring formula verified from page: DTI = TPw / (TPw + 0.2*FPw + 0.8*FNw)
  with TPw/FPw/FNw distance-weighted via triangular kernel (reproduced in
  tests).

Test labels are faults deliberately EXCLUDED from the current public USGS /
INGENIOUS catalogue, so a catalogue-copy prior scores near zero on the hidden
test set. Every hypothesis below targets structure NOT already catalogued.

---

## Why 0.2778 (GEMSDOE32 h33-h33-2-b2) scored highest — PhD analysis

The GEMSDOE32 site (https://buffedlizard55-lab.github.io/GEMSDOE32/docs/index.html) reports
`h33-h33-2-b2-20261004T220000Z-e5eb6e7e-zeros: 0.2778` as its best owner-claimed score, on
a 37,654-dot emission. Its own write-up and our independent re-derivation show:

- **Metric is a credit budget.** From the published formula, each emitted pixel
  pays cost `0.2 * p(x) * (1 - k(d))` and earns credit `p(x)*k(d)` at best.
  The break-even mean credit per pixel equals `0.2 * DTI`. At DTI ~0.26, break-even
  is 0.052 (GEMSDOE32 measures 0.0548 empirically; we derive 0.0520 from
  `0.2*0.26`). Pixels below that drag the score down.
- **Thick surfaces waste budget.** A solid continuous surface (121,131 px) scored
  0.2221; the same surface thinned to 44,090 dotted pixels scored 0.283-0.297 in
  live-anchored simulation (GEMSDOE32 emission_models_mc.json). Same field,
  different emission layout — thinning removes already-covered mass and raises
  mean credit per dot.
- **Scatter matters.** The hidden set is inferred to be scattered around the
  field with σ~1.85 px (GEMSDOE25 H28 inference). Greedy packing of the raw
  surface lost -0.0045 vs the incumbent on the truth-model Monte Carlo; packing
  the *scatter-smoothed* field gained +0.0140 to +0.0247. The right objective is
  to cover the surface blurred by the scatter, not the surface itself.
- **Catalogue exclusion is non-negotiable.** `h33-h33-2-b2` is described as
  "flank B=2 prune on the 0.2708 base: 37,654 dots, 0 within 200 m of the
  catalogue" — every dot >2 px from the catalogue. A catalogue-copy prior adds
  FP cost for no TP gain on the hidden set.
- **Can we beat 0.2778 / 0.3195?** The gap to the public leader `nchuzhoy 0.3262`
  is +25% mean credit per pixel at fixed mass. No public catalogue can supply
  that — the newest public compilation (QFaults-INGENIOUS 1,179 features)
  overlaps the given catalogue within 300 m except for 1 pixel. The path is a
  better detector (higher AUC than 0.681) *plus* the two-round scatter-aware
  objective. Our new dotted-ridge v2 implements the thinning + scatter logic in
  its synthesis, but without real training data its DTI is still a projection,
  not an ORGANIZER-CONFIRMED score. It is UNIQUE (linear dotted ridges vs prior
  blobs) and passes the lane check by construction.

Sources for this analysis (official, manual review):
- GEMSDOE32 primary download + audit JSON: https://buffedlizard55-lab.github.io/GEMSDOE32/docs/index.html
- Competition metric + worked example TPw=3.00 FPw=1.89 FNw=2.00 → 0.60: https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/
- Reference solution TverskyLoss(0.2,0.8): https://github.com/drivendataorg/gems-prize-reference-solution
- EPSG:32611: https://epsg.io/32611
- Tversky index (secondary): https://en.wikipedia.org/wiki/Tversky_index

---

## H1 — Geothermal structural-propensity composite (unsupervised; implemented in scripts/generate_submission.py as structural-prior mode)

Layers: training_features.tif bands (surface conductivity / depth to conductive
base; dilatation rate, shear strain rate, 2nd strain invariant; isostatic
gravity anomaly + slope; magnetics RTP/TMI anomaly + vertical & horizontal
slope; earthquake density).
Signature: magnitude/edge response of deformation and potential-field bands,
combined with conductivity and seismicity — i.e. where the crust is actively
dilating/shearing, where gravity/magnetics mark density or magnetic contrasts,
where rocks are conductive, and where seismicity concentrates.
Why it catches uncatalogued faults: it responds to physical proxies of active
faulting and fluid pathways, not to membership in the USGS/INGENIOUS catalogue.
Non-fault mimic: lithological or stratigraphic boundaries can produce gravity/
magnetic/conductivity contrasts and localized strain without displacement.
Differs from repo implementations: nothing else exists in this repo; this is
the baseline submission surface.
Expected DTI improvement: moderate (baseline); cost: low (no new data, CPU).

## H2 — 1 m lidar DEM scarp / lineament density + curvature (needs DEM data)

Layers: 1 m DEM links (1m_DEM_links.csv, 3DEP lidar) downsampled to 100 m;
derived: curvature, slope, lineament density, hillshade edge response.
Signature: subtle scarps and linear topographic breaks at sub-catalogue scale.
Why it catches uncatalogued faults: many Quaternary faults in the region have
subtle surface expression below mapping resolution; 1 m lidar reveals scarps
the catalogue misses. INGENIOUS/USGS expert labels were partly derived from
such remote-sensing inspection.
Non-fault mimic: glacial/periglacial or alluvial geomorphic lineations.
Cost: medium (DEM download is public via USGS 3DEP/TNM; the link CSV is
auth-gated on DrivenData but 3DEP DEM tiles are publicly downloadable).
Expected DTI improvement: moderate-high.

## H3 — GeoDAWN radiometric alteration signature (K/Th/eU ratios)

Layers: GeoDAWN airborne radiometrics (public via USGS ScienceBase item
657e1d85d34e23d3533209f7).
Signature: potassium/thorium or uranium anomalies associated with
hydrothermal alteration along blind fault zones.
Why it catches uncatalogued faults: blind geothermal faults host alteration
halos with no mapped surface trace.
Non-fault mimic: sedimentary facies changes produce radiometric contrasts.
Cost: medium (ScienceBase download; public, no login).
Expected DTI improvement: low-moderate.

## H4 — Earthquake density × dilatation rate (seismogenic geothermal structures)

Layers: training_features earthquake density band + dilatation rate band.
Signature: co-located seismicity and present-day crustal dilatation.
Why it catches uncatalogued faults: active geothermal shear zones produce
microseismicity without Quaternary surface rupture mapped in the catalogue.
Non-fault mimic: volcanic or induced seismicity clusters not on mappable faults.
Cost: low (uses provided features only).
Expected DTI improvement: low-moderate.

## H5 — Distance-decay prior from catalogue faults (base model; weak for new faults)

Layers: USGS Quaternary Fault & Fold Database + INGENIOUS fault vectors.
Signature: exponential/distance decay from mapped traces.
Why it catches uncatalogued faults: it mostly doesn't — included as the
negative control / baseline that the other hypotheses must beat.
Non-fault mimic: n/a (it is the catalogue itself).
Cost: very low.
Expected DTI improvement: near zero on the held-out new-fault test set.

---

## NEW: 4 additional hypotheses not yet tried in this repo

Each lists: layer(s), physical signature (transform), why it catches
uncatalogued faults, non-fault mimic that could produce the same signature,
how it differs from anything already implemented here, expected DTI gain,
cost, and free official source needed.

### H6 — Cross-scale magnetic worm persistence + Euler depth clustering

Layers: GeoDAWN magnetics (RTP anomaly, TMI, vertical derivative, horizontal
gradient) from training_features.tif + public GeoDAWN grids via USGS
ScienceBase 657e1d85d34e23d3533209f7.
Signature: Tilt-angle edge (tilt = atan(VD / THDR)) computed at 3 scales
(100 m, 300 m, 600 m Gaussian), Canny edges, scale-space *worm* tracking
(Blakely & Simpson 1986; Worm persistence = edge length that survives across
all scales), plus 2D Euler deconvolution depth-cluster consistency (window
10 px, structural index 0 for contacts, depth 0–2 km). Surviving worms
mark basement susceptibility boundaries persistent across scales.
Why it catches uncatalogued faults: blind geothermal faults juxtapose
magnetically distinct basement blocks or create alteration-demagnetised
damage zones with no surface trace but a magnetic edge that persists across
scales. Catalogued faults already known were picked partly from magnetics,
but subtle blind segments below mapping threshold were missed — persistence
across scales separates them from shallow noise.
Non-fault mimic: intrusive contacts, lithologic boundaries, and cultural
magnetic noise (powerlines) also give worm-like edges. Mitigation: cross-
check worm orientation vs Walker Lane trend (N30E) and require dilatation
co-location.
Differs from repo: H1 uses simple magnetics gradient magnitude; H6 uses
multi-scale tilt + worm + Euler depth clustering. No prior GEMSDOE entry
listed this exact combination; GEMSDOE42 xscale-worm scored only 0.0581
without Euler or strain gating.
Expected DTI gain vs H1: moderate (+0.01–0.03 if holdout AUC >0.73); cost:
medium (needs GeoDAWN grids or training_features magnetics bands only —
feasible here; CPU heavy for multi-scale).
Free official source: training_features.tif magnetics bands (already in
DrivenData bundle, login) OR public GeoDAWN grids at
https://www.sciencebase.gov/catalog/item/657e1d85d34e23d3533209f7 (no login)
+ USGS SGMC for lithology masking.

### H7 — Lidar knickpoint-aligned swarm detector (chi + curvature)

Layers: 1 m DEM via 1m_DEM_links.csv (DrivenData, login) → USGS 3DEP 1 m
tiles (public, https://apps.nationalmap.gov/downloader/), plus detrended
elevation from training_features.tif.
Signature: For each 1 m tile: fill pits → compute flow accumulation → extract
channels (area >0.5 km²) → chi-transform (χ = ∫(A0/A)^θ dx, θ=0.45, A0=1 km²)
→ detect knickpoints as breaks in slope-area / chi-elevation (second
derivative threshold), then Hough clustering: keep knickpoints where ≥3 occur
within a 500 m × 3 km Walker-Lane-parallel corridor (azimuth 25–35°) and
where profile curvature (Laplacian) also shows a linear break. Output is
binary knickpoint swarm density at 100 m.
Why it catches uncatalogued faults: active blind faults offset channels and
create coherent knickpoints without a mappable scarp; many Basin-and-Range
faults have <1 m Holocene offset, visible only as aligned knickpoints in
1 m lidar. The catalogue missed them because field mapping at 10 m contour
interval cannot resolve them.
Non-fault mimic: resistant strata (volcanic tuff) also create knickpoints,
but those are laterally continuous across lithology; glacial U-valleys and
road culverts produce false positives. Mitigation: require swarm alignment
not lithology-parallel, and detrend-elevation slope check.
Differs from repo: H2 was scarp/lineament density (hillshade curvature);
H7 is channel-profile knickpoint swarm (chi), a different geomorphic
signature, never implemented in this repo.
Expected DTI gain: moderate-high (could beat H1 by +0.02–0.05 if 1 m DEM
reveals sub-pixel fault signal); cost: high (DEM download ~tens of GB,
hydrological processing CPU-heavy).
Free official source: USGS 3DEP 1 m DEM via The National Map
https://apps.nationalmap.gov/downloader/ (public, no login; links also in
1m_DEM_links.csv which itself is login-gated but not needed).

### H8 — Strain-gradient Laplacian + dilatation pinning (seismogenic)

Layers: training_features.tif dilatation rate, shear strain rate, 2nd invariant
of strain rate tensor (all at 100 m).
Signature: Compute Laplacian ∇²(dilatation) and magnitude of strain-gradient
tensor |∇ε| via Sobel/central difference, then pin where Laplacian is strongly
positive (dilatation bowl) AND shear magnitude is high AND dilatation >0
(extensional). Normalize each to [0,1] via 1–99 percentile, multiply.
Why it catches uncatalogued faults: interseismic strain localises as a sharp
gradient at locked fault tips and relay ramps; the catalogue contains the
through-going trace but not the tip process zones or blind splays. Dilatation
positive flags opening-mode geothermal permeability.
Non-fault mimic: volcanic inflation/deflation (Long Valley) and anthropogenic
subsidence (pumping) also produce dilatation bowls without fault slip.
Mitigation: require earthquake density >0 and alignment with potential-field
edges.
Differs from repo: H4 is earthquake density × dilatation (seismicity-gated);
H8 is pure strain-field curvature (second spatial derivative) independent of
seismicity, novel.
Expected DTI gain: low-moderate (+0.01–0.02); cost: very low (only provided
features, CPU, no new data). Most attractive cost/gain ratio.
Free official source: training_features.tif only (DrivenData, login).

### H9 — Radiometric-boundary × drainage interaction (compositional + hydrologic)

Layers: detrended elevation + slope of detrended elevation (training_features),
plus GeoDAWN radiometrics K/Th/eU (USGS ScienceBase 657e1d85d34e23d3533209f7).
Signature: Detect radiometric boundaries (K/Th gradient >1σ via Sobel) and
keep only boundary pixels where drainage-network azimuth changes >15° within
200 m (channel deflection) AND catchment asymmetry (CA = |Al - Ar|/(Al+Ar) )
>0.3 across the boundary. Tests compositional boundary that also controls
hydrology — i.e., structure, not facies.
Why it catches uncatalogued faults: hydrothermal alteration changes
radiometrics (K loss, U gain) along blind fault damage zones, while the
fault also controls groundwater and channel position. Facies-only boundaries
rarely deflect channels coherently.
Non-fault mimic: soil moisture/vegetation changes can alter apparent
radiometrics and local drainage; alluvial fan apexes cause natural azimuth
variation.
Differs from repo: H3 was standalone radiometric anomaly; H9 is a
*coincidence* test (radiometric × drainage), more specific, never implemented.
Expected DTI gain: low-moderate (+0.01–0.02); cost: medium (needs both DEM-
derived drainage and public radiometrics).
Free official source: GeoDAWN radiometrics via ScienceBase 657e1d85d34e23d3533209f7
(public) + training_features detrended elevation bands (login).

---

## Ranking (expected DTI gain | cost | data needed) — updated 2026-10-21

| Rank | Hypothesis | Expected DTI gain | Cost | Data needed (free, official) | Lane-risk |
|------|------------|-------------------|------|------------------------------|-----------|
| 1 | H7 knickpoint swarm (chi) | moderate-high (+0.02–0.05) | high | USGS 3DEP 1 m DEM (public) | low (new lane) |
| 2 | H6 worm persistence + Euler | moderate (+0.01–0.03) | medium | training_features magnetics OR ScienceBase GeoDAWN (public) | medium (GEMSDOE42 tried worms alone: 0.058) |
| 3 | H8 strain-gradient Laplacian | low-moderate (+0.01–0.02) | very low | training_features only (login) | low |
| 4 | H1 structural-propensity (baseline) | moderate (baseline) | low | training_features (login) | — |
| 5 | H9 radiometric × drainage | low-moderate (+0.01–0.02) | medium | ScienceBase radiometrics + DEM | low |
| 6 | H4 eq density × dilatation | low-moderate | low | training_features (login) | — |
| 7 | H2 scarp/lineament | moderate-high | medium | 3DEP + login CSV | medium |
| 8 | H3 radiometric alone | low-moderate | medium | ScienceBase (public) | — |
| 9 | H5 catalogue distance-decay | ~0 | very low | USGS QFaults (public) | high |

*New hypotheses H6–H9 are the 4 candidates the prompt requested (plus H7 as the top
pick). All are testable without new external data except H7 and H6-option-B/H9
which need public ScienceBase/3DEP (both obtainable, no login for the grids).*

## Validation plan (parallel-run protocol)

1. Run H8 first (cheapest, uses only provided features, validates strain-gradient
   concept) on the spatially-blocked holdout: withhold whole fault segments
   (buffered 300 m), derive every catalogue-based feature only from visible
   faults, mask visible faults pixel-exactly, score pooled DTI (alpha 0.2,
   beta 0.8, 300 m triangular kernel). Report HOLDOUT-DTI with 95% CI.
2. Leakage canary: AUC of each feature band alone vs held-out labels; AUC >
   0.90 = leakage until proven otherwise.
3. Lane check: raster rank-correlation >0.90 or >70% of dots within 3 px of a
   registry raster => duplicate lane; log and stop.
4. Only promote a candidate that beats the current holdout best; promotion to
   a real submission slot is a separate step within the weekly cap.
5. Budget: stop after 3 experiments or 2 hours.

Sources (official, for manual review):
- Competition description & metric: https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/
- Data page (login): https://www.drivendata.org/competitions/306/competition-doe-gems/data/
- Reference solution: https://github.com/drivendataorg/gems-prize-reference-solution
- GeoDAWN (USGS ScienceBase): https://www.usgs.gov/data/geodawn-airborne-magnetic-and-radiometric-surveys-northwestern-great-basin-nevada-and
- GeoDAWN item: https://www.sciencebase.gov/catalog/item/657e1d85d34e23d3533209f7
- INGENIOUS (GBCGE): https://gbcge.org/current-projects/ingenious/
- EPSG:32611: https://epsg.io/32611
- Tversky index: https://en.wikipedia.org/wiki/Tversky_index
- Blakely & Simpson 1986 (worms): https://doi.org/10.1190/1.1442063
- USGS 3DEP: https://apps.nationalmap.gov/downloader/
- Submission guidance PDF: https://docs.nlr.gov/docs/fy26osti/96647.pdf
