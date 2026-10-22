# GEMSDOE80 — unique GEMS Prize GeoTIFF

**Re-read this README and `docs/SESSION_PROMPT.md` at the start of every session.**

> **OK TO DOWNLOAD AND SUBMIT (format).**
> Primary file: [`docs/downloads/gemsdoe80-h81a-relay-20261022T010259Z-zeros.tif`](docs/downloads/gemsdoe80-h81a-relay-20261022T010259Z-zeros.tif)
>
> Single-band float32, EPSG:32611, 3730×3292, transform (100, 0, 243350, 0, −100, 4508550),
> **all 12,279,160 cells finite and in [0, 1]**, no nodata tag, 0 NaN, 7415 off-catalogue unit dots.
> This encoding avoids DrivenData's `Predicted values must be in range [0, 1]`.
>
> **Not organizer-scored.** HOLDOUT-DTI = **0.0607** (evaluator `GEMSDOE80-evaluate_holdout.py-v1`,
> 18,313 withheld positives; 95% CI [0.0541, 0.0636] bootstrap). Verdict **promote** as a
> packaged candidate. No weekly slot spent. Unique vs H33/H19-5 (Spearman 0.02, 3 px overlap 0.22).

GitHub Pages: executive summary and the download button are at the top of
[`docs/index.html`](docs/index.html) · how-to-submit: [`docs/executive-summary.html`](docs/executive-summary.html)

**Unique name:** `gemsdoe80-h81a-relay-20261022T010259Z-zeros`

**Note (85/140):** `H81A relay tip-link d8-35 NMS 7415; off-cat buf4; HOLDOUT-DTI only; unique vs H33/H80`

SHA-256 zeros: `701ed7fed8611861f2fe2ec6437c850aede0a6a792b5fc3dbbca10191953bbbd`
SHA-256 zip: `5fae213ad545cc41c2ec2c1199496ad9f54f0151a963b48989976c1fbf1c4420`

Previous packaged candidate (archived, not primary): `gemsdoe80-h80c-strain-accom-resid-20261021T233457Z-zeros` (HOLDOUT 0.0507, 28,807 dots).

## 0. What this is

DrivenData competition [#306 GEMS Prize](https://www.drivendata.org/competitions/306/competition-doe-gems/):
pixel-wise probability of geothermal-indicative faults in the GeoDAWN region
(NW Great Basin, Walker Lane). Metric: distance-weighted Tversky index, α=0.2, β=0.8,
triangular kernel R=300 m
([problem description](https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/)).
Test labels are faults **not** in the public USGS / INGENIOUS catalogue.
Public leaderboard 2026-10-21: **#1 xiaofanhu 0.4413** —
higher scores exist and the goal is to beat the 0.2778 lane.

## 1. Why H33-2-B2 was cited at 0.2778, and whether we can beat it

H33-2-B2 (GEMSDOE32) is a **metric-geometry** file: take a 0.2708 dotted ridge,
delete every dot within 2 px of the catalogue, keep 37,654 unit pixels, write
all-finite zeros on the official 3730×3292 grid. DTI is a covering budget —
`sparse off-catalogue mass beats a thick surface`. The break-even credit per
dot is `0.2·DTI` (≈0.052 at the leader's level; GEMSDOE32 measured 0.0548).
Dots below that bar hurt the score. H33-2-B2 kept only dots above the bar by
flank-pruning — it optimized **geometry**, not detection.

The 0.2778 figure is **not ORGANIZER-CONFIRMED in this checkout**
(GEMSDOE32 labelled that file UNSCORED / projected 0.2747). Official public
leaderboard 2026-10-21: **#1 xiaofanhu 0.4413** (26 submissions). Higher scores
exist. Beating 0.2778 on the live set requires a *different detector* that
still hits hidden faults with higher credit density. Copying H33 is forbidden
(lane protocol). H81-A is unique vs that family (Spearman ≈ 0.02, 3 px overlap
0.22) and HOLDOUT-DTI 0.0607 > H80-C 0.0507. Live DTI is unknown until a receipt
exists.

### Scientific why — geothermal structural settings

Faulds et al. (2021, OSTI 1724082) inventoried 250 Great Basin geothermal
systems: **step-overs / relay ramps in normal fault zones = 32%**, fault
terminations = 25%, fault intersections = 22%. Giddens & Faulds (2025, SGW)
measured 450+ step-overs: mean relay-ramp width 2.8 km (3.3 km for high-T
producing systems), overlapping hard-linked geometries most productive. These
bonanza zones sit **between** mapped traces — the breaching faults are
off-catalogue by definition. The test set (hidden faults not in USGS) is
enriched in exactly those relay-linkage segments. A detector that extrapolates
between fault tips should have higher credit density than a generic ridge.

## 2. Lane H81-A (experiment 3 of 3, new best of 6 total)

**Hypothesis.** Unmapped geothermal-indicative faults occupy relay ramps /
step-overs between overlapping normal-fault tips — the most common geothermal
structural setting (32% step-overs + 25% terminations). A geometric linkage
between nearby visible tips (gap 0.8–3.5 km, mean 2.8 km) predicts the
breaching fault that hosts permeability, by construction off the catalogue and
in a dilational quadrant.

**Mechanism.** Endpoint detection (fault pixel with exactly one 8-neighbor) on
the visible skeleton; cKDTree linkage of endpoint pairs gap 8–35 px (0.8–3.5 km,
Giddens mean 2.8 km); create midpoints, dilate 2 px, distance-to-midpoint
`exp(-d/10)`; gate on positive dilatation (band 8, volumetric strain); exclude
4 px catalogue collar; NMS local maxima, sep=3 px, cap 34k.

**Named non-fault mimic.** Fluvial braid terminations, dune alignments,
cadastral grids can mimic tip-like terminations; geodetic interpolation seams
can mimic dilatation.

**Leakage canary.** Max per-feature AUC 0.736 (linkage_score) — under the 0.90
leakage bar. Dilatation alone 0.568. No leakage flag.

**Registry.** vs H80-C and H33 family: not a duplicate (Spearman ≤ 0.02,
3 px overlap ≤ 0.23 < 0.70 threshold).

**Comparison on same hide-and-recover split:**
- H81-A: HOLDOUT 0.0607 (5903 dots holdout, 7415 submit), beats random 0.0123
- H81-C: HOLDOUT 0.0563 (27275 dots) — second best
- H81-B: HOLDOUT 0.0233 — not selected
- Previous H80-C: HOLDOUT 0.0507 (28807 dots) — old best

Full card: [`research/run-card.json`](research/run-card.json) · hypotheses:
[`docs/hypotheses.md`](docs/hypotheses.md) · 5 ranked candidates before this run

## 3. Reproduce / validate

```bash
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
# competition rasters in data/ (gitignored): sample_submission.tif, labels.tif, training_features.tif
python scripts/prepare_data.py
python scripts/run_h81_experiments.py          # 3-experiment budget, H81-A/B/C
python scripts/validate_submission.py docs/downloads/gemsdoe80-h81a-relay-20261022T010259Z-zeros.tif --template data/sample_submission.tif
python scripts/evaluate_holdout.py docs/downloads/gemsdoe80-h81a-relay-20261022T010259Z-zeros.tif
```

Holdout = hide-and-recover: withhold whole fault segments + buffer, derive
catalogue-based features from visible faults only, mask visible pixels, score
pooled DTI (triangular kernel, 300 m).

## 4. Repository layout

```
docs/index.html              GitHub Pages — download + go/no-go first (H81 primary)
docs/executive-summary.html  how to submit, field-by-field
docs/downloads/              H81 zeros + zip + nan twin + audit JSON (H80 archived)
docs/SESSION_PROMPT.md       session charter (re-read every time)
scripts/evaluate_holdout.py  shared holdout
scripts/submission_writer.py shared writer (zeros + nan)
scripts/run_h81_experiments.py  current lane (H81-A/B/C)
scripts/run_h80_experiments.py  previous lane (H80-A/B/C)
scripts/validate_submission.py  format validator (CRS/shape/transform/values)
src/gems80/                  metric, I/O, emission, uniqueness
research/run-card.json       machine-readable card (this run)
research/evidence/           manifests + winner audit + geothermal literature
```

Large rasters stay in `data/` (gitignored). They were restored from the
owner-mirror pins in GEMSDOE24 / GEMSDOE (sha256 verified against GEMSDOE32's
manifest). Not organizer-authenticated.

## 5. Irregularities flagged

- No DrivenData login in this environment; rasters came from GitHub owner mirrors (pins verified, authenticity unverified — see `docs/sources.md` B).
- `sample_submission.tif` mirror contains 60,988 ones on catalogue traces (not
  "total fault absence"); used only as the grid/footprint template.
- Component-resample bootstrap CI was biased (unique() collisions) — now bootstrap with fixed seed, CI [0.0541, 0.0636] for H81.
- Shipped TIF excludes the *full* catalogue, so it scores ~0 against withheld
  catalogue pixels by construction. HOLDOUT-DTI 0.0607 used a visible-only mask.
- User-cited live highs 0.2778 / 0.3195 / 0.4413 are not mutually dated; 0.4413
  was read from the official leaderboard on 2026-10-21 (see `docs/sources.md` C).

## 6. Official links

- https://www.drivendata.org/competitions/306/competition-doe-gems/
- https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/  (metric + format)
- https://www.drivendata.org/competitions/306/competition-doe-gems/leaderboard/
- https://github.com/drivendataorg/gems-prize-reference-solution
- https://www.usgs.gov/data/geodawn-airborne-magnetic-and-radiometric-surveys-northwestern-great-basin-nevada-and  (GeoDAWN DOI 10.5066/P93LGLVQ)
- https://gbcge.org/current-projects/ingenious/  (INGENIOUS DOE DE-EE0009254)
- https://epsg.io/32611
- https://en.wikipedia.org/wiki/Tversky_index
- https://www.osti.gov/servlets/purl/1724082  (Faulds favorable settings 32% step-over)
- https://pangea.stanford.edu/ERE/pdf/IGAstandard/SGW/2025/Giddens.pdf  (Giddens step-over analysis, 2.8 km mean)

Core values: **Maximize P(Win)** · **Own the Outcome**. Negative results are deliverables.

---

## APPENDIX — session prompt (charter)

> Generate a unique TIF for DrivenData #306. Do not copy previous GEMSDOE
> submissions. Make it obvious whether the TIF is OK to download and submit.
> Validate in-project (no human inspector). Follow the parallel-run protocol
> (lane uniqueness 0.90 / 70% @ 3 px; reuse holdout + writer; label HOLDOUT-DTI
> or ORGANIZER-CONFIRMED; leakage canary AUC 0.90; run card; budget 3 experiments
> or 2 hours). Study GEMSDOE32 H33-2-B2 (user-cited 0.2778) and aim higher
> without copying it. Official leaderboard high 0.4413 (2026-10-21). Build GitHub
> Pages with executive summary first and a one-click TIF. Open a PR and merge.
> No hallucinations. Flag irregularities. Full operational text:
> [`docs/SESSION_PROMPT.md`](docs/SESSION_PROMPT.md).

Every session re-reads this README and the prompt before changing code.

---

## APPENDIX — deep research store (geothermal vent science)

This section is the start of the requested "full list that follow our
requirements" — a living evidence store for everyday use, so we do not have to
manually re-check. It will be expanded in future sessions; it is auditable with
official links.

**Geothermal fault permeability literature (Great Basin, Walker Lane):**

| Source | Claim (what it does NOT prove) |
|--------|-------------------------------|
| Faulds & Hinz 2015 (OSTI 1724082) — Favorable Tectonic and Structural Settings, GRC | 32% step-overs / relay ramps, 25% terminations, 22% intersections among 250 Great Basin systems. Host more than main fault segments. |
| Faulds et al. 2021 — Assessment of Favorable Structural Settings (GRC) | Same inventory, adds that producing systems in step-overs prefer overlapping hard-linked (73%) with relay width 3.4 km. |
| Giddens & Faulds 2025 (SGW 2025) — Analysis of Step-Overs in Normal Fault Zones | 450+ step-overs mapped; mean relay width 2.8 km, high-T 3.3 km; Walker Lane subset: right-stepping 74%, hard-linked 59%, width 2.0 km. LiDAR + 10 m DEM used. |
| Glen & Earney 2024 (USGS) — GeoDAWN DOI 10.5066/P93LGLVQ | Magnetic flight lines 200 m (Area 1) / 400 m (Area 2), drape 100–150 m; geotiffs public. Area 2 coarser than 300 m kernel — edge detection limited there. |
| Glen et al. 2022 (GRC 1034630) — Geophysical Characterization, Granite Springs Valley | Derivative/filtering (tilt, HG) delineate buried faults; joint gravity-magnetic modeling reveals hidden geothermal system. |
| Kreemer & Young 2022 — Crustal strain rates (SRL) | Geodetic strain correlates with geothermal loci; dilatation/shear bands in training_features are from this class of model. |
| Siler 2022 (USGS DOI 10.5066/P9YL58W6) — Slip/dilation tendency for Quaternary faults | Faults optimally oriented to dilation in current stress field maintain permeability. |
| DOE GDR 1391 (INGENIOUS compilation, DOI 10.15121/1881483) | Regional dataset: well/spring T & chemistry, Q-fault traces in footprint (CC BY 4.0) — free official source, not yet used in H81 but available for H82. |

**External data that could be used (free, official, shareable):**
- USGS 3DEP 1 m DEM — https://apps.nationalmap.gov/lidar-explorer/ / S3 `usgs-lidar` — 706/716 tiles for the prize footprint (OCR from PDF). Free, no license restriction. Needed for H81-E lidar curvature.
- USGS GeoDAWN grids — https://doi.org/10.5066/P93LGLVQ — CC0, high-res mag & radiometric (K/Th/U).
- GDR 1391 — https://gdr.openei.org/submissions/1391 — well/spring chemistry + fault traces — CC BY 4.0.
- USGS SGMC geologic map faults — free, used for SGMC-off diagnostic.

**Current feed idea (not yet built):** A nightly fetch of the DrivenData leaderboard page (public HTML) + USGS 3DEP availability JSON to keep the proxy frame
up to date without manual visits. Blocked today by egress to drivendata.org via
fetch tool only; shell curl fails. Future work: add a GitHub Action that runs
on a runner with full egress and commits the snapshot.

---

## APPENDIX — limitations and what is still needed

- **Data placement**: `data/` is gitignored. Run `bash scripts/download_competition_data.sh` on an unrestricted machine, then `python scripts/prepare_data.py`. This env restored via `scripts/restore_mirrors.py --core` (pins verified, not organizer-authenticated).
- **1 m DEM**: not yet downloaded (~400 GB for full footprint). `docs/irregularities.md` IR-01 tracks licensing.
- **Organizer receipt**: no live submission has been made; no ORGANIZER-CONFIRMED score exists.
- **H81-E**: requires 1 m DEM curvature orthogonal to extension — check 3DEP obtainability before marking viable (it is obtainable — see above).
- **Training**: GPU not needed for H81 (CPU morphology + distance transform); reference UNet would need GPU but is not on the critical path for the dotted lane.
