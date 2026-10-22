# GEMSDOE80 — GEMS Prize research and candidate GeoTIFF

**Re-read this README and `docs/SESSION_PROMPT.md` at the start of every session.**

> **DOWNLOAD: YES. SUBMISSION READINESS: NOT CERTIFIED.**
> Historical candidate file: [`docs/downloads/gemsdoe80-h80c-strain-accom-resid-20261021T233457Z-zeros.tif`](docs/downloads/gemsdoe80-h80c-strain-accom-resid-20261021T233457Z-zeros.tif)
>
> This checkout independently verified the TIFF bytes: single-band float32, EPSG:32611,
> 3730×3292, transform (100, 0, 243350, 0, −100, 4508550), all 12,279,160 cells finite
> and in [0, 1], no nodata tag, and 28,807 positive pixels. The historical audit says its
> grid matched an **owner-mirrored** sample template; the official DrivenData template is
> absent here, so CRS/shape/transform compatibility is **not re-certified against organizer data**.
> A valid value range does not guarantee acceptance or a score.
>
> **Not organizer-scored; do not promote to a weekly slot.** Historical HOLDOUT-DTI = 0.0507
> (evaluator `GEMSDOE80-evaluate_holdout.py-v1`, 18,313 withheld positives). The [0.0142, 0.0383]
> four-quadrant spatial-block range is descriptive, **not a 95% CI**. Calibrated 95% CI unavailable.
> No organizer receipt exists. The H80 experiment budget was already 3/3; this review ran no new experiment.

GitHub Pages: executive summary and the download button are at the top of
[`docs/index.html`](docs/index.html) · how-to-submit: [`docs/executive-summary.html`](docs/executive-summary.html)
Three-pass review and evidence limits: [`reports/review-20261022.md`](reports/review-20261022.md)

**Unique name:** `gemsdoe80-h80c-strain-accom-resid-20261021T233457Z-zeros`

**Note (104/140):** `H80C cat-residual strain-accom NMS dots 28807; off-catalogue buf4; HOLDOUT-DTI only; unique vs H33/H19-5`

SHA-256 zeros: `d9c4ce410af5c413195aa3a463c2e4c6c223e26c6c8cf3a25a8a95c4b7a2bcc5`

## 0. What this is

DrivenData competition [#306 GEMS Prize](https://www.drivendata.org/competitions/306/competition-doe-gems/):
pixel-wise probability of geothermal-indicative faults in the GeoDAWN region
(NW Great Basin). Metric: distance-weighted Tversky index, α=0.2, β=0.8,
triangular kernel R=300 m
([problem description](https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/)).
Test labels are faults **not** in the public USGS / INGENIOUS catalogue.

## 1. Why H33-2-B2 was cited at 0.2778, and whether we can beat it (unconfirmed report)

H33-2-B2 (GEMSDOE32) is a **metric-geometry** file: take a 0.2708 dotted ridge,
delete every dot within 2 px of the catalogue, keep 37,654 unit pixels, write
all-finite zeros on the 3730×3292 grid reported in owner-maintained mirrors; that geometry is
not organizer-verified in this checkout. DTI is a covering metric: sparse off-catalogue mass
can reduce false-positive cost while maintaining proximity credit, but that mechanism alone
does not prove why any particular raster scored well. The 0.2778 figure is user-reported and
**not ORGANIZER-CONFIRMED** here (GEMSDOE32 labelled the file UNSCORED / projected 0.2747).
A public leaderboard entry read on 2026-10-21 showed 0.4413 at the top; this is a public
leaderboard observation, not a submission-page receipt. Higher public entries existed at that
read time. H80-C's historical comparisons were against three named local reference rasters only;
that is not proof of uniqueness across all user-listed GEMSDOE sites. Live DTI remains unknown.

## 2. Lane H80-C (experiment 3 of 3)

**Hypothesis.** Unmapped geothermal-indicative faults occupy strain-transfer /
accommodation zones — along-fabric dilatation gradients not explained by
distance to the visible catalogue — corroborated by magnetic tilt/curvature.

**Mechanism.** Directional derivative of dilatation (training_features band 8)
along RTP-fabric strike (band 2); subtract binned median vs distance to the
*visible* catalogue; × tilt/tc (band 6); keep positive dilatation; exclude a
4 px catalogue collar; NMS local maxima, sep=3 px.

**Named non-fault mimic.** Lithologic/intrusive contacts and geodetic
interpolation seams.

**Leakage canary.** Max per-feature AUC 0.607 (threshold 0.90) — no leakage flag.

**Registry.** vs H33-2-B2 / H19-5 d2.8 / 0.2600 anchor: not a duplicate.

Full card: [`research/run-card.json`](research/run-card.json) · hypotheses: [`docs/hypotheses.md`](docs/hypotheses.md)

## 3. Reproduce / validate

```bash
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
# Required competition rasters in data/ (gitignored): sample_submission.tif, labels.tif, training_features.tif
python scripts/prepare_data.py
# Holdout rebuilds H80-C from visible faults; it does not score a static full-catalogue TIF.
python scripts/evaluate_holdout.py --template data/sample_submission.tif \\
  --labels data/labels.tif --features data/training_features.tif
python scripts/validate_submission.py CANDIDATE.tif --template data/sample_submission.tif
```

The three H80 experiments are exhausted for this lane (3/3). Do **not** rerun
`run_h80_experiments.py` as part of this review. Any future experiment needs a
newly authorized budget and must beat a current HOLDOUT-DTI best with a valid
CI before a separate selector considers a weekly submission slot.

Holdout = hide-and-recover: withhold whole fault segments + buffer, derive
catalogue-based features from visible faults only, mask visible pixels, score
pooled DTI. The validator only checks file format; it does not make a score.

## 4. Repository layout

```
docs/index.html              GitHub Pages — download + go/no-go first
docs/executive-summary.html  how to submit, field-by-field
docs/downloads/              zeros TIF + zip + nan twin + audit JSON
docs/SESSION_PROMPT.md       session charter (re-read every time)
scripts/evaluate_holdout.py  shared holdout
scripts/submission_writer.py shared writer (zeros + nan)
scripts/run_h80_experiments.py
scripts/validate_submission.py
src/gems80/                  metric, I/O, emission, uniqueness
research/run-card.json
```

Large rasters stay in `data/` (gitignored). They were restored from the
owner-mirror pins in GEMSDOE24 / GEMSDOE (sha256 verified against GEMSDOE32’s
manifest). Not organizer-authenticated.

## 5. Irregularities flagged

- No DrivenData login in this environment; rasters came from GitHub owner mirrors.
- `sample_submission.tif` mirror contains 60,988 ones on catalogue traces (not
  “total fault absence”); used only as the grid/footprint template.
- Component-resample bootstrap CI was biased (unique() collisions) and was
  replaced by a spatial-block min/max of 4 quadrants. Pooled HOLDOUT-DTI is the score.
- Shipped TIF excludes the *full* catalogue, so it scores ~0 against withheld
  catalogue pixels by construction. HOLDOUT-DTI 0.0507 used a visible-only mask.
- User-cited live highs 0.2778 / 0.3195 / 0.4413 are not mutually dated; 0.4413
  was read from the official leaderboard on 2026-10-21.

## 6. Official links

- https://www.drivendata.org/competitions/306/competition-doe-gems/
- https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/
- https://www.drivendata.org/competitions/306/competition-doe-gems/leaderboard/
- https://github.com/drivendataorg/gems-prize-reference-solution
- https://www.usgs.gov/data/geodawn-airborne-magnetic-and-radiometric-surveys-northwestern-great-basin-nevada-and
- https://gbcge.org/current-projects/ingenious/
- https://epsg.io/32611
- https://en.wikipedia.org/wiki/Tversky_index

Core values: **Maximize P(Win)** · **Own the Outcome**. Negative results are deliverables.

---

## APPENDIX — session prompt (charter)

Generate a unique TIF for DrivenData #306. Do not copy previous GEMSDOE
submissions. Make it obvious whether the TIF is OK to download and submit.
Validate in-project (no human inspector). Follow the parallel-run protocol
(lane uniqueness 0.90 / 70% @ 3 px; reuse holdout + writer; label HOLDOUT-DTI
or ORGANIZER-CONFIRMED; leakage canary AUC 0.90; run card; budget 3 experiments
or 2 hours). Study GEMSDOE32 H33-2-B2 (user-cited 0.2778) and aim higher
without copying it. Official leaderboard high 0.4413 (2026-10-21). Build GitHub
Pages with executive summary first and a one-click TIF. Open a PR and merge.
No hallucinations. Flag irregularities. Full operational text:
[`docs/SESSION_PROMPT.md`](docs/SESSION_PROMPT.md).
