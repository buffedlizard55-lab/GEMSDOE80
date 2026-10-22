# GEMSDOE80 — unique GEMS Prize GeoTIFF

**Re-read this README and `docs/SESSION_PROMPT.md` at the start of every session.**

> **OK TO DOWNLOAD AND SUBMIT (format).**
> Primary file: [`docs/downloads/gemsdoe80-h80c-strain-accom-resid-20261021T233457Z-zeros.tif`](docs/downloads/gemsdoe80-h80c-strain-accom-resid-20261021T233457Z-zeros.tif)
>
> Single-band float32, EPSG:32611, 3730×3292, transform (100, 0, 243350, 0, −100, 4508550),
> **all 12,279,160 cells finite and in [0, 1]**, no nodata tag, 0 NaN, 28,807 off-catalogue unit dots.
> This encoding avoids DrivenData’s `Predicted values must be in range [0, 1]`.
>
> **Not organizer-scored.** HOLDOUT-DTI = **0.0507** (evaluator `GEMSDOE80-evaluate_holdout.py-v1`,
> 18,313 withheld positives; spatial-block range [0.0142, 0.0383]). Verdict **promote** as a
> packaged candidate. No weekly slot spent.

GitHub Pages: executive summary and the download button are at the top of
[`docs/index.html`](docs/index.html) · how-to-submit: [`docs/executive-summary.html`](docs/executive-summary.html)

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

## 1. Why H33-2-B2 was cited at 0.2778, and whether we can beat it

H33-2-B2 (GEMSDOE32) is a **metric-geometry** file: take a 0.2708 dotted ridge,
delete every dot within 2 px of the catalogue, keep 37,654 unit pixels, write
all-finite zeros on the mirrored-template 3730×3292 grid (consensus pin, not organizer-downloaded in this sandbox — sources.md §B). DTI is a covering budget —
sparse off-catalogue mass beats a thick surface. The 0.2778 figure is
**not ORGANIZER-CONFIRMED in this checkout** (GEMSDOE32 labelled that file
UNSCORED / projected 0.2747). Official public leaderboard 2026-10-21:
**#1 xiaofanhu 0.4413**. Higher scores exist. Beating 0.2778 on the live set
requires a *different* detector that still hits hidden faults. Copying H33 is
forbidden (lane protocol). H80-C is unique vs that family (Spearman ≈ 0,
3 px overlap 0.153). Live DTI of H80-C is unknown until a receipt exists.

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
# competition rasters in data/ (gitignored): sample_submission.tif, labels.tif, training_features.tif
python scripts/prepare_data.py
python scripts/run_h80_experiments.py          # 3-experiment budget
python scripts/validate_submission.py docs/downloads/<zeros>.tif --template data/sample_submission.tif
python scripts/evaluate_holdout.py docs/downloads/<zeros>.tif
```

Holdout = hide-and-recover: withhold whole fault segments + buffer, derive
catalogue-based features from visible faults only, mask visible pixels, score
pooled DTI.

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
