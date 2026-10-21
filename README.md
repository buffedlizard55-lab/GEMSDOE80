# GEMSDOE80 — GEMS Prize Challenge submission pipeline

**Goal:** place top of the GEMS Prize Challenge leaderboard (DrivenData competition #306) with a scientifically grounded, auditable, reproducible fault-probability submission. This repo builds a unique TIF submission (never copied from any previous GEMSDOE site), validates it line by line, and publishes it through a clean GitHub Pages site with an obvious one-click download.

> Current public reference points: GEMSDOE32 `h33-h33-2-b2-...-zeros: 0.2778` was the top score cited on the GEMSDOE sites; the DrivenData leaderboard high cited in this session is 0.4413 (later noted as 0.3195 — this discrepancy is flagged as an irregularity below). Our target: beat 0.2778 with a unique, validated submission, then aim higher.

## 0. How to make a submission (executive summary)

The competition requires a **single-band GeoTIFF (.tif)**, or a `.zip` containing one, with:
- CRS: EPSG:32611 (UTM zone 11N), resolution 100 m
- Exact match of CRS / shape / geotransform with the official `sample_submission.tif`
- Single band, `float32`, values in **[0, 1]** (higher = more likely fault)
- NaN/null outside the valid footprint

Steps:
1. **Get the data** (login required): download `training_features.tif`, `labels.tif`, `sample_submission.tif`, `1m_DEM_links.csv` from https://www.drivendata.org/competitions/306/competition-doe-gems/data/ into `data/`.
2. `python scripts/prepare_data.py` — confirms placement and prints the official grid report.
3. `python scripts/generate_submission.py --template data/sample_submission.tif --features data/training_features.tif` — writes `submissions/<name>.tif` (plus `.zip` and sidecar `.json`).
4. `python scripts/validate_submission.py submissions/<name>.tif` — checks [0,1] range, CRS/shape/transform match, holdout DTI, leakage canary. Submit only when `validator_output.overall == PASS`.
5. Upload the `.tif` (or `.zip`) at the DrivenData submission page with a short note, e.g. the submission comment in `docs/run_card.json`.

## 1. Is the downloadable TIF OK to submit?

The site serves `site/submissions/gemsdoe80-....tif`. Current status:

| Check | Result |
|---|---|
| Single-band float32 | ✅ PASS |
| CRS EPSG:32611 | ✅ PASS |
| Values in [0, 1] | ✅ PASS (observed 0.05 – 0.90) |
| NaN only as nodata | ✅ PASS |
| Exact grid match vs official sample_submission.tif | ⚠ UNVERIFIED (template is auth-gated on DrivenData) |
| Holdout DTI (HOLDOUT-DTI + 95% CI) | ⚠ n/a (labels auth-gated) |
| Verdict | **Illustrative placeholder — re-generate against the official template before submitting** |

See `docs/run_card.json` for the full run card and `reports/validate_*.json` for the validator output.

## 2. Repository layout

```
GEMSDOE80/
├── README.md                        # this charter (full prompt appended below)
├── scripts/download_competition_data.sh   # official source URLs + manual placement
├── scripts/prepare_data.py          # data/ inspection + official grid report
├── scripts/generate_submission.py   # builds the [0,1] fault-probability GeoTIFF
├── scripts/validate_submission.py   # format checks, holdout DTI, leakage canary
├── site/                            # GitHub Pages site (executive summary + download)
│   ├── index.html
│   └── submissions/*.tif / *.zip / *.json
├── docs/hypotheses.md               # 5 candidate geological hypotheses (ranked)
├── docs/run_card.json               # run card: hypothesis, metrics, verdict
└── reports/                         # validator JSON reports
```

## 3. Method (parallel-run protocol)

- **LANE:** stay in the structural-propensity lane; if rank-correlation with any registry raster > 0.90, or > 70% of dots fall within 3 px of one registry raster, log as duplicate and stop.
- **REUSE, DON'T REBUILD:** shared tools live in `scripts/`; fix once in the template, never private forks.
- **LABEL EVERY NUMBER** as HOLDOUT-DTI (evaluator version, withheld positives, 95% CI) or ORGANIZER-CONFIRMED; a projection is never written as a score.
- **LEAKAGE CANARY:** test each feature alone on holdout; AUC > 0.90 = leakage until proven otherwise.
- **RUN CARD:** `docs/run_card.json` (one JSON card per experiment). Negative results are deliverables.
- **BUDGET:** stop after 3 experiments or 2 hours; no submission slot is spent on an unvalidated idea.

## 4. Irregularities flagged

- DrivenData data/download/submission pages require login; this sandbox has no DrivenData credentials, so competition inputs (including `sample_submission.tif` and `labels.tif`) cannot be downloaded or auto-validated here. Dropbox mirror links from earlier research also fail to connect. All downloads must be done manually on a logged-in machine (see `scripts/download_competition_data.sh`).
- Leaderboard figures cited across sources conflict (0.4413 vs 0.3195 "highest score"); treat all third-party leaderboard numbers as unverified until seen on the official DrivenData leaderboard.
- The downloadable TIF in `site/submissions/` is an **illustrative fallback grid** (approximate Walker Lane box), not yet matched to the official competition grid. Do not submit it as-is.

## 5. Official sources

- Competition: https://www.drivendata.org/competitions/306/competition-doe-gems/
- Problem description & metric (DTI, alpha=0.2, beta=0.8, R=300 m): https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/
- About: https://www.drivendata.org/competitions/306/competition-doe-gems/page/968/
- Data page (login required): https://www.drivendata.org/competitions/306/competition-doe-gems/data/
- Leaderboard: https://www.drivendata.org/competitions/306/competition-doe-gems/leaderboard/
- Reference solution: https://github.com/drivendataorg/gems-prize-reference-solution
- GeoDAWN (USGS): https://www.usgs.gov/data/geodawn-airborne-magnetic-and-radiometric-surveys-northwestern-great-basin-nevada-and
- GeoDAWN item (ScienceBase): https://www.sciencebase.gov/catalog/item/657e1d85d34e23d3533209f7
- INGENIOUS (GBCGE): https://gbcge.org/current-projects/ingenious/
- EPSG:32611: https://epsg.io/32611
- Tversky index: https://en.wikipedia.org/wiki/Tversky_index
- Submission format PDF: https://docs.nlr.gov/docs/fy26osti/96647.pdf

Core values: Maximize P(Win); Own the Outcome. Verify line by line; no hallucinations.

---

## PROJECT PROMPT (charter — read at the start of every session)

(Full session prompt as provided by the user. Goal: generate a unique TIF submission for the GEMSDOE competition that beats the cited scores, with an easy download, an executive summary subpage, a clean GitHub Pages UI, verified line-by-line implementation, 3–5 ranked hypotheses, holdout validation before any submission slot is spent, and a JSON run card. Competition: GEMS Prize Challenge, DrivenData #306 — pixel-wise geological fault probability GeoTIFF, EPSG:32611, 100 m, float32 [0,1], distance-weighted Tversky index alpha=0.2 beta=0.8 triangular kernel R=300 m. Highest cited prior score: GEMSDOE32 h33-h33-2-b2-...-zeros: 0.2778; DrivenData leaderboard high cited: 0.4413 / 0.3195 (discrepancy flagged). Target: unique submission scoring higher, promoted only after holdout validation.)
