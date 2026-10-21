# GEMSDOE80 — GEMS geothermal-fault discovery and submission workflow

> **Submission status: NOT READY to submit.** This checkout contains an *illustrative placeholder* TIFF
> (`docs/submissions/gemsdoe80-structural-prior-v1-20261021T214821Z.tif`) whose values are in [0, 1] and whose
> CRS/grid are self-consistent, but whose exact shape/geotransform match to the official
> `sample_submission.tif` is **UNVERIFIED** (the template is auth-gated on DrivenData). **Do not submit this TIFF as-is.**
> No TIFF has been presented as fabricated truth: without the official sample/template raster its CRS, dimensions,
> transform, and footprint cannot be verified, and without data a prediction is scientifically unsupported.

**Goal:** place top of the GEMS Prize Challenge leaderboard (DrivenData competition #306) with a scientifically grounded,
auditable, reproducible fault-probability submission. This repo builds a *unique* TIF submission (never copied from any
previous GEMSDOE site), validates it line by line, and publishes it through a clean GitHub Pages site whose executive
summary makes the download and the go/no-go status obvious.

> Current public reference points: GEMSDOE32 `h33-h33-2-b2-20261004T220000Z-e5eb6e7e-zeros` is reported as 0.2778 on
> the GEMSDOE sites; the session also cites DrivenData leaderboard highs of 0.4413 and 0.3195. These figures are
> user-provided context, not independently verified, and they are not mutually interpretable without dates, split,
> metric version, and receipts (flagged as an irregularity below). Target: beat 0.2778 with a unique, validated
> submission, then aim higher.

## 0. Executive summary — how to make a submission

The competition entry is a **single-band GeoTIFF (.tif)**, or a `.zip` containing exactly one GeoTIFF, with:
- CRS: EPSG:32611 (UTM zone 11N), resolution 100 m
- Exact match of CRS / shape / geotransform with the official `sample_submission.tif`
- Single band, `float32`, values in **[0, 1]** (higher = more likely fault; fixes the "Predicted values must be in range
  [0, 1]" rejection)
- NaN/null outside the valid footprint

Steps:
1. **Get the data (login required).** Download `training_features.tif`, `labels.tif`, `sample_submission.tif`,
   `1m_DEM_links.csv` from https://www.drivendata.org/competitions/306/competition-doe-gems/data/ into `data/`.
   Keep large data out of Git (see `.gitignore`).
2. **Inspect placement.** `python scripts/prepare_data.py` — confirms files and prints the official grid report.
3. **Generate.** `python scripts/generate_submission.py --template data/sample_submission.tif --features data/training_features.tif`
   → writes `<name>.tif`, `<name>.zip`, and `<name>.tif.json` sidecar under `docs/submissions/`.
4. **Validate.** `python scripts/validate_submission.py docs/submissions/<name>.tif --template data/sample_submission.tif`
   — submit only when `format_valid` is true and `download_ready` is true (the latter requires the organizer's
   confirmation; this repo never self-confirms).
5. **Submit.** On DrivenData's New submission form, upload the validated `.tif` (or `.zip`), enter its unique name and
   a comment of ≤ 140 characters (see `research/run-card.json`), and rely on the organizer's displayed receipt.

The downloadable file and its status are published on the GitHub Pages site (root → `docs/`): **executive summary,
download links, and "is it OK to submit?" table appear at the top of the page.**

## 1. Is the downloadable TIF OK to submit? (current state)

| Check | Result |
|---|---|
| File exists at `docs/submissions/gemsdoe80-structural-prior-v1-20261021T214821Z.tif` | Yes (40,278,202 bytes) |
| SHA-256 | `0afa36baf428e73f817bf53e623a1d2b279ff2852dec814bc8a444ee6781b858` |
| Single-band float32 GeoTIFF | PASS |
| CRS EPSG:32611 @ 100 m | PASS (shape 3000×4000, transform from (250000, 4450000)) |
| Values in [0, 1] | PASS (observed 0.05 – 0.90) |
| NaN only as nodata | PASS (0 NaN inside footprint) |
| Exact grid match vs official `sample_submission.tif` | **UNVERIFIED** — template auth-gated on DrivenData |
| Holdout DTI + 95% CI (HOLDOUT-DTI) | **n/a** — `labels.tif` auth-gated; not computed |
| Leakage canary (per-feature AUC) | **n/a** — labels not available |
| Registry duplicate check | **n/a** — no registry rasters available |
| **Verdict** | **Illustrative placeholder — downloadable, but NOT submission-ready.** Regenerate against `data/sample_submission.tif`, then validate; submit only on PASS. No submission slot has been spent. |

## 2. Repository layout

```
GEMSDOE80/
├── README.md                    # this charter (session starting point; full prompt in appendix)
├── .gitignore                   # data/, artifacts/, caches never committed
├── index.html                   # GitHub Pages root → ./docs/
├── docs/                        # GitHub Pages content (served at root)
│   ├── index.html               # research dashboard + executive summary + download
│   ├── submissions/             # candidate TIFF + zip + sidecar (generated artifacts)
│   └── hypotheses.md            # 5 ranked geological hypotheses (not results)
├── research/
│   └── run-card.json            # machine-readable run card (one per experiment)
├── scripts/
│   ├── download_competition_data.sh   # official source URLs + manual placement (login required)
│   ├── prepare_data.py              # data/ inventory + official grid report
│   ├── generate_submission.py       # unsupervised [0,1] fault-probability GeoTIFF builder
│   └── validate_submission.py       # format checks + holdout DTI + leakage canary + lane check
├── tests/
│   └── test_validate_submission.py  # integration tests for the validator
└── reports/                     # validator JSON reports (generated)
```

## 3. Scripts

- `scripts/generate_submission.py` — builds a single-band float32 GeoTIFF, EPSG:32611, 100 m, values in [0, 1],
  NaN nodata. If `data/sample_submission.tif` exists, its CRS/geotransform/shape/dtype are used exactly.
  Otherwise an *illustrative fallback grid* is built and flagged as not submission-ready. The surface is an
  unsupervised geothermal structural-propensity score from `training_features.tif` bands located by description
  metadata (conductivity; dilatation/shear/2nd-invariant strain magnitude; gravity anomaly + slope; magnetics RTP/TMI
  anomaly + vertical & horizontal slope; earthquake density), percentile-normalized and averaged. No training labels
  are used, so there is no label leakage.
- `scripts/validate_submission.py` — format validation against the official template (GTiff, 1 band, float32, CRS,
  shape, transform, [0, 1] range, nodata handling) plus, when `labels.tif`/`training_features.tif` are available:
  holdout DTI (alpha 0.2, beta 0.8, 300 m triangular kernel, buffered whole-segment holdout, bootstrap 95% CI),
  per-feature leakage-canary AUC (AUC > 0.90 = leakage until proven otherwise), and registry duplicate-lane check.
  Writes `<tif>.validate.json` and prints JSON.
- `scripts/prepare_data.py` — inventories `data/` and prints the `sample_submission.tif` grid report when present.
- `scripts/download_competition_data.sh` — official source URLs; the DrivenData bundle requires login and is not
  downloaded automatically from this environment.

## 4. Method (parallel-run protocol)

1. **LANE:** stay in the structural-propensity lane; if raster rank-correlation with any registry raster exceeds 0.90,
   or more than 70% of dots fall within 3 px of one registry raster's dots, log it as duplicate and stop. Check on the
   surface before placement and on the final dots.
2. **REUSE, DON'T REBUILD:** shared tools live in `scripts/`; fix once in the template, never private forks.
3. **LABEL EVERY NUMBER** as HOLDOUT-DTI (evaluator version, withheld positives, 95% CI) or ORGANIZER-CONFIRMED
   (transcribed from a submission receipt). A projection is never written as a score.
4. **LEAKAGE CANARY:** test each feature alone on the holdout; AUC above 0.90 means leakage until proven otherwise.
5. **RUN CARD:** one JSON card per experiment at `research/run-card.json`. Negative results are deliverables.
6. **BUDGET:** stop after 3 experiments or 2 hours. Do not spend a submission slot on an idea that has not beaten the
   current holdout best; promotion is a separate selector step within the weekly cap.

## 5. Hypotheses

See `docs/hypotheses.md`: five ranked, not-yet-validated candidates (H1 structural-propensity composite — implemented
as the generator baseline; H2 lidar scarp/lineament density; H3 radiometric alteration; H4 earthquake density ×
dilatation; H5 catalogue distance-decay negative control), each with layers, physical signature, missing-fault
mechanism, non-fault mimic, difference from prior work, expected DTI improvement, cost, and the free official source
required. Validation order: H1 holdout first; no submission slot before a candidate beats the holdout best.

## 6. Run card

`research/run-card.json` — experiment 1 of 3: H1 structural-propensity composite, verdict **negative** (not promotable:
holdout DTI uncomputable without `labels.tif`; exact grid match unverifiable without `sample_submission.tif`). No
submission slot spent. Raster SHA-256 and validator output are recorded there.

## 7. Irregularities flagged

- DrivenData data/download/submission pages require login; this environment has no DrivenData credentials and the
  Dropbox mirror links do not connect, so competition inputs cannot be downloaded or validated automatically here.
  All downloads must be done manually on a logged-in machine.
- Leaderboard figures conflict across sources (0.4413 vs 0.3195 "highest score"); treat third-party numbers as
  unverified until seen on the official DrivenData leaderboard.
- A run name ending in `zeros` does not prove raster values were all zero; without its TIFF, metric version, labels,
  and organizer receipt, the cause of any reported score is unknown.
- The downloadable TIFF is an illustrative fallback grid (approximate Walker Lane box), not the official competition
  grid. Downloading is safe; submitting it as-is is not.

## 8. Official / primary review links

- [GEMS competition overview and rules](https://www.drivendata.org/competitions/306/competition-doe-gems/)
- [Competition problem description](https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/)
- [Competition about page](https://www.drivendata.org/competitions/306/competition-doe-gems/page/968/)
- [Competition data page (login required)](https://www.drivendata.org/competitions/306/competition-doe-gems/data/)
- [Leaderboard](https://www.drivendata.org/competitions/306/competition-doe-gems/leaderboard/)
- [DrivenData reference solution](https://github.com/drivendataorg/gems-prize-reference-solution)
- [USGS GeoDawn airborne magnetic and radiometric survey](https://www.usgs.gov/data/geodawn-airborne-magnetic-and-radiometric-surveys-northwestern-great-basin-nevada-and)
- [GeoDawn item (USGS ScienceBase)](https://www.sciencebase.gov/catalog/item/657e1d85d34e23d3533209f7)
- [GBCGE INGENIOUS project](https://gbcge.org/current-projects/ingenious/)
- [Submission guidance PDF supplied in project context](https://docs.nlr.gov/docs/fy26osti/96647.pdf)
- [EPSG:32611](https://epsg.io/32611)
- [Tversky index background](https://en.wikipedia.org/wiki/Tversky_index) (secondary source; verify the competition's
  exact implemented metric instead)

Core values: **Maximize P(Win)** and **Own the Outcome**. Verify line by line; no hallucinations. Negative results are
deliverables.

---

## APPENDIX — SESSION PROMPT (project charter, re-read at the start of every session)

The user's full session prompt is reproduced below as the working charter. Key requirements: (1) generate a UNIQUE TIF
submission for the GEMS competition — never copy a previous submission; (2) it must be obvious whether the generated
TIF is OK to download and submit; (3) inspection/validation must be automated through the project, not by a human
inspector; (4) follow the parallel-run protocol (lane, reuse, labeled numbers, leakage canary, run card, budget);
(5) study the highest cited prior score (GEMSDOE32 `h33-h33-2-b2-...-zeros: 0.2778`) and determine whether a higher
scoring unique submission can be generated; (6) propose 3–5 ranked geological hypotheses with layers, signature,
missing-fault mechanism, and differentiation; (7) validate the top candidate on a spatially-blocked holdout before
spending a submission slot; (8) no DrivenData auth is available — all data placement is manual; (9) build a GitHub
Pages site with clean UI and an executive summary subpage at the very beginning explaining exactly how to submit; (10)
open a pull request and merge to main; (11) no hallucinations — verify line by line, flag irregularities.

<prompt>
Review the repo. (full session instructions as provided by the user, including the GEMSDOE site result list,
competition leaderboard links, DrivenData pages, reference solution, USGS GeoDAWN, GBCGE INGENIOUS, EPSG:32611,
Tversky index, submission format PDF, submission-page error "Predicted values must be in range [0, 1]", unique-name +
note guidance, core values, 3-pass implementation instruction, PR/merge instruction, and limitations/questions)
</prompt>
