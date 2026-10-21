# GEMSDOE80 — GEMS geothermal-fault discovery and submission workflow

> **Submission status: NOT READY to submit (by design).** This checkout contains *illustrative placeholder* TIFFs
> whose values are verified in **[0, 1]** (fixing the portal error `Predicted values must be in range [0, 1]`)
> and whose grids are self-consistent (EPSG:32611 @ 100 m), but whose exact shape/geotransform match to the
> official `sample_submission.tif` is **UNVERIFIED** (template is auth-gated on DrivenData). **Do not submit any
> TIFF below as-is without re-gridding to the official template.** Values are checked line by line; grid alignment
> must be checked on a logged-in machine.

**Goal:** place top of the GEMS Prize Challenge leaderboard (DrivenData competition #306) with a scientifically grounded,
auditable, reproducible fault-probability submission. This repo builds a *unique* TIF submission (never copied from any
previous GEMSDOE site — the primary is linear dotted-ridge, not prior blobs), validates it line by line without human
intervention, and publishes it through a clean GitHub Pages site whose executive summary + go/no-go status are obvious
at the top of the page.

> Current public reference points: GEMSDOE32 `h33-h33-2-b2-20261004T220000Z-e5eb6e7e-zeros: 0.2778` (verified on the
> GEMSDOE32 site; 37,654 dots, 0 within 200 m of catalogue) is the highest among the 80+ GEMSDOE sites the prompt
> lists; the prompt also cites DrivenData leaderboard highs of 0.4413 / 0.3195 / 0.3262 (nchuzhoy 0.3262 on 2026-10-04
> per GEMSDOE32's verified leaderboard read) — these need reconciliation with dates/splits/receipts and are flagged
> below. Target: beat 0.2778 with a unique, credit-dense dotted synthesis, then aim higher. PhD-level analysis of why
> 0.2778 won and how to beat it is in `docs/hypotheses.md` and the executive summary.

## 0. Executive summary — how to make a submission

The competition entry is a **single-band GeoTIFF (.tif)**, or a `.zip` containing exactly one GeoTIFF, with:
- CRS: EPSG:32611 (UTM zone 11N), resolution 100 m
- Exact match of CRS / shape / geotransform with the official `sample_submission.tif`
- Single band, `float32`, values in **[0, 1]** (higher = more likely fault; fixes the "Predicted values must be in range [0, 1]" rejection — caused by sentinel `-3.4e38` or NaN leaking outside [0,1])
- All-finite, `nodata=None` (GEMSDOE32-verified safe format: every cell finite, in [0,1], no sentinel)

Steps:
1. **Get the data (login required).** Download `training_features.tif`, `labels.tif`, `sample_submission.tif`,
   `1m_DEM_links.csv` from https://www.drivendata.org/competitions/306/competition-doe-gems/data/ into `data/`.
   Keep large data out of Git (see `.gitignore`).
2. **Inspect placement.** `python scripts/prepare_data.py` — confirms files and prints the official grid report.
3. **Generate (unique, not copied).** `python scripts/generate_submission.py --template data/sample_submission.tif --features data/training_features.tif --mode dotted-ridge --dots 37654`
   → writes `<name>.tif`, `<name>.zip`, and `<name>.tif.json` sidecar under `docs/submissions/`.
4. **Validate (no human needed).** `python scripts/validate_submission.py docs/submissions/<name>.tif --template data/sample_submission.tif`
   — submit only when `format_valid` is true and `illustrative_fallback` is false (the latter requires the organizer's confirmation; this repo never self-confirms).
5. **Submit.** On DrivenData's New submission form, upload the validated `.tif` (or `.zip`), enter its unique name and
   a comment of ≤ 140 characters (see `research/run-card.json` — e.g. `Dotted Walker-Lane ridges 37654 dots, 2% scatter, credit-dense thinning, strain+mag prior; label-free, unique.`), and rely on the organizer's displayed receipt.

The downloadable files and their status are published on the GitHub Pages site (root → `docs/`): **executive summary
subpage at the very top (see below), download buttons, and "is it OK to submit?" tables appear immediately.**

- **Subpage:** `docs/executive-summary.html` — the single page to read before you download. It explains exactly what to put in *File to submit* and *Note (optional)*, how the `[0,1]` error was fixed, and when the file becomes submission-ready. Linked from the top of `docs/index.html`.
- **Dashboard:** `docs/index.html` — executive summary (above) + downloadable artifacts + readiness table + why 0.2778 won + hypotheses queue + run card + official links.

## 1. Is the downloadable TIF OK to submit? (current state)

| Check | Result — PRIMARY ★ dotted-ridge v2 (recommended) | Result — secondary structural-prior v2 | Legacy v1 |
|---|---|---|---|
| File | `docs/submissions/gemsdoe80-dotted-ridge-v2-20261021T231505Z.tif` (622 KB) | `gemsdoe80-structural-prior-v2-20261021T231523Z.tif` (39 MB) | `gemsdoe80-structural-prior-v1-20261021T214821Z.tif` (39 MB) |
| SHA-256 | `dc43a85e55f9687b8f5d4797eed7ae1097e5d9c6b8be659a3811f057c1212219` | `a7a476f00630b0c784c740211798e04e09d71ddd4e1a4c8f19938ae1991f713a` | `0afa36baf428e73f817bf53e623a1d2b279ff2852dec814bc8a444ee6781b858` |
| Single-band float32 GeoTIFF | PASS | PASS | PASS |
| CRS EPSG:32611 @ 100 m | PASS (3000×4000, nodata None) | PASS (3000×4000, nodata None) | PASS (3000×4000, nodata NaN) |
| Values in [0, 1] | PASS (0.02 – 1.0, **0 cells outside**) | PASS (0.05 – 0.90, **0 outside**) | PASS (0.05 – 0.90, 0 outside) |
| No NaN/Inf/sentinel -3.4e38 inside | PASS (0 NaN, 0 Inf, 0 sentinel, nodata None) | PASS (0 NaN, 0 Inf, 0 sentinel) | PASS (0 NaN, self-validates) |
| Self-validation format_valid | PASS | PASS | PASS |
| Exact grid vs official `sample_submission.tif` | **UNVERIFIED** — template auth-gated; fallback is 3000×4000 vs official ~3292×3730 per GEMSDOE32 | same | same |
| Holdout DTI + 95% CI (HOLDOUT-DTI) | **n/a** — `labels.tif` auth-gated | n/a | n/a |
| Leakage canary (per-feature AUC) | **n/a** — labels not available (synthesis uses no labels, so no leakage) | n/a | n/a |
| Registry lane check | **n/a** — no registry raster, but linear dotted vs blob = expected <<0.90 correlation, <70% dot overlap | n/a | n/a |
| **Verdict** | **Unique, downloadable, [0,1]-fixed — but illustr. fallback, NOT submission-ready until re-gridded** | Same | Same |
| Unique name (≤140-char note) | `gemsdoe80-dotted-ridge-v2-20261021T231505Z` / `Dotted Walker-Lane ridges 37654 dots, 2% scatter, credit-dense thinning, strain+mag prior; label-free, unique.` (137 chars) | structural-prior note (135 chars) in run-card | see sidecar |

**Download:** buttons at the top of `docs/index.html` and on `docs/executive-summary.html`. Also `docs/submissions/<name>.zip` for portal upload. No submission slot has been spent.

## 2. Repository layout

```
GEMSDOE80/
├── README.md                    # this charter (re-read at every session)
├── .gitignore                   # data/, artifacts/, caches never committed
├── index.html                   # GitHub Pages root → ./docs/
├── docs/                        # GitHub Pages content (served at root)
│   ├── index.html               # research dashboard + executive summary + download (top)
│   ├── executive-summary.html   # SUBPAGE: exact submission steps (File + Note)
│   ├── submissions/             # candidate TIFFs + zips + sidecars (generated, tracked)
│   │   ├── gemsdoe80-dotted-ridge-v2-20261021T231505Z.tif  ★ primary (622 KB, 37654 dots)
│   │   ├── gemsdoe80-structural-prior-v2-20261021T231523Z.tif
│   │   └── gemsdoe80-structural-prior-v1-20261021T214821Z.tif  (legacy)
│   └── hypotheses.md            # 9 ranked geological hypotheses (H1–H9, 4 new in 2026-10-21)
├── research/
│   └── run-card.json            # machine-readable run card (experiment 2 of 3)
├── scripts/
│   ├── download_competition_data.sh   # official URLs + manual placement (login required)
│   ├── prepare_data.py              # data/ inventory + official grid report
│   ├── generate_submission.py       # unsupervised [0,1] GeoTIFF builder (structural-prior + dotted-ridge)
│   └── validate_submission.py       # format checks + holdout DTI + leakage canary + lane check
├── tests/
│   └── test_validate_submission.py  # integration tests for the validator
└── reports/                     # validator JSON reports (generated)
```

## 3. Scripts

- `scripts/generate_submission.py` — builds a single-band float32 GeoTIFF, EPSG:32611, 100 m, values strictly in [0, 1], all-finite `nodata=None` (GEMSDOE32-safe, fixes sentinel bug). If `data/sample_submission.tif` exists, its CRS/geotransform/shape/dtype are used exactly. Otherwise an illustrative fallback is built and flagged as not submission-ready. Modes: `structural-prior` = continuous geothermal structural-propensity from `training_features.tif` bands (conductivity; dilatation/shear/2nd-invariant; gravity anomaly+slope; magnetics RTP/TMI + vertical/horizontal slope; earthquake density), percentile-normalized and averaged; `dotted-ridge` = unique sparse synthesis of 12 Walker-Lane-parallel linear ridges + relay ramps + dotted sampling (≈37654 dots, 2% scatter) that raises mean credit per dot. No training labels are used — no label leakage.
- `scripts/validate_submission.py` — format validation against the official template (GTiff, 1 band, float32, CRS, shape, transform, [0,1] range, sentinel/NaN/Inf checks) plus, when `labels.tif`/`training_features.tif` are available: holdout DTI (alpha 0.2, beta 0.8, 300 m triangular kernel, buffered whole-segment holdout, bootstrap 95% CI), per-feature leakage-canary AUC (AUC > 0.90 = leakage until proven otherwise), and registry duplicate-lane check (rank-correlation >0.90 or >70% dots within 3 px = duplicate). Writes `<tif>.validate.json` and prints JSON; exit 0 = format OK.
- `scripts/prepare_data.py` — inventories `data/` and prints the `sample_submission.tif` grid report when present.
- `scripts/download_competition_data.sh` — official source URLs; the DrivenData bundle requires login and is not downloaded automatically from this environment.

## 4. Method (parallel-run protocol)

1. **LANE:** stay in the dotted-ridge lane (linear Walker-Lane ridges + thinning + scatter); if raster rank-correlation with any registry raster exceeds 0.90, or >70% of dots fall within 3 px of one registry raster's dots, log duplicate and stop. Check on the surface before placement AND on the final dots.
2. **REUSE, DON'T REBUILD:** shared tools live in `scripts/` (cached feature stack, `evaluate_holdout.py` analogue in `validate_submission.py`, `submission_writer.py` analogue in `generate_submission.py`); fix once in the template, never private forks.
3. **LABEL EVERY NUMBER** as HOLDOUT-DTI (evaluator version, withheld positives, 95% CI) or ORGANIZER-CONFIRMED (transcribed from a submission receipt). A projection is never written as a score.
4. **LEAKAGE CANARY:** test each feature alone on the holdout; AUC above 0.90 means leakage until proven otherwise.
5. **RUN CARD:** one JSON card per experiment at `research/run-card.json`. Negative results are deliverables. Current: experiment 2 of 3 (dotted-ridge v2), verdict negative (not promotable: holdout + grid unverifiable).
6. **BUDGET:** stop after 3 experiments or 2 hours. Do not spend a submission slot on an idea that has not beaten the current holdout best; promotion is a separate selector step within the weekly cap shown on the submission page.

## 5. Hypotheses — why 0.2778 won, and 4 new contrarian detectors

See `docs/hypotheses.md`: nine ranked candidates. H1 is the baseline; H2–H5 are prior baselines; **H6–H9 are new (added 2026-10-21)** to satisfy the prompt's "3–5 not yet tried" request, each with layers, physical signature (transform), why it catches uncatalogued faults, non-fault mimic, how it differs from anything already implemented in this repo, expected DTI gain, cost, and the free official source required:

- **H6 worm persistence + Euler depth clustering** — multi-scale tilt-angle worms (100→600 m) that survive across scales + 2D Euler depth clustering (0–2 km); catches blind basement faults. Differs from H1's simple gradient magnitude. Needs GeoDAWN magnetics (training_features or ScienceBase 657e1d85d34e23d3533209f7, public).
- **H7 knickpoint-aligned swarm (chi, TOP NEW PICK)** — chi-transformed channel profiles → coherent knickpoints within 500 m × 3 km Walker-Lane corridors; catches blind faults offsetting channels with no scarp. Needs 3DEP 1 m DEM (public via https://apps.nationalmap.gov/downloader/).
- **H8 strain-gradient Laplacian (cheapest)** — ∇²(dilatation) × |∇ε| pinning of extensional bowls; catches locked fault tips/relay ramps. Needs only training_features.tif (login), very low cost.
- **H9 radiometric × drainage coincidence** — K/Th boundary only where drainage azimuth changes >15° & catchment asymmetry >0.3; separates structural alteration from facies. Needs ScienceBase radiometrics + DEM.

Ranking: H7 (moderate-high, +0.02–0.05) > H6 (moderate) > H8 (low-moderate, cheapest) > H9 > others. Validation order: H8 first (cheapest, uses provided features) on the spatially-blocked holdout; no submission slot before a candidate beats the holdout best. Updated ranking table and PhD analysis of the 0.2778 credit-density mechanism are in `docs/hypotheses.md`.

The dotted-ridge v2 synthesis already demonstrates the **emission** side of the 0.2778 mechanism (thinning + scatter + catalogue exclusion), but the **detector** side still needs real holdout validation once data are placed.

## 6. Run card

`research/run-card.json` — experiment 2 of 3: dotted-ridge v2 synthesis (37,654 linear dots, unique vs prior blobs), verdict **negative** (not promotable: holdout DTI uncomputable without `labels.tif`; exact grid match unverifiable without `sample_submission.tif`; lane/leakage also not runnable without registry/labels). No submission slot spent. Prior experiment 1 (structural-prior v1 blob) is recorded inside the run card as `prior_experiment`. Raster SHA-256 and full validator output are recorded there and in each sidecar.

## 7. Irregularities flagged

- DrivenData data/download/submission pages require login; this environment has no DrivenData credentials and the Dropbox mirror links are not on the allowed-hosts list, so competition inputs cannot be downloaded or validated automatically here. All downloads must be done manually on a logged-in machine — the scripts then validate without human inspection.
- Leaderboard figures conflict across sources (0.4413 vs 0.3195 vs 0.3262 vs 0.2778 cited as "highest"); treat third-party numbers as unverified until seen on a dated official leaderboard or an ORGANIZER-CONFIRMED receipt. Latest verified read we can reference is GEMSDOE32's snapshot: #1 nchuzhoy 0.3262, #2 DARD 0.3195 on 2026-10-04.
- A run name ending in `zeros` does not prove raster values were all zero — GEMSDOE32's `zeros` file is all-finite with values in [0,1] and no nodata tag; the file the prompt cites (`h33-h33-2-b2-...-zeros: 0.2778`) is reported as 37,654 dots, not zeros. Without its TIFF, metric version, labels, and receipt, the cause of any reported score is unknown but the `zeros` suffix is not a literal description.
- The downloadable TIFFs are illustrative fallback grids (approximate Walker Lane box, 3000×4000) not the official competition grid (3292×3730 per GEMSDOE32). Downloading is safe; submitting any of them as-is would fail the grid-match check (but not the [0,1] check — that error is already fixed). Re-grid to the official template before a real submission.
- Illustrative synthesis uses no training labels, so no holdout DTI can be claimed — any "projected 0.27xx" would be a model, not a score.

## 8. Official / primary review links

- [GEMS competition overview and rules](https://www.drivendata.org/competitions/306/competition-doe-gems/)
- [Competition problem description](https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/)
- [Competition about page](https://www.drivendata.org/competitions/306/competition-doe-gems/page/968/)
- [Competition data page (login required)](https://www.drivendata.org/competitions/306/competition-doe-gems/data/)
- [Leaderboard](https://www.drivendata.org/competitions/306/competition-doe-gems/leaderboard/)
- [DrivenData reference solution](https://github.com/drivendataorg/gems-prize-reference-solution)
- [USGS GeoDawn airborne magnetic and radiometric surveys](https://www.usgs.gov/data/geodawn-airborne-magnetic-and-radiometric-surveys-northwestern-great-basin-nevada-and)
- [GeoDawn item (USGS ScienceBase, 657e1d85d34e23d3533209f7)](https://www.sciencebase.gov/catalog/item/657e1d85d34e23d3533209f7)
- [USGS 3DEP / The National Map elevation](https://apps.nationalmap.gov/downloader/)
- [Blakely & Simpson 1986 — worm persistence](https://doi.org/10.1190/1.1442063)
- [GBCGE INGENIOUS project](https://gbcge.org/current-projects/ingenious/)
- [Submission guidance PDF (NLR, 96647)](https://docs.nlr.gov/docs/fy26osti/96647.pdf)
- [EPSG:32611](https://epsg.io/32611)
- [Tversky index background](https://en.wikipedia.org/wiki/Tversky_index) (secondary source; verify the competition's exact metric instead — see Problem description's formula)
- [GEMSDOE32 docs — analysis of the 0.2778 scorer](https://buffedlizard55-lab.github.io/GEMSDOE32/docs/index.html)

Core values: **Maximize P(Win)** and **Own the Outcome**. Verify line by line; no hallucinations. Negative results are deliverables. The executive-summary subpage makes the download and the go/no-go obvious without a human inspector.

---

## APPENDIX — SESSION PROMPT (project charter, re-read at the start of every session)

The user's full session prompt is reproduced below as the working charter. Key requirements: (1) generate a UNIQUE TIF submission for the GEMS competition — never copy a previous submission; (2) it must be obvious whether the generated TIF is OK to download and submit; (3) inspection/validation must be automated through the project, not by a human inspector; (4) follow the parallel-run protocol (lane, reuse, labeled numbers, leakage canary, run card, budget); (5) study the highest cited prior score (GEMSDOE32 `h33-h33-2-b2-...-zeros: 0.2778`) and determine whether a higher scoring unique submission can be generated; (6) propose 3–5 ranked geological hypotheses with layers, signature, missing-fault mechanism, and differentiation; (7) validate the top candidate on a spatially-blocked holdout before spending a submission slot; (8) no DrivenData auth is available — all data placement is manual; (9) build a GitHub Pages site with clean UI and an executive summary subpage at the very beginning explaining exactly how to submit; (10) open a pull request and merge to main; (11) no hallucinations — verify line by line, flag irregularities.

<prompt>
Review the repo. (full session instructions as provided by the user, including the GEMSDOE site result list, competition leaderboard links, DrivenData pages, reference solution, USGS GeoDAWN, GBCGE INGENIOUS, EPSG:32611, Tversky index, submission format PDF, submission-page error "Predicted values must be in range [0, 1]", unique-name + note guidance, core values, 3-pass implementation instruction, PR/merge instruction, and limitations/questions)
</prompt>
