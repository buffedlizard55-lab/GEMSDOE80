# Source register — what was checked, how, and what it does NOT prove

All checks: 2026-10-21 (UTC), from this sandbox. Egress is limited to github.com, codeload,
api.github.com, registry.npmjs.org, pypi.org, files.pythonhosted.org. Official pages were read
through the page-fetch tool. Raw DrivenData data, Dropbox, USGS file servers and epsg.io were
**not** downloadable from the shell (curl returned no connection).

Labels used below: **OFFICIAL-PAGE** = read from an official page today; **THIRD-PARTY** = read
from another team's public repository (not authenticated, not an organizer source); **NOT
VERIFIED** = stated somewhere but not checked.

## A. Official competition sources (OFFICIAL-PAGE)

| # | Claim | Source (verifiable) | Status |
|---|---|---|---|
| A1 | Metric: distance-weighted Tversky, α=0.2, β=0.8, triangular kernel R=300 m (3 px at 100 m); p(x)∈[0,1]; TP_w / FP_w / FN_w definitions | https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/ (section "Performance metric") | OFFICIAL-PAGE |
| A2 | Submission: single layer float32, values 0–1, EPSG 32611, 100 m, same bounds as training data, outside bounds null/NaN | same page, "Submission format" | OFFICIAL-PAGE |
| A3 | `sample_submission.tif` "predicts total fault absence" (i.e. all-zero) | same page | OFFICIAL-PAGE |
| A4 | Features: `training_features.tif` (EPSG 32611, 100 m) with bands listed in the page; `1m_DEM_links.csv`; labels = USGS Quaternary faults + INGENIOUS | same page, "Provided features", "Labels" | OFFICIAL-PAGE |
| A5 | Test labels = expert-identified faults **not** in the USGS database; public/private split; Phase 1 private, Phase 2 expanded labels | https://www.drivendata.org/competitions/306/competition-doe-gems/ and page 967 | OFFICIAL-PAGE |
| A6 | Competition ends 2026-12-03 11:59 p.m. UTC (homepage) | https://www.drivendata.org/competitions/306/competition-doe-gems/ | OFFICIAL-PAGE (a third-party irregularity note disputes the cutoff; see irregularities) |
| A7 | Data page requires login | https://www.drivendata.org/competitions/306/competition-doe-gems/data/ → redirects to login | OFFICIAL-PAGE |
| A8 | Rules: external data allowed if licensed and shareable with organizers | homepage "Use of external data"; rules document linked from it | OFFICIAL-PAGE (rules page itself redirects to a HeroX PDF; full rules PDF below) |
| A9 | Rules PDF (Sept 2026): $300k prize, eligibility, Phase 1 / Phase 2 structure | https://docs.nlr.gov/docs/fy26osti/96647.pdf | OFFICIAL-PAGE (sections 1.1–1.3 read) |
| A10 | Public leaderboard top rows (see README §4) | https://www.drivendata.org/competitions/306/competition-doe-gems/leaderboard/ | OFFICIAL-PAGE, "best public DW-Tversky" — not a submission receipt |
| A11 | GeoDAWN: CC0 data release, DOI 10.5066/P93LGLVQ; magnetic flight lines 200 m (Area 1) and 400 m (Area 2); geoTIFF grids included | https://www.usgs.gov/data/geodawn-airborne-magnetic-and-radiometric-surveys-northwestern-great-basin-nevada-and and https://www.sciencebase.gov/catalog/item/657e1d85d34e23d3533209f7 | OFFICIAL-PAGE (ScienceBase file downloads not tested from sandbox) |
| A12 | DrivenData reference solution: reads `data/numeric_features.tif` and `data/labels.tif`; binarises labels `y<1 → 0`; writes output with the labels' CRS and transform; uses Tversky loss α=0.2, β=0.8 | https://github.com/drivendataorg/gems-prize-reference-solution (commit `aebe92f7`, notebook cells 2, 6, 19) | THIRD-PARTY-HOSTED but the organizer's own repository |

Not yet verified from an official source: the exact grid (shape, geotransform) of
`sample_submission.tif`, the band names/descriptions inside `training_features.tif`, the weekly
submission cap, and the organizer's exact deadline timezone convention.

## B. Third-party claims about the grid and data (THIRD-PARTY, NOT VERIFIED)

Several public GEMSDOE repositories state the same grid: EPSG:32611, 3730 × 3292 cells, geotransform
`(100, 0, 243350, 0, -100, 4508550)`, and the SHA-256 `2176d08e…54d35cbc` for `sample_submission.tif`.

| Repository (commit read) | Where | Note |
|---|---|---|
| buffedlizard55-lab/GEMSDOE32 (`f62038b5`) | `data/prepared_manifest.json`, `docs/downloads/*-audit.json` | `data/restore_receipt.json` lists a local storage root (`/home/user/GEMSDOE32/.cache/gems_data`) and file hashes. It does not document where the competition files came from. GEMSDOE30's log (IR-30-024) says comparable files came from owner mirrors. |
| buffedlizard55-lab/GEMSDOE70 (`158aec19`) | `docs/downloads/*.json` | repeats the grid |
| buffedlizard55-lab/GEMSDOE77 (`6c363745`) | `data/restore_receipt.json`, docs | repeats the grid |
| buffedlizard55-lab/GEMSDOE49 (`8ab48f22`) | `data/bridge/*.json` | repeats the grid |
| buffedlizard55-lab/GEMSDOE66 (`2734d5d4`) | `artifacts/*` | repeats the grid |

Consensus among owner-published repos is **not** organizer verification.  (The hashes quoted are as printed in those repos; they were not compared with any organizer file.) These repos also appear to
redistribute competition-derived inputs; that is a licensing question for the user (see irregularities
IR-01). This project does **not** download or use those copies.

The earlier placeholder in this repo used a different grid (3000 × 4000 at `(250000, 4450000)`) that
is not supported by any source checked here, and it was removed.

## C. Owner-reported scores (not organizer receipts)

| Claim | Where | Verified? |
|---|---|---|
| `h33-h33-2-b2-20261004T220000Z-e5eb6e7e-zeros` = 0.2778 | user brief; GEMSDOE32 site | GEMSDOE32's own `docs/downloads/submissions_manifest.json` lists the primary file with `"receipt": null`, and its audit note says `UNSCORED`; the same entry gives a catalogue-proxy DTI of 0.0067. **Not verified.** |
| 0.2778 appears on the public leaderboard | official leaderboard page, rank #48 (extradr19) | The leaderboard exposes no file hash or name, so this row cannot be linked to the GEMSDOE32 file. **Coincidence or not is unknown.** |
| 0.3195 is the highest score | user brief | **Contradicted** by the official public leaderboard, where the top score is 0.4413 and 0.3195 does not appear in the top 50 rows read. |

## D. Method references used only as background (NOT VERIFIED)

- Tversky index (Wikipedia): background only. The competition's own formula (A1) governs.
- The "GEMSDOE30 proxy frame" analysis in GEMSDOE30's irregularity log (IR-30-039, IR-30-046) argues
  that catalogue-proxy scores rank the owner's files in the opposite order to the public leaderboard.
  This is a third-party argument. It is plausible given that A5 says test faults are *not* in the
  catalogue, but it is **not verified** here.
