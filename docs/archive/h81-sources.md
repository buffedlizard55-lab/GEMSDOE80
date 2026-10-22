# Source register — what was checked, how, and what it does NOT prove

All checks: 2026-10-22 (UTC), from this sandbox. Egress is limited to github.com, codeload,
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
| A4 | Features: `training_features.tif` (EPSG 32611, 100 m) with 19 bands listed in the page; `1m_DEM_links.csv`; labels = USGS Quaternary faults + INGENIOUS | same page, "Provided features", "Labels" | OFFICIAL-PAGE |
| A5 | Test labels = expert-identified faults **not** in the USGS database; public/private split; Phase 1 private, Phase 2 expanded labels | https://www.drivendata.org/competitions/306/competition-doe-gems/ and page 967 | OFFICIAL-PAGE |
| A6 | Competition ends 2026-12-03 11:59 p.m. UTC (homepage) | https://www.drivendata.org/competitions/306/competition-doe-gems/ | OFFICIAL-PAGE (a third-party irregularity note disputes the cutoff; see irregularities IR-11) |
| A7 | Data page requires login | https://www.drivendata.org/competitions/306/competition-doe-gems/data/ → redirects to login | OFFICIAL-PAGE |
| A8 | Rules: external data allowed if licensed and shareable with organizers | homepage "Use of external data"; rules document linked from it | OFFICIAL-PAGE (rules page itself redirects to a HeroX PDF; full rules PDF below) |
| A9 | Rules PDF (Sept 2026): $300k prize, eligibility, Phase 1 / Phase 2 structure | https://docs.nlr.gov/docs/fy26osti/96647.pdf | OFFICIAL-PAGE (sections 1.1–1.3 read) |
| A10 | Public leaderboard top rows (2026-10-21 — xiaofanhu 0.4413 #1, extradr19 0.2778 #48) | https://www.drivendata.org/competitions/306/competition-doe-gems/leaderboard/ | OFFICIAL-PAGE, "best public DW-Tversky" — not a submission receipt; no file hash shown |
| A11 | GeoDAWN: CC0 data release, DOI 10.5066/P93LGLVQ; magnetic flight lines 200 m (Area 1) and 400 m (Area 2); geoTIFF grids included | https://www.usgs.gov/data/geodawn-airborne-magnetic-and-radiometric-surveys-northwestern-great-basin-nevada-and and https://www.sciencebase.gov/catalog/item/657e1d85d34e23d3533209f7 | OFFICIAL-PAGE (ScienceBase file downloads not tested from sandbox) |
| A12 | DrivenData reference solution: reads `data/numeric_features.tif` and `data/labels.tif`; binarises labels `y<1 → 0`; writes output with the labels' CRS and transform; uses Tversky loss α=0.2, β=0.8 | https://github.com/drivendataorg/gems-prize-reference-solution (commit `aebe92f7`, notebook cells 2, 6, 19) | THIRD-PARTY-HOSTED but the organizer's own repository |
| A13 | INGENIOUS project: DOE DE-EE0009254, 2021-2025, Great Basin regional dataset compilation (DOI 10.15121/1881483), heat flow maps (DOI 10.5066/P9BZPVUC) | https://gbcge.org/current-projects/ingenious/ | OFFICIAL-PAGE |
| A14 | Walker's Lane and western Great Basin strain: GeoDAWN area 149,030 line-km, 51,857 km²; Area 1 (Clayton Valley, 200 m) Li-clay focus, Area 2 (geothermal) 400 m | https://www.usgs.gov/data/geodawn-airborne-magnetic-and-radiometric-surveys-northwestern-great-basin-nevada-and | OFFICIAL-PAGE |
| A15 | EPSG:32611 = UTM Zone 11N, WGS84 | https://epsg.io/32611 | OFFICIAL-PAGE (not fetched via shell; via page-fetch tool is allowed, but grid numbers come from template file, not epsg.io) |
| A16 | Tversky index formula (background) | https://en.wikipedia.org/wiki/Tversky_index | OFFICIAL-PAGE (Wikipedia, used only for intuition; competition's formula A1 governs) |

Not yet verified from an official source: the exact grid (shape, geotransform) of
`sample_submission.tif` — those numbers (3730×3292, transform 100,0,243350,0,-100,4508550)
were read from the restored file at `data/sample_submission.tif` (SHA 2176d08e…),
pins prove mirror consistency, not organizer authenticity. The weekly submission cap
(3/week) is from rules PDF §3.2, not re-checked on the submission page (needs login).

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

Consensus among owner-published repos is **not** organizer verification. (The hashes quoted are as printed in those repos; they were not compared with any organizer file.) These repos also appear to
redistribute competition-derived inputs; that is a licensing question for the user (see irregularities
IR-01). This project verified the restored files' hashes against `research/evidence/gems32-data-manifest.json` (pinned at `data/bridge/labels.tif` SHA 7ba308cc…); they match.

## C. Owner-reported scores (not organizer receipts)

| Claim | Where | Verified? |
|---|---|---|
| `h33-h33-2-b2-20261004T220000Z-e5eb6e7e-zeros` = 0.2778 | user brief; GEMSDOE32 site | GEMSDOE32's own `docs/downloads/submissions_manifest.json` lists the primary file with `"receipt": null`, and its audit note says `UNSCORED`; the same entry gives a catalogue-proxy DTI of 0.0067. **Not verified.** |
| 0.2778 appears on the public leaderboard | official leaderboard page, rank #48 (extradr19) on 2026-10-21 | The leaderboard exposes no file hash or name, so this row cannot be linked to the GEMSDOE32 file. **Coincidence or not is unknown.** |
| 0.3195 is the highest score | user brief | **Contradicted** by the official public leaderboard on 2026-10-21, where the top score is 0.4413 and 0.3195 does not appear in the top 50 rows read. |
| 0.4413 is the highest score | official leaderboard page, 2026-10-21, #1 xiaofanhu 26 submissions | **OFFICIAL-PAGE** — best public DW-Tversky as displayed. Not a receipt. |
| This repo's HOLDOUT-DTI 0.0607 = live DTI 0.0607 | nowhere — we label HOLDOUT-DTI, not ORGANIZER-CONFIRMED | **A projection is never written as a score** (protocol rule 3). |

## D. Method references — scientific basis for H81 (OFFICIAL-PAGE where possible)

| Ref | Claim | Source |
|-----|-------|--------|
| D1 | Step-overs / relay ramps host 32% of Great Basin geothermal systems; terminations 25%, intersections 22% | Faulds & Hinz 2015 (OSTI 1724082) — https://www.osti.gov/servlets/purl/1724082 |
| D2 | Same inventory adds producing systems prefer overlapping hard-linked (73%) with average relay width 3.4 km | Faulds et al. 2021 — Assessment of Favorable Structural Settings (GRC, OSTI 1724082) |
| D3 | 450+ step-overs mapped; mean relay width 2.8 km, high-T 3.3 km; Walker Lane: right-stepping 74%, hard-linked 59%, width 2.0 km | Giddens & Faulds 2025 (SGW) — https://pangea.stanford.edu/ERE/pdf/IGAstandard/SGW/2025/Giddens.pdf (official GRC-adjacent) |
| D4 | GeoDAWN magnetic derivative filtering (tilt, HG) delineates buried faults; joint gravity-magnetic modeling reveals blind geothermal system | Glen et al. 2022 (GRC 1034630) — https://publications.mygeoenergynow.org/grc/1034630.pdf / https://doi.org/10.5066/P93LGLVQ |
| D5 | Quaternary fault slip/dilation tendency for Great Basin (open USGS release) | Siler 2022 DOI 10.5066/P9YL58W6 — https://doi.org/10.5066/P9YL58W6 |
| D6 | 3DEP 1 m lidar used to map Quaternary faults; hillshade reveals extent | Giddens et al. 2025 SGW paper § Methods (official) |
| D7 | INGENIOUS regional dataset compilation — well/spring T & chemistry, Q-fault traces (CC BY 4.0) | https://doi.org/10.15121/1881483 and https://gdr.openei.org/submissions/1391 |
| D8 | Strain rates correlate with geothermal loci in Great Basin | Kreemer & Young 2022 (SRL) — https://doi.org/10.1785/0220220153 |

## E. What was NOT verified

- The exact grid numbers (3730×3292 etc.) from an organizer file — pins prove mirror consistency only.
- The live leaderboard scores below rank 50 — page truncated.
- The weekly submission cap on the portal itself — needs login.
- The organizer's deadline timezone — homepage says 2026-12-03 11:59 p.m. UTC; a third-party log says PDF Appendix A.1 says 5 p.m. ET (IR-11).
- That the 3DEP 1 m DEM tiles are already cloud-optimized for this footprint — only the tile list (706/716) from the competition PDF (OCR) has been seen.
