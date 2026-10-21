# GEMSDOE80 — GEMS geothermal-fault discovery and submission workflow

> **Submission status: NOT READY.** There is currently no competition raster, feature stack, label raster, holdout evaluator, submission writer, or validated output in this checkout. **Do not download or submit a TIFF from this repository yet.** No TIFF has been fabricated: without the official sample/template raster its CRS, dimensions, transform, and footprint cannot be verified, and without data a prediction would be scientifically unsupported.

## Project charter (read before every project session)

The objective is to improve the probability of identifying previously unmapped faults relevant to geothermal systems in the GEMS competition, while avoiding leakage, duplicate strategies, unsupported claims, and invalid submissions. Own the full outcome: research, validation, reproducibility, and a clear go/no-go decision. Treat failures and negative results as useful evidence. The Arena core values are **Maximize P(Win)** and **Own the Outcome**.

### Non-negotiable operating rules

1. Work only on the assigned Arena branch. Respect the active session lane; compare a candidate raster with every available registry raster **before placement and on final dots**. Stop and mark it duplicate if rank correlation exceeds 0.90 or over 70% of candidate dots lie within 3 pixels of one registry raster's dots.
2. Reuse the shared cached feature stack, `evaluate_holdout.py`, and `submission_writer.py` if they exist. They do not exist in this initial checkout. Do not create a private fork of shared tools. A holdout must hide whole fault segments with a buffer, derive catalogue-based features only from visible faults, mask visible faults pixel-exactly, and evaluate pooled DTI at alpha 0.2, beta 0.8, with a 300 m triangular kernel. Test each feature alone first; AUC > 0.90 is a leakage alarm until explained.
3. Label every performance number as **HOLDOUT-DTI** (evaluator version, withheld-positive count, 95% CI) or **ORGANIZER-CONFIRMED** (transcribed from a submission receipt). Never present a projection as a score.
4. Limit each investigation to three experiments or two hours, whichever comes first. Do not use a weekly submission slot or select/promote a submission before a candidate beats the current spatial holdout best. Promotion is a separate selector step governed by the competition's current weekly cap.
5. Make the downloadable TIFF unmistakable: expose it only after checks confirm one-band GeoTIFF, exact sample CRS/shape/transform, valid footprint, finite values, and predictions in [0,1]. Give each candidate a unique name and a comment of at most 140 characters. A file ending in `-zeros` is not evidence of its contents.
6. Research from official, free sources; record source URLs, access date, what each source supports, and any limitations. Never claim external availability or scientific validation without checking it. No manual inspection is a substitute for an automated, reproducible holdout.
7. Keep this README as the session starting point. Re-read the charter, inspect the repository and data manifests, review the run card, and verify official source links before each project session.

## Current repository review (2026-10-21)

This branch starts from a single initial commit containing only an 11-byte `README.md`. It has no application, documentation site, tests, scripts, dependencies, data, submission raster, cached feature stack, registry rasters, evaluator, writer, or data manifest. Therefore no line-by-line execution of the requested training/holdout pipeline is possible yet. `data/` placement is the genuine blocker; no external data was downloaded in this session.

The reported result lists are user-provided context, not independently verified submissions or receipts. In that context, `h33-h33-2-b2-20261004T220000Z-e5eb6e7e-zeros` is reported as 0.2778, but the supplied prompt also gives leaderboard highs of 0.4413 and 0.3195. These numbers are not mutually interpretable without dates, leaderboard/split, metric version, and receipts. The filename suffix does not establish that the raster is all zero. Thus we cannot establish why that run scored as reported, nor claim a reproducible improvement. A zero-valued prediction would have no true positives under a conventional positive-label DTI; if that is what the file contains, its nonzero reported score would require the precise competition evaluator and label/footprint semantics to explain. No such conclusion is made here.

### Current experimental decision

**Negative / not validated.** No experiment was run; no performance number is available. There is no defensible score comparison, no raster-overlap calculation, no leakage-canary AUC, and no validator output to report. The run card in [`research/run-card.json`](research/run-card.json) records these unknowns as null rather than inventing them. The validator in [`scripts/validate_submission.py`](scripts/validate_submission.py) can check a candidate against an official sample template once data and `rasterio` are available; it does not generate predictions or certify scientific merit.

## Candidate hypotheses to evaluate (not results)

These are research candidates, not claims that a result has been tested or that no other GEMSDOE branch has tried it. This checkout has no implementation to compare against; the supplied results list mentions many terrain, magnetic, radiometric, structural, hydrothermal, and lineament variants, so global uniqueness cannot be established from names alone. All four candidates can be derived from provided DEM/geophysical layers if the corresponding bands are present; confirm band names, units, masks, and grid before implementation. None should enter a slot without the required blocked holdout and registry-duplicate checks.

| Rank | Candidate and layers | Target signature and missing-fault mechanism | Main non-fault mimic | Difference / validation cost |
|---|---|---|---|---|
| 1 (tentative) | **Drainage-network lateral deflection / catchment imbalance** from the provided DEM | Extract flow paths and compare opposing catchment geometry and channel azimuth changes across a corridor. A concealed, permeability-controlling structure could deflect drainage or organize aligned channels where the fault catalogue has no mapped trace. | Lithologic contrasts, landslides, road/ditch routing, and erosional capture. | Unlike a direct scarp/curvature score, tests network-scale hydrologic organization. Low-to-medium compute; requires the 1 m DEM/flow-processing stack and buffered-segment holdout. Expected DTI benefit is uncertain; ranked first only for testability and mechanistic separation, not a promised gain. Hydro-conduit approaches are named in the supplied results, so duplicate check is mandatory. |
| 2 | **DEM channel knickpoint and longitudinal-profile breaks** from the provided DEM | Detect spatially coherent channel-profile slope changes on multiple nearby channels; fault-controlled uplift or permeability contrasts may generate aligned breaks beyond mapped traces. | Resistant bedrock boundaries, glacial/episodic erosion, dams, and DEM artifacts. | Distinct from map-plane scarp ridges; needs robust channel extraction and segment-level holdout. Medium compute. Could overlap prior topographic/hydro ideas. |
| 3 | **Cross-scale magnetic edge persistence** from the provided magnetic channels (band identities and units must be confirmed) | Compare edge orientation/position across scales and retain only stable boundaries that do not coincide with mapped faults. Buried structures may generate magnetic susceptibility contrasts without surface scarps. | Intrusive contacts, alteration, cultural magnetic sources, and survey seams. | Tests scale-persistent edge geometry rather than a single magnetic tilt/derivative or generic fusion. Low-to-medium compute, but magnetic/radiometric/tilt strategies appear in the supplied prior-results list; high duplicate risk. |
| 4 | **Radiometric boundary–drainage interaction** from available radiometric channels plus DEM | Test whether compositional boundaries align with channel deflection or drainage-density transitions away from catalogued traces; alteration/permeability architecture could affect both. | Lithologic contacts, soil/vegetation moisture, and acquisition/processing boundaries. | A specific interaction term rather than a standalone radiometric score; medium compute. Similar to listed radiometric-fusion/hydro methods, so low novelty confidence and strict overlap checks required. |

**Ranking caveat:** this is a pre-data prioritization, not an empirical DTI ranking. No official source availability is asserted beyond the user-referenced competition data and USGS survey pages; their contents, licenses, and accessibility have not been independently checked from this network. No external-data-dependent candidate is represented as viable. The immediate next task is obtaining the competition's authorized data bundle and recovering the shared evaluator/writer, not trying these ideas now.

## Submission workflow

1. Obtain the competition data through the authorized DrivenData account and place it under `data/` (the competition download redirects to login per the project context). Keep large data out of Git.
2. Inventory actual files and metadata. Confirm the official sample-submission GeoTIFF and training features/labels are present. The official reference solution instructs users to download competition data into `data/`; see [DrivenData's reference-solution repository](https://github.com/drivendataorg/gems-prize-reference-solution) and the [competition data page](https://www.drivendata.org/competitions/306/competition-doe-gems/data/).
3. Restore or implement the shared cached feature stack, segment-buffered evaluator, and writer in the shared template (none exists here). Run feature-alone leakage canaries and no more than three experiments / two hours.
4. Compare holdout DTI and its 95% CI against the current best using identical evaluator version, withheld-positive set, and 300 m kernel. Run the surface/final raster rank-correlation and 3-pixel dot-overlap duplicate tests against available registry rasters.
5. Write a unique candidate GeoTIFF and run `python scripts/validate_submission.py OUTPUT.tif --template OFFICIAL_SAMPLE.tif`. Only a successful report justifies labelling it **READY TO DOWNLOAD**. This checkout currently has no output TIFF, so there is no download link and no authorized submission recommendation.
6. On the competition's New submission form, upload that validated single-band `.tif` (or a zip containing exactly one GeoTIFF), enter the unique name and a comment <=140 characters, and rely on the organizer's displayed validation/receipt for confirmation. Never claim a score before the receipt exists.

### Official / primary review links

- [GEMS competition overview and rules](https://www.drivendata.org/competitions/306/competition-doe-gems/)
- [Competition problem description](https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/)
- [Competition about page](https://www.drivendata.org/competitions/306/competition-doe-gems/page/968/)
- [Competition data page (login required per project context)](https://www.drivendata.org/competitions/306/competition-doe-gems/data/)
- [DrivenData reference solution](https://github.com/drivendataorg/gems-prize-reference-solution)
- [USGS GeoDawn airborne magnetic and radiometric survey](https://www.usgs.gov/data/geodawn-airborne-magnetic-and-radiometric-surveys-northwestern-great-basin-nevada-and)
- [GBCGE INGENIOUS project](https://gbcge.org/current-projects/ingenious/)
- [Submission guidance PDF supplied by the project](https://docs.nlr.gov/docs/fy26osti/96647.pdf)
- [EPSG:32611](https://epsg.io/32611)
- [Tversky index background](https://en.wikipedia.org/wiki/Tversky_index) (secondary source; verify the competition's exact implemented metric instead)

External pages and data are not asserted as checked unless explicitly stated above. This environment's outbound access is restricted; unavailable sources must be checked later from an authorized environment.

## Limitations / next actions

- **Blocking:** request/obtain the authorized competition data and sample raster; restore the project template's cache, `evaluate_holdout.py`, and `submission_writer.py`. Do not bypass access controls or synthesize replacement truth data.
- Confirm the actual evaluator and metric from competition materials; reconcile the 0.2778, 0.4413, and 0.3195 figures by date, split, metric, and receipt.
- Add tests with synthetic rasters for CRS/shape/transform mismatch, nodata/NaN inside the footprint, values outside [0,1], and valid output. The validator here is a safety check, not an evaluator.
- Build the requested public-facing status/dashboard only after its source-of-truth manifests, metrics, and artifact generation are implemented. A static page must never imply an unvalidated file is ready.
- No competition login is available to this process; no automatic submission, receipt, current leaderboard scrape, or slot selection was performed.
