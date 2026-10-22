# Lane H80-C review addendum — 2026-10-22 UTC

## Scope and parallel-lane boundary

This note covers only the H80-C strain-accommodation lane from the active Arena branch. While PR #12 was open, `main` advanced through PR #9 with a different H80-DPVC detector and repository reorganization. That new lane is preserved unchanged; this note does not compare, select, promote, or validate its detector. The shared site/run-card currently on `main` names H80-DPVC as primary. Promotion between lanes belongs to the separate selector, not this lane review.

## H80-C result and file status

- H80-C is a historical candidate already in this checkout before the review. No new raster was generated in this session.
- Before the concurrent merge, the TIFF bytes were independently checked: one float32 band; EPSG:32611; 3730×3292; transform `(100, 0, 243350, 0, -100, 4508550)`; all values finite in `[0,1]`; no nodata tag; 28,807 positive pixels; SHA-256 `d9c4ce410af5c413195aa3a463c2e4c6c223e26c6c8cf3a25a8a95c4b7a2bcc5`.
- The official DrivenData sample template was not available. The remembered grid pin and prior comparison are owner-mirror-derived, not organizer-authenticated. The concurrent main branch added a mirrored sample raster and a synthetic-template revalidation; neither authenticates organizer provenance. That merge also removed the H80-C TIFF from the current Pages downloads. Therefore the historical raster is **not currently available from this branch's download link and is not certified for submission**. No receipt or live score exists.
- The prior H80-C HOLDOUT-DTI is 0.0507 under `GEMSDOE80-evaluate_holdout.py-v1`, with 18,313 withheld positives. Its four-quadrant range `[0.0142, 0.0383]` is descriptive, not a calibrated 95% CI. The run-card CI is unavailable. No H80-C experiment was rerun: its budget was already 3/3 and current competition inputs/registry were not authenticated in this branch session.
- Historical registry uniqueness evidence was limited to three named local raster references. It does not prove uniqueness against all GEMSDOE submissions. No current surface/final-dot check was possible.

## Review findings from the pre-reorganization branch

The branch's user-facing `OK TO SUBMIT` wording overstated what its historic mirror-based report proved. Its optional validator also mixed format checks with static-raster holdout scoring and a positive-only AUC computation. I changed those files and tests in commit `6425b3f` and ran **26 passed, 2 skipped** on that pre-reorganization tree. Because `main` subsequently replaced that implementation, those code edits are not asserted as changes to the current H80-DPVC implementation. The test result is for the exact pre-merge commit only.

The current repository's H80-DPVC validator/evaluator, added by the other lane, is outside this review's scope. Its own source/provenance and claims should be assessed separately; mirror hashes establish integrity to the mirror, not authenticity from DrivenData.

## Verdict

**Negative for H80-C slot promotion.** No current authenticated template, no re-computed holdout with valid CI, no full registry, and no organizer receipt. Do not spend a weekly slot based on H80-C. No weekly slot was used by this review.

Manual reference links: [competition rules](https://www.drivendata.org/competitions/306/competition-doe-gems/), [problem description](https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/), [competition data page (login-gated)](https://www.drivendata.org/competitions/306/competition-doe-gems/data/), [organizer-hosted reference solution](https://github.com/drivendataorg/gems-prize-reference-solution), [USGS GeoDAWN](https://www.usgs.gov/data/geodawn-airborne-magnetic-and-radiometric-surveys-northwestern-great-basin-nevada-and).
