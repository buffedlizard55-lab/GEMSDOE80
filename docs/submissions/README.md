# Submission files (validated only)

**Status: one validated candidate is present. Use the zeros file only.**

| File | Role | SHA-256 |
|---|---|---|
| `gemsdoe80-h80c-strain-accom-resid-20261021T233457Z-zeros.tif` | **Primary — download and submit this** (all 12,279,160 cells finite, in [0,1], no nodata tag) | `d9c4ce410af5c413195aa3a463c2e4c6c223e26c6c8cf3a25a8a95c4b7a2bcc5` |
| `gemsdoe80-h80c-strain-accom-resid-20261021T233457Z-zeros.zip` | Same GeoTIFF in a single-TIF zip (portal also accepts this) | see `docs/run-card.json` |
| `gemsdoe80-h80c-strain-accom-resid-20261021T233457Z-nan.tif` | NaN-outside twin — **DO NOT SUBMIT** (NaN outside footprint can trigger `Predicted values must be in range [0, 1]`) | `1a5949ad343d6cf2842bb7081a3cee7e72a4831850a48f68d1af4d9a215e8887` |

Unique name: `gemsdoe80-h80c-strain-accom-resid-20261021T233457Z-zeros`
Note (104/140): `H80C cat-residual strain-accom NMS dots 28807; off-catalogue buf4; HOLDOUT-DTI only; unique vs H33/H19-5`

Format was validated in-project by re-reading the written bytes (`scripts/validate_submission.py`,
`tests/test_shipped_tif.py`). **No organizer score exists** — HOLDOUT-DTI 0.0507 is a
catalogue hide-and-recover proxy, not a leaderboard number. No weekly slot was spent by this repo;
spending one is a separate human selector decision. Grid provenance: mirrored template pin
(consensus 3730x3292 EPSG:32611 grid shared by owner repos), not an organizer download made
inside this sandbox — see `docs/sources.md` section B and IR-04/IR-27+.
