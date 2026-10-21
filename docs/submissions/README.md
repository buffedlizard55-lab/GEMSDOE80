# Submission files

**Status: no submission-ready TIF exists in this repository.**

The earlier placeholder raster (`gemsdoe80-structural-prior-v1-20261021T214821Z`) was removed on
2026-10-21. It was built from invented Gaussian "structural zones" on a 3000×4000 grid that does not
match any verified official grid, so it must not be downloaded or submitted. It remains in git
history if an audit is needed.

A file appears here only when `scripts/generate_submission.py` succeeds against the official
`data/sample_submission.tif` and `data/training_features.tif`, and `scripts/validate_submission.py`
passes. Even then `download_ready` stays `false` until an organizer receipt is recorded.
