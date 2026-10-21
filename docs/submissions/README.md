# Submission files (validated only)

**Status: empty. No file here is approved for download or submission.**

A file is copied into this folder only after all of these pass, and never by the generator:

1. `scripts/generate_submission.py` writes a candidate to `artifacts/candidates/` (git-ignored).
2. `scripts/validate_submission.py` reports `format_valid: true` against the official
   `sample_submission.tif` (grid, CRS, transform, values in [0,1], NaN outside the footprint).
3. The candidate beats the current holdout best on the internal HOLDOUT-DTI protocol, with a CI,
   and the leakage canary is at or below AUC 0.90.
4. A person chooses the weekly slot (three per week, official rules §3.2).

The earlier placeholder rasters (`gemsdoe80-structural-prior-v1/v2`, `gemsdoe80-dotted-ridge-v2`) were
removed from this folder. Their validation reports are kept under `reports/`. The rasters themselves are
in git history and fail the official grid check.
