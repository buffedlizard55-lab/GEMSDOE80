# Submission files (validated)

**Status: PRIMARY is `gemsdoe80-h81a-relay-20261022T010259Z-zeros.tif` (H81-A relay-ramp, 7,415 dots, HOLDOUT 0.0607). Archived: `gemsdoe80-h80c-strain-accom-resid-20261021T233457Z-zeros.tif` (H80-C, 28,807 dots, HOLDOUT 0.0507). Both are format-valid and downloadable from the GitHub Pages site.**

A file is copied into this folder only after all of these pass:

1. `scripts/run_h81_experiments.py` (or `run_h80_experiments.py`) writes a candidate to `docs/downloads/` via `scripts/submission_writer.py`.
2. `scripts/validate_submission.py` reports `format_valid: true` against the official `sample_submission.tif` (grid, CRS, transform, values in [0,1], NaN outside the footprint — re-read from bytes).
3. The candidate beats the current holdout best on the internal HOLDOUT-DTI protocol (hide-and-recover whole segments + buffer, mask visible, pooled DTI α=0.2 β=0.8 R=300 m), with a CI, and the leakage canary is at or below AUC 0.90.
4. Lane uniqueness vs GEMSDOE family: Spearman |r| ≤0.90 and 3 px overlap ≤0.70.
5. A person chooses the weekly slot (three per week, official rules §3.2).

No weekly slot has been spent. Promotion here means "package as a candidate."

Archived H80 remains valid for re-use; its reports are in `reports/` and `research/run-card.json` history. The earlier placeholder rasters (`gemsdoe80-structural-prior-v1/v2`, `gemsdoe80-dotted-ridge-v2`) were removed from this folder because their revalidation reports show they fail the official grid; they are kept in git history.

Current primary audit:
- zeros SHA256 `701ed7fed8611861f2fe2ec6437c850aede0a6a792b5fc3dbbca10191953bbbd`
- zip SHA256 `5fae213ad545cc41c2ec2c1199496ad9f54f0151a963b48989976c1fbf1c4420`
- 7,415 positive pixels, 0 on-catalogue, 0 NaN inside footprint, 12,279,160 finite in [0,1]
- HOLDOUT-DTI 0.0607, CI [0.0541, 0.0636], beats H80 0.0507 and random 0.0123 on same split
- Duplicate check: vs H80 Spearman 0.02 / overlap 0.22 → not duplicate
