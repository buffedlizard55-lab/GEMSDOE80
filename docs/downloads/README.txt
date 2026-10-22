PRIMARY SUBMISSION FILE (use this):
  gemsdoe80-h81a-relay-20261022T010259Z-zeros.tif

OK TO UPLOAD: yes (all 12,279,160 cells finite, values in [0, 1],
official 3730x3292 EPSG:32611 grid). Organizer score: none.
Previous candidate archived:
  gemsdoe80-h80c-strain-accom-resid-20261021T233457Z-zeros.tif

Name: gemsdoe80-h81a-relay-20261022T010259Z-zeros
Note: H81A relay tip-link d8-35 NMS 7415; off-cat buf4; HOLDOUT-DTI only; unique vs H33/H80 (85 chars)

sha256 zeros: 701ed7fed8611861f2fe2ec6437c850aede0a6a792b5fc3dbbca10191953bbbd
sha256 zip:   5fae213ad545cc41c2ec2c1199496ad9f54f0151a963b48989976c1fbf1c4420
sha256 nan:   2eaf644564a1db30935f72add58e5ebdcd62f4f885260507b503555731aa410b

Format: single-band float32, EPSG:32611, 100 m, 3730x3292,
transform (100, 0, 243350, 0, -100, 4508550), no nodata tag, 0 NaN inside footprint,
0 on-catalogue pixels, 7415 positive dots.

HOLDOUT-DTI 0.0607 (95% CI [0.0541, 0.0636], 18,313 withheld positives, 5903 dots on holdout)
beats random 0.0123 and beats H80 0.0507 on same split.
Not organizer-scored. See docs/run-card.json and research/run-card.json.

How to submit:
1. Download docs/downloads/gemsdoe80-h81a-relay-20261022T010259Z-zeros.tif (or .zip)
2. Open https://www.drivendata.org/competitions/306/competition-doe-gems/submissions/
3. File to submit = that .tif/.zip
4. Note = H81A relay tip-link d8-35 NMS 7415; off-cat buf4; HOLDOUT-DTI only; unique vs H33/H80

Validate locally:
  python scripts/validate_submission.py docs/downloads/gemsdoe80-h81a-relay-20261022T010259Z-zeros.tif --template data/sample_submission.tif
  python -m pytest tests/test_shipped_tif.py -v
