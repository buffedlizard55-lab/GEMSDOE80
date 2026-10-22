# Historical candidate copies (not currently certified)

This folder contains byte-identical copies of the H80-C historical candidate also exposed under
`docs/downloads/`. They are retained for audit/history; they are **not** approved for weekly-slot
promotion. The official competition template and data are absent from the current checkout.
A previous format audit matched an owner-mirrored template, which is not organizer-authenticated.

The current local byte check confirms single-band float32, finite [0,1] values, and the file's
stored CRS/shape/transform. This does not independently establish that those metadata match the
official `sample_submission.tif`, nor guarantee portal acceptance. See the status banner at
[`docs/index.html`](../index.html), the run card, and `reports/review-20261022.md`.

No new experiments or weekly submission slots were used in the 2026-10-22 review. The H80
experiment budget was already exhausted (3/3).