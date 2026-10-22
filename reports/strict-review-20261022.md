# H80-DPVC strict protocol review — 2026-10-22

## Pass 1: inspect and independently verify bytes
Read the standing README, prior-session reviews, preregistration, shared writer and
evaluator entrypoints, generator, evidence and tests. The shipped raster is NOT a
new experiment of this session. Independently re-opened TIFF/ZIP: single float32
band, all finite values in [0,1], grid consistent with mirrored pin, hash matches
run card; ZIP contains exactly those bytes. See machine release audit.

## Pass 2: bugs and assumptions
- **Blocking:** the inherited density exemption violates the user's all-registry
  overlap rule. Recorded literal overlap is 1.0; label protocol duplicate and stop.
  This does not prove copying. A dense positive raster makes this rule very
  restrictive, but that is not authority to silently weaken it.
- Earlier tests asserted promotion rather than testing the requested gate. New
  regression cases cover dense-raster rejection, absent evidence, NaN, infinity,
  invalid ranges and threshold boundaries. Generator now consults the strict gate.
- Removed submission approval from packaging; a separate selector must choose slots.
- Authentication/provenance: mirrors prove integrity against owner pins, not official
  origin. Competition inputs and the full raster registry are absent in this checkout.
  Holdout and canary results are inherited JSON evidence, not independently re-run.
- Research note's section 1.1 dilates both masks and multiplies the dilations. That
  is **not the same metric** as the truth-cell/dot-cell definition in the shared
  metric code. Its constant-TP argument is a hypothesis, not causal proof from
  unavailable hidden labels. Do not use its numeric extrapolations for promotion.
- Neither user-reported 0.2778 nor conflicting maximum-score claims have a receipt.
  Historical prose claiming ORGANIZER-CONFIRMED owner-reported numbers is invalid.
- An improved fault map is not a direct prediction of geothermal vents; permeability,
  temperature and fluids are additional requirements. No vent discovery is verified.

## Pass 3: requirements reconciliation
The one-click TIFF remains available **for research only**; the executive summary
and landing page now lead with DO NOT SUBMIT. No new TIFF was generated: doing so
would violate the explicit duplicate-stop and exhausted experiment budget. No slot
was consumed. Negative results are the deliverable required by the protocol.

The five ranked hypotheses already recorded in research/h80_hypotheses.md and
registry/preregistration.json are historical, not five new viable ideas. Their
expected improvements are unvalidated projections, not scores. No new idea can be
proposed as obtainable solely from a source URL: the external official hosts are
outside this sandbox's outbound allowlist. Future work must first authenticate
inputs, reconcile the metric prose with code, audit fold-local feature derivation
and CI resampling, and compare against the **current best**, not just H33. The
existing artifact's inherited HOLDOUT-DTI is 0.16199 (gems57-pooled-hide-v2,
22,663 positive appearances across two possibly overlapping seed draws; reported
95% CI [0.14569, 0.17759]). Independence of that pooled CI remains unverified.

## Source ledger and limits
- Organizer repository obtainable via GitHub API this session, HEAD
  `aebe92f7c8a990f0e3443451b7a825d9afd6336b`:
  https://github.com/drivendataorg/gems-prize-reference-solution
- Competition problem / format (manual-review reference, not fetched this session):
  https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/
- Competition data (login-gated per user; not independently fetched this session):
  https://www.drivendata.org/competitions/306/competition-doe-gems/data/
- Official GeoDAWN source (manual-review reference; not obtainable from allowed hosts):
  https://www.usgs.gov/data/geodawn-airborne-magnetic-and-radiometric-surveys-northwestern-great-basin-nevada-and
- Reference H33 site (owner source, not organizer receipt):
  https://buffedlizard55-lab.github.io/GEMSDOE32/docs/index.html

GitHub Pages API confirmed a built docs site in this session. Other sessions have
open PRs; no other lane or PR is merged by this review. No live leaderboard claims
are made. Full scientific validation and a new promoted unique artifact remain
unfulfilled, explicitly blocked rather than fabricated.
