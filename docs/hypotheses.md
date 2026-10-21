# Candidate geological hypotheses — GEMS fault prediction

**Status (2026-10-21): none validated.** The holdout cannot be computed: `labels.tif` and
`training_features.tif` are behind the DrivenData login. Ranks below are *qualitative priors* for
testability and cost. **No expected DTI gain is given as a number, because no number can be
justified before a holdout run.** Any figure shown later must be labelled HOLDOUT-DTI with its
protocol, count and CI.

Scoring facts (official page, "Performance metric"): DTI = TP_w / (TP_w + 0.2·FP_w + 0.8·FN_w),
triangular kernel R = 300 m = 3 px. Test faults are expert-labelled faults that are **not** in the
USGS database, so a prior that copies the catalogue is expected to score poorly on the test set.
(A third-party log argues the same thing with its own data. That claim is not verified here; see `docs/irregularities.md` IR-22.)

## Candidates

### C1 — Linear ridge continuity from magnetic edges (highest expected discovery value)
- **Layers:** magnetics RTP and TMI, and their vertical and horizontal slope bands
  (`training_features.tif`; exact band names unverified, see IR-08).
- **Signature:** long, continuous ridges in the horizontal-gradient magnitude or Hessian ridge
  response, kept only if they connect over ≥ 2 km. Targets concealed fault traces that show as
  linear magnetic contrasts before they reach the surface.
- **Why it can find uncatalogued faults:** the USGS catalogue is built from mapped surface traces.
  A linear magnetic discontinuity that is not on the map is exactly the missing-fault class.
- **Named non-fault mimic:** intrusive dikes and lithologic contacts (linear magnetic contrast, no
  displacement); survey-line artefacts. GeoDAWN Area 2 flight lines are 400 m apart (USGS page), so
  fine linear artefacts are under-sampled there.
- **Difference from the repo:** H1 (`scripts/generate_submission.py`) averages |slope| magnitudes
  with no linearity or continuity test.
- **Cost:** low–medium (CPU, features already provided).
- **Validation:** holdout after data placement (run E2 below).

### C2 — Multi-physics coincidence gating (cheapest test of the H1 idea)
- **Layers:** magnetics edge response (from C1), isostatic-gravity horizontal gradient, conductivity
  or depth-to-conductive-base contrast, strain-rate magnitude, earthquake density.
- **Signature:** pixels where at least three *independent* property contrasts coincide within a
  300 m radius. Built as a min or geometric mean across terms, not H1's arithmetic mean.
- **Why it can find uncatalogued faults:** a fault zone is expected to perturb several properties at
  once. Single-layer anomalies are cheap to produce, and a coincidence rule suppresses them.
- **Named non-fault mimic:** lithologic contacts that change density, magnetisation and conductivity
  together.
- **Difference from the repo:** H1 is an arithmetic mean over terms. C2 is a coincidence (AND-type)
  rule with a spatial radius.
- **Cost:** low.
- **Validation:** holdout after data placement (run E3 below).

### C3 — Metric-aware emission (scoring-side, not a new fault signal)
- **Layers:** whichever probability surface is chosen (H1, C1 or C2).
- **Mechanism:** the FP term is Σ p(x)·(1 − k), so low-confidence mass costs 0.2 per unit, and the
  FN term costs 0.8 per missed label. Emitting dots at a spacing matched to the 3 px kernel, and
  thresholding p on low-confidence pixels, trades a small FP cost for recall. The effect is
  measured on the holdout only.
- **Named mimic:** none; this optimises the scorer, not geology. Any apparent gain can come from
  dot density alone, so it must be reported against a density-matched control.
- **Difference from prior GEMSDOE sites:** they thin dots using catalogue proximity. A third-party
  audit records a catalogue-proxy DTI of 0.0067 for the file behind the 0.2778 claim (IR-05). That
  proxy is itself disputed (IR-19), so it cannot rank candidates.
- **Cost:** very low.

### C4 — 1 m lidar scarp and lineament detection (high value, blocked here)
- **Layers:** 1 m DEM tiles listed in `1m_DEM_links.csv` (DrivenData login) or public USGS 3DEP
  tiles (not reachable from this sandbox; obtainability **not verified** in this session).
- **Signature:** sub-metre scarps and linear topographic breaks, detected with curvature and
  hillshade-edge filters, then resampled to 100 m.
- **Mimic:** glacial and alluvial lineations, roads, canal banks.
- **Status:** cannot be validated without the DEM tiles. Listed for the user to decide whether the
  download is worth the effort.

### C5 — Catalogue distance-decay negative control (baseline)
- **Purpose:** a control, not a candidate. If a catalogue-only prior scores well on the internal
  holdout, the holdout is leaking. Any candidate must beat this control on the same holdout.

## Ranking (qualitative)

| Rank | Candidate | Expected discovery value | Cost | Blocked by |
|---|---|---|---|---|
| 1 | C1 ridge continuity | high (if magnetic lineaments coincide with new faults) | low–medium | labels + features download |
| 2 | C2 coincidence gating | medium | low | labels + features download |
| 3 | C3 metric-aware emission | none on its own; affects score mechanics | very low | holdout |
| 4 | C4 lidar scarps | high | high | DEM download, not reachable here |
| 5 | C5 control | none (sanity check) | very low | labels |

## Validation order (budget: 3 experiments or 2 h; no submission slot)

E1 is tooling and does not count against the experiment budget. E2–E4 are the three experiments.

1. **E1 (done, tooling only):** repair the validator DTI and canary; unit tests pass. Not a
   competition experiment.
2. **E2:** H1 baseline on the holdout (sets the reference "holdout best").
3. **E3:** C1 on the same holdout. Promote only if it beats E2 with a CI that does not overlap.
4. **E4:** C2 or C3 layered on the winner. This is the third and last experiment in the budget.

A submission slot is spent only after a candidate beats the current holdout best and the validator
passes. Weekly cap: three submissions (official rules §3.2).
