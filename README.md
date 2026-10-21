# GEMSDOE80 — auditable fault-discovery research

## Submission status: **DO NOT SUBMIT the existing TIFF**

[Executive summary](docs/executive-summary.html) · [Website](docs/index.html) ·
[Three-pass audit](reports/review-20261021.md) · [Run card](research/run-card.json)

The inherited downloadable TIFF is synthetic. Automated validation against a hash-pinned
owner-mirrored template now **fails shape and transform**. A range in [0,1] alone does
not make a competition submission valid. No new submission, experiment, or slot was used
in this review. Scientific promotion is **negative / blocked**, not a low measured score.

## Start every session here — project charter

User requirements (condensed, **not a verbatim copy** of the lengthy original prompt):

* Build a unique scientific GeoTIFF, never copy another submission; prominently state
  whether it is safe to submit and offer an easy download only with truthful status.
* Fault discovery supports geothermal exploration; do not conflate faults with vents.
* Assign one geological lane; stop if surface rank correlation exceeds 0.90 or more than
  70% of final dots lie within 3 px of any registry entry. Check before placement and after dots.
* Reuse the template cached stack, `evaluate_holdout.py`, `submission_writer.py`; no private
  evaluator fork. Whole-segment buffered hide-and-recover; visible-only catalogue features;
  exact visible-fault mask; pooled DTI (requested alpha 0.2, beta 0.8, 300 m triangular kernel).
* Per-feature holdout AUC >0.90 is a leakage canary until explained. Use negatives too.
* Score labels: HOLDOUT-DTI with evaluator version, withheld positives and 95% CI, or
  ORGANIZER-CONFIRMED with a submission-page receipt. User reports are not verified scores.
* Compare 3–5 geological hypotheses before implementation; name layers, transform,
  missing-fault mechanism, non-fault mimic, novelty limits and cost.
* At most three experiments or two hours. Promotion is a separate selector step and must
  beat the current holdout best within the verified weekly cap. Never auto-submit here.
* Automate inspections, use free official data with auditable provenance, flag unknowns;
  no fabricated evidence, scores, leaderboard feeds, readiness, or uniqueness claims.
* Three passes: implement/check; review/fix; reconcile with requirements. Maintain a clean
  Pages site and executive subpage, unique name and comment ≤140 characters when warranted.
* Open and merge a pull request when checks pass. Preserve negative results and next steps.
* Core values: **Maximize P(Win)** through honest risk/evidence tradeoffs; **Own the Outcome**
  end-to-end, including fixing unsafe inherited behavior instead of concealing blockers.

The full original prompt remains in the session conversation. The previous README falsely
claimed to reproduce it in full but contained a placeholder; that claim has been removed.

## Reproduce the audit

```bash
python -m venv .venv
.venv/bin/pip install -r requirements.txt
.venv/bin/python scripts/restore_mirrors.py  # small template; integrity ≠ organizer authenticity
.venv/bin/python -m pytest -q
.venv/bin/python scripts/validate_submission.py \
  docs/submissions/gemsdoe80-structural-prior-v1-20261021T214821Z.tif \
  --template data/sample_submission.tif  # expected exit 1: wrong grid
```

`restore_mirrors.py --core` can restore owner-mirrored features/labels (~420 MB),
using immutable GitHub refs and size/hash checks. Data stay ignored. The historical
`download_competition_data.sh` is only guidance, not an authenticated downloader.
`generate_submission.py` now exits 2 without writing anything until the shared writer,
assigned lane and evaluation gates are restored. It must not silently invent geology.
The private baseline and synthetic feature-generation helpers have been retired.

Format validation and promotion are separate. `format_valid` is a local format check;
`download_ready` stays false because scientific gates are unavailable. Organizer acceptance
comes **after upload** and is not a prerequisite for local format validation.

## Research findings and limits

[Ranked hypotheses](docs/hypotheses.md) and [source-grounded audit](reports/review-20261021.md)
explain the cited GEMSDOE32 artifact. It prunes catalogue-near dots from an existing sparse
base; `zeros` describes outside-footprint encoding, not an all-zero predictor. This may
reduce false-positive mass, but a higher leaderboard score is not proven or promised.
The cited 0.2778 and conflicting leaderboard highs are user reports, not receipts.

Shared cache/evaluator/writer, complete registry, current best and weekly cap are absent.
No single method paragraph designates this session's lane. That blocks an honest new
experiment before placement. External websites outside GitHub/package hosts are not
accessible from this sandbox; official URLs are review leads, not verified downloads.
CPU training is supported by the official reference; GPU is helpful, not mandatory.

Official reference, inspected at pinned revision:
https://github.com/drivendataorg/gems-prize-reference-solution/tree/aebe92f7c8a990f0e3443451b7a825d9afd6336b

Other official sources (not independently fetched here):
* https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/
* https://www.drivendata.org/competitions/306/competition-doe-gems/data/
* https://www.drivendata.org/competitions/306/competition-doe-gems/leaderboard/
* https://doi.org/10.5066/P93LGLVQ
* https://docs.nlr.gov/docs/fy26osti/96647.pdf

Next session must restore the shared infrastructure and assigned lane before pursuing
hypothesis 1. No manual visual sign-off can substitute for the automated gates.
