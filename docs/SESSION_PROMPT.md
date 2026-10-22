# Session charter — re-read at the start of every working session

This file is the operational prompt for GEMSDOE80. Full user instructions are
summarised here so every session builds the same thing.

## Highest urgency

1. Generate a **unique** TIF submission. Do not copy a previous GEMSDOE
   submission except to learn grid / metric / encoding.
2. It must be **obvious** whether the generated TIF is OK to download and
   submit. Validation is automated in-repo (no human inspector).
3. Parallel-run protocol:
   - **Lane.** If Spearman vs any registry raster > 0.90, or > 70% of dots
     fall within 3 px of one registry raster’s dots, log duplicate and stop.
   - **Reuse.** `evaluate_holdout.py`, `submission_writer.py`, cached features.
   - **Label every number** HOLDOUT-DTI (evaluator version, withheld positives,
     95% CI / spatial-block range) or ORGANIZER-CONFIRMED. Never write a
     projection as a score.
   - **Leakage canary.** Per-feature AUC > 0.90 = leakage until proven otherwise.
   - **Run card.** One JSON: hypothesis; mechanism; named non-fault mimic;
     holdout DTI + CI; correlation/overlap vs registry; raster sha256;
     validator output; submission name + ≤140 char note; verdict promote/negative.
   - **Budget.** Stop after 3 experiments or 2 hours. Do not pick the weekly
     slot; promotion to a real slot is a separate selector.
4. Study GEMSDOE32 `h33-h33-2-b2-…-zeros` (user-cited 0.2778) and ask whether
   a unique file can score higher. Official public leaderboard #1 (2026-10-21):
   xiaofanhu **0.4413**.
5. GitHub Pages with executive summary first, one-click TIF download.
6. Open a pull request and merge to main.
7. No hallucinations. Verify line by line. Flag irregularities.

## Official links

- https://www.drivendata.org/competitions/306/competition-doe-gems/
- https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/
- https://www.drivendata.org/competitions/306/competition-doe-gems/page/968/
- https://www.drivendata.org/competitions/306/competition-doe-gems/data/
- https://www.drivendata.org/competitions/306/competition-doe-gems/leaderboard/
- https://github.com/drivendataorg/gems-prize-reference-solution
- https://www.usgs.gov/data/geodawn-airborne-magnetic-and-radiometric-surveys-northwestern-great-basin-nevada-and
- https://gbcge.org/current-projects/ingenious/
- https://epsg.io/32611
- https://en.wikipedia.org/wiki/Tversky_index
- https://docs.nlr.gov/docs/fy26osti/96647.pdf

## Format that avoids `Predicted values must be in range [0, 1]`

Single-band float32, EPSG:32611, 3730×3292, transform
(100, 0, 243350, 0, −100, 4508550), all cells finite, values in [0, 1],
**no nodata tag** (zeros outside the footprint). NaN-outside is the
competition’s written format but has tripped the same portal error.

## Core values

Maximize P(Win). Own the Outcome.
