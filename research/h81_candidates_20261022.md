# H81 candidates — 2026-10-22 session (ranked, line-by-line vs repo state)

Current shipped primary: H80-DPVC hero, HOLDOUT-DTI (gems57-pooled-hide-v2,
K=22,663 withheld positives) = 0.16199; its own unshipped E80-2 ablation =
0.18937 (same evaluator/K). Promotion gate for this session: a new arm must
beat BOTH pooled controls re-evaluated on identical folds in
`scripts/run_h81_experiment.py`. No number here is ORGANIZER-CONFIRMED;
0.2778 (H33) and 0.4413 (leaderboard top) are owner/user-reported in the
standing prompt, not receipts in this repo.

| Rank | ID | Layers (repo channels) | Physical signature | Why it catches catalogue-missing faults | Non-fault mimic | Differs from implemented? | Cost |
|---|---|---|---|---|---|---|---|
| 1 | H81-1 conduit-gated meta | static cols 18,19,20 (MT clay-cap feeder, cond gradient/worm), 22,23 (magnetite-destruction demag), 48 (transtension strain) gating OOF meta surface | Multiplicative gate 0.60+0.40·mean-rank(conduit cols) on self-masked fracture-transfer surface; skeleton-submod 300 m, collar 2.40 px, 37,654 dots | Hydrothermal fluid upflow along young blind splays creates conductive smectite clay cap + demagnetized conduit before a Quaternary scarp exists; USGS/INGENIOUS is scarp-led | Playa evaporite / lacustrine clay conductivity — gate only re-weights meta dots (fracture support required) | DPVC used these channels as HGB inputs and a 28% blind sub-vent allocator tranche; never as a multiplicative gate on the meta surface with E80-2 emission | Low (tested this session, Exp. 1) |
| 2 | H81-2 pure meta transfer | OOF meta surface (59-feature self-masked stack) | E80-2 emission transferred to full catalogue: skeleton-submod, collar 2.40, 37,654 dots, no 28% vent-tranche dilution | Same fracture-transfer physics as DPVC core; removes allocator dilution that holdout penalized (0.18937 vs 0.16199) | Same as DPVC (Lahontan benches, TIGER cuts) | Emission-only transfer of an already-validated arm — not a new detector; shipped only if H81-1 fails and it re-beats hero | Low (control arm this session) |
| 3 | H81-3 radiometric sinter dipole gate | cols 25 (sinter leach edge), 27 (K/Th adularia), 28 (U/Th redox grad) | Cross-strike dipole gate on meta surface | Silica sinter/adularia halos flank blind upflow zones with no relief | Rhyolitic tuff K enrichment | Channels were HGB inputs only, never an explicit dipole gate | Low — NOT RUN (budget) |
| 4 | H81-4 updip-relay amplification | cols 38–42 (updip60 steps/tilt), 58 (vis_updip_relay_coupler) | Exponentiate meta surface by coupler rank in d_cat>10 px band only | Relay ramps between en-échelon range-front segments over deep fill | Lithologic basement contacts | Coupler was one of 59 inputs; never a distance-banded amplifier | Low — NOT RUN (budget) |
| 5 | H81-5 volcanic-vent conduit gate | cols 29,30 (GDR vent aniso near/far), 32 (volc vent conduit) | Gate meta dots within strike-anisotropic vent plumes | Quaternary vent alignments mark active extensional conduits absent from fault catalogue | Vent fields without through-going faults | Used inside DPVC sub-vent tranche only | Low — NOT RUN (budget) |

Validation (H81-1, this session): leakage canary on the conduit composite
alone (symmetric AUC ≤ 0.90 required), dual-seed (20/30) hide-and-recover
with whole-segment withholding + 3 px buffer, catalogue features derived
from visible faults only, visible faults masked pixel-exactly, pooled DTI
α=0.2 β=0.8 300 m triangular kernel, paired spatial-block CIs vs controls,
pre/post registry uniqueness gates (|ρ|≤0.90, ≤70% dots within 3 px),
15-check submission validator. Results: `evidence/h81_experiment.json`.

No external data beyond the SHA-256-pinned manifest was needed; all layers
above already exist in the cached stack. No weekly DrivenData slot spent.
