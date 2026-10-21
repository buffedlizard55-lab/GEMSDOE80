# Research queue — hypotheses, not performance claims

Reviewed against this checkout's broad structural composite and GEMSDOE32 H33 source
at `f62038b54dcb7767d4d997f80b033cdc3d42bef9`. Novelty across every GEMSDOE repository
has **not** been established. No method paragraph designates a single lane in the current
prompt: these are a queue only, not permission to run several lanes.

Ranking is qualitative scientific judgment about test priority / cost, NOT a DTI projection.
No candidate is viable for promotion until its data and the shared tools are available.

| Priority | Candidate / specific layers | Transform and missing-fault mechanism | Named non-fault mimic | Difference from reviewed implementation | Cost / status |
|---|---|---|---|---|---|
| 1 | Compensated gravity step under cover: `iso_grav_anom`, `depth_to_base_surf` | Paired opposite-sign second derivatives, with a basement-depth step oriented normal to the inferred boundary; test displacement beneath basin cover | Intrusive density contact | Requires signed cross-field step geometry, not average absolute anomalies or H33 straight-boundary coherence | Medium; core mirror listed, shared cache missing |
| 2 | Asymmetric incision across a candidate line: raw DEM and derived catchments | Compare drainage-area-normalized channel steepness on opposing flanks; persistent asymmetry could indicate concealed differential uplift | Resistant lithologic boundary | Not H33-B's detrended-valley collinearity; requires real flow-routing elevations, NOT detrended elevations | High; blocked on verified raw DEM acquisition |
| 3 | Conductive damage-zone double edge: `cond_surf`, `depth_to_base_surf` | Two parallel, opposite-gradient conductivity edges enclosing a narrow zone, corroborated by a depth step; may indicate unmapped fractured rock | Clay-rich sediment ribbon | Detects a bounded zone rather than H33-D single coherent gradient or scalar conductivity | Medium; shared core cache missing |
| 4 | Radiometric alteration asymmetry: GeoDAWN K/Th/U and `rtp` | Compare cross-strike K/Th residuals after lithology stratification with a magnetic low; test alteration along buried fluid pathways | Weathering/soil redistribution | Signed cross-strike contrast instead of repo scalar average or H33 multi-physics orientation consensus | High; ratios/lithology provenance unverified |

None is described as proven to catch faults: those are falsifiable mechanisms, not findings.
Potential-field inversion is nonunique; confidence must come from spatial transfer tests and
negative controls, not a visually plausible anomaly. Predicting fault traces is not direct
prediction of geothermal vents or temperature.

## External sources and obtainability

* Raw DEM: USGS 3DEP / The National Map, https://apps.nationalmap.gov/downloader/ .
* Radiometrics: USGS GeoDAWN release https://doi.org/10.5066/P93LGLVQ .
* Lithology: USGS State Geologic Map Compilation https://doi.org/10.3133/ds1052 .

These official source URLs are leads, **not accessed/availability-verified datasets in this
session**. Network access is limited to GitHub and package registries. External-dependent
candidates remain blocked rather than being proposed as immediately viable.

## Top-candidate preflight

No experiment run. Before generating a surface: identify the assigned lane, restore the
shared cached stack, `evaluate_holdout.py`, `submission_writer.py`, registry and current best;
withhold entire fault segments and buffers, rebuild catalogue features from visible faults,
mask visible faults pixel-exactly, run every feature's AUC canary with positives AND negatives,
and evaluate pooled HOLDOUT-DTI with segment/block-bootstrap CI. Compare the surface before
placement, then compare final dots. Stop on correlation >0.90 OR overlap >70% within 3 px.
Do not substitute a synthetic trial, pixel split, mean of fold DTIs, or one-value bootstrap.
