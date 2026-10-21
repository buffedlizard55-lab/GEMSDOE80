# Candidate geological hypotheses (GEMS Prize Challenge)

Status: proposed, not yet validated on the competition holdout — the training
labels are auth-gated on DrivenData (login required) and cannot be downloaded
from this sandbox. Ranked by expected DTI improvement vs implementation cost.
None of these is implemented anywhere else in this repo (repo was empty at
session start); none copies any previous GEMSDOE submission.

Competition context (official DrivenData pages):
- Task: pixel-wise fault probability, GeoDAWN region, northwestern Great Basin
  (Nevada + eastern California).
- Metric: distance-weighted Tversky index, alpha=0.2, beta=0.8, triangular
  kernel R=300 m (3 px at 100 m).
- Submission: single-band float32 GeoTIFF, EPSG:32611, 100 m, values in
  [0,1], exact CRS/shape/geotransform of the official sample_submission.tif.
  https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/

Test labels are faults deliberately EXCLUDED from the current public USGS /
INGENIOUS catalogue, so a catalogue-copy prior scores near zero on the hidden
test set. Every hypothesis below targets structure NOT already catalogued.

## H1 — Geothermal structural-propensity composite (unsupervised; implemented in scripts/generate_submission.py)

Layers: training_features.tif bands (surface conductivity / depth to conductive
base; dilatation rate, shear strain rate, 2nd strain invariant; isostatic
gravity anomaly + slope; magnetics RTP/TMI anomaly + vertical & horizontal
slope; earthquake density).
Signature: magnitude/edge response of deformation and potential-field bands,
combined with conductivity and seismicity — i.e. where the crust is actively
dilating/shearing, where gravity/magnetics mark density or magnetic contrasts,
where rocks are conductive, and where seismicity concentrates.
Why it catches uncatalogued faults: it responds to physical proxies of active
faulting and fluid pathways, not to membership in the USGS/INGENIOUS catalogue.
Non-fault mimic: lithological or stratigraphic boundaries can produce gravity/
magnetic/conductivity contrasts and localized strain without displacement.
Differs from repo implementations: nothing else exists in this repo; this is
the baseline submission surface.
Expected DTI improvement: moderate (baseline); cost: low (no new data, CPU).

## H2 — 1 m lidar DEM scarp / lineament density + curvature (needs DEM data)

Layers: 1 m DEM links (1m_DEM_links.csv, 3DEP lidar) downsampled to 100 m;
derived: curvature, slope, lineament density, hillshade edge response.
Signature: subtle scarps and linear topographic breaks at sub-catalogue scale.
Why it catches uncatalogued faults: many Quaternary faults in the region have
subtle surface expression below mapping resolution; 1 m lidar reveals scarps
the catalogue misses. INGENIOUS/USGS expert labels were partly derived from
such remote-sensing inspection.
Non-fault mimic: glacial/periglacial or alluvial geomorphic lineations.
Cost: medium (DEM download is public via USGS 3DEP/TNM; the link CSV is
auth-gated on DrivenData but 3DEP DEM tiles are publicly downloadable).
Expected DTI improvement: moderate-high.

## H3 — GeoDAWN radiometric alteration signature (K/Th/eU ratios)

Layers: GeoDAWN airborne radiometrics (public via USGS ScienceBase item
657e1d85d34e23d3533209f7).
Signature: potassium/thorium or uranium anomalies associated with
hydrothermal alteration along blind fault zones.
Why it catches uncatalogued faults: blind geothermal faults host alteration
halos with no mapped surface trace.
Non-fault mimic: sedimentary facies changes produce radiometric contrasts.
Cost: medium (ScienceBase download; public, no login).
Expected DTI improvement: low-moderate.

## H4 — Earthquake density × dilatation rate (seismogenic geothermal structures)

Layers: training_features earthquake density band + dilatation rate band.
Signature: co-located seismicity and present-day crustal dilatation.
Why it catches uncatalogued faults: active geothermal shear zones produce
microseismicity without Quaternary surface rupture mapped in the catalogue.
Non-fault mimic: volcanic or induced seismicity clusters not on mappable faults.
Cost: low (uses provided features only).
Expected DTI improvement: low-moderate.

## H5 — Distance-decay prior from catalogue faults (base model; weak for new faults)

Layers: USGS Quaternary Fault & Fold Database + INGENIOUS fault vectors.
Signature: exponential/distance decay from mapped traces.
Why it catches uncatalogued faults: it mostly doesn't — included as the
negative control / baseline that the other hypotheses must beat.
Non-fault mimic: n/a (it is the catalogue itself).
Cost: very low.
Expected DTI improvement: near zero on the held-out new-fault test set.

## Ranking (expected DTI gain | cost | data needed)

| Rank | Hypothesis | Expected DTI gain | Cost | Data needed (free, official) |
|------|------------|-------------------|------|------------------------------|
| 1 | H1 structural-propensity composite | moderate | low | training_features.tif (DrivenData, login) |
| 2 | H4 earthquake density × dilatation | low-moderate | low | training_features.tif (DrivenData, login) |
| 3 | H2 lidar scarp/lineament density | moderate-high | medium | 1m_DEM_links.csv (DrivenData login) + USGS 3DEP tiles (public) |
| 4 | H3 radiometric alteration | low-moderate | medium | USGS ScienceBase GeoDAWN (public) |
| 5 | H5 catalogue distance-decay | ~0 | very low | USGS Quaternary Faults (public) |

## Validation plan (parallel-run protocol)

1. Run H1 on the spatially-blocked holdout: withhold whole fault segments
   (buffered), derive every feature only from visible faults, score pooled DTI
   (alpha 0.2, beta 0.8, 300 m triangular kernel). Report HOLDOUT-DTI with
   95% CI.
2. Leakage canary: AUC of each feature band alone vs held-out labels; AUC >
   0.90 = leakage until proven otherwise.
3. Lane check: raster rank-correlation > 0.90 or > 70% of dots within 3 px of a
   registry raster => duplicate lane; log and stop.
4. Only promote a candidate that beats the current holdout best; promotion to
   a real submission slot is a separate step within the weekly cap.

Sources (official, for manual review):
- Competition description & metric: https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/
- Data page (login): https://www.drivendata.org/competitions/306/competition-doe-gems/data/
- Reference solution: https://github.com/drivendataorg/gems-prize-reference-solution
- GeoDAWN (USGS ScienceBase): https://www.usgs.gov/data/geodawn-airborne-magnetic-and-radiometric-surveys-northwestern-great-basin-nevada-and
- GeoDAWN item: https://www.sciencebase.gov/catalog/item/657e1d85d34e23d3533209f7
- INGENIOUS (GBCGE): https://gbcge.org/current-projects/ingenious/
- EPSG:32611: https://epsg.io/32611
- Tversky index: https://en.wikipedia.org/wiki/Tversky_index
