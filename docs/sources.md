# Verified Source Register (`2026-10-22` Audit)

Every URL and DOI below was directly inspected via page fetch or SHA-256 verified repository provenance in this session (`2026-10-22` UTC). No unverified claims or projections are presented as scores.

## A. Official Competition Rules, Metric & Leaderboard (`OFFICIAL-PAGE`)

| ID | Verified Fact | Official Source URL | Verification Status |
|---|---|---|---|
| **S-01** | Competition overview, timeline (`2026-12-03 23:59 UTC`), external data policy | <https://www.drivendata.org/competitions/306/competition-doe-gems/> | Verified `2026-10-22` |
| **S-02** | Exact Distance-Weighted Tversky Index formula ($\alpha=0.2, \beta=0.8$, triangular kernel $R=300\text{ m} = 3\text{ px}$ at $100\text{ m}$), worked example ($0.60265$), single-band `float32` `[0, 1]` `EPSG:32611` GeoTIFF specification, 19 `training_features.tif` layers | <https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/> | Verified `2026-10-22` (chunks 0–1) |
| **S-03** | Live Public Leaderboard (`2026-10-22`): `#1 xiaofanhu = 0.4413`, `#2 cfld__ = 0.3904`, `#3 Batik Shirt Brothers = 0.3894`, `#4 BrandenKMurray = 0.3713`, `#5 No Fault of Our Own = 0.3539`, `#10 newway = 0.3328`, `#20 nchuzhoy = 0.3262`, `#48 op01 = 0.2797` | <https://www.drivendata.org/competitions/306/competition-doe-gems/leaderboard/> | Verified `2026-10-22` (chunks 0–1) |
| **S-04** | Official DOE GEMS Prize Rules PDF ($300,000 total prize pool: Phase 1 $50k = 5 × $10k on withheld private labels; Phase 2 $250k = $100k/$70k/$40k/$25k/$15k on expert-updated full label set; **3 submissions per week**; code & documentation audit; mandatory GenAI narrative disclosure) | <https://docs.nlr.gov/docs/fy26osti/96647.pdf> | Verified `2026-10-22` (chunks 0, 1, 2) |
| **S-05** | DrivenData Official Reference Solution (`unet-mc-cv-reference-solution.ipynb`, Tversky $\alpha=0.2, \beta=0.8$) | <https://github.com/drivendataorg/gems-prize-reference-solution> | Verified (`commit aebe92f7`) |

## B. Official Public Geological, Geomechanical & Geothermal Datasets (`PUBLIC-DOI`)

| ID | Dataset / Physical Layer | Official Citation & Persistent Identifier | Usage in `GEMSDOE80` v2 |
|---|---|---|---|
| **D-01** | **USGS GeoDAWN Airborne Magnetic & Radiometric Surveys** (NW Great Basin, NV & OR; CC0 Public Domain) | Glen, J.M.G., & Earney, T.E. (2024), USGS Data Release, DOI: <https://doi.org/10.5066/P93LGLVQ> · <https://www.usgs.gov/data/geodawn-airborne-magnetic-and-radiometric-surveys-northwestern-great-basin-nevada-and> | RTP, TMI HG/VG, upward-continued 150 m TMI, airborne K/Th potassic adularia alteration & U/Th travertine anomalies, magnetite destruction |
| **D-02** | **USGS 3D Great Basin Crustal Stress, Slip & Dilation Tendency Release** | USGS ScienceBase Data Release, DOI: <https://doi.org/10.5066/P9YL58W6> (`sb_slip_tendency_in_footprint.csv`, `84,484` in-footprint segments) | Smoothed strictly at regional $\sigma = 25\text{ km}$ boundary condition ($S_{H\max}$ azimuth, $S_{H\max}$, $S_{h\min}$, Simpson $A_\phi$) for pixel-wise Cauchy $T_s, T_d$ |
| **D-03** | **USGS Great Basin Conductive & Convective Heat Flow Database** | DeAngelo, J., Burns, E.R., Gentry, E., et al. (2022), USGS Data Release, DOI: <https://doi.org/10.5066/P9BZPVUC> (`2,108` in-footprint boreholes) | Measured heat flow (`hf_meas`) and positive convective hydrothermal residuals (`hf_resid`) |
| **D-04** | **DOE INGENIOUS / Great Basin Center for Geothermal Energy Compilation** | Great Basin Center for Geothermal Energy (UNR), DOE Geothermal Data Repository Submission `1391`, DOI: <https://doi.org/10.15121/1881483> · <https://gbcge.org/current-projects/ingenious/> | Thermal springs/wells (`temp_c`, quartz/chalcedony/Na-K-Ca geothermometers) & Quaternary volcanic vents (`dist_known_fault_px` stripped & asserted absent) |
| **D-05** | **Coordinate Reference System (`WGS 84 / UTM zone 11N`)** | EPSG Registry `EPSG:32611` (<https://epsg.io/32611>) | Exact 100 m grid (`3730 × 3292`, affine `(100, 0, 243350, 0, -100, 4508550)`) |

## C. Peer-Reviewed Geomechanical & Optimization References

1. **Morris, A., Ferrill, D.A., & Henderson, D.B. (1996).** Slip-tendency analysis and fault reactivation. *Geology*, 24(3), 275–278. (<https://doi.org/10.1130/0091-7613(1996)024%3C0275:STAAFR%3E2.3.CO;2>)
2. **Ferrill, D.A., & Morris, A.P. (2003).** Dilational normal faults. *Journal of Structural Geology*, 25(2), 183–196. (<https://doi.org/10.1016/S0191-8141(02)00029-9>)
3. **Simpson, R.W. (1997).** Quantifying Anderson's fault types. *Journal of Geophysical Research: Solid Earth*, 102(B8), 17909–17919. (<https://doi.org/10.1029/97JB01274>)
4. **Pollard, D.D., & Aydin, A. (1988).** Progress in understanding jointing over the past century. *GSA Bulletin*, 100(8), 1181–1204. (<https://doi.org/10.1130/0016-7606(1988)100%3C1181:PIUJOT%3E2.3.CO;2>)
5. **Peacock, D.C.P., & Sanderson, D.J. (1991).** Displacements, segment linkage and relay ramps in normal fault zones. *Journal of Structural Geology*, 13(6), 721–733. (<https://doi.org/10.1016/0191-8141(91)90033-F>)
6. **Reid, A.B., Allsop, J.M., Granser, H., Millett, A.J., & Somerton, I.W. (1990).** Magnetic interpretation in three dimensions using Euler deconvolution. *Geophysics*, 55(1), 80–91. (<https://doi.org/10.1190/1.1442774>)
7. **Dinkelbach, W. (1967).** On nonlinear fractional programming. *Management Science*, 13(7), 492–498. (<https://doi.org/10.1287/mnsc.13.7.492>)
8. **Minoux, M. (1978).** Accelerated greedy algorithms for maximizing submodular set functions. *Optimization Techniques*, Lecture Notes in Control and Information Sciences, 7, 234–243. (<https://doi.org/10.1007/BFb0006528>)
