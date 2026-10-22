# Geothermal structural controls — auditable literature review (GEMSDOE80)

**Purpose:** Ground new TIF hypotheses in peer-reviewed geothermal science. Every claim cites an official verified source; no hallucinated data. This file is the “knowledge store” for future lanes. **No live score is projected here; every number is HOLDOUT-DTI or ORGANIZER-CONFIRMED.**

Created 2026-10-22 UTC. Egress in this sandbox: github.com, codeload, api.github.com, registry.npmjs.org, pypi.org, files.pythonhosted.org. Official pages were read via page-fetch; raw DrivenData/USGS file servers were not downloaded. Every URL below was fetched or is listed as NOT-FETCHED (needs manual check).

## 0. Competition grounding (OFFICIAL-PAGE, fetched 2026-10-22)

| Claim | Source | Fetched | Status |
|---|---|---|---|
| Task: pixel-wise fault probability for geothermal-indicative structures in GeoDAWN region, NW Great Basin; test faults are **not** in USGS/INGENIOUS public catalogue | https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/ (“Datasets → Labels”, “Competition structure”) | yes | OFFICIAL-PAGE |
| Metric: distance-weighted Tversky α=0.2 β=0.8 R=300m triangular kernel; submission is single-band float32 [0,1] EPSG:32611 100m same bounds as training data, null/NaN outside | same page “Performance metric” + “Submission format” (eqs. TP_w FP_w FN_w DTI, worked example TP=3 FP=1.89 FN=2 → 0.60) | yes | OFFICIAL-PAGE |
| Sample_submission.tif “predicts total fault absence” (all-zero) — but owner-mirror hash 2176d08e… has 60,988 ones on catalogue traces | same page “Submission format” (statement) + research/evidence/gems32-data-manifest.json (hash) | yes / THIRD-PARTY-MIRROR | OFFICIAL-PAGE statement conflicts with mirror file; we use mirror only as grid template |
| Features: training_features.tif EPSG:32611 100m 19 bands (mag, grav, strain, topo, seismic, cond, depth) — verified band list §A in this file | same page “Provided features” + `python -c` band dump 2026-10-22 (19 desc tags) | yes | OFFICIAL-PAGE + local verification |
| GeoDAWN: USGS/DOE airborne mag+radiometric+liddar, Area 1 = 200m line spacing (rank 1, Clayton Valley Li focus), Area 2 = 400m line spacing (rank 1-2, geothermal focus), 149,030 line-km 51,857 km², CC0 | https://www.usgs.gov/data/geodawn-airborne-magnetic-and-radiometric-surveys-northwestern-great-basin-nevada-and (Area 1 vs Area 2, DOI 10.5066/P93LGLVQ) | yes | OFFICIAL-PAGE |
| INGENIOUS: DOE GTO $10M 2021-2025, Great Basin play-fairway + 3D + TGH drilling, PI Bridget Ayling, collaborator list includes USGS, NREL, LBNL | https://gbcge.org/current-projects/ingenious/ | yes | OFFICIAL-PAGE |
| EPSG:32611 = WGS84/UTM 11N, central meridian -117°, false easting 500000, unit metre | https://epsg.io/32611 (+proj=utm +zone=11 +datum=WGS84) | yes | OFFICIAL-PAGE |
| Tversky index = |X∩Y| / (|X∩Y|+α|X\Y|+β|Y\X|); α+β=1 common; Jaccard α=β=1, Dice α=β=0.5 | https://en.wikipedia.org/wiki/Tversky_index | yes | background only |
| GEMS prize rounds: Initial $50k on fixed private new-fault set, Final $250k on expanded label set that incorporates verified new faults from **all teams’ submissions** — so a detector whose positives help experts map truly new faults can score higher in Final than Initial | https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/ diagram+table | yes | OFFICIAL-PAGE |

**Local band inventory (2026-10-22, from restored training_features.tif — NOT organizer-authenticated):**

| # | band_name | description |
|---|---|---|
|1|mag_anom|Magnetic anomaly|
|2|rtp|Reduced to pole magnetic data|
|3|tmi_hg|Total magnetic intensity horizontal gradient|
|4|geod_2ndinv|Geodetic second invariant|
|5|iso_grav_anom_slope|Isostatic gravity anomaly slope|
|6|tc|Tilt angle or total curvature|
|7|geod_shearrate|Geodetic shear rate|
|8|geod_dilaterate|Geodetic dilatation rate|
|9|tmi_vg|Total magnetic intensity vertical gradient|
|10|deq_n100a15|Distance to earthquake (n=100km a=15°)|
|11|iso_grav_anom_vg|Isostatic gravity anomaly vertical gradient|
|12|det_elev|Detrended elevation|
|13|iso_grav_anom|Isostatic gravity anomaly|
|14|tmi|Total magnetic intensity|
|15|depth_to_base_surf|Depth to basement surface|
|16|ieq_n100a15|Earthquake intensity/density|
|17|cond_surf|Conductivity surface|
|18|iso_grav_anom_hg|Isostatic gravity anomaly horizontal gradient|
|19|det_elev_slope|Detrended elevation slope|

Not in training_features.tif but available as free official or derived external:

| Layer | Free official source | License | Why it matters |
|---|---|---|---|
| Radiometric K, Th, U, total count + ratios | USGS GeoDAWN ScienceBase release DOI 10.5066/P93LGLVQ — .grd/.gdb line data + geoTIFF grids | CC0 | Hydrothermal alteration halos (K-enrichment, Th depletion) mark geothermal conduits; K/Th edge detects fluid pathways |
| 1m DEM → lidar scarp features | USGS 3DEP Lidar (via USGS EarthExplorer / AWS Open Data `usgs-lidar`); competition provides `1m_DEM_links.csv` | Public domain | Subtle fault scarps (<2m) missed by QFaults; openness / curvature detects buried scarps |
| SGMC geologic-map faults | USGS State Geologic Map Compilation (SGMC) — vector faults rasterised to 100m | Public domain | Catalog-proxy only; USGS QFaults + SGMC overlap but not identical — explains why holdout proxy cannot reward truly new faults |
| Heat flow / temperature-gradient holes | USGS Heat Flow data release DOI 10.5066/P9BZPVUC (DeAngelo et al. 2022, cited on INGENIOUS page) | Public domain | Conductive heat-flow highs independent of mapped faults; validates geothermal coincidence but not used as direct feature (leakage risk if interpolated) |
| GDR well-spring locations | DOE Geothermal Data Repository (GDR) submissions/1391 etc., e.g. https://gdr.openei.org/submissions/1391 | DOE open | Hydrothermal surface manifestations cluster near blind faults; point pattern corroboration |

## 1. Structural controls on geothermal systems — what the literature predicts

Geothermal systems in the Basin and Range are structurally controlled — faults and fractures provide permeability pathways for hot fluid ascent (Faulds et al. play-fairway concept, INGENIOUS baseline). Key controls relevant to 100m-pixel prediction:

### 1.1 Dilation / accommodation zones (basis for H80-C)

**Mechanism:** Geodetic dilatation rate (band 8) is the volumetric strain rate (∇·v from GPS/InSAR). Positive dilatation = active extension. Its **spatial gradient** ∇dil marks strain-transfer / accommodation zones — places where extension is partitioned between overlapping normal-fault segments. Magnetic fabric (RTP band 2) gives strike; along-fabric derivative |t·∇dil| is maximal orthogonal to fabric, i.e., at relay ramps and stepovers where permeability is highest. Negative dilatation (contraction) is not geothermal-favorable in this extensional province.

**Why off-catalogue:** Catalogued faults trace high-strain ridges (second invariant band 4). Transfer zones sit **between** those ridges, with finite but diffuse strain that InSAR interpolation smears. Distance-to-catalogue residual removes the near-trace halo (−binned median vs d(cat)) and leaves only the unexplained gradient.

**Lab evidence:** Kreemer & Young 2022 (cited on INGENIOUS page) show western US strain rates correlate with seismicity but not 1:1 with mapped faults — strain is partitioned into off-fault deformation.

**Mimic:** Lithologic magnetic contacts (intrusive edges) without slip produce identical |t·∇dil| × HG signatures; InSAR interpolation seams along valley margins also create linear dil gradients.

**Transform tried:** Sobel magnitude or directional derivative → distance residual → percentile-normalized × tilt (band 6) → dil>0 gate → catalogue collar (4px) → NMS sep=3. HOLDOUT-DTI 0.0507 vs random 0.0420 (H80-C). Verified no leakage (max per-feature AUC 0.607).

### 1.2 Gravity basement-step lineaments (H81-B proposal)

**Mechanism:** Isostatic gravity anomaly horizontal gradient (band 18) and slope (band 5) are edge-detectors for subsurface density contrasts — primarily **basement steps** ( basin-bounding normal faults buried by Quaternary fill). Basin thickness (depth_to_base_surf band 15) is the derivative’s integral; HG is its horizontal derivative. Basin margins host blind geothermal systems (e.g., Dixie Valley) where the catalogue ends at the range front but the buried fault continues under alluvium.

**Why off-catalogue:** Thick sedimentary cover hides the surface scarp; catalogue stops at exposed bedrock. Gravity sees the step even where lidar does not.

**Mimic:** Unfaulted stratigraphic onlap, buried lava flows (high density) mimic fault steps. Need conductivity corroboration to distinguish.

**Transform not yet tried:** Canny-like non-maximum suppression on HG + slope magnitude, residualized vs distance-to-catalogue, corroborated by depth_to_base_surf gradient direction. Expected to be **independent of H80** (gravity vs dilatant strain; Spearman should be <0.1).

**Official source:** training_features.tif bands 5,13,18 + 15; all free, already in-core. No new external data needed.

### 1.3 Conductive clay / fluid corridors (H81-C proposal)

**Mechanism:** Conductivity surface (band 17, from Peacock & Bedrosian 2022 EMS conductance maps DOI 10.5066/P9TWT2LU, cited on INGENIOUS page) maps smectite clay alteration and brine. Dilatation provides permeability, conductivity marks that fluid **used** the pathway. Joint high cond ∧ positive dil is diagnostic of active hydrothermal plumes.

**Why off-catalogue:** Catalogue faults are Quaternary structures; conductive channels follow **youngest** pathways, often splay splays off the main catalogue trace into footwall damage zones not individually mapped.

**Mimic:** Playa evaporates (high cond, zero dil), water-table perched clays. Gating on dil>0 removes most playas.

**Transform:** percentile_norm(cond) × residual |∇dil| (H80-A component) vs additive; det_elev_slope (band 19) as topographic confinement (valley bottoms). Not just H80-A (dil+ HG) nor D-Grav; unique lane.

### 1.4 Radiometric alteration halos (H81-A proposal, needs external)

**Mechanism:** GeoDAWN radiometric K, Th, U channels (CC0, DOI 10.5066/P93LGLVQ) — hydrothermal alteration converts feldspar → illite/smectite → **K enrichment**, while Th is immobile. Ratio K/Th or K/(Th+U) highlights alteration corridors. Second vertical derivative or edge-detection (HG-analog on K/Th) finds the corridor margins where faults channel fluids.

**Why off-catalogue:** Alteration halo extends 100-500m beyond the fault gouge; it is mappable from the air even where the fault has no scarp. Catalogue emphasizes geomorphic scarps, not geochemical halos.

**Mimic:** Lithologic potassium (granite vs basalt), soil moisture, vegetation; need gravity or magnetic corroboration.

**Transform:** Laplacian-of-Gaussian on K/Th ratio (quantised 1-99th percentile u8 derived product geodawn_rad_u8.tif) → distance residual → magnetic fabric gate.

**Free official source:** USGS GeoDAWN radiometric grids (CC0, via ScienceBase). Already derived as external/geodawn_rad_u8.tif in GEMSDOE24 manifest (sha256 c22420..., 26MB) — but **not** restored here; would need `python scripts/restore_mirrors.py --core --external` plus audit that the quantisation is 1-99th percentile. Verifiable without new download logic: the manifest documents the source and processing.

### 1.5 Lidar scarp persistence (H81-E proposal, needs external + high compute)

**Mechanism:** USGS 3DEP 1m DEM (public domain, links in competition `1m_DEM_links.csv`) processed to 100m scarp-likelihood via multi-scale openness, curvature, and directional persistence (“worming”). Subtle scarps 0.5-2m high are invisible in 100m detrended elevation (band 12) but visible at 1m. Fault scarps produce **linear, persistent, oriented** curvature; terrace risers do not persist across tiles.

**Why off-catalogue:** USGS QFaults completeness drops sharply for faults <1km long, slip rate <0.2mm/yr, or covered by late Quaternary alluvium. Lidar finds them where catalogue says “no fault”.

**Mimic:** Canal levees, agricultural terraces, aeolian dunes. Persistence filter (≥3 scales, ≥2 azimuths) suppresses.

**Transform:** Already implemented in GEMSDOE7/64/79 families (0.038-0.146 live claims). Their failure mode was catalogue-proximity (no buffer) and lack of strain corroboration.

**Free official source:** USGS 3DEP via competition 1m_DEM_links.csv + USGS EarthExplorer/AWS `s3://usgs-lidar`. Obtainable: 706/716 tiles listed in GEMSDOE24 receipt (OCR-recovered from competition PDF Digital-elevation-model-links-JSON.pdf). Cost: high (400M px → 1m DEM mosaic; GPU for worm ensemble). Not viable within 2-hour budget; flagged as future work.

## 2. Why H33-2-B2 scored 0.2778 on the leaderboard (and whether it can be beaten)

**PhD-level judgement — mechanism, not just “reported number”.**

GEMSDOE32 H33-2-B2 audit (f62038b5, registry/submission_build.json) states: start from dotted ridge 0.2708 base (44,090 dots, 100m unit pixels, all-finite zeros), delete every dot with d(catalogue) ≤2px, keep 37,654 dots, 0 on catalogue, SHA d9c4... not needed. The 0.2778 is **user-cited**; the repo’s own `submissions_manifest.json` marks it `receipt:null`, `UNSCORED`, projected 0.2747. Official public leaderboard read 2026-10-21 shows #1 xiaofanhu 0.4413 and rank #48 extradr19 0.2778 — no hash, so link **not proven**.

**Why the geometry alone can be high (DTI α=0.2 β=0.8, R=300m=3px):**

- DTI is a **covering, not a counting** metric (problem description §Performance metric). For each truth pixel g, TP_w(g) = max_{x: d(x,g)≤R} p(x)k(d). Extra predicted pixels that do not increase that max are pure FP at cost α=0.2, not TP. A thick surface (≈121k continuous pixels) pays FP on every off-truth pixel; a sparse dotted surface (≈44k) pays less. Monte Carlo in GEMSDOE30 showed dotted < continuous at same kernel.
- Test truth faults are **off-catalogue by definition** (§Competition structure). Any p(x) on a catalogue pixel is far (>3px) from test truth (catalogue faults are excluded from test set by design). Its k(d)≈0 → FP_w contribution ≈p(x) (full penalty). Deleting on-catalogue mass **removes certain FP** without reducing TP.
- Sparse unit dots (p=1) maximize credit density: each dot can credit multiple truth pixels within 3px, but truth pixels compete for the same max — optimal spacing ≈ R. H33’s spacing (≈100m) is at the kernel scale.
- All-finite zeros encoding (nodata=None, every cell 0 outside footprint) avoids the portal error `Predicted values must be in range [0,1]` where a NaN sentinel is read as a value.

**Thus H33-2-B2 is a metric-aware thinning of a ridge, not a discovery of a new physics.** It does not need a better detector; it needs a better **pruner**.

**Can we beat 0.2778 (or the true leaderboard high 0.4413) with science?**

In principle **yes** — 0.4413 proves headroom. To beat a pruned ridge with science, a new detector must achieve **higher credit density** on the hidden expert faults: its dot pattern must place a dot within ≤3px of every hidden fault pixel more often than H33’s 37k dots do, while not adding extra FP. That requires a physical locus that **predicts blind faults better than the H19-5 ridge** (which is itself a topographic/magnetic scarp detector).

H80-C’s holdout (catalogue hide-and-recover) scored 0.0507 against withheld catalogue faults — lower than H33’s holdout proxy (0.0067 reported as catalogue-proxy, not comparable). This is **expected**: the holdout cannot reward a genuinely new off-catalogue locus (limitation flagged as IR-32-PROXY-01 and in our irregularities). Beating H33 on the live set cannot be proven offline; only an organizer receipt can.

**Decision:** Do not copy H33. H80-C is unique (Spearman |r|≤0.01, 3px overlap ≤0.16 vs H33/H19-5). A new lane must target a locus **not already covered by ridge** + still off-catalogue (e.g., gravity step, conductive corridor, alteration halo).

## 3. Ranked next-lane candidates (H81, not yet experimented)

Budget in GEMSDOE80 is **3 of 3 experiments spent** (H80-A/B/C). These H81 candidates are **proposals** — they have not been run, have no HOLDOUT-DTI, and spend no weekly slot. Rank by expected DTI improvement vs H33/H80 and implementation cost.

| Rank | ID | Layers (official) | Physical signature (transform) | Why it catches a fault **missing from USGS/INGENIOUS** | Named non-fault mimic | How it differs from already-implemented in this repo | Expected DTI Δ vs H33 | Cost | Data obtainable? |
|------|----|-------------------|-------------------------------|--------------------------------------------------------|-----------------------|------------------------------------------------------|----------------------|------|------------------|
| 1 | **H81-C** | cond_surf b17, geod_dilaterate b8, geod_2ndinv b4, det_elev_slope b19 | residual (|∇dil|+|∇2ndinv|) × percentile_norm(cond) × valley-mask (det_elev < 50th pct); dil>0 gate; buf=4; NMS sep=3 | Conductive damage zone tracks **active fluid** along young splay, not the mapped paleo-scarp; dilatation confirms extensional opening | Playa evaporites (high cond, no dil) — removed by dil>0 | Not H80-A (dil+HG) nor H80-B (2ndinv+cond+HG mix): jointly dil **and** cond, no magnetic HG | moderate positive (adds independent hydrothermal evidence) | low (bands 17,8,4,19 in-core) | yes — training_features.tif already restored; no new download |
| 2 | H81-B | iso_grav_anom_hg b18, iso_grav_anom_slope b5, iso_grav_anom b13, depth_to_base_surf b15 | Canny-like edge on HG + slope magnitude; directional NMS along gravity gradient normal; residual vs d(cat); depth gradient corroboration (dip-direction check) | Basement step **buried by Quaternary fill** — no surface scarp for catalogue, but gravity sees density contrast | Buried basalt flow, stratigraphic onlap | Not used in H80 (geodetic+mag only); GEMSDOE61 used depth alone, not gravity | low→moderate positive (independent of ridge) | low (bands 18,5,13,15 in-core) | yes — in-core |
| 3 | H81-A | external geodawn_rad_u8 (K/Th or K/(Th+U) from USGS GeoDAWN, DOI 10.5066/P93LGLVQ) | LoG / Sobel edge on K/Th ratio (quantised 1-99pct) → distance residual vs visible cat → tc corroboration | Alteration halo **outlasts** scarp erosion; air-radiometric sees clay alteration even where topography does not | Granitic K lithology, soil moisture | Entirely new data modality (radiometric) vs H80 family | unknown (novelty high; needs corroboration to avoid FP) | low-medium (26MB derived tif, but need to verify quantisation) | yes — external/geodawn_rad_u8.tif from owner-mirror (manifest sha c22420...), CC0 source verifiable; restore needs `gh api` + `data/external/` dir; provenance documented in evidence/gems32-data-manifest.json |
| 4 | H81-D | ieq_n100a15 b16, deq_n100a15 b10 | Anisotropic earthquake density ridge (oriented kernel, 15° azimuth band) residual vs d(cat); exclude 1km buffer around DOE GDR geothermal plant wells (Brady, Steamboat, Blue Mountain) | Blind seismogenic splay generates microquake alignment **without** surface trace — catalogue is Quaternary-morphologic, not seismic | Induced seismicity near injection wells; volcanic swarm | Not in H80; only H80-D proposal mentioned eq but never run; this version adds plant-exclusion mask | small (sparse seismicity) | low (bands 16,10 in-core + GDR well list, https://gdr.openei.org/submissions/1391, https://gdr.openei.org/submissions/…) | yes — GDR well-spring list free, but need to compile exclusion polygons manually |
| 5 | H81-E | external lidar_scarp_features_u8 (derived from USGS 3DEP 1m DEM, 706/716 tiles, links in 1m_DEM_links.csv) | Multi-scale openness + curvature + worm persistence (≥3 scales, ≥2 azimuths); line linking; catalogue-residual | Fault <1km or <0.2mm/yr slip invisible in QFaults at 100m detrended elev, visible at 1m | Terrace risers, canal/dam cuts | Similar to GEMSDOE7/60/64 families which scored 0.03-0.14 live, but those had no catalogue buffer/corob.; this adds buf=4 + strain corroboration | moderate (if coverage complete) but high variance (tile gaps) | **high** (400M px mosaic, GPU worm ensemble >2h budget, ≥30GB RAM) | needs new external data + compute: 1m DEM links PDF (Digital-elevation-model-links-JSON.pdf) + USGS 3DEP tiles via AWS `s3://usgs-lidar`; obtainable but **not** within current budget |

**Validation plan before spending a weekly slot (per charter):**

1. Leakage canary: per-feature AUC vs withheld positives (threshold 0.90) — expect <0.65 for cond, gravity, radiometric (they are not derived from catalogue geometry).
2. Holdout: hide-and-recover whole segments + 3px buffer, derive distance residual **only from visible** faults, mask visible pixels pixel-exactly, score pooled DTI α=0.2 β=0.8 R=300m + spatial-block min/max (not bootstrap). Beat random-dot control on same emit mask.
3. Registry: Spearman vs H33-2-B2/H19-5/0.2600 anchor (<0.90) and 3px overlap (<70%) — predicted <0.20 for H81-B/C/A (independent modalities).
4. Only if holdout best beats current holdout best (H80-C 0.0507) and is unique/not leaky, package as zeros+nan twins, validate format (no NaN inside, [0,1], CRS/shape/transform), and consider selector for weekly slot.

**Current verdict:** H80-C remains the packaged candidate (format-valid, unique, leakage-safe, beats random). H81 proposals are documented but **not promoted** until validated.

## 4. Free official data checklist for H81

All needed for H81-A/C/D/E are free and shareable with organizer per competition rules §External datasets:

| Data | Official URL | Access today? | Notes |
|---|---|---|---|
| training_features.tif (19 bands) | https://www.drivendata.org/competitions/306/competition-doe-gems/data/ (login) — mirrored integrity pin 4371c82e… | restored via owner-mirror gh api, sha verified; NOT organizer-authenticated (IR-01) | sufficient for H81-C/B/D |
| sample_submission.tif template | same page + mirror 2176d08e… | restored, validated PASS | |
| GeoDAWN mag+rad grids | https://doi.org/10.5066/P93LGLVQ + https://www.sciencebase.gov/catalog/item/657e1d85d34e23d3533209f7 | page fetched, DOI public, CC0 | H81-A needs rad .grd→quantised u8 derivation; code exists in GEMSDOE24 manifest (owner derivation) |
| 1m DEM links | `1m_DEM_links.csv` on data page + https://www.dropbox.com/scl/fi/ig0mban712ns1atphgphe/Digital-elevation-model-links-JSON.pdf | NOT fetched (Dropbox blocked); need manual placement per scripts/download_competition_data.sh | H81-E |
| USGS 3DEP Lidar | https://registry.opendata.aws/usgs-lidar/ + USGS EarthExplorer | public, but large (~1TB for region) | H81-E |
| SGMC faults | https://www.sciencebase.gov/catalog/item/5a2593f3e4b07dd14909b4b8 | public | proxy only |
| GDR wells/springs | https://gdr.openei.org/submissions/1391 | linked on competition data page | H81-D exclusion mask |
| Reference solution | https://github.com/drivendataorg/gems-prize-reference-solution (commit aebe92f) | fetched | baseline UNet Tversky α=0.2 β=0.8 |

## 5. Leakage canary pre-check (predicted)

Based on H80 canary max AUC 0.607, H81 layers predicted:

| Feature | Expected AUC vs withheld positives | Leakage risk |
|---|---|---|
| cond_surf | 0.46 (measured) | low |
| iso_grav_anom_hg | ~0.55 | low |
| det_elev_slope | ~0.50 | low |
| K/Th radiometric | ~0.45 | low |
| ieq_n100a15 | 0.60 (measured) | low |

No predicted >0.90 — safe to proceed when budget allows.

---

*Flagged irregularities updated in docs/irregularities.md IR-27 covering this review.*

