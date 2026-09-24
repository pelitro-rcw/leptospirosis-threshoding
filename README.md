# Analytical Pipeline - Leptospirosis
**Clean-session pipe safety:** all executable scripts use the native R `|>` pipe (R >= 4.1). No stage depends on `%>%` being attached by `dplyr` or `magrittr`.

## Purpose
This archive is a reproducible NCR and Philippine regional leptospirosis surveillance pipeline adapted from the supplied dengue analytical pipeline. The leptospirosis outcome is **`LC_DOH` on the `Regional Data` sheet**. The workbook's `QC Data` outcome is dengue and is not used as the NCR leptospirosis series.

## Main deliverables

### Stage 1: NCR eta and STA/LTA derivation
`Stage1_eta_threshold_derivation_NCR.R` derives the Constant TA activation/deactivation thresholds and selects Constant/Continuous STA/LTA windows. Rolling windows are calendar-contiguity safe and cannot bridge the missing 2020-2021 interval.

### Stage 2: NCR outbreak-threshold drift QC
`Stage2_outbreak_threshold_drift_NCR_Lepto.R` generates **one single-panel figure restricted to 2023-2025** plus supporting threshold tables. Earlier observed years remain available only as legitimate donor years for threshold estimation. Leptospirosis observations are absent for 2020-2021, so a literal pandemic-included counterfactual is not estimated. Instead the QC compares observed donors inside each calendar y-5..y-1 window with a gap-aware last-five-observed-years specification. No missing pandemic values are imputed. The displayed weekly leptospirosis y-axis uses the original data-driven scale with automatic headroom so peaks and threshold curves are not cropped. Figure A1 draws only **Leptospirosis Cases**, **Outbreak Threshold (mean + 2SD)**, and **Alarm Threshold (mean + SD)** in its legend.

### Stage A: NCR six-detector analysis
`StageA_NCR_Lepto_Analysis.R` evaluates exactly these six dashboard detectors:
1. Constant Transmission Acceleration
2. Continuous Transmission Acceleration
3. Outbreak Threshold
4. Alarm Threshold
5. WHO 75th Percentile Threshold
6. WHO 90th Percentile Threshold

Key exports include:
- `FigureA2_panel_a_DominanceMatrix.{pdf,png}`
- `FigureA2_panel_b_Scatter.{pdf,png}`
- `FigureA3_Comparison_2019_2024_Constant_TA.{pdf,png}`
- `FigureA3A_Constant_Transmission_Acceleration_<year>.{pdf,png}` for each evaluable year
- `FigureA3_HeadToHead_panel_a_BurdenAccuracy.{pdf,png}`
- `FigureA3_HeadToHead_panel_b_Timeliness.{pdf,png}`
- `FigureA3_HeadToHead_panel_c_FalseAlarms.{pdf,png}`

The 2019 and 2024 peaks are derived from NCR `LC_DOH`; dengue's hard-coded 2024 W34/W47 override is not carried into leptospirosis. All weekly NCR timing plots use the original data-driven weekly-case scale with headroom. The secondary **R(t)** axis is also uncapped and data-driven, but uses only about 4-5 well-spaced tick labels so large ratios do not compress the right-axis numbers. Wider left/right plot margins are reserved to prevent either y-axis from being cropped.

### Stage B: Regional analysis
`StageB_Regional_Lepto_Analysis.R` evaluates the same six detectors across **all 17 canonical Philippine regions**. Regional eligibility is now coverage-based rather than incidence-based: a region-year is evaluable with at least 40 observed leptospirosis weeks, and a region requires at least five evaluable years. The source workbook satisfies this rule for all 17 regions. Low-incidence regions are therefore retained rather than silently excluded by a peak-case cutoff.

Requested map/dominance exports include:
- `FigureB5_panel_a_DetectorMap.{pdf,png}`
- `FigureB5_panel_b_DominanceProbability.{pdf,png}`
- `FigureB5_Combined_RegionalDominance.{pdf,png}`

Dominance probability is calculated over all six requested detectors. The detector map displays the observed composite winner among Constant TA, Continuous TA, Outbreak Threshold, Alarm Threshold, WHO 75th, and WHO 90th. Before export, Stage B hard-checks that the analytical table, dissolved geometry, map join, and dominance-probability table contain exactly the same 17 canonical regions. The bundled `PackedSpatVector` is explicitly unpacked before `sf` conversion. The earlier `consensus_pass`/`sig_star_dominance` schema failure is also removed from the executable map path and the required regional schema is validated explicitly.

## Four rainfall threshold products
`StageC_RF_Threshold_Derivation.R` derives **four separate environmental products** from `RF_HDX`:

1. **RF Product 1:** antecedent cumulative rainfall in **mm**, preceding the first Constant TA activation.
2. **RF Product 2:** antecedent cumulative rainfall in **mm**, preceding the first Outbreak Threshold crossing.
3. **RF Product 3:** antecedent **rainfall STA/LTA R(t)**, preceding the first Constant TA activation.
4. **RF Product 4:** antecedent **rainfall STA/LTA R(t)**, preceding the first Outbreak Threshold crossing.

For Products 1 and 2, the environmental feature is cumulative rainfall during **t-4..t-1** at the disease-trigger week. For Products 3 and 4, rainfall is converted to its own STA/LTA ratio series and the maximum antecedent rainfall R(t) during **t-4..t-1** is used. The disease-trigger week is excluded. A1/A2 remains the independent epidemiological truth standard: an RF-supported disease trigger is a **true alarm iff the trigger week lies in A1 OR A2**; an RF-supported trigger outside A1 union A2 is a **false alarm**.

The rainfall R(t) is **empirically derived from rainfall itself**. Candidate rainfall STA/LTA windows (2/13, 3/20, 4/26; two-week guard) are evaluated independently for Products 3 and 4; the rainfall windows and thresholds are not copied from the leptospirosis case detector. All four operational RF products are now defined prospectively in the **upper-tail direction (`RF >= threshold`)**. Low rainfall can therefore never be labelled a stronger RF trigger merely because it falls farther below an ROC cutoff.

The fixed point estimate is derived with an A1/A2-consistent false-alarm-aware threshold search. Candidate upper-tail thresholds qualify the first Constant TA or Outbreak Threshold anchor using the antecedent RF feature and are evaluated against the established A1/A2 true/false classification for sensitivity, specificity, balanced accuracy, PPV, trigger rate, false-alarm rate, and false alarms per evaluable unit. Thresholds within 0.01 of the best balanced accuracy are ranked by lower false-alarm burden, lower false-alarm rate, higher PPV, sensitivity and specificity, then by the more conservative higher threshold. `RF_Operational_Threshold_Candidate_Grid.csv` exports every candidate and identifies the selected point estimate for auditability.

Thresholds are **scope-specific but fixed over time**. NCR Products 1-4 are derived only from NCR history and then frozen across NCR evaluation years; their empirical stability bounds come from leave-one-year-out (LOYO) re-derivation. Regional Products 1-4 are derived from the pooled 17-region history and frozen across all 17 regional reruns; their bounds come from leave-one-region-out (LORO) re-derivation; for RF STA/LTA products the rainfall STA/LTA window is re-selected inside every held-out fold. Thus an NCR plot is evaluated against an NCR-derived cutoff rather than a pooled Regional cutoff.

For every product the ordered triplet is `lower <= point <= upper`. These are empirical **threshold-stability bounds**, not separately re-estimated annual thresholds. Operational interpretation is fixed and monotonic: **RF < lower = no activation; lower <= RF < point = precautionary activation; point <= RF < upper = formal RF trigger; RF >= upper = strong RF trigger**. The adopted point estimate remains the primary threshold for headline RF-gated performance. Observed RF values may vary by year and region; the relevant fixed threshold triplet does not.

`StageD_RF_Product1_NCR_Regional.R`, `StageE_RF_Product2_NCR_Regional.R`, `StageF_RF_Product3_NCR_Regional.R`, and `StageG_RF_Product4_NCR_Regional.R` each rerun the **full NCR and full 17-region analysis separately**, with output folders suffixed `_RF_Product1` through `_RF_Product4`. Disease-only outputs are preserved.

### RF temporal overlays in NCR and Regional products
Every RF-product rerun makes the environmental precursor visible in the exported figures rather than applying it only as an internal gate. At disease-trigger week **t**, the RF feature is computed strictly from **t-4..t-1** and the trigger week itself is excluded. Products 1/3 pair the RF overlay with the first RF-qualified **Constant TA** onset; Products 2/4 pair it with the first RF-qualified **Outbreak Threshold** onset. If an earlier ungated onset does not meet the RF point threshold, the overlay proceeds to the next qualifying onset. A1/A2, not the rainfall timing rule, determines whether the resulting RF-supported disease trigger is true or false.

In NCR RF-product folders, the existing per-year Constant TA panels and `FigureA3_Comparison_2019_2024_Constant_TA` inherit this overlay. `RF_Trigger_Overlay_NCR.csv` records the displayed RF week, disease-anchor week, threshold, operator, RF feature value, and antecedent window. In Regional RF-product folders, `RF_Trigger_Overlay_Regional.csv` provides the same audit trail by region-year, and `FigureB6_RF_Trigger_Overlay_Regional_<YEAR>.{pdf,png}` provides 17-region small-multiple overlays for every evaluable year.

### Stage C direct RF threshold head-to-head
`StageC_RF_Threshold_Derivation.R` now performs an immediate derivation-stage comparison before any RF-gated detector rerun. Product 1 is paired with Product 3 for the Constant TA anchor, and Product 2 is paired with Product 4 for the Outbreak Threshold anchor. The comparison uses identical REGION/YR/WN rows, paired DeLong tests for ROC AUC after direction-standardizing each predictor, exact paired sign-flip tests across all 17 regions for AUC and point-threshold operating metrics, and leave-one-region-out AUC stability. NCR and pooled 17-region discrimination are both exported. The Stage C conclusion is explicitly derivation-stage evidence only; the later Stage H comparison remains the downstream operational-performance decision.

### RF mm versus RF STA/LTA head-to-head
`StageH_RF_mm_vs_Rt_HeadToHead.R` performs two matched comparisons:
- Product 1 versus Product 3: both preceding **Constant TA**.
- Product 2 versus Product 4: both preceding **Outbreak Threshold**.

For NCR and Regional analyses, the comparison is organized into the same three reporting domains: **Burden Accuracy, Timeliness, and False Alarms**. Formal inference is anchor-focused: P1 versus P3 is tested on Constant TA performance, while P2 versus P4 is tested on Outbreak Threshold performance. NCR is paired by evaluable year and uses exact paired sign-flip tests; Regional inference treats each of the 17 regions as a cluster and applies the same exact paired sign-flip test to region-level paired summaries. Paired bootstrap 95% confidence intervals and Bonferroni-adjusted p-values are reported across the eight prespecified metrics within each anchor and spatial scale. Each anchor produces separate panel A, B, and C files plus a combined A-C figure.

Stage H also reports a bootstrap **dominance probability**, compares fixed-threshold stability using the Stage C LORO lower/point/upper framework, and evaluates balanced accuracy at the lower, point, and upper threshold bounds. `RF_HeadToHead_OverallDecision.csv` returns one of three conclusions for NCR and Regional separately: **RF mm preferred**, **RF STA/LTA R(t) preferred**, or **No clear empirical superiority**. A winner is not forced: the prespecified rule requires majority directional superiority, no Bonferroni-significant deterioration, inferential or dominance support, and no conflicting majority across the three threshold-bound sensitivity checks. The adopted point estimate remains the primary operational threshold for headline performance.

## Reproducible execution
1. Extract the archive.
2. Open R/RStudio and set the working directory to this project folder.
3. Install dependencies once with `source("R/install_dependencies.R")`.
4. Run `source("run_all.R")`.

`run_all.R` executes a parser/schema/**geometry-conversion** preflight first, then Stage 1, Stage 2, disease-only NCR, disease-only Regional, the four-product RF derivation, RF Products 1-4 NCR/Regional reruns, the RF-mm versus RF-R(t) head-to-head stage, and finally `scripts/99_validate_outputs.R`. The validator stops the run if a required figure/table is missing or empty, Figure A1 does not contain exactly 2023-2025, the NCR or Regional results omit any requested detector, the regional map/dominance exports do not contain all 17 canonical regions, any of the four RF thresholds is non-finite, the empirical rainfall R(t) windows/operators are missing, either RF representation head-to-head is incomplete, or an RF-product overlay fails the fixed antecedent t-4..t-1 RF contract.

### Start-to-end runtime safeguards
The pipeline no longer uses schema-fragile `do.call(rbind, ...)` aggregation. Heterogeneous detector, year, region, bootstrap, and significance-table branches are combined with `dplyr::bind_rows()`, and the NCR yearly-summary empty-data branch explicitly returns the same RF-overlay columns as observed years. The preflight rejects any reintroduced schema-fragile base row-binding call.

`run_all.R` executes every stage through a stage-aware warning/error wrapper. Warnings are deduplicated and written to `outputs/_run_logs/Warnings.csv` instead of accumulating into an opaque `50 or more warnings` message. If a stage errors, `outputs/_run_logs/Errors.csv` records the exact script and message before execution stops. A successful run still requires the final output validator to pass.

## Local geometry
A pinned Philippine GADM level-1 RDS file is bundled under `data/geometry/gadm/` to support regional mapping. The bundled file may deserialize as a `PackedSpatVector`; the preflight and Stage B explicitly unpack it with `terra` before `sf` conversion, preventing the `st_as_sf(PackedSpatVector)` failure. `terra` is a required dependency but is not attached to the search path, avoiding its masking of common plotting/data functions.

## Interpretation
The RF thresholds are exploratory environmental gating thresholds derived from the supplied dataset. They should not be interpreted as externally validated intervention thresholds until tested prospectively and against independent rainfall sources.

### RF operational direction and comparator robustness
The current RF products use **only the upper-tail `>=` operational rule**. Stage C, the RF handoff, NCR/Regional execution, and final validation reject any Product 1-4 whose operator is not `>=`. This prevents low rainfall from being misclassified as a strong rainfall activation.



## RF-product figure display policy

For RF Products 1-4, NCR and Regional time-series figures are rainfall-centred. Weekly leptospirosis cases are displayed as light-blue bars on the left axis. Products 1 and 2 use an antecedent rainfall **mm** right axis and never display a leptospirosis R(t) axis. Products 3 and 4 use the empirically derived rainfall STA/LTA **R(t)** right axis. The adopted RF point estimate and its fixed lower/upper empirical stability bounds are drawn in the same environmental units. The first qualifying RF event is annotated with the observed RF value, fixed point estimate, lower/upper bounds, and its operational status (precautionary, formal, or strong). The corresponding Constant TA or Outbreak Threshold disease trigger is marked after the RF antecedent signal. Disease-only leptospirosis figures retain their established TA/threshold structure; the lower/point/upper RF additions apply only to the four RF-product reruns. Exported figures suppress subtitles/captions and use padded, collected legends to prevent overlap or cropping.

### RF-product figure display rule
RF Product 1/3 figures emphasize rainfall preceding Constant TA and retain only the Constant TA disease anchor; RF Product 2/4 figures emphasize rainfall preceding the Outbreak Threshold and retain only the Outbreak Threshold disease anchor. Products 1/2 use a rainfall-mm right axis; Products 3/4 use rainfall STA/LTA R(t). The no-RF leptospirosis figures retain their prior TA/threshold structure. Legends are bottom-aligned and horizontal, using multiple rows where needed.

## Legend layout revision
All publication legends now use bottom placement, horizontal reading order, compact key/text spacing, and multi-row wrapping for long legends. The NCR drift and regional detector-map legends wrap to two rows; RF-product legends use compact horizontal environmental keys. This is a display-only revision and does not alter analytical calculations.

## Stage C primary qualifying-trigger rule

Stage C RF threshold derivation uses one primary positive disease anchor per region-year: the first detector-active Constant TA or Outbreak Threshold week that falls within A1 or A2. Detector-active weeks outside A1/A2 remain false-alarm comparators; later true active weeks are audit-only. P1/P3 share the Constant TA truth inventory and P2/P4 share the Outbreak Threshold truth inventory. Stage C NCR and Regional figures are exported separately to preserve readable axes and legends.

## Independent per-region branch

The pipeline now distinguishes three analytical scales:

1. **NCR**: the original city-scale/NCR-style analysis.
2. **Regional pooled**: the existing pooled 17-region analysis used for regional comparison, maps, dominance probabilities, and pooled RF thresholds.
3. **Independent per-region replication**: Stage I applies the NCR-style pipeline separately to each of the 17 regions. Each region receives a disease-only run and four RF-product runs. The RF runs use thresholds derived independently for that region in Stage C, with region-specific LOYO stability bounds.

Independent outputs are written under:

`outputs/StageI_PerRegion_Independent/<REGION>/DiseaseOnly/`

and

`outputs/StageI_PerRegion_Independent/<REGION>/RF_Product1/` through `RF_Product4/`.

For RF products in the pooled Regional branch, the primary environmental outputs are the RF-to-anchor files and figures. Products 1 and 3 map rainfall preceding Constant TA. Products 2 and 4 map rainfall preceding Outbreak Threshold. The six-detector Figure B5 outputs are retained only as diagnostic case-detector context.
