# Validation Report

This revision uses four validation layers.

1. **Static source/data-contract QA before packaging.** Required scripts, detector names, requested figure stems, four RF product runners, empirical rainfall STA/LTA derivation, operator-aware gating, 17-region coverage, single-panel Figure A1 with a one-row horizontal three-item legend, dynamic weekly-case/R(t) axes, regional geometry handling, and RF mm versus RF R(t) comparisons are checked. R lexical delimiters/quotes are checked for balance and the ZIP container is tested after creation.
2. **R-native preflight.** `scripts/00_self_check.R` parses every R file, verifies the required workbook sheet/columns and exactly 17 canonical regions, deserializes the bundled GADM RDS, unpacks a `PackedSpatVector` with `terra`, and dry-runs conversion to `sf` before analysis.
3. **Attachment-order independence.** All executable data pipelines use R's native `|>` operator rather than `%>%`. The preflight scans the project and stops if the legacy magrittr pipe is reintroduced, eliminating the reported `could not find function "%>%"` failure from clean-session execution.

4. **R-native final output contract.** `scripts/99_validate_outputs.R` runs last and stops the pipeline if any required PDF/PNG/CSV is absent or empty; if Figure A1 years are not exactly 2023-2025; if NCR or Regional outputs omit any of the six detectors; if the regional map or dominance table omits any of the 17 canonical regions; if any of Products 1-4 lacks a finite threshold; if Products 3/4 lack empirical rainfall STA/LTA windows or valid comparison operators; if any of the four NCR/Regional reruns is incomplete; or if the RF mm versus RF R(t) head-to-head outputs are incomplete.

The authoring container does not provide an R runtime, so it is not possible to truthfully claim a numerical end-to-end R execution inside this container. The archive therefore includes fail-fast R-native parse, schema, geometry, regional-coverage, RF-handoff, and final-output validation so the target R environment cannot silently report success after an incomplete run.

## RF operator runtime regression correction
The current RF threshold handoff enforces `>=` for all four operational products. Point estimates are false-alarm-aware and scope-specific; NCR uses NCR-derived thresholds and Regional uses pooled-17-region thresholds.
A prior implementation compared a named scalar operator with `identical()`, which returns `FALSE`
even when the printed value is the same. The pipeline now strips names/whitespace with
`normalize_rf_operator()` and applies all gates through `apply_rf_gate()`. Preflight and final
validation explicitly tests the upper-tail operational comparator and the ordered no/precautionary/formal/strong classification. Non-upper-tail RF operators are rejected.


## 2026-08-18 row-binding/runtime hardening

The reported runtime failure `Error in rbind(deparse.level, ...) : numbers of columns of arguments do not match` was traced to schema-fragile base row binding, including the NCR yearly-summary builder after RF-overlay fields were added. The empty-year branch had the pre-overlay schema while observed-year branches carried `First_RF_Trigger_Week`, `RF_Anchor_Detector`, and `RF_Scale`.

Corrections in this build:

- All executable `do.call(rbind, ...)` calls were replaced with `dplyr::bind_rows()`.
- The NCR yearly-summary empty-data branch now returns the same RF-overlay schema as observed-year branches.
- Stage A asserts the exact yearly-summary column contract and 13-row 2013-2025 skeleton before plotting/export.
- The preflight rejects future reintroduction of schema-fragile `do.call(rbind, ...)` calls.
- `surveillance` is included in the required dependency installer and in the full-pipeline preflight package audit.
- `run_all.R` records deduplicated stage warnings in `outputs/_run_logs/Warnings.csv` and stage errors in `outputs/_run_logs/Errors.csv`. This prevents an accumulated `50 or more warnings` message from obscuring the actual failing stage while retaining a reproducible warning audit.

A successful run remains conditional on `scripts/99_validate_outputs.R` reporting `FINAL OUTPUT VALIDATION: PASS`.

## RF figure scope correction

The RF visualization layer is explicitly separated from the disease-only visualization layer. Disease-only NCR timing panels retain the previously validated Constant/Continuous TA R(t), activation threshold, Outbreak Threshold, and trigger-marker structure. RF Product 1/3 timing figures retain only the Constant TA disease anchor and suppress Outbreak Threshold visual elements; RF Product 2/4 timing figures retain only the Outbreak Threshold disease anchor and suppress Constant TA visual elements. Products 1/2 use an antecedent rainfall-mm secondary axis; Products 3/4 use the empirically derived rainfall STA/LTA R(t) secondary axis. Constant-TA eta_on is reported as anchor text in P1/P3 rather than plotted against the environmental axis, avoiding a dimensionally invalid case-R(t)-versus-rainfall axis. RF subtitles/captions are suppressed, while non-RF figures retain their prior structure. Legend boxes are bottom-aligned and horizontal, with multi-row guides where needed.

A source/data-contract audit confirmed balanced R delimiters/quotes, zero legacy `%>%` use, zero vertical legend-box/direction declarations, distinct mm and rainfall-R(t) axis contracts, anchor-specific RF output naming, and complete 17-region / 102 region-year workbook coverage for the six evaluable years. This environment does not provide an R runtime, so numerical end-to-end execution must still be performed in the target R installation; the package's R-native preflight and final validator remain authoritative for runtime completion.

## RF fixed lower/point/upper operational thresholds

Stage C now exports one fixed lower bound, adopted point estimate, and upper bound for each of Products 1-4 separately for NCR and Regional analysis. NCR uses NCR-derived point estimates with leave-one-year-out stability bounds; Regional uses pooled-17-region point estimates with leave-one-region-out stability bounds. Products 3/4 re-select the rainfall STA/LTA window inside each stability fold. The point estimate remains the primary threshold for headline metrics. The lower/point/upper interpretation is monotonic in rainfall: below lower = no activation, lower-to-point = precautionary, point-to-upper = formal, and at/above upper = strong. RF figures display the scope-appropriate fixed triplet on the correct environmental axis.
