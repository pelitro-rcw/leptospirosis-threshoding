# Implementation Crosswalk

| Requested revision | Implementation | Main output / validation |
|---|---|---|
| NCR outbreak-threshold drift QC | `Stage2_outbreak_threshold_drift_NCR_Lepto.R` | One single-panel `FigureA1_NCR_OutbreakThreshold_Drift` restricted to 2023-2025 + CSVs |
| Figure A1 legend | Stage 2 | Exactly: Leptospirosis Cases; Outbreak Threshold (mean + 2SD); Alarm Threshold (mean + SD) |
| Weekly leptospirosis axes | Stages 2/A | Data-driven scale with headroom; no forced 800 ceiling; expanded margins prevent clipping |
| Readable R(t) secondary axis | Stage A | Uncapped dynamic R(t) with about 4-5 well-spaced labels and anti-crop right margin |
| NCR six-detector dominance/scatter | Stage A | Constant TA, Continuous TA, Outbreak Threshold, Alarm Threshold, WHO 75th, WHO 90th |
| 2019 vs 2024 Constant TA + per-year figures | Stage A | `FigureA3_Comparison_2019_2024_Constant_TA` + per-year FigureA3A exports |
| NCR A-C head-to-head | Stage A | Burden Accuracy, Timeliness, False Alarms panels |
| Regional map must include all 17 regions | Stage B | Coverage-based eligibility; hard assertions on analytical table, geometry dissolve/join, map table, and dominance table |
| Fix PackedSpatVector map crash | Stage B + preflight | `terra::unwrap()` before `sf::st_as_sf()`; conversion dry-run in `00_self_check.R` |
| Fix `consensus_pass` regional failure | Stage B | Required schema created explicitly; obsolete unused plotting/transmute failure path removed |
| RF Product 1 | Stage C + D | RF mm preceding first Constant TA; full NCR + 17-region rerun |
| RF Product 2 | Stage C + E | RF mm preceding first Outbreak Threshold; full NCR + 17-region rerun |
| RF Product 3 | Stage C + F | Empirical rainfall STA/LTA R(t) preceding first Constant TA; full NCR + 17-region rerun |
| RF Product 4 | Stage C + G | Empirical rainfall STA/LTA R(t) preceding first Outbreak Threshold; full NCR + 17-region rerun |
| Rainfall R(t) must be empirically derived | Stage C | Rainfall-specific candidate STA/LTA windows and multi-method activation thresholds; no borrowing of leptospirosis case eta |
| Empirical RF trigger direction | Stage C/config/A/B | All operational RF products use the prespecified upper-tail `>=` direction; NCR and Regional use scope-specific fixed thresholds written to the handoff |
| RF mm vs RF R(t) head-to-head | Stage H | P1 vs P3 and P2 vs P4, each in Burden Accuracy, Timeliness, False Alarms; NCR + paired 17-region Regional outputs |
| RF temporal trigger overlay | Stage A/B RF reruns | For the first ungated target Constant-TA (P1/P3) or Outbreak-Threshold (P2/P4) anchor, the figures search the fixed 1-4-week warning horizon and show the actual first formal upper-tail RF threshold crossing if one exists. If the point threshold is not reached, no RF-trigger line is drawn; the maximum antecedent RF can still be reported as no/precautionary activation. Regional RF reruns add 17-region yearly Figure B6 overlays and audit CSVs. |
| Avoid year-gap bridging | Shared continuity helpers | STA/LTA and antecedent RF windows require exact seven-day continuity |
| Fail on incomplete output | `scripts/99_validate_outputs.R` | Verifies four RF products, empirical RF R(t) fields/operators, all 17 regions, six detectors, and all head-to-head outputs |
| Fixed RF lower/point/upper bounds | Stage C + shared config + RF reruns | Scope-specific fixed bounds: NCR uses LOYO stability around NCR-derived P1-P4 thresholds; Regional uses LORO stability around pooled-17-region thresholds. The point estimate remains primary; lower/upper provide monotonic precautionary/strong upper-tail activation levels |
| Direct RF mm vs RF STA/LTA derivation head-to-head | Stage C | P1 vs P3 and P2 vs P4 on common weeks; NCR and pooled AUC, paired DeLong, 17-region exact sign-flip tests, LORO AUC stability, and A-C head-to-head figures |

| Prospective RF false-alarm window | Stage C + shared config + NCR/Regional RF reruns | A1/A2 remains the independent truth standard: an RF-supported disease anchor is true iff its trigger week lies in A1 OR A2; otherwise it is a false alarm. Overlays use the actual first qualifying RF trigger. |
