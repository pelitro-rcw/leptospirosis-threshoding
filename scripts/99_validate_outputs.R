# =============================================================================
# FINAL OUTPUT VALIDATOR
# =============================================================================
# Called last by run_all.R. A successful pipeline run is reported only when all
# required analytical products exist, are non-empty, and contain the requested
# detector/year content. This converts silent/partial figure failures into a
# hard error with an explicit missing-output list.
# =============================================================================
root <- Sys.getenv("TA_PROJECT_ROOT", unset = "")
if (!nzchar(root)) root <- normalizePath(getwd(), winslash = "/")
source(file.path(root, "R", "00_config.R"))

# Operational RF products must use the upper-tail >= rule. Validate both the
# comparator and the lower/point/upper activation hierarchy before reading outputs.
z <- apply_rf_gate(c(0, 1, 2), 1, ">=")
if (!identical(z, c(FALSE, TRUE, TRUE))) {
  stop("Final RF upper-tail comparator validation failed.", call. = FALSE)
}
cls <- classify_rf_activation(c(1, 2, 3, 4), point = 3, lower = 2, upper = 4,
                              operator = ">=")
if (!identical(cls, c("No RF activation", "Precautionary RF trigger",
                      "Formal RF trigger", "Strong RF trigger"))) {
  stop("Final RF operational-bound classification validation failed.", call. = FALSE)
}

assert_files <- function(paths, label) {
  info <- file.info(paths)
  bad <- paths[!file.exists(paths) | is.na(info$size) | info$size <= 0]
  if (length(bad)) {
    stop(label, " validation failed. Missing/empty output(s):\n  ",
         paste(bad, collapse = "\n  "), call. = FALSE)
  }
  invisible(TRUE)
}

pair <- function(dir, stem) file.path(dir, paste0(stem, c(".pdf", ".png")))

stage1 <- file.path(DIR_OUTPUT, "Stage1_eta_thresholds_NCR")
assert_files(c(
  pair(stage1, "Fig_eta_threshold_derivation_QC"),
  file.path(stage1, "eta_thresholds_derived.R"),
  file.path(stage1, "sta_lta_windows_derived.R"),
  file.path(stage1, "Table_FINAL_eta_decision.csv")
), "Stage 1")

stage2 <- file.path(DIR_OUTPUT, "Stage2_outbreak_threshold_drift_NCR_Lepto")
assert_files(c(
  pair(stage2, "FigureA1_NCR_OutbreakThreshold_Drift"),
  file.path(stage2, "FigureA1_NCR_OutbreakThreshold_Drift_2023_2025.csv"),
  file.path(stage2, "FigureA1_threshold_drift_comparison.csv")
), "Stage 2")

# Figure A1 must contain exactly the requested displayed years.
a1 <- utils::read.csv(file.path(stage2, "FigureA1_NCR_OutbreakThreshold_Drift_2023_2025.csv"),
                      stringsAsFactors = FALSE)
if (!"YR" %in% names(a1) || !identical(sort(unique(as.integer(a1$YR))), 2023:2025)) {
  stop("FigureA1 validation failed: displayed years must be exactly 2023, 2024, 2025.",
       call. = FALSE)
}
a1_required_cols <- c("Cases", "Outbreak_Calendar", "Alarm_Calendar")
missing_a1_cols <- setdiff(a1_required_cols, names(a1))
if (length(missing_a1_cols)) {
  stop("FigureA1 validation failed: missing plotted series column(s): ",
       paste(missing_a1_cols, collapse = ", "), call. = FALSE)
}

required_methods <- c(
  "Constant Transmission Acceleration",
  "Continuous Transmission Acceleration",
  "Outbreak Threshold",
  "Alarm Threshold",
  "WHO 75th Percentile Threshold",
  "WHO 90th Percentile Threshold"
)

years_required <- c(2018L, 2019L, 2022L, 2023L, 2024L, 2025L)

validate_stage_a <- function(suffix = "") {
  d <- file.path(DIR_OUTPUT, paste0("StageA_NCR_Lepto_Analysis", suffix))
  is_rf <- nzchar(suffix)
  product_id <- if (is_rf) suppressWarnings(as.integer(sub(".*RF_Product", "", suffix))) else NA_integer_
  anchor_is_cta <- is_rf && product_id %in% c(1L, 3L)
  anchor_is_ot  <- is_rf && product_id %in% c(2L, 4L)

  stems <- c(
    "FigureA2_panel_a_DominanceMatrix",
    "FigureA2_panel_b_Scatter",
    "FigureA3_HeadToHead_panel_a_BurdenAccuracy",
    "FigureA3_HeadToHead_panel_b_Timeliness",
    "FigureA3_HeadToHead_panel_c_FalseAlarms"
  )
  if (!is_rf || anchor_is_cta) {
    stems <- c(stems,
      "FigureA3_Comparison_2019_2024_Constant_TA",
      "FigureA3A_Constant_Transmission_Acceleration_Multipanel")
    yearly_stem <- "FigureA3A_Constant_Transmission_Acceleration_"
  } else if (anchor_is_ot) {
    stems <- c(stems,
      "FigureA3_Comparison_2019_2024_Outbreak_Threshold",
      "FigureA3A_Outbreak_Threshold_Multipanel")
    yearly_stem <- "FigureA3A_Outbreak_Threshold_"
  } else {
    stop("Stage A validator could not infer RF product from suffix: ", suffix, call. = FALSE)
  }

  files <- unlist(lapply(stems, function(s) pair(d, s)), use.names = FALSE)
  files <- c(files,
             unlist(lapply(years_required, function(y)
               pair(d, paste0(yearly_stem, y))), use.names = FALSE),
             file.path(d, "FigureA2_method_summary_primary.csv"),
             file.path(d, "FigureA2_dominance_matrix_long.csv"),
             file.path(d, "FigureA3A_Constant_Transmission_Acceleration_Yearly_Summary.csv"),
             file.path(d, "evaluation_framework", "Table2_Per_Year_Detail.csv"))
  assert_files(files, paste0("Stage A", suffix))

  sm <- utils::read.csv(file.path(d, "FigureA2_method_summary_primary.csv"),
                        stringsAsFactors = FALSE)
  mcol <- intersect(c("Method", "Detector"), names(sm))
  if (!length(mcol)) stop("Stage A", suffix, ": method summary has no Method/Detector column.", call. = FALSE)
  missing_methods <- setdiff(required_methods, unique(as.character(sm[[mcol[1]]])))
  if (length(missing_methods)) {
    stop("Stage A", suffix, " is missing requested detector(s): ",
         paste(missing_methods, collapse = ", "), call. = FALSE)
  }
}

validate_stage_b <- function(suffix = "") {
  d <- file.path(DIR_OUTPUT, paste0("StageB_Regional_Lepto_Analysis", suffix))
  files <- c(
    pair(d, "FigureB5_Combined_RegionalDominance"),
    pair(d, "FigureB5_panel_a_DetectorMap"),
    pair(d, "FigureB5_panel_b_DominanceProbability"),
    pair(d, "FigureB4_Method_Summary"),
    file.path(d, "StageB_Regional_8Metric_Summary.csv"),
    file.path(d, "StageB_Regional_Dominance_Probabilities.csv"),
    file.path(d, "StageB_Regional_Dominance_Matrix.csv"),
    file.path(d, "StageB_Detector_Map_RegionTable.csv"),
    file.path(d, "StageB_Detector_Map_JoinAudit.csv"),
    file.path(d, "StageB_Regional_Bootstrap_Replicates.csv")
  )
  assert_files(files, paste0("Stage B", suffix))

  sm <- utils::read.csv(file.path(d, "StageB_Regional_8Metric_Summary.csv"),
                        stringsAsFactors = FALSE)
  if (!"Method" %in% names(sm)) stop("Stage B", suffix, ": summary has no Method column.", call. = FALSE)
  missing_methods <- setdiff(required_methods, unique(as.character(sm$Method)))
  if (length(missing_methods)) {
    stop("Stage B", suffix, " is missing requested detector(s): ",
         paste(missing_methods, collapse = ", "), call. = FALSE)
  }

  # Bootstrap schema/runtime contract. Every one of the 17 regions must have
  # BOOT_N_CI replicates for each of the six requested methods, even when an RF
  # gate produces zero qualifying triggers. Zero-trigger replicates are real
  # observations and must never disappear as NULL rows.
  boot_rep <- utils::read.csv(
    file.path(d, "StageB_Regional_Bootstrap_Replicates.csv"),
    stringsAsFactors = FALSE
  )
  required_boot_cols <- c(
    "REGION", ".rep", "Method", "TAM", "N_True_Alarms", "Sensitivity",
    "Mean_Lead_Time", "WP", "PPV", "ALY", "N_False_Alarms"
  )
  missing_boot_cols <- setdiff(required_boot_cols, names(boot_rep))
  if (length(missing_boot_cols)) {
    stop("Stage B", suffix, " bootstrap replicate table is missing column(s): ",
         paste(missing_boot_cols, collapse = ", "), call. = FALSE)
  }
  if (any(is.na(boot_rep$Method)) || any(!nzchar(as.character(boot_rep$Method)))) {
    stop("Stage B", suffix, " bootstrap replicate table contains missing Method values.",
         call. = FALSE)
  }
  boot_region_counts <- table(as.character(boot_rep$REGION))
  expected_per_region <- REGIONAL_BOOT_N_CI * length(required_methods)
  if (length(boot_region_counts) != 17L ||
      any(as.integer(boot_region_counts) != expected_per_region)) {
    stop("Stage B", suffix,
         " bootstrap replicate table must contain exactly ", expected_per_region,
         " rows per region (", REGIONAL_BOOT_N_CI, " replicates x ", length(required_methods),
         " methods).", call. = FALSE)
  }

  # Regional Figure B5 schema contract. This specifically guards the previously
  # observed consensus_pass / sig_star_dominance failure and confirms that the
  # map result table was fully constructed before the stage is accepted.
  region_table <- utils::read.csv(
    file.path(d, "StageB_Detector_Map_RegionTable.csv"),
    stringsAsFactors = FALSE
  )
  required_region_cols <- c(
    "REGION", "Observed_Winner", "Dominance_Probability",
    "consensus_R1_pass", "consensus_pass", "consensus_tier",
    "consensus_winner", "p_consensus", "sig_star_dominance",
    "Leader_Label_short"
  )
  missing_region_cols <- setdiff(required_region_cols, names(region_table))
  if (length(missing_region_cols)) {
    stop("Stage B", suffix, " regional map table is missing required column(s): ",
         paste(missing_region_cols, collapse = ", "), call. = FALSE)
  }

  canonical_17 <- c(
    "BARMM", "CAR", "MIMAROPA", "NCR",
    paste("REGION", c("I", "II", "III", "IV-A", "V", "VI", "VII", "VIII",
                      "IX", "X", "XI", "XII", "XIII"))
  )
  observed_regions <- unique(as.character(region_table$REGION))
  unknown_regions <- setdiff(observed_regions, canonical_17)
  if (length(unknown_regions)) {
    stop("Stage B", suffix, " regional map table contains non-canonical region(s): ",
         paste(unknown_regions, collapse = ", "), call. = FALSE)
  }
  missing_regions <- setdiff(canonical_17, observed_regions)
  if (!setequal(observed_regions, canonical_17) || length(observed_regions) != 17L ||
      anyDuplicated(region_table$REGION)) {
    stop("Stage B", suffix,
         ": regional map table must contain exactly all 17 canonical Philippine regions. Missing: ",
         paste(missing_regions, collapse = ", "), call. = FALSE)
  }

  dom <- utils::read.csv(
    file.path(d, "StageB_Regional_Dominance_Probabilities.csv"),
    stringsAsFactors = FALSE
  )
  if (!"REGION" %in% names(dom)) {
    stop("Stage B", suffix, ": dominance-probability table has no REGION column.",
         call. = FALSE)
  }
  prob_cols <- c("P_ConstantTA", "P_ContinuousTA", "P_OutbreakThreshold",
                 "P_AlarmThreshold", "P_WHO75", "P_WHO90")
  missing_prob_cols <- setdiff(prob_cols, names(dom))
  if (length(missing_prob_cols)) {
    stop("Stage B", suffix, " dominance-probability table is missing six-detector column(s): ",
         paste(missing_prob_cols, collapse = ", "), call. = FALSE)
  }
  dom_regions <- unique(as.character(dom$REGION))
  if (!setequal(dom_regions, canonical_17) || length(dom_regions) != 17L) {
    stop("Stage B", suffix, ": dominance-probability table must contain all 17 regions. Missing: ",
         paste(setdiff(canonical_17, dom_regions), collapse = ", "), call. = FALSE)
  }

  join_audit <- utils::read.csv(
    file.path(d, "StageB_Detector_Map_JoinAudit.csv"), stringsAsFactors = FALSE)
  if (!all(c("REGION", "status") %in% names(join_audit))) {
    stop("Stage B", suffix, ": map join audit has an unexpected schema.", call. = FALSE)
  }
  canonical_join <- join_audit[join_audit$REGION %in% canonical_17, , drop = FALSE]
  if (nrow(canonical_join) != 17L || anyDuplicated(canonical_join$REGION) ||
      any(as.character(canonical_join$status) != "matched")) {
    bad <- canonical_join$REGION[as.character(canonical_join$status) != "matched"]
    stop("Stage B", suffix, ": all 17 canonical regions must be matched in the map geometry join. Problem region(s): ",
         paste(unique(bad), collapse = ", "), call. = FALSE)
  }


  if (nzchar(suffix)) {
    product_id <- suppressWarnings(as.integer(sub(".*RF_Product", "", suffix)))
    if (!product_id %in% 1:4) stop("Stage B RF validator cannot infer product from ", suffix, call. = FALSE)
    anchor_stem <- if (product_id %in% c(1L, 3L)) "ConstantTA" else "OutbreakThreshold"
    rf_files <- c(
      file.path(d, "RF_to_Anchor_Trigger_Detail.csv"),
      file.path(d, "RF_to_Anchor_Regional_Summary.csv"),
      file.path(d, "RF_Product_Output_Guide.txt"),
      pair(d, paste0("FigureB7_RF_to_", anchor_stem, "_Map"))
    )
    assert_files(rf_files, paste0("Stage B RF-to-anchor", suffix))
    rs <- utils::read.csv(file.path(d, "RF_to_Anchor_Regional_Summary.csv"), stringsAsFactors = FALSE)
    need_rs <- c("REGION", "RF_Product", "RF_Representation", "Disease_Anchor",
                 "N_Anchor_Trigger_Weeks", "N_RF_Supported_Anchor_Weeks",
                 "N_RF_Supported_True_Alarms", "N_RF_Supported_False_Alarms",
                 "RF_Supported_PPV", "RF_Qualification_Rate")
    miss_rs <- setdiff(need_rs, names(rs))
    if (length(miss_rs) || nrow(rs) != 17L || anyDuplicated(rs$REGION) ||
        !setequal(as.character(rs$REGION), canonical_17)) {
      stop("Stage B", suffix, " RF-to-anchor summary must contain a complete 17-region RF->anchor mapping. Missing columns: ",
           paste(miss_rs, collapse=", "), call. = FALSE)
    }
    expected_anchor <- if (product_id %in% c(1L,3L)) "Constant TA" else "Outbreak Threshold"
    if (any(as.integer(rs$RF_Product) != product_id) ||
        any(as.character(rs$Disease_Anchor) != expected_anchor)) {
      stop("Stage B", suffix, " RF-to-anchor mapping has the wrong product/anchor relation.", call. = FALSE)
    }
  }
}

# Disease-only products.
validate_stage_a("")
validate_stage_b("")

stagec <- file.path(DIR_OUTPUT, "StageC_RF_Threshold_Derivation")
assert_files(c(
  file.path(stagec, "Product1_First_Constant_TA_Events.csv"),
  file.path(stagec, "Product2_First_Outbreak_Threshold_Events.csv"),
  file.path(stagec, "Product1_First_A1A2_True_Constant_TA_Anchor.csv"),
  file.path(stagec, "Product2_First_A1A2_True_Outbreak_Threshold_Anchor.csv"),
  file.path(stagec, "Product1_StageC_Derivation_Inventory_A1A2.csv"),
  file.path(stagec, "Product2_StageC_Derivation_Inventory_A1A2.csv"),
  file.path(stagec, "RF_Threshold_Derivation_Support.csv"),
  file.path(stagec, "Product1_All_Constant_TA_Onsets_A1A2.csv"),
  file.path(stagec, "Product2_All_Outbreak_Threshold_Onsets_A1A2.csv"),
  file.path(stagec, "Product1_All_Constant_TA_TriggerWeeks_A1A2.csv"),
  file.path(stagec, "Product2_All_Outbreak_Threshold_TriggerWeeks_A1A2.csv"),
  file.path(stagec, "RF_Product_Anchor_Truth_Mapping.csv"),
  file.path(stagec, "RF_mm_ROC_Input.csv"),
  file.path(stagec, "RF_Rt_ROC_Input.csv"),
  file.path(stagec, "Product3_First_Constant_TA_RF_Rt_Events.csv"),
  file.path(stagec, "Product4_First_Outbreak_Threshold_RF_Rt_Events.csv"),
  file.path(stagec, "RF_Threshold_Derivation_Summary.csv"),
  file.path(stagec, "RF_Threshold_Derivation_Summary_NCR.csv"),
  file.path(stagec, "RF_Threshold_Derivation_Summary_AllScopes.csv"),
  file.path(stagec, "RF_Threshold_Derivation_Summary_PerRegion.csv"),
  file.path(stagec, "RF_Threshold_LOYO_Stability_PerRegion.csv"),
  file.path(stagec, "RF_STA_LTA_Window_Candidates_PerRegion.csv"),
  file.path(stagec, "RF_Threshold_Derivation_Support_PerRegion.csv"),
  file.path(stagec, "RF_Threshold_Operational_Bounds.csv"),
  file.path(stagec, "RF_Threshold_Operational_Bounds_NCR.csv"),
  file.path(stagec, "RF_Threshold_Bound_Sensitivity.csv"),
  file.path(stagec, "RF_Threshold_Bound_Sensitivity_NCR.csv"),
  file.path(stagec, "RF_Threshold_LOYO_Stability_NCR.csv"),
  file.path(stagec, "RF_Threshold_HeadToHead_CommonRows.csv"),
  file.path(stagec, "RF_Threshold_HeadToHead_Overall.csv"),
  file.path(stagec, "RF_Threshold_HeadToHead_ByRegion.csv"),
  file.path(stagec, "RF_Threshold_HeadToHead_ExactPairedTests.csv"),
  file.path(stagec, "RF_Threshold_HeadToHead_LORO.csv"),
  file.path(stagec, "RF_Threshold_HeadToHead_LORO_Tests.csv"),
  file.path(stagec, "RF_Threshold_HeadToHead_Decision.csv"),
  file.path(stagec, "RF_STA_LTA_Window_Candidates.csv"),
  file.path(stagec, "RF_Operational_Threshold_Candidate_Grid.csv"),
  file.path(stagec, "RF_Rt_Threshold_Empirical_Methods.csv"),
  file.path(stagec, "RF_Threshold_LORO_Stability.csv"),
  file.path(stagec, "rf_thresholds_derived.R"),
  file.path(stagec, "StageC_Methods_Note.txt"),
  pair(stagec, "FigureC_NCR_RF_Threshold_Products"),
  pair(stagec, "FigureC_Regional_RF_Threshold_Products"),
  pair(stagec, "FigureC_NCR_RF_STA_LTA_Window_Selection"),
  pair(stagec, "FigureC_Regional_RF_STA_LTA_Window_Selection"),
  pair(stagec, "FigureC_NCR_RF_Operational_Threshold_Bounds"),
  pair(stagec, "FigureC_Regional_RF_Operational_Threshold_Bounds"),
  pair(stagec, "FigureC_NCR_RF_HeadToHead_ConstantTA"),
  pair(stagec, "FigureC_Regional_RF_HeadToHead_ConstantTA"),
  pair(stagec, "FigureC_NCR_RF_HeadToHead_OutbreakThreshold"),
  pair(stagec, "FigureC_Regional_RF_HeadToHead_OutbreakThreshold")
), "Stage C")

# Stage C truth-contract validation: RF threshold derivation must use the same
# A1/A2 true/false alarm definition as Stage A/B, not a separate future-anchor
# definition. Event must be identical to IsTrue in both RF representations.
for (f in c("RF_mm_ROC_Input.csv", "RF_Rt_ROC_Input.csv")) {
  zz <- utils::read.csv(file.path(stagec, f), stringsAsFactors = FALSE)
  need_truth <- c("REGION", "YR", "WN", "InA1", "InA2", "IsTrue", "Alarm_Class", "Event")
  miss_truth <- setdiff(need_truth, names(zz))
  if (length(miss_truth)) {
    stop("Stage C A1/A2 truth validation failed for ", f, ": missing ",
         paste(miss_truth, collapse = ", "), call. = FALSE)
  }
  event <- as.integer(zz$Event)
  istrue <- as.integer(as.logical(zz$IsTrue))
  ina1 <- as.logical(zz$InA1)
  ina2 <- as.logical(zz$InA2)
  if (any(event != istrue, na.rm = TRUE) ||
      any(istrue != as.integer(ina1 | ina2), na.rm = TRUE)) {
    stop("Stage C A1/A2 truth validation failed for ", f,
         ": Event must equal IsTrue = InA1 OR InA2.", call. = FALSE)
  }
}

# A primary positive is unique per product/scope/region-year: it is the first
# detector-active trigger week within A1 OR A2. Later true weeks are audit-only.
for (f in c("RF_mm_ROC_Input.csv", "RF_Rt_ROC_Input.csv")) {
  zz <- utils::read.csv(file.path(stagec, f), stringsAsFactors = FALSE)
  grp <- intersect(c("Scope", "Product_ID", "REGION", "YR"), names(zz))
  pos <- zz[as.integer(zz$Event) == 1L, , drop = FALSE]
  if (nrow(pos)) {
    key <- do.call(paste, c(pos[grp], sep = "|"))
    if (any(table(key) > 1L)) {
      stop("Stage C primary-anchor validation failed for ", f,
           ": more than one positive derivation anchor exists in a product/scope/region-year.",
           call. = FALSE)
    }
  }
}

for (f in c("Product1_First_A1A2_True_Constant_TA_Anchor.csv",
            "Product2_First_A1A2_True_Outbreak_Threshold_Anchor.csv")) {
  aa <- utils::read.csv(file.path(stagec, f), stringsAsFactors = FALSE)
  if (nrow(aa) && (!all(as.logical(aa$IsTrue)) ||
                   anyDuplicated(paste(aa$REGION, aa$YR, sep = "|")))) {
    stop("Stage C first-true-anchor validation failed for ", f,
         ": anchors must be A1/A2 true and unique per region-year.", call. = FALSE)
  }
}

support <- utils::read.csv(file.path(stagec, "RF_Threshold_Derivation_Support.csv"),
                           stringsAsFactors = FALSE)
need_support <- c("Scope", "Product_ID", "Anchor", "N_Evaluable_Units",
                  "N_Units_With_Any_Detector_Trigger",
                  "N_Units_With_A1A2_True_Trigger", "N_Qualifying_Years",
                  "N_Primary_True_Anchors", "N_False_Comparator_Triggers",
                  "N_Later_True_Triggers_AuditOnly")
if (!all(need_support %in% names(support)) || nrow(support) != 8L ||
    !identical(sort(unique(as.integer(support$Product_ID))), 1:4)) {
  stop("Stage C derivation-support audit is incomplete or has an unexpected schema.",
       call. = FALSE)
}

per_region_rf <- utils::read.csv(
  file.path(stagec, "RF_Threshold_Derivation_Summary_PerRegion.csv"),
  stringsAsFactors = FALSE, check.names = FALSE)
need_pr <- c("REGION", "Product_ID", "Anchor", "RF_Scale", "Adopted_Threshold",
             "Threshold_Lower", "Threshold_Upper", "Comparison_Operator",
             "RF_STA", "RF_LTA", "RF_Guard")
miss_pr <- setdiff(need_pr, names(per_region_rf))
if (length(miss_pr) || nrow(per_region_rf) != 68L ||
    length(unique(per_region_rf$REGION)) != 17L || any(table(per_region_rf$REGION) != 4L)) {
  stop("Stage C independent per-region threshold table must contain exactly Products 1-4 for all 17 regions. Missing columns: ",
       paste(miss_pr, collapse=", "), call. = FALSE)
}
if (any(!is.finite(as.numeric(per_region_rf$Adopted_Threshold))) ||
    any(!is.finite(as.numeric(per_region_rf$Threshold_Lower))) ||
    any(!is.finite(as.numeric(per_region_rf$Threshold_Upper))) ||
    any(as.character(per_region_rf$Comparison_Operator) != ">=")) {
  stop("Stage C independent per-region RF thresholds contain invalid/non-upper-tail values.", call. = FALSE)
}
for (rr in unique(as.character(per_region_rf$REGION))) {
  z <- per_region_rf[per_region_rf$REGION == rr, , drop=FALSE]
  z <- z[order(as.integer(z$Product_ID)), , drop=FALSE]
  if (!identical(as.integer(z$Product_ID), 1:4) ||
      !identical(as.character(z$Anchor), c("Constant TA", "Outbreak Threshold", "Constant TA", "Outbreak Threshold"))) {
    stop("Stage C per-region RF product-to-anchor mapping failed for ", rr, ".", call. = FALSE)
  }
}

anchor_map <- utils::read.csv(file.path(stagec, "RF_Product_Anchor_Truth_Mapping.csv"),
                              stringsAsFactors = FALSE)
if (!all(c("RF_Product", "Anchor_Truth_Product_ID", "Disease_Anchor") %in% names(anchor_map)) ||
    nrow(anchor_map) != 4L ||
    !identical(as.integer(anchor_map$RF_Product), 1:4) ||
    !identical(as.integer(anchor_map$Anchor_Truth_Product_ID), c(1L, 2L, 1L, 2L)) ||
    !identical(as.character(anchor_map$Disease_Anchor),
               c("Constant TA", "Outbreak Threshold", "Constant TA", "Outbreak Threshold"))) {
  stop("Stage C anchor mapping validation failed: expected P1/P3 -> Constant TA and P2/P4 -> Outbreak Threshold.",
       call. = FALSE)
}

rf <- utils::read.csv(file.path(stagec, "RF_Threshold_Derivation_Summary.csv"),
                      stringsAsFactors = FALSE)
rf_ncr <- utils::read.csv(file.path(stagec, "RF_Threshold_Derivation_Summary_NCR.csv"),
                          stringsAsFactors = FALSE)
required_rf_cols <- c("Scope", "Product_ID", "Anchor", "RF_Scale", "Adopted_Threshold",
                      "Threshold_Lower", "Threshold_Upper",
                      "Precautionary_Activation_Bound", "Strong_Activation_Bound",
                      "Bounds_Method", "N_Stability_Thresholds",
                      "Units", "Comparison_Operator", "ROC_Direction",
                      "RF_STA", "RF_LTA", "RF_Guard",
                      "False_Alarm_Rate", "False_Alarm_Weeks",
                      "False_Alarm_Weeks_per_Evaluable_Unit", "PPV", "Balanced_Accuracy",
                      "Operational_Direction", "N_Evaluable_Units",
                      "N_RF_Eligible_Qualifying_Units", "N_RF_Eligible_Qualifying_Years",
                      "N_Units_With_Any_Detector_Trigger",
                      "N_Units_With_A1A2_True_Trigger", "N_Qualifying_Years",
                      "N_Primary_True_Anchors", "N_False_Comparator_Triggers")
validate_rf_summary <- function(z, scope_label) {
  miss <- setdiff(required_rf_cols, names(z))
  if (length(miss)) {
    stop("Stage C ", scope_label, " RF summary missing column(s): ",
         paste(miss, collapse = ", "), call. = FALSE)
  }
  if (!identical(sort(unique(as.integer(z$Product_ID))), 1:4) || nrow(z) != 4L ||
      any(!is.finite(as.numeric(z$Adopted_Threshold)))) {
    stop("Stage C ", scope_label, " must contain one finite threshold for Products 1-4.",
         call. = FALSE)
  }
  lo <- suppressWarnings(as.numeric(z$Threshold_Lower))
  pt <- suppressWarnings(as.numeric(z$Adopted_Threshold))
  hi <- suppressWarnings(as.numeric(z$Threshold_Upper))
  if (any(!is.finite(lo)) || any(!is.finite(hi)) || any(lo > pt) || any(hi < pt)) {
    stop("Stage C ", scope_label, " RF bounds must satisfy lower <= point <= upper.",
         call. = FALSE)
  }
  if (any(as.character(z$Comparison_Operator) != ">=")) {
    stop("Stage C ", scope_label,
         " must use the upper-tail RF >= threshold concept for all four products.",
         call. = FALSE)
  }
  if (!identical(as.character(z$RF_Scale[match(1:4, z$Product_ID)]),
                 c("mm", "mm", "R(t)", "R(t)"))) {
    stop("Stage C ", scope_label, " RF product scales must be mm, mm, R(t), R(t).",
         call. = FALSE)
  }
  rt_rows <- z$Product_ID %in% c(3L, 4L)
  if (any(!is.finite(as.numeric(z$RF_STA[rt_rows]))) ||
      any(!is.finite(as.numeric(z$RF_LTA[rt_rows]))) ||
      any(!is.finite(as.numeric(z$RF_Guard[rt_rows])))) {
    stop("Stage C ", scope_label,
         " Products 3/4 require empirically selected rainfall STA/LTA windows.",
         call. = FALSE)
  }
  invisible(TRUE)
}
validate_rf_summary(rf, "Regional")
validate_rf_summary(rf_ncr, "NCR")
for (z in list(rf, rf_ncr)) {
  cta <- z[z$Product_ID %in% c(1L, 3L), , drop = FALSE]
  ot <- z[z$Product_ID %in% c(2L, 4L), , drop = FALSE]
  support_cols <- c("N_Evaluable_Units", "N_Units_With_Any_Detector_Trigger",
                    "N_Units_With_A1A2_True_Trigger", "N_Qualifying_Years",
                    "N_Primary_True_Anchors", "N_False_Comparator_Triggers")
  if (nrow(cta) != 2L || nrow(ot) != 2L ||
      any(vapply(support_cols, function(cc) length(unique(cta[[cc]])) != 1L, logical(1))) ||
      any(vapply(support_cols, function(cc) length(unique(ot[[cc]])) != 1L, logical(1)))) {
    stop("Stage C support parity failed: P1/P3 must share Constant-TA support and P2/P4 must share Outbreak-Threshold support.",
         call. = FALSE)
  }
}

# The full candidate grid must demonstrate that the adopted point estimate was
# selected from an upper-tail (>=), false-alarm-aware search for all four
# products in both scopes.
candidate_grid <- utils::read.csv(
  file.path(stagec, "RF_Operational_Threshold_Candidate_Grid.csv"),
  stringsAsFactors = FALSE
)
need_grid <- c("Scope", "Product_ID", "Anchor", "RF_Scale", "Comparison_Operator",
               "Threshold", "Youden_J", "PPV", "False_Alarm_Rate",
               "False_Alarm_Weeks_per_Evaluable_Unit", "Is_Selected_Point_Estimate")
miss_grid <- setdiff(need_grid, names(candidate_grid))
if (length(miss_grid) || !setequal(unique(as.integer(candidate_grid$Product_ID)), 1:4) ||
    !setequal(unique(as.character(candidate_grid$Scope)), c("NCR", "Regional (pooled 17 regions)")) ||
    any(as.character(candidate_grid$Comparison_Operator) != ">=")) {
  stop("Stage C threshold candidate-grid validation failed: both scopes and Products 1-4 must use the upper-tail >= rule.",
       call. = FALSE)
}
sel_grid <- candidate_grid[as.logical(candidate_grid$Is_Selected_Point_Estimate), , drop = FALSE]
if (nrow(sel_grid) != 8L ||
    any(table(paste(sel_grid$Scope, sel_grid$Product_ID, sep = "|")) != 1L)) {
  stop("Stage C threshold candidate-grid validation failed: exactly one selected point estimate is required for each scope/product.",
       call. = FALSE)
}
if (!all(c("TP", "FP", "FN", "TN", "Sensitivity", "Specificity") %in% names(candidate_grid)) ||
    any(!is.finite(as.numeric(sel_grid$TP))) || any(as.numeric(sel_grid$TP) < 1L) ||
    any(!is.finite(as.numeric(sel_grid$Sensitivity))) || any(as.numeric(sel_grid$Sensitivity) <= 0)) {
  stop("Stage C point-threshold validation failed: every adopted threshold must produce at least one true prospective RF warning.",
       call. = FALSE)
}

# Point-estimate values in the full candidate grid must exactly agree with the
# scope-specific four-product summaries. This prevents a later reporting table
# from silently using a different threshold than the operational handoff.
for (scope_label in c("NCR", "Regional (pooled 17 regions)")) {
  sm <- if (scope_label == "NCR") rf_ncr else rf
  for (pid in 1:4) {
    gg <- sel_grid[sel_grid$Scope == scope_label & sel_grid$Product_ID == pid, , drop = FALSE]
    ss <- sm[sm$Product_ID == pid, , drop = FALSE]
    if (nrow(gg) != 1L || nrow(ss) != 1L ||
        !isTRUE(all.equal(as.numeric(gg$Threshold[1]), as.numeric(ss$Adopted_Threshold[1]), tolerance = 1e-9))) {
      stop("Stage C selected threshold differs between candidate grid and derivation summary for ",
           scope_label, " Product ", pid, ".", call. = FALSE)
    }
  }
}

# Stability bounds must come from complete scope-appropriate re-derivations.
loro <- utils::read.csv(file.path(stagec, "RF_Threshold_LORO_Stability.csv"),
                        stringsAsFactors = FALSE)
loyo <- utils::read.csv(file.path(stagec, "RF_Threshold_LOYO_Stability_NCR.csv"),
                        stringsAsFactors = FALSE)
need_stab <- c("Product_ID", "Threshold", "AUC", "False_Alarm_Rate", "Comparison_Operator")
if (!all(c(need_stab, "Held_Out_Region", "RF_STA", "RF_LTA") %in% names(loro)) ||
    !all(c(need_stab, "Held_Out_Year", "RF_STA", "RF_LTA") %in% names(loyo)) ||
    nrow(loro) != 68L || nrow(loyo) != 24L ||
    any(as.character(loro$Comparison_Operator) != ">=") ||
    any(as.character(loyo$Comparison_Operator) != ">=") ||
    any(!is.finite(as.numeric(loro$Threshold))) ||
    any(!is.finite(as.numeric(loyo$Threshold)))) {
  stop("Stage C stability validation failed: expected 17 LORO and 6 NCR LOYO upper-tail threshold re-derivations for each of four products.",
       call. = FALSE)
}
for (stab in list(loro, loyo)) {
  rt <- stab$Product_ID %in% c(3L, 4L)
  if (any(!is.finite(as.numeric(stab$RF_STA[rt]))) ||
      any(!is.finite(as.numeric(stab$RF_LTA[rt])))) {
    stop("Stage C RF STA/LTA stability validation failed: Products 3/4 must re-select a finite rainfall window in every held-out fold.",
         call. = FALSE)
  }
}

validate_bounds_table <- function(path, scope_label) {
  z <- utils::read.csv(path, stringsAsFactors = FALSE)
  need <- c("Scope", "Product_ID", "Threshold_Lower", "Adopted_Threshold",
            "Threshold_Upper", "Precautionary_Activation_Bound",
            "Strong_Activation_Bound", "Comparison_Operator", "Bounds_Method")
  miss <- setdiff(need, names(z))
  if (length(miss) || nrow(z) != 4L || !setequal(as.integer(z$Product_ID), 1:4) ||
      any(as.character(z$Comparison_Operator) != ">=")) {
    stop("Stage C ", scope_label, " operational RF bounds table is incomplete or not upper-tail.",
         call. = FALSE)
  }
  invisible(TRUE)
}
validate_bounds_table(file.path(stagec, "RF_Threshold_Operational_Bounds.csv"), "Regional")
validate_bounds_table(file.path(stagec, "RF_Threshold_Operational_Bounds_NCR.csv"), "NCR")

validate_bound_perf <- function(path, scope_label) {
  z <- utils::read.csv(path, stringsAsFactors = FALSE)
  need <- c("Scope", "Product_ID", "Bound", "Operational_Role", "Threshold",
            "Sensitivity", "Specificity", "Balanced_Accuracy", "PPV",
            "False_Alarm_Rate", "Trigger_Rate")
  miss <- setdiff(need, names(z))
  if (length(miss) || nrow(z) != 12L ||
      !setequal(unique(as.character(z$Bound)), c("Lower", "Point", "Upper")) ||
      !setequal(unique(as.integer(z$Product_ID)), 1:4)) {
    stop("Stage C ", scope_label,
         " bound sensitivity must contain Lower/Point/Upper rows for all four products, including false-alarm rate.",
         call. = FALSE)
  }
  invisible(TRUE)
}
validate_bound_perf(file.path(stagec, "RF_Threshold_Bound_Sensitivity.csv"), "Regional")
validate_bound_perf(file.path(stagec, "RF_Threshold_Bound_Sensitivity_NCR.csv"), "NCR")

# Direct Stage-C RF mm versus RF STA/LTA head-to-head validation.
h2h_overall <- utils::read.csv(file.path(stagec, "RF_Threshold_HeadToHead_Overall.csv"),
                               stringsAsFactors = FALSE)
need_h2h_overall <- c(
  "Scope", "Anchor", "RF_mm_Product", "RF_Rt_Product", "N_Common_Weeks", "N_Events",
  "N_Positive_Warning_Weeks",
  "AUC_RF_mm", "AUC_RF_Rt", "Delta_AUC_Rt_minus_mm", "p_Paired_DeLong",
  "p_Paired_DeLong_Bonferroni", "AUC_Significant_005",
  "RF_mm_Adopted_Threshold", "RF_Rt_Adopted_Threshold",
  "Sensitivity_RF_mm", "Sensitivity_RF_Rt", "Specificity_RF_mm", "Specificity_RF_Rt",
  "Balanced_Accuracy_RF_mm", "Balanced_Accuracy_RF_Rt",
  "PPV_RF_mm", "PPV_RF_Rt", "False_Alarm_Rate_RF_mm", "False_Alarm_Rate_RF_Rt"
)
miss_h2h_overall <- setdiff(need_h2h_overall, names(h2h_overall))
if (length(miss_h2h_overall) || nrow(h2h_overall) != 4L ||
    !setequal(unique(as.character(h2h_overall$Scope)), c("Pooled 17 regions", "NCR")) ||
    !setequal(unique(as.character(h2h_overall$Anchor)), c("Constant TA", "Outbreak Threshold"))) {
  stop("Stage C validation failed: direct RF head-to-head overall table must contain NCR and pooled 17-region comparisons for both anchors.",
       call. = FALSE)
}
pooled_h2h <- h2h_overall[h2h_overall$Scope == "Pooled 17 regions", , drop = FALSE]
if (nrow(pooled_h2h) != 2L ||
    any(!is.finite(suppressWarnings(as.numeric(pooled_h2h$AUC_RF_mm)))) ||
    any(!is.finite(suppressWarnings(as.numeric(pooled_h2h$AUC_RF_Rt))))) {
  stop("Stage C validation failed: pooled 17-region direct RF head-to-head AUCs must be finite for both anchors and representations.",
       call. = FALSE)
}

h2h_region <- utils::read.csv(file.path(stagec, "RF_Threshold_HeadToHead_ByRegion.csv"),
                              stringsAsFactors = FALSE)
if (!all(c("Anchor", "REGION", "AUC_RF_mm", "AUC_RF_Rt", "Balanced_Accuracy_RF_mm",
           "Balanced_Accuracy_RF_Rt") %in% names(h2h_region)) ||
    nrow(h2h_region) != 34L ||
    !setequal(unique(as.character(h2h_region$REGION)), canonical_17)) {
  stop("Stage C validation failed: direct RF head-to-head regional table must contain all 17 regions for both anchors.",
       call. = FALSE)
}

h2h_tests <- utils::read.csv(file.path(stagec, "RF_Threshold_HeadToHead_ExactPairedTests.csv"),
                             stringsAsFactors = FALSE)
expected_h2h_metrics <- c("AUC", "Sensitivity", "Specificity", "Balanced_Accuracy", "PPV", "False_Alarm_Rate")
if (!all(c("Anchor", "Metric", "N_Paired_Regions", "Mean_Difference_Rt_minus_mm",
           "p_exact_signflip", "p_bonferroni", "Preferred") %in% names(h2h_tests)) ||
    nrow(h2h_tests) != 12L ||
    !setequal(unique(as.character(h2h_tests$Metric)), expected_h2h_metrics) ||
    any(as.integer(h2h_tests$N_Paired_Regions) < 0L) ||
    any(as.integer(h2h_tests$N_Paired_Regions) > 17L)) {
  stop("Stage C validation failed: exact paired RF derivation tests must contain six regional metrics for both anchors, with paired-region counts between 0 and 17 according to event availability.",
       call. = FALSE)
}

h2h_loro <- utils::read.csv(file.path(stagec, "RF_Threshold_HeadToHead_LORO.csv"),
                            stringsAsFactors = FALSE)
if (!all(c("Anchor", "Held_Out_Region", "RF_mm_AUC", "RF_Rt_AUC",
           "AUC_Difference_Rt_minus_mm") %in% names(h2h_loro)) ||
    nrow(h2h_loro) != 34L ||
    !setequal(unique(as.character(h2h_loro$Held_Out_Region)), canonical_17)) {
  stop("Stage C validation failed: RF derivation LORO head-to-head must contain all 17 held-out regions for both anchors.",
       call. = FALSE)
}

h2h_decision <- utils::read.csv(file.path(stagec, "RF_Threshold_HeadToHead_Decision.csv"),
                                stringsAsFactors = FALSE)
if (!all(c("Anchor", "Derivation_Stage_Conclusion", "Significant_Evidence", "Note") %in% names(h2h_decision)) ||
    nrow(h2h_decision) != 2L ||
    !setequal(unique(as.character(h2h_decision$Anchor)), c("Constant TA", "Outbreak Threshold"))) {
  stop("Stage C validation failed: RF derivation head-to-head decision table must contain one conclusion per anchor.",
       call. = FALSE)
}

rt_rows <- rf$Product_ID %in% c(3L, 4L)
if (any(!is.finite(as.numeric(rf$RF_STA[rt_rows]))) ||
    any(!is.finite(as.numeric(rf$RF_LTA[rt_rows]))) ||
    any(!is.finite(as.numeric(rf$RF_Guard[rt_rows])))) {
  stop("Stage C validation failed: Products 3/4 require empirically selected RF STA/LTA windows.",
       call. = FALSE)
}

# Full NCR + Regional reruns for all four environmental products.
# Independent NCR-style replication for all 17 regions: disease-only + P1-P4.
stagei <- file.path(DIR_OUTPUT, "StageI_PerRegion_Independent")
manifest_i <- file.path(stagei, "PerRegion_Independent_Run_Manifest.csv")
assert_files(manifest_i, "Stage I independent per-region analysis")
mi <- utils::read.csv(manifest_i, stringsAsFactors = FALSE)
if (!all(c("REGION", "Product", "Analysis", "Output_Directory", "Completed") %in% names(mi)) ||
    nrow(mi) != 85L || length(unique(mi$REGION)) != 17L ||
    any(table(mi$REGION) != 5L) || any(!as.logical(mi$Completed))) {
  stop("Stage I must contain 17 regions x 5 independent analyses (DiseaseOnly + P1-P4).", call. = FALSE)
}
for (rr in unique(as.character(mi$REGION))) {
  safe_rr <- gsub("[^A-Za-z0-9]+", "_", rr); safe_rr <- gsub("^_+|_+$", "", safe_rr)
  for (product in 0:4) {
    label <- if (product == 0L) "DiseaseOnly" else paste0("RF_Product", product)
    d_i <- file.path(stagei, safe_rr, label)
    common_i <- c(
      pair(d_i, "FigureA2_panel_a_DominanceMatrix"),
      pair(d_i, "FigureA2_panel_b_Scatter"),
      file.path(d_i, "FigureA2_method_summary_primary.csv"),
      file.path(d_i, "evaluation_framework", "Table2_Per_Year_Detail.csv")
    )
    if (product == 0L || product %in% c(1L,3L)) {
      common_i <- c(common_i, pair(d_i, "FigureA3_Comparison_2019_2024_Constant_TA"))
    } else {
      common_i <- c(common_i, pair(d_i, "FigureA3_Comparison_2019_2024_Outbreak_Threshold"))
    }
    if (product > 0L) {
      common_i <- c(common_i, file.path(d_i, "RF_Gate_Metadata.csv"),
                    file.path(d_i, paste0("RF_Operational_Weekly_", safe_rr, ".csv")))
    }
    assert_files(common_i, paste0("Stage I ", rr, " / ", label))
  }
}

for (i in 1:4) {
  validate_stage_a(paste0("_RF_Product", i))
  validate_stage_b(paste0("_RF_Product", i))
  md_a <- utils::read.csv(file.path(DIR_OUTPUT,
    paste0("StageA_NCR_Lepto_Analysis_RF_Product", i), "RF_Gate_Metadata.csv"),
    stringsAsFactors = FALSE)
  md_b <- utils::read.csv(file.path(DIR_OUTPUT,
    paste0("StageB_Regional_Lepto_Analysis_RF_Product", i), "RF_Gate_Metadata.csv"),
    stringsAsFactors = FALSE)
  required_md <- c("RF_Product", "RF_Scale", "Threshold_Scope", "Threshold", "Lower_Bound", "Upper_Bound",
                   "Precautionary_Activation_Bound", "Strong_Activation_Bound",
                   "Primary_Metrics_Threshold", "Operator", "Units")
  if (!all(required_md %in% names(md_a)) || !all(required_md %in% names(md_b))) {
    stop("RF Product ", i, " metadata schema is incomplete.", call. = FALSE)
  }
  if (as.integer(md_a$RF_Product[1]) != i || as.integer(md_b$RF_Product[1]) != i) {
    stop("RF Product ", i, " rerun metadata has the wrong product identifier.", call. = FALSE)
  }
  if (as.character(md_a$Operator[1]) != ">=" ||
      as.character(md_b$Operator[1]) != ">=") {
    stop("RF Product ", i, " rerun metadata must use the upper-tail >= operator.", call. = FALSE)
  }
  expected_op_a <- ">="
  expected_op_b <- ">="
  expected_thr_a <- as.numeric(rf_ncr$Adopted_Threshold[rf_ncr$Product_ID == i])
  expected_thr_b <- as.numeric(rf$Adopted_Threshold[rf$Product_ID == i])
  if (normalize_rf_operator(md_a$Operator[1]) != expected_op_a ||
      normalize_rf_operator(md_b$Operator[1]) != expected_op_b) {
    stop("RF Product ", i, " rerun must use the upper-tail >= operator.", call. = FALSE)
  }
  if (as.character(md_a$Threshold_Scope[1]) != "NCR" ||
      as.character(md_b$Threshold_Scope[1]) != "Regional") {
    stop("RF Product ", i, " rerun threshold scope is incorrect.", call. = FALSE)
  }
  if (!isTRUE(all.equal(as.numeric(md_a$Threshold[1]), expected_thr_a, tolerance = 1e-9)) ||
      !isTRUE(all.equal(as.numeric(md_b$Threshold[1]), expected_thr_b, tolerance = 1e-9))) {
    stop("RF Product ", i,
         " rerun threshold does not match the correct scope-specific Stage C threshold.",
         call. = FALSE)
  }
  expected_lo_a <- as.numeric(rf_ncr$Threshold_Lower[rf_ncr$Product_ID == i])
  expected_hi_a <- as.numeric(rf_ncr$Threshold_Upper[rf_ncr$Product_ID == i])
  expected_prec_a <- as.numeric(rf_ncr$Precautionary_Activation_Bound[rf_ncr$Product_ID == i])
  expected_strong_a <- as.numeric(rf_ncr$Strong_Activation_Bound[rf_ncr$Product_ID == i])
  expected_lo_b <- as.numeric(rf$Threshold_Lower[rf$Product_ID == i])
  expected_hi_b <- as.numeric(rf$Threshold_Upper[rf$Product_ID == i])
  expected_prec_b <- as.numeric(rf$Precautionary_Activation_Bound[rf$Product_ID == i])
  expected_strong_b <- as.numeric(rf$Strong_Activation_Bound[rf$Product_ID == i])
  checks <- list(
    list(md = md_a, lo = expected_lo_a, hi = expected_hi_a,
         prec = expected_prec_a, strong = expected_strong_a, scope = "NCR"),
    list(md = md_b, lo = expected_lo_b, hi = expected_hi_b,
         prec = expected_prec_b, strong = expected_strong_b, scope = "Regional")
  )
  for (cc in checks) {
    md <- cc$md
    if (!isTRUE(all.equal(as.numeric(md$Lower_Bound[1]), cc$lo, tolerance = 1e-9)) ||
        !isTRUE(all.equal(as.numeric(md$Upper_Bound[1]), cc$hi, tolerance = 1e-9)) ||
        !isTRUE(all.equal(as.numeric(md$Precautionary_Activation_Bound[1]), cc$prec, tolerance = 1e-9)) ||
        !isTRUE(all.equal(as.numeric(md$Strong_Activation_Bound[1]), cc$strong, tolerance = 1e-9))) {
      stop("RF Product ", i, " ", cc$scope,
           " uncertainty bounds do not match Stage C.", call. = FALSE)
    }
    if (as.character(md$Primary_Metrics_Threshold[1]) != "Point estimate") {
      stop("RF Product ", i, " must retain the point estimate as the primary headline threshold.",
           call. = FALSE)
    }
  }
  expected_scale <- if (i <= 2L) "mm" else "rt"
  if (tolower(as.character(md_a$RF_Scale[1])) != expected_scale ||
      tolower(as.character(md_b$RF_Scale[1])) != expected_scale) {
    stop("RF Product ", i, " rerun metadata has the wrong RF scale.", call. = FALSE)
  }

  # RF temporal-overlay contract: the RF feature at disease week t is calculated
  # strictly from t-4..t-1 and compared with the fixed point estimate. A1/A2
  # remains the independent true/false alarm reference in downstream metrics.
  a_dir <- file.path(DIR_OUTPUT, paste0("StageA_NCR_Lepto_Analysis_RF_Product", i))
  b_dir <- file.path(DIR_OUTPUT, paste0("StageB_Regional_Lepto_Analysis_RF_Product", i))
  a_overlay_path <- file.path(a_dir, "RF_Trigger_Overlay_NCR.csv")
  b_overlay_path <- file.path(b_dir, "RF_Trigger_Overlay_Regional.csv")
  a_operational_path <- file.path(a_dir, "RF_Operational_Weekly_NCR.csv")
  b_operational_path <- file.path(b_dir, "RF_Operational_Weekly_Regional.csv")
  assert_files(c(a_overlay_path, b_overlay_path, a_operational_path, b_operational_path),
               paste0("RF Product ", i, " temporal/operational outputs"))
  assert_files(unlist(lapply(years_required, function(y)
    pair(b_dir, paste0("FigureB6_RF_Trigger_Overlay_Regional_", y))),
    use.names = FALSE), paste0("RF Product ", i, " regional overlay figures"))
  rf_overlay_regions <- c(
    "BARMM", "CAR", "MIMAROPA", "NCR",
    paste("REGION", c("I", "II", "III", "IV-A", "V", "VI", "VII", "VIII",
                      "IX", "X", "XI", "XII", "XIII"))
  )
  rf_region_stems <- paste0("FigureB6A_RF_Overlay_",
                            gsub("[^A-Za-z0-9]+", "_", rf_overlay_regions))
  assert_files(unlist(lapply(rf_region_stems, function(stem) pair(b_dir, stem)),
                      use.names = FALSE),
               paste0("RF Product ", i, " 17-region RF overlay figures"))

  check_overlay <- function(path, label) {
    z <- utils::read.csv(path, stringsAsFactors = FALSE)
    need <- c("RF_Product", "RF_Scale", "Anchor", "RF_Trigger_WN",
              "Anchor_Trigger_WN", "RF_Warning_Horizon_Start_WN", "RF_Warning_Horizon_End_WN",
              "RF_Window_Start_WN", "RF_Window_End_WN",
              "RF_Feature_At_RF_Trigger", "RF_Max_Preceding_Feature", "RF_Feature_At_Anchor", "RF_Threshold", "RF_Threshold_Lower",
              "RF_Threshold_Upper", "RF_Precautionary_Bound", "RF_Strong_Bound",
              "RF_Activation_Status", "RF_Operator", "RF_Lead_Weeks")
    miss <- setdiff(need, names(z))
    if (length(miss)) {
      stop(label, " is missing overlay column(s): ", paste(miss, collapse = ", "),
           call. = FALSE)
    }
    if (any(as.character(z$RF_Operator[!is.na(z$RF_Operator) & nzchar(as.character(z$RF_Operator))]) != ">=")) {
      stop(label, " contains a non-upper-tail RF operator.", call. = FALSE)
    }
    expected_anchor <- if (i %in% c(1L, 3L)) "Constant TA" else "Outbreak Threshold"
    observed_anchor <- unique(as.character(z$Anchor[!is.na(z$Anchor) & nzchar(as.character(z$Anchor))]))
    if (length(observed_anchor) && any(observed_anchor != expected_anchor)) {
      stop(label, " has the wrong disease anchor. Expected ", expected_anchor,
           "; observed: ", paste(observed_anchor, collapse = ", "), call. = FALSE)
    }
    zlo <- suppressWarnings(as.numeric(z$RF_Threshold_Lower))
    zpt <- suppressWarnings(as.numeric(z$RF_Threshold))
    zhi <- suppressWarnings(as.numeric(z$RF_Threshold_Upper))
    z_has_threshold <- is.finite(zpt)
    if (any(z_has_threshold & (!is.finite(zlo) | !is.finite(zhi) | zlo > zpt | zhi < zpt))) {
      stop(label, " has invalid RF lower/point/upper bounds.", call. = FALSE)
    }
    anchor_finite <- is.finite(suppressWarnings(as.numeric(z$Anchor_Trigger_WN)))
    if (any(anchor_finite)) {
      aw_all <- as.integer(z$Anchor_Trigger_WN[anchor_finite])
      whs <- as.integer(z$RF_Warning_Horizon_Start_WN[anchor_finite])
      whe <- as.integer(z$RF_Warning_Horizon_End_WN[anchor_finite])
      if (any(whs != pmax(1L, aw_all - RF_PRECEDING_WEEKS)) ||
          any(whe != pmax(1L, aw_all - 1L))) {
        stop(label, " has an invalid t-4..t-1 antecedent RF window before the disease anchor.",
             call. = FALSE)
      }
      no_formal <- anchor_finite & !is.finite(suppressWarnings(as.numeric(z$RF_Trigger_WN)))
      if (any(no_formal)) {
        feat_anchor <- suppressWarnings(as.numeric(z$RF_Feature_At_Anchor[no_formal]))
        point <- suppressWarnings(as.numeric(z$RF_Threshold[no_formal]))
        if (any(is.finite(feat_anchor) & is.finite(point) & feat_anchor >= point)) {
          stop(label, " reports no formal RF trigger even though the antecedent t-4..t-1 RF feature at the disease anchor reaches the point threshold.",
               call. = FALSE)
        }
      }
    }
    finite <- anchor_finite &
              is.finite(suppressWarnings(as.numeric(z$RF_Trigger_WN)))
    if (any(finite)) {
      aw <- as.integer(z$Anchor_Trigger_WN[finite])
      rw <- as.integer(z$RF_Trigger_WN[finite])
      ws <- as.integer(z$RF_Window_Start_WN[finite])
      we <- as.integer(z$RF_Window_End_WN[finite])
      lead <- as.integer(z$RF_Lead_Weeks[finite])
      if (any(!is.finite(lead)) || any(lead < 1L | lead > RF_PRECEDING_WEEKS) ||
          any(rw != aw - lead) || any(we != rw) ||
          any(ws != pmax(1L, rw - RF_PRECEDING_WEEKS + 1L))) {
        stop(label, " violates the t-4..t-1 antecedent RF overlay definition.",
             call. = FALSE)
      }
      feat <- suppressWarnings(as.numeric(z$RF_Feature_At_RF_Trigger[finite]))
      thrz <- suppressWarnings(as.numeric(z$RF_Threshold[finite]))
      opz <- as.character(z$RF_Operator[finite])
      ok <- mapply(function(f, th, op) {
        if (!is.finite(f) || !is.finite(th)) return(FALSE)
        isTRUE(apply_rf_gate(f, th, op))
      }, feat, thrz, opz)
      if (any(!ok)) {
        stop(label, " contains an anchor whose preceding RF signal does not pass its adopted upper-tail threshold.",
             call. = FALSE)
      }
      st <- as.character(z$RF_Activation_Status[finite])
      if (any(!is.na(st) & nzchar(st) & !st %in% c("Formal RF trigger", "Strong RF trigger"))) {
        stop(label, " contains a point-gated anchor not classified as a formal/strong RF trigger.",
             call. = FALSE)
      }
    }
    invisible(TRUE)
  }
  check_operational <- function(path, label) {
    z <- utils::read.csv(path, stringsAsFactors = FALSE)
    need <- c("RF_Product", "RF_Scale", "RF_Feature", "Threshold_Lower",
              "Threshold_Point", "Threshold_Upper", "Comparison_Operator",
              "Precautionary_RF_Signal", "Formal_RF_Signal", "Strong_RF_Signal",
              "Antecedent_t4_to_t1_Meets_Point_Threshold", "RF_Activation_Status")
    miss <- setdiff(need, names(z))
    if (length(miss)) {
      stop(label, " is missing operational-bound column(s): ",
           paste(miss, collapse = ", "), call. = FALSE)
    }
    lo <- suppressWarnings(as.numeric(z$Threshold_Lower))
    pt <- suppressWarnings(as.numeric(z$Threshold_Point))
    hi <- suppressWarnings(as.numeric(z$Threshold_Upper))
    if (any(is.finite(pt) & (!is.finite(lo) | !is.finite(hi) | lo > pt | hi < pt))) {
      stop(label, " has invalid lower/point/upper RF bounds.", call. = FALSE)
    }
    valid_status <- c("No RF activation", "Precautionary RF trigger",
                      "Formal RF trigger", "Strong RF trigger")
    st <- as.character(z$RF_Activation_Status)
    st <- st[!is.na(st) & nzchar(st)]
    if (length(st) && any(!st %in% valid_status)) {
      stop(label, " contains an unexpected RF activation status.", call. = FALSE)
    }
    invisible(TRUE)
  }
  check_operational(a_operational_path, paste0("RF Product ", i, " NCR operational table"))
  check_operational(b_operational_path, paste0("RF Product ", i, " Regional operational table"))
  check_overlay(a_overlay_path, paste0("RF Product ", i, " NCR overlay"))
  check_overlay(b_overlay_path, paste0("RF Product ", i, " Regional overlay"))
}

# RF mm versus RF STA/LTA R(t) head-to-head outputs.
stageh <- file.path(DIR_OUTPUT, "StageH_RF_mm_vs_Rt_HeadToHead")
assert_files(c(
  file.path(stageh, "RF_HeadToHead_NCR_Long.csv"),
  file.path(stageh, "RF_HeadToHead_Regional_Long.csv"),
  file.path(stageh, "RF_HeadToHead_NCR_AnchorPaired_ByYear.csv"),
  file.path(stageh, "RF_HeadToHead_NCR_ExactPairedTests.csv"),
  file.path(stageh, "RF_HeadToHead_Regional_AnchorPaired_ByRegion.csv"),
  file.path(stageh, "RF_HeadToHead_Regional_ClusterPairedTests.csv"),
  file.path(stageh, "RF_HeadToHead_Regional_PairedTests.csv"),
  file.path(stageh, "RF_HeadToHead_DominanceProbability.csv"),
  file.path(stageh, "RF_HeadToHead_ThresholdStability.csv"),
  file.path(stageh, "RF_HeadToHead_BoundSensitivity.csv"),
  file.path(stageh, "RF_HeadToHead_OverallDecision.csv"),
  file.path(stageh, "RF_HeadToHead_Summary.csv"),
  pair(stageh, "FigureH_NCR_RFmm_vs_RFRt_Preceding_ConstantTA"),
  pair(stageh, "FigureH_NCR_RFmm_vs_RFRt_Preceding_OutbreakThreshold"),
  pair(stageh, "FigureH_Regional_RFmm_vs_RFRt_Preceding_ConstantTA"),
  pair(stageh, "FigureH_Regional_RFmm_vs_RFRt_Preceding_OutbreakThreshold"),
  pair(stageh, "FigureH_NCR_Preceding_ConstantTA_panel_a_BurdenAccuracy"),
  pair(stageh, "FigureH_NCR_Preceding_ConstantTA_panel_b_Timeliness"),
  pair(stageh, "FigureH_NCR_Preceding_ConstantTA_panel_c_FalseAlarms"),
  pair(stageh, "FigureH_NCR_Preceding_OutbreakThreshold_panel_a_BurdenAccuracy"),
  pair(stageh, "FigureH_NCR_Preceding_OutbreakThreshold_panel_b_Timeliness"),
  pair(stageh, "FigureH_NCR_Preceding_OutbreakThreshold_panel_c_FalseAlarms"),
  pair(stageh, "FigureH_Regional_Preceding_ConstantTA_panel_a_BurdenAccuracy"),
  pair(stageh, "FigureH_Regional_Preceding_ConstantTA_panel_b_Timeliness"),
  pair(stageh, "FigureH_Regional_Preceding_ConstantTA_panel_c_FalseAlarms"),
  pair(stageh, "FigureH_Regional_Preceding_OutbreakThreshold_panel_a_BurdenAccuracy"),
  pair(stageh, "FigureH_Regional_Preceding_OutbreakThreshold_panel_b_Timeliness"),
  pair(stageh, "FigureH_Regional_Preceding_OutbreakThreshold_panel_c_FalseAlarms")
), "Stage H")

hh_reg <- utils::read.csv(file.path(stageh, "RF_HeadToHead_Regional_Long.csv"),
                          stringsAsFactors = FALSE)
hh_ncr <- utils::read.csv(file.path(stageh, "RF_HeadToHead_NCR_Long.csv"),
                          stringsAsFactors = FALSE)
expected_anchors <- c("Preceding Constant TA", "Preceding Outbreak Threshold")
if (!"REGION" %in% names(hh_reg)) stop("Stage H regional comparison lacks REGION.", call. = FALSE)
if (!setequal(unique(as.character(hh_reg$Anchor)), expected_anchors) ||
    !setequal(unique(as.character(hh_ncr$Anchor)), expected_anchors)) {
  stop("Stage H must contain both RF representation comparison anchors.", call. = FALSE)
}
for (a in expected_anchors) {
  rr <- unique(as.character(hh_reg$REGION[hh_reg$Anchor == a]))
  if (length(rr) != 17L) stop("Stage H must compare 17 regions for ", a, ".", call. = FALSE)
}

# Formal anchor-focused RF representation inference must exist for both scales.
hh_ncr_test <- utils::read.csv(file.path(stageh, "RF_HeadToHead_NCR_ExactPairedTests.csv"),
                               stringsAsFactors = FALSE)
hh_reg_test <- utils::read.csv(file.path(stageh, "RF_HeadToHead_Regional_ClusterPairedTests.csv"),
                               stringsAsFactors = FALSE)
need_test_cols <- c("Scope", "Anchor", "Anchor_Method", "Metric", "Domain",
                    "N_Paired_Units", "Mean_RF_mm", "Mean_RF_Rt",
                    "Direction_Adjusted_Improvement", "Difference_CI_Lower",
                    "Difference_CI_Upper", "p_exact_signflip", "p_bonferroni",
                    "Significant_005", "Preferred", "Significant_Favored")
for (obj in list(hh_ncr_test, hh_reg_test)) {
  miss <- setdiff(need_test_cols, names(obj))
  if (length(miss)) {
    stop("Stage H formal paired-inference table is missing column(s): ",
         paste(miss, collapse = ", "), call. = FALSE)
  }
  if (!setequal(unique(as.character(obj$Anchor)), expected_anchors)) {
    stop("Stage H formal paired inference must contain both RF anchors.", call. = FALSE)
  }
  if (any(!is.na(obj$p_bonferroni) & (obj$p_bonferroni < 0 | obj$p_bonferroni > 1))) {
    stop("Stage H adjusted p-values must lie in [0,1].", call. = FALSE)
  }
  for (a in expected_anchors) {
    if (sum(obj$Anchor == a) != 8L) {
      stop("Stage H must test exactly eight prespecified metrics for ", a, ".", call. = FALSE)
    }
  }
}
hh_reg_pairs <- utils::read.csv(file.path(stageh, "RF_HeadToHead_Regional_AnchorPaired_ByRegion.csv"),
                                stringsAsFactors = FALSE)
hh_ncr_pairs <- utils::read.csv(file.path(stageh, "RF_HeadToHead_NCR_AnchorPaired_ByYear.csv"),
                                stringsAsFactors = FALSE)
for (a in expected_anchors) {
  rr <- unique(as.character(hh_reg_pairs$REGION[hh_reg_pairs$Anchor == a]))
  if (length(rr) != 17L) {
    stop("Stage H Regional paired source data must retain all 17 regional clusters for ", a, ".",
         call. = FALSE)
  }
  yy <- unique(hh_ncr_pairs$Unit[hh_ncr_pairs$Anchor == a])
  if (length(yy) < 2L) {
    stop("Stage H NCR paired source data must contain at least two evaluable years for ", a, ".",
         call. = FALSE)
  }
}
if (any(hh_reg_test$N_Paired_Units < 0L | hh_reg_test$N_Paired_Units > 17L, na.rm = TRUE)) {
  stop("Stage H Regional paired-unit counts are outside the valid 0-17 range.", call. = FALSE)
}

hh_dom <- utils::read.csv(file.path(stageh, "RF_HeadToHead_DominanceProbability.csv"),
                          stringsAsFactors = FALSE)
if (nrow(hh_dom) != 4L ||
    !setequal(unique(as.character(hh_dom$Scope)), c("NCR", "Regional")) ||
    !setequal(unique(as.character(hh_dom$Anchor)), expected_anchors) ||
    any(!is.na(hh_dom$Pr_RFRt_Dominates) & (hh_dom$Pr_RFRt_Dominates < 0 | hh_dom$Pr_RFRt_Dominates > 1)) ||
    any(!is.na(hh_dom$Pr_RFmm_Dominates) & (hh_dom$Pr_RFmm_Dominates < 0 | hh_dom$Pr_RFmm_Dominates > 1))) {
  stop("Stage H dominance-probability output is incomplete or invalid.", call. = FALSE)
}

hh_stab <- utils::read.csv(file.path(stageh, "RF_HeadToHead_ThresholdStability.csv"),
                           stringsAsFactors = FALSE)
if (nrow(hh_stab) != 8L ||
    !setequal(unique(as.character(hh_stab$Scope)), c("NCR", "Regional")) ||
    !setequal(unique(as.integer(hh_stab$Product_ID)), 1:4) ||
    any(!is.finite(hh_stab$Point_Estimate)) ||
    any(hh_stab$Lower_Bound > hh_stab$Point_Estimate) ||
    any(hh_stab$Upper_Bound < hh_stab$Point_Estimate) ||
    any(as.character(hh_stab$Comparison_Operator) != ">=")) {
  stop("Stage H threshold-stability comparison must contain valid scope-specific Products 1-4 under the upper-tail >= rule.", call. = FALSE)
}

hh_bound <- utils::read.csv(file.path(stageh, "RF_HeadToHead_BoundSensitivity.csv"),
                            stringsAsFactors = FALSE)
if (nrow(hh_bound) != 12L ||
    !setequal(unique(as.character(hh_bound$Scope)), c("NCR", "Regional")) ||
    !setequal(unique(as.character(hh_bound$Bound)), c("Lower", "Point", "Upper")) ||
    !setequal(unique(as.character(hh_bound$Anchor)), expected_anchors) ||
    !all(c("RF_mm_False_Alarm_Rate", "RF_Rt_False_Alarm_Rate",
           "Preferred_Lower_False_Alarm") %in% names(hh_bound))) {
  stop("Stage H bound-sensitivity comparison must contain scope-specific Lower/Point/Upper results for both anchors, including false-alarm rate.",
       call. = FALSE)
}

hh_decision <- utils::read.csv(file.path(stageh, "RF_HeadToHead_OverallDecision.csv"),
                               stringsAsFactors = FALSE)
allowed_decisions <- c("RF mm preferred", "RF STA/LTA R(t) preferred",
                       "No clear empirical superiority")
if (nrow(hh_decision) != 4L ||
    !setequal(unique(as.character(hh_decision$Scope)), c("NCR", "Regional")) ||
    !setequal(unique(as.character(hh_decision$Anchor)), expected_anchors) ||
    any(!as.character(hh_decision$Overall_Conclusion) %in% allowed_decisions)) {
  stop("Stage H overall empirical-decision table is incomplete or has an invalid conclusion.",
       call. = FALSE)
}

# run_all.R maintains a stage-aware error log from the current clean run.
# Reaching this validator with any recorded stage error would indicate an
# inconsistent execution state and must never be reported as success.
error_log <- file.path(DIR_OUTPUT, "_run_logs", "Errors.csv")
if (!file.exists(error_log)) {
  stop("Run-level error log is missing; run_all.R must orchestrate the full pipeline.",
       call. = FALSE)
}
error_rows <- tryCatch(utils::read.csv(error_log, stringsAsFactors = FALSE),
                       error = function(e) NULL)
if (is.null(error_rows) || !all(c("Stage", "Error") %in% names(error_rows))) {
  stop("Run-level error log has an invalid schema.", call. = FALSE)
}
if (nrow(error_rows) > 0L) {
  stop("Run-level error log is not empty; pipeline cannot be declared successful.",
       call. = FALSE)
}

report_dir <- file.path(DIR_OUTPUT, "_validation")
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(data.frame(
  Check = c("Required outputs", "FigureA1 display years", "NCR six-detector suite",
            "Regional six-detector suite and 17-region map",
            "RF Product 1 mm", "RF Product 2 mm",
            "RF Product 3 empirical R(t)", "RF Product 4 empirical R(t)",
            "RF mm versus R(t) head-to-head"),
  Status = rep("PASS", 9L), stringsAsFactors = FALSE
), file.path(report_dir, "Final_Output_Validation.csv"), row.names = FALSE)
cat("\nFINAL OUTPUT VALIDATION: PASS\n")
