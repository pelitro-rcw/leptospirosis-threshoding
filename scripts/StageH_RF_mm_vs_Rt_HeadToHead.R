# =============================================================================
# STAGE H - FORMAL HEAD-TO-HEAD: RF mm VS EMPIRICAL RF STA/LTA R(t)
# =============================================================================
# Anchor-focused inferential comparisons:
#   P1 (RF mm)  vs P3 (RF R(t)) -> Constant TA
#   P2 (RF mm)  vs P4 (RF R(t)) -> Outbreak Threshold
#
# NCR inference is paired by evaluable year using exact paired sign-flip tests.
# Regional inference is paired by region (17 regional clusters), using exact
# paired sign-flip tests on region-level summaries. Bootstrap percentile CIs
# preserve the paired unit. Bonferroni adjustment is performed separately
# within each anchor and spatial scale across the eight prespecified metrics.
#
# The point estimate remains the primary RF threshold. Lower/upper empirical
# stability bounds are evaluated as a threshold-sensitivity diagnostic and are
# NOT substituted for the primary point-estimate head-to-head metrics.
# =============================================================================

REQUIRED_PACKAGES <- c("dplyr", "tidyr", "readr", "ggplot2", "patchwork")
root <- Sys.getenv("TA_PROJECT_ROOT", unset = "")
if (!nzchar(root)) root <- normalizePath(getwd(), winslash = "/")
source(file.path(root, "R", "00_config.R"))
require_packages(REQUIRED_PACKAGES, "Stage H RF head-to-head")
invisible(lapply(REQUIRED_PACKAGES, function(p) {
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}))
source(file.path(DIR_R, "01_publication_theme.R"), local = TRUE)

OUT_DIR <- file.path(DIR_OUTPUT, "StageH_RF_mm_vs_Rt_HeadToHead")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

METHOD_ORDER <- c(
  "Constant Transmission Acceleration",
  "Continuous Transmission Acceleration",
  "Outbreak Threshold",
  "Alarm Threshold",
  "WHO 75th Percentile Threshold",
  "WHO 90th Percentile Threshold"
)
METHOD_SHORT <- c(
  "Constant Transmission Acceleration" = "Constant TA",
  "Continuous Transmission Acceleration" = "Continuous TA",
  "Outbreak Threshold" = "Outbreak",
  "Alarm Threshold" = "Alarm",
  "WHO 75th Percentile Threshold" = "WHO 75th",
  "WHO 90th Percentile Threshold" = "WHO 90th"
)

METRIC_SPECS <- data.frame(
  Metric = c("TAM", "N_True_Alarms", "PPV", "Sensitivity",
             "Mean_Lead_Time", "WP", "ALY", "N_False_Alarms"),
  Domain = c(rep("Burden Accuracy", 4), rep("Timeliness", 3), "False Alarms"),
  Metric_Label = c("TAM", "True alarms", "PPV", "Sensitivity",
                   "Mean lead time", "Warning persistence", "ALY", "False alarms"),
  Higher_Is_Better = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE),
  NCR_Column = c("TAM_yr", "n_True", "PPV_yr", "Sensitivity_yr",
                 "MLT_yr_wks", "WP_yr_wks", "ALY_yr", "n_False"),
  Regional_Column = c("TAM", "N_True_Alarms", "PPV", "Sensitivity",
                      "Mean_Lead_Time", "WP", "ALY", "N_False_Alarms"),
  stringsAsFactors = FALSE
)

PAIR_CONFIG <- list(
  list(
    anchor = "Preceding Constant TA",
    anchor_short = "Constant TA",
    disease_method = "Constant Transmission Acceleration",
    ncr_detector = "Constant_TA",
    mm_suffix = "_RF_Product1", rt_suffix = "_RF_Product3",
    mm_product = 1L, rt_product = 3L,
    stem = "Preceding_ConstantTA"
  ),
  list(
    anchor = "Preceding Outbreak Threshold",
    anchor_short = "Outbreak Threshold",
    disease_method = "Outbreak Threshold",
    ncr_detector = "Outbreak_Threshold",
    mm_suffix = "_RF_Product2", rt_suffix = "_RF_Product4",
    mm_product = 2L, rt_product = 4L,
    stem = "Preceding_OutbreakThreshold"
  )
)

# -----------------------------------------------------------------------------
# Utilities
# -----------------------------------------------------------------------------
safe_mean <- function(x) {
  x <- suppressWarnings(as.numeric(x)); x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else mean(x)
}

safe_median <- function(x) {
  x <- suppressWarnings(as.numeric(x)); x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else stats::median(x)
}

seed_from_text <- function(x, offset = 0L) {
  z <- utf8ToInt(enc2utf8(paste(x, collapse = "|")))
  as.integer((sum(z * seq_along(z)) + as.integer(offset)) %% 2000000000L + 1001L)
}

paired_boot_ci <- function(d, B = 5000L, conf = 0.95, seed = 1L) {
  d <- suppressWarnings(as.numeric(d)); d <- d[is.finite(d)]
  n <- length(d)
  if (!n) return(c(NA_real_, NA_real_))
  if (n == 1L) return(c(d, d))
  set.seed(seed)
  reps <- replicate(B, mean(d[sample.int(n, n, replace = TRUE)]))
  alpha <- (1 - conf) / 2
  as.numeric(stats::quantile(reps, probs = c(alpha, 1 - alpha),
                             na.rm = TRUE, names = FALSE, type = 8))
}

exact_signflip_p <- function(d, exact_max_n = 20L, B = 100000L, seed = 1L) {
  d <- suppressWarnings(as.numeric(d)); d <- d[is.finite(d)]
  n <- length(d)
  if (!n) return(NA_real_)
  obs <- abs(mean(d))
  if (all(abs(d) <= sqrt(.Machine$double.eps))) return(1)
  tol <- sqrt(.Machine$double.eps)
  if (n <= exact_max_n) {
    # 2^17 = 131072 at the Regional scale, which is small enough for an exact
    # paired randomization distribution while treating region as the cluster.
    idx <- 0:(2^n - 1)
    perm_stat <- numeric(length(idx))
    for (j in seq_len(n)) {
      signs_j <- ifelse(bitwAnd(idx, bitwShiftL(1L, j - 1L)) == 0L, -1, 1)
      perm_stat <- perm_stat + signs_j * d[j]
    }
    perm_stat <- abs(perm_stat / n)
    return(mean(perm_stat >= (obs - tol)))
  }
  set.seed(seed)
  perm_stat <- replicate(B, abs(mean(d * sample(c(-1, 1), n, replace = TRUE))))
  (1 + sum(perm_stat >= (obs - tol))) / (B + 1)
}

paired_test_table <- function(dat, scope) {
  # dat: one row per paired unit x metric with RF_mm / RF_Rt.
  out <- dat |>
    dplyr::group_by(Anchor, Anchor_Method, Metric, Domain, Metric_Label, Higher_Is_Better) |>
    dplyr::group_modify(~ {
      z <- .x |>
        dplyr::filter(is.finite(RF_mm), is.finite(RF_Rt))
      raw_d <- z$RF_Rt - z$RF_mm
      adj_d <- if (isTRUE(.y$Higher_Is_Better[1])) raw_d else -raw_d
      seed <- seed_from_text(c(scope, .y$Anchor, .y$Metric), 700L)
      ci_raw <- paired_boot_ci(raw_d, B = 5000L, seed = seed)
      ci_adj <- if (isTRUE(.y$Higher_Is_Better[1])) ci_raw else -rev(ci_raw)
      p <- exact_signflip_p(raw_d, seed = seed + 17L)
      sd_d <- stats::sd(adj_d, na.rm = TRUE)
      std_eff <- if (is.finite(sd_d) && sd_d > 0) mean(adj_d) / sd_d else NA_real_
      data.frame(
        Scope = scope,
        N_Paired_Units = nrow(z),
        Mean_RF_mm = safe_mean(z$RF_mm),
        Mean_RF_Rt = safe_mean(z$RF_Rt),
        Mean_Difference_Rt_minus_mm = safe_mean(raw_d),
        Difference_CI_Lower = ci_raw[1],
        Difference_CI_Upper = ci_raw[2],
        Direction_Adjusted_Improvement = safe_mean(adj_d),
        Direction_Adjusted_CI_Lower = ci_adj[1],
        Direction_Adjusted_CI_Upper = ci_adj[2],
        Paired_Standardized_Effect = std_eff,
        p_exact_signflip = p,
        stringsAsFactors = FALSE
      )
    }) |>
    dplyr::ungroup() |>
    dplyr::group_by(Anchor) |>
    dplyr::mutate(
      p_bonferroni = pmin(1, p_exact_signflip * sum(is.finite(p_exact_signflip))),
      Significant_005 = is.finite(p_bonferroni) & p_bonferroni < 0.05,
      Preferred = dplyr::case_when(
        !is.finite(Direction_Adjusted_Improvement) ~ NA_character_,
        Direction_Adjusted_Improvement > 0 ~ "RF STA/LTA R(t)",
        Direction_Adjusted_Improvement < 0 ~ "RF mm",
        TRUE ~ "Tie"
      ),
      Significant_Favored = dplyr::case_when(
        !Significant_005 ~ "Not significant",
        Direction_Adjusted_Improvement > 0 ~ "RF STA/LTA R(t)",
        Direction_Adjusted_Improvement < 0 ~ "RF mm",
        TRUE ~ "Tie"
      )
    ) |>
    dplyr::ungroup()
  out
}

# -----------------------------------------------------------------------------
# Existing six-detector descriptive product contrasts (retained for context)
# -----------------------------------------------------------------------------
read_ncr_summary <- function(suffix) {
  path <- file.path(DIR_OUTPUT, paste0("StageA_NCR_Lepto_Analysis", suffix),
                    "FigureA2_method_summary_primary.csv")
  if (!file.exists(path)) stop("Missing NCR RF product summary: ", path, call. = FALSE)
  x <- utils::read.csv(path, stringsAsFactors = FALSE)
  expected <- c("Method", "TAM", "N_True_Alarms_yr", "PPV", "Sensitivity",
                "Mean_Lead_Time", "WP", "ALY", "N_False_Alarms_yr")
  if (!all(expected %in% names(x))) {
    stop("NCR RF product summary has an unexpected schema: ", path, call. = FALSE)
  }
  x
}

read_regional_summary <- function(suffix) {
  path <- file.path(DIR_OUTPUT, paste0("StageB_Regional_Lepto_Analysis", suffix),
                    "StageB_Regional_8Metric_Summary.csv")
  if (!file.exists(path)) stop("Missing Regional RF product summary: ", path, call. = FALSE)
  x <- utils::read.csv(path, stringsAsFactors = FALSE)
  expected <- c("REGION", "Method", METRIC_SPECS$Regional_Column)
  if (!all(expected %in% names(x))) {
    stop("Regional RF product summary has an unexpected schema: ", path, call. = FALSE)
  }
  x
}

make_ncr_descriptive_long <- function(cfg) {
  cols <- c("TAM", "N_True_Alarms_yr", "PPV", "Sensitivity",
            "Mean_Lead_Time", "WP", "ALY", "N_False_Alarms_yr")
  mm <- read_ncr_summary(cfg$mm_suffix) |>
    dplyr::filter(Method %in% METHOD_ORDER) |>
    dplyr::select(Method, dplyr::all_of(cols)) |>
    tidyr::pivot_longer(-Method, names_to = "Metric_Source", values_to = "RF_mm")
  rt <- read_ncr_summary(cfg$rt_suffix) |>
    dplyr::filter(Method %in% METHOD_ORDER) |>
    dplyr::select(Method, dplyr::all_of(cols)) |>
    tidyr::pivot_longer(-Method, names_to = "Metric_Source", values_to = "RF_Rt")
  metric_lookup <- data.frame(
    Metric_Source = cols,
    Metric = METRIC_SPECS$Metric,
    Domain = METRIC_SPECS$Domain,
    Metric_Label = METRIC_SPECS$Metric_Label,
    Higher_Is_Better = METRIC_SPECS$Higher_Is_Better,
    stringsAsFactors = FALSE)
  mm |>
    dplyr::inner_join(rt, by = c("Method", "Metric_Source")) |>
    dplyr::left_join(metric_lookup, by = "Metric_Source") |>
    dplyr::mutate(
      Anchor = cfg$anchor,
      RF_mm_Product = cfg$mm_product,
      RF_Rt_Product = cfg$rt_product,
      Difference_Rt_minus_mm = RF_Rt - RF_mm,
      Direction_Adjusted_Improvement = ifelse(Higher_Is_Better,
                                               Difference_Rt_minus_mm,
                                               -Difference_Rt_minus_mm),
      Preferred = dplyr::case_when(
        Direction_Adjusted_Improvement > 0 ~ "RF STA/LTA R(t)",
        Direction_Adjusted_Improvement < 0 ~ "RF mm",
        TRUE ~ "Tie"))
}

make_regional_descriptive_long <- function(cfg) {
  mm <- read_regional_summary(cfg$mm_suffix) |>
    dplyr::filter(Method %in% METHOD_ORDER) |>
    dplyr::select(REGION, Method, dplyr::all_of(METRIC_SPECS$Regional_Column)) |>
    tidyr::pivot_longer(-c(REGION, Method), names_to = "Regional_Column", values_to = "RF_mm")
  rt <- read_regional_summary(cfg$rt_suffix) |>
    dplyr::filter(Method %in% METHOD_ORDER) |>
    dplyr::select(REGION, Method, dplyr::all_of(METRIC_SPECS$Regional_Column)) |>
    tidyr::pivot_longer(-c(REGION, Method), names_to = "Regional_Column", values_to = "RF_Rt")
  lookup <- METRIC_SPECS |>
    dplyr::select(Regional_Column, Metric, Domain, Metric_Label, Higher_Is_Better)
  mm |>
    dplyr::inner_join(rt, by = c("REGION", "Method", "Regional_Column")) |>
    dplyr::left_join(lookup, by = "Regional_Column") |>
    dplyr::mutate(
      Anchor = cfg$anchor,
      RF_mm_Product = cfg$mm_product,
      RF_Rt_Product = cfg$rt_product,
      Difference_Rt_minus_mm = RF_Rt - RF_mm,
      Direction_Adjusted_Improvement = ifelse(Higher_Is_Better,
                                               Difference_Rt_minus_mm,
                                               -Difference_Rt_minus_mm),
      Preferred = dplyr::case_when(
        Direction_Adjusted_Improvement > 0 ~ "RF STA/LTA R(t)",
        Direction_Adjusted_Improvement < 0 ~ "RF mm",
        TRUE ~ "Tie"))
}

ncr_long <- dplyr::bind_rows(lapply(PAIR_CONFIG, make_ncr_descriptive_long))
reg_long <- dplyr::bind_rows(lapply(PAIR_CONFIG, make_regional_descriptive_long))
readr::write_csv(ncr_long, file.path(OUT_DIR, "RF_HeadToHead_NCR_Long.csv"), na = "")
readr::write_csv(reg_long, file.path(OUT_DIR, "RF_HeadToHead_Regional_Long.csv"), na = "")

for (a in unique(reg_long$Anchor)) {
  rr <- unique(as.character(reg_long$REGION[reg_long$Anchor == a]))
  if (length(rr) != 17L) {
    stop("RF regional head-to-head must contain 17 regions for ", a,
         "; found ", length(rr), ".", call. = FALSE)
  }
}

# -----------------------------------------------------------------------------
# NCR anchor-focused paired YEAR data and exact inference
# -----------------------------------------------------------------------------
read_ncr_yearly <- function(suffix, detector_code) {
  path <- file.path(DIR_OUTPUT, paste0("StageA_NCR_Lepto_Analysis", suffix),
                    "evaluation_framework", "Table2A_Per_Year_Detail_With_Compartments.csv")
  if (!file.exists(path)) stop("Missing NCR per-year RF product table: ", path, call. = FALSE)
  x <- utils::read.csv(path, stringsAsFactors = FALSE)
  required <- c("Year", "Detector", "TAM_yr", "n_True", "PPV_yr",
                "First_A1_True_Week", "MLT_yr_wks", "WP_yr_wks", "ALY_yr", "n_False")
  if (!all(required %in% names(x))) {
    stop("NCR per-year RF product table has an unexpected schema: ", path, call. = FALSE)
  }
  x |>
    dplyr::filter(Detector == detector_code) |>
    dplyr::mutate(Sensitivity_yr = as.numeric(is.finite(First_A1_True_Week)))
}

make_ncr_anchor_paired <- function(cfg) {
  mm0 <- read_ncr_yearly(cfg$mm_suffix, cfg$ncr_detector)
  rt0 <- read_ncr_yearly(cfg$rt_suffix, cfg$ncr_detector)
  rows <- lapply(seq_len(nrow(METRIC_SPECS)), function(i) {
    sp <- METRIC_SPECS[i, ]
    mm <- mm0 |>
      dplyr::transmute(Unit = as.integer(Year), RF_mm = suppressWarnings(as.numeric(.data[[sp$NCR_Column]])))
    rt <- rt0 |>
      dplyr::transmute(Unit = as.integer(Year), RF_Rt = suppressWarnings(as.numeric(.data[[sp$NCR_Column]])))
    dplyr::inner_join(mm, rt, by = "Unit") |>
      dplyr::mutate(
        Anchor = cfg$anchor,
        Anchor_Method = cfg$disease_method,
        Metric = sp$Metric,
        Domain = sp$Domain,
        Metric_Label = sp$Metric_Label,
        Higher_Is_Better = sp$Higher_Is_Better,
        RF_mm_Product = cfg$mm_product,
        RF_Rt_Product = cfg$rt_product,
        Difference_Rt_minus_mm = RF_Rt - RF_mm,
        Direction_Adjusted_Improvement = ifelse(Higher_Is_Better,
                                                 Difference_Rt_minus_mm,
                                                 -Difference_Rt_minus_mm)
      )
  })
  dplyr::bind_rows(rows)
}

ncr_anchor_paired <- dplyr::bind_rows(lapply(PAIR_CONFIG, make_ncr_anchor_paired))
ncr_tests <- paired_test_table(ncr_anchor_paired, "NCR")
readr::write_csv(ncr_anchor_paired,
                 file.path(OUT_DIR, "RF_HeadToHead_NCR_AnchorPaired_ByYear.csv"), na = "")
readr::write_csv(ncr_tests,
                 file.path(OUT_DIR, "RF_HeadToHead_NCR_ExactPairedTests.csv"), na = "")

# -----------------------------------------------------------------------------
# Regional anchor-focused paired REGION data and cluster-level exact inference
# -----------------------------------------------------------------------------
make_regional_anchor_paired <- function(cfg) {
  mm0 <- read_regional_summary(cfg$mm_suffix) |>
    dplyr::filter(Method == cfg$disease_method)
  rt0 <- read_regional_summary(cfg$rt_suffix) |>
    dplyr::filter(Method == cfg$disease_method)
  rows <- lapply(seq_len(nrow(METRIC_SPECS)), function(i) {
    sp <- METRIC_SPECS[i, ]
    mm <- mm0 |>
      dplyr::transmute(Unit = as.character(REGION), RF_mm = suppressWarnings(as.numeric(.data[[sp$Regional_Column]])))
    rt <- rt0 |>
      dplyr::transmute(Unit = as.character(REGION), RF_Rt = suppressWarnings(as.numeric(.data[[sp$Regional_Column]])))
    dplyr::inner_join(mm, rt, by = "Unit") |>
      dplyr::mutate(
        REGION = Unit,
        Anchor = cfg$anchor,
        Anchor_Method = cfg$disease_method,
        Metric = sp$Metric,
        Domain = sp$Domain,
        Metric_Label = sp$Metric_Label,
        Higher_Is_Better = sp$Higher_Is_Better,
        RF_mm_Product = cfg$mm_product,
        RF_Rt_Product = cfg$rt_product,
        Difference_Rt_minus_mm = RF_Rt - RF_mm,
        Direction_Adjusted_Improvement = ifelse(Higher_Is_Better,
                                                 Difference_Rt_minus_mm,
                                                 -Difference_Rt_minus_mm)
      )
  })
  dplyr::bind_rows(rows)
}

regional_anchor_paired <- dplyr::bind_rows(lapply(PAIR_CONFIG, make_regional_anchor_paired))
for (a in unique(regional_anchor_paired$Anchor)) {
  rr <- unique(regional_anchor_paired$REGION[regional_anchor_paired$Anchor == a])
  if (length(rr) != 17L) {
    stop("Anchor-focused Regional RF comparison must contain all 17 regions for ", a,
         "; found ", length(rr), ".", call. = FALSE)
  }
}
regional_tests <- paired_test_table(regional_anchor_paired, "Regional")
readr::write_csv(regional_anchor_paired,
                 file.path(OUT_DIR, "RF_HeadToHead_Regional_AnchorPaired_ByRegion.csv"), na = "")
readr::write_csv(regional_tests,
                 file.path(OUT_DIR, "RF_HeadToHead_Regional_ClusterPairedTests.csv"), na = "")
# Backward-compatible filename retained, now containing the stronger exact
# region-cluster paired inference rather than the old asymptotic Wilcoxon test.
readr::write_csv(regional_tests,
                 file.path(OUT_DIR, "RF_HeadToHead_Regional_PairedTests.csv"), na = "")

# -----------------------------------------------------------------------------
# Dominance probability across the eight prespecified metrics
# -----------------------------------------------------------------------------
compute_dominance <- function(dat, scope, B = 5000L) {
  dplyr::bind_rows(lapply(unique(dat$Anchor), function(anchor_i) {
    z <- dat |>
      dplyr::filter(Anchor == anchor_i, is.finite(Direction_Adjusted_Improvement)) |>
      dplyr::group_by(Metric) |>
      dplyr::mutate(
        .scale = stats::sd(Direction_Adjusted_Improvement, na.rm = TRUE),
        .scale = ifelse(!is.finite(.scale) | .scale <= sqrt(.Machine$double.eps),
                        max(abs(Direction_Adjusted_Improvement), na.rm = TRUE), .scale),
        .scale = ifelse(!is.finite(.scale) | .scale <= sqrt(.Machine$double.eps), 1, .scale),
        Standardized_Improvement = Direction_Adjusted_Improvement / .scale
      ) |>
      dplyr::ungroup() |>
      dplyr::group_by(Unit) |>
      dplyr::summarise(Unit_Score = mean(Standardized_Improvement, na.rm = TRUE), .groups = "drop") |>
      dplyr::filter(is.finite(Unit_Score))
    n <- nrow(z)
    if (!n) {
      return(data.frame(Scope = scope, Anchor = anchor_i, N_Units = 0L,
                        Observed_Dominance_Score = NA_real_,
                        Pr_RFRt_Dominates = NA_real_, Pr_RFmm_Dominates = NA_real_,
                        Bootstrap_CI_Lower = NA_real_, Bootstrap_CI_Upper = NA_real_,
                        stringsAsFactors = FALSE))
    }
    set.seed(seed_from_text(c(scope, anchor_i, "dominance"), 910L))
    boot <- replicate(B, mean(z$Unit_Score[sample.int(n, n, replace = TRUE)]))
    ci <- as.numeric(stats::quantile(boot, c(.025, .975), na.rm = TRUE, names = FALSE, type = 8))
    data.frame(
      Scope = scope, Anchor = anchor_i, N_Units = n,
      Observed_Dominance_Score = mean(z$Unit_Score),
      Pr_RFRt_Dominates = mean(boot > 0),
      Pr_RFmm_Dominates = mean(boot < 0),
      Bootstrap_CI_Lower = ci[1], Bootstrap_CI_Upper = ci[2],
      stringsAsFactors = FALSE)
  }))
}

dominance_tbl <- dplyr::bind_rows(
  compute_dominance(ncr_anchor_paired, "NCR"),
  compute_dominance(regional_anchor_paired, "Regional")
)
readr::write_csv(dominance_tbl,
                 file.path(OUT_DIR, "RF_HeadToHead_DominanceProbability.csv"), na = "")

# -----------------------------------------------------------------------------
# Fixed-threshold stability and lower/point/upper sensitivity diagnostics
# -----------------------------------------------------------------------------
stagec <- file.path(DIR_OUTPUT, "StageC_RF_Threshold_Derivation")
summary_reg_path <- file.path(stagec, "RF_Threshold_Derivation_Summary.csv")
summary_ncr_path <- file.path(stagec, "RF_Threshold_Derivation_Summary_NCR.csv")
loro_path <- file.path(stagec, "RF_Threshold_LORO_Stability.csv")
loyo_path <- file.path(stagec, "RF_Threshold_LOYO_Stability_NCR.csv")
bound_reg_path <- file.path(stagec, "RF_Threshold_Bound_Sensitivity.csv")
bound_ncr_path <- file.path(stagec, "RF_Threshold_Bound_Sensitivity_NCR.csv")
need_stagec <- c(summary_reg_path, summary_ncr_path, loro_path, loyo_path,
                 bound_reg_path, bound_ncr_path)
if (!all(file.exists(need_stagec))) {
  stop("Stage H requires scope-specific Stage C threshold summaries, stability files, and bound-sensitivity files.",
       call. = FALSE)
}
rf_summary_reg <- utils::read.csv(summary_reg_path, stringsAsFactors = FALSE)
rf_summary_ncr <- utils::read.csv(summary_ncr_path, stringsAsFactors = FALSE)
loro <- utils::read.csv(loro_path, stringsAsFactors = FALSE)
loyo <- utils::read.csv(loyo_path, stringsAsFactors = FALSE)
bound_reg <- utils::read.csv(bound_reg_path, stringsAsFactors = FALSE)
bound_ncr <- utils::read.csv(bound_ncr_path, stringsAsFactors = FALSE)
need_sum <- c("Product_ID", "RF_Scale", "Adopted_Threshold", "Threshold_Lower",
              "Threshold_Upper", "Comparison_Operator", "N_Stability_Thresholds")
for (ss in list(rf_summary_reg, rf_summary_ncr)) {
  if (!all(need_sum %in% names(ss))) {
    stop("Stage C RF threshold summary lacks columns required for head-to-head stability analysis.",
         call. = FALSE)
  }
  if (any(as.character(ss$Comparison_Operator) != ">=")) {
    stop("Stage H requires the current upper-tail RF threshold framework (>=).", call. = FALSE)
  }
}

stability_scope <- function(summary_tbl, stability_tbl, scope, stability_type) {
  dplyr::bind_rows(lapply(1:4, function(pid) {
    ss <- summary_tbl |> dplyr::filter(Product_ID == pid) |> dplyr::slice(1)
    zz <- stability_tbl |> dplyr::filter(Product_ID == pid, is.finite(Threshold))
    point <- as.numeric(ss$Adopted_Threshold[1])
    lo <- as.numeric(ss$Threshold_Lower[1])
    hi <- as.numeric(ss$Threshold_Upper[1])
    data.frame(
      Scope = scope,
      Product_ID = pid,
      RF_Scale = as.character(ss$RF_Scale[1]),
      Point_Estimate = point,
      Lower_Bound = lo,
      Upper_Bound = hi,
      Absolute_Bound_Width = hi - lo,
      Relative_Bound_Width = (hi - lo) / max(abs(point), sqrt(.Machine$double.eps)),
      Comparison_Operator = ">=",
      Stability_Type = stability_type,
      N_Stability = nrow(zz),
      Stability_Direction_Consistency = 1,
      Stability_Threshold_Median = safe_median(zz$Threshold),
      Stability_Threshold_IQR = if (nrow(zz) >= 2L) stats::IQR(zz$Threshold, na.rm = TRUE, type = 8) else NA_real_,
      # Legacy names retained for downstream compatibility.
      N_LORO = nrow(zz),
      LORO_Direction_Stability = 1,
      LORO_Threshold_Median = safe_median(zz$Threshold),
      LORO_Threshold_IQR = if (nrow(zz) >= 2L) stats::IQR(zz$Threshold, na.rm = TRUE, type = 8) else NA_real_,
      stringsAsFactors = FALSE)
  }))
}

stability_by_product <- dplyr::bind_rows(
  stability_scope(rf_summary_ncr, loyo, "NCR", "LOYO by year"),
  stability_scope(rf_summary_reg, loro, "Regional", "LORO by region")
)
threshold_stability <- dplyr::bind_rows(lapply(PAIR_CONFIG, function(cfg) {
  stability_by_product |>
    dplyr::filter(Product_ID %in% c(cfg$mm_product, cfg$rt_product)) |>
    dplyr::mutate(
      Anchor = cfg$anchor,
      Representation = ifelse(Product_ID == cfg$mm_product, "RF mm", "RF STA/LTA R(t)")
    )
}))
readr::write_csv(threshold_stability,
                 file.path(OUT_DIR, "RF_HeadToHead_ThresholdStability.csv"), na = "")

bound_all <- dplyr::bind_rows(
  bound_ncr |> dplyr::mutate(Scope = "NCR"),
  bound_reg |> dplyr::mutate(Scope = "Regional")
)
bound_sensitivity <- dplyr::bind_rows(lapply(PAIR_CONFIG, function(cfg) {
  dplyr::bind_rows(lapply(c("NCR", "Regional"), function(scope_i) {
    bb <- bound_all |>
      dplyr::filter(Scope == scope_i,
                    Product_ID %in% c(cfg$mm_product, cfg$rt_product)) |>
      dplyr::mutate(
        Anchor = cfg$anchor,
        Representation = ifelse(Product_ID == cfg$mm_product, "RF mm", "RF STA/LTA R(t)"),
        Balanced_Accuracy = ifelse(is.finite(Balanced_Accuracy), Balanced_Accuracy,
                                   rowMeans(cbind(Sensitivity, Specificity), na.rm = FALSE))
      )
    mm <- bb |>
      dplyr::filter(Product_ID == cfg$mm_product) |>
      dplyr::select(Bound, RF_mm_Balanced_Accuracy = Balanced_Accuracy,
                    RF_mm_Sensitivity = Sensitivity, RF_mm_Specificity = Specificity,
                    RF_mm_False_Alarm_Rate = False_Alarm_Rate,
                    RF_mm_Trigger_Rate = Trigger_Rate)
    rt <- bb |>
      dplyr::filter(Product_ID == cfg$rt_product) |>
      dplyr::select(Bound, RF_Rt_Balanced_Accuracy = Balanced_Accuracy,
                    RF_Rt_Sensitivity = Sensitivity, RF_Rt_Specificity = Specificity,
                    RF_Rt_False_Alarm_Rate = False_Alarm_Rate,
                    RF_Rt_Trigger_Rate = Trigger_Rate)
    dplyr::inner_join(mm, rt, by = "Bound") |>
      dplyr::mutate(
        Scope = scope_i,
        Anchor = cfg$anchor,
        Preferred_Balanced_Accuracy = dplyr::case_when(
          RF_Rt_Balanced_Accuracy > RF_mm_Balanced_Accuracy ~ "RF STA/LTA R(t)",
          RF_Rt_Balanced_Accuracy < RF_mm_Balanced_Accuracy ~ "RF mm",
          TRUE ~ "Tie"),
        Preferred_Lower_False_Alarm = dplyr::case_when(
          RF_Rt_False_Alarm_Rate < RF_mm_False_Alarm_Rate ~ "RF STA/LTA R(t)",
          RF_Rt_False_Alarm_Rate > RF_mm_False_Alarm_Rate ~ "RF mm",
          TRUE ~ "Tie")
      )
  }))
}))
readr::write_csv(bound_sensitivity,
                 file.path(OUT_DIR, "RF_HeadToHead_BoundSensitivity.csv"), na = "")

# -----------------------------------------------------------------------------
# Overall empirical decision: winner, no forced winner when evidence is mixed.
# -----------------------------------------------------------------------------
all_tests <- dplyr::bind_rows(ncr_tests, regional_tests)
make_decision <- function(scope, anchor) {
  tt <- all_tests |> dplyr::filter(Scope == scope, Anchor == anchor)
  dd <- dominance_tbl |> dplyr::filter(Scope == scope, Anchor == anchor) |> dplyr::slice(1)
  bs <- bound_sensitivity |> dplyr::filter(Scope == scope, Anchor == anchor)
  cfg <- PAIR_CONFIG[[which(vapply(PAIR_CONFIG, function(x) identical(x$anchor, anchor), logical(1)))[1]]]
  st <- threshold_stability |> dplyr::filter(Scope == scope, Anchor == anchor)
  mm_st <- st |> dplyr::filter(Product_ID == cfg$mm_product) |> dplyr::slice(1)
  rt_st <- st |> dplyr::filter(Product_ID == cfg$rt_product) |> dplyr::slice(1)

  wins_rt <- sum(tt$Preferred == "RF STA/LTA R(t)", na.rm = TRUE)
  wins_mm <- sum(tt$Preferred == "RF mm", na.rm = TRUE)
  sig_rt <- sum(tt$Significant_Favored == "RF STA/LTA R(t)", na.rm = TRUE)
  sig_mm <- sum(tt$Significant_Favored == "RF mm", na.rm = TRUE)
  bound_rt <- sum(bs$Preferred_Balanced_Accuracy == "RF STA/LTA R(t)", na.rm = TRUE)
  bound_mm <- sum(bs$Preferred_Balanced_Accuracy == "RF mm", na.rm = TRUE)
  pr_rt <- if (nrow(dd)) dd$Pr_RFRt_Dominates[1] else NA_real_
  pr_mm <- if (nrow(dd)) dd$Pr_RFmm_Dominates[1] else NA_real_

  # Pre-specified, transparent rule. Statistical significance is not required
  # in every metric; a winner must show majority directional superiority, no
  # Bonferroni-significant deterioration, and either inferential support or a
  # dominance probability >= 0.75. Bound sensitivity must not favor the other
  # representation in at least two of the three Lower/Point/Upper checks.
  rt_preferred <- wins_rt >= 5L && sig_mm == 0L &&
    ((sig_rt >= 1L) || (is.finite(pr_rt) && pr_rt >= 0.75)) && bound_mm < 2L
  mm_preferred <- wins_mm >= 5L && sig_rt == 0L &&
    ((sig_mm >= 1L) || (is.finite(pr_mm) && pr_mm >= 0.75)) && bound_rt < 2L
  decision <- if (rt_preferred && !mm_preferred) "RF STA/LTA R(t) preferred" else
    if (mm_preferred && !rt_preferred) "RF mm preferred" else
      "No clear empirical superiority"

  data.frame(
    Scope = scope,
    Anchor = anchor,
    RF_mm_Product = cfg$mm_product,
    RF_Rt_Product = cfg$rt_product,
    Metrics_Favoring_RFmm = wins_mm,
    Metrics_Favoring_RFRt = wins_rt,
    Significant_Metrics_Favoring_RFmm = sig_mm,
    Significant_Metrics_Favoring_RFRt = sig_rt,
    Pr_RFmm_Dominates = pr_mm,
    Pr_RFRt_Dominates = pr_rt,
    Bounds_Favoring_RFmm = bound_mm,
    Bounds_Favoring_RFRt = bound_rt,
    RFmm_Relative_Bound_Width = if (nrow(mm_st)) mm_st$Relative_Bound_Width[1] else NA_real_,
    RFRt_Relative_Bound_Width = if (nrow(rt_st)) rt_st$Relative_Bound_Width[1] else NA_real_,
    RFmm_LORO_Direction_Stability = if (nrow(mm_st)) mm_st$LORO_Direction_Stability[1] else NA_real_,
    RFRt_LORO_Direction_Stability = if (nrow(rt_st)) rt_st$LORO_Direction_Stability[1] else NA_real_,
    Overall_Conclusion = decision,
    stringsAsFactors = FALSE)
}

overall_decision <- dplyr::bind_rows(lapply(c("NCR", "Regional"), function(scope_i) {
  dplyr::bind_rows(lapply(vapply(PAIR_CONFIG, function(x) x$anchor, character(1)), function(anchor_i) {
    make_decision(scope_i, anchor_i)
  }))
}))
readr::write_csv(overall_decision,
                 file.path(OUT_DIR, "RF_HeadToHead_OverallDecision.csv"), na = "")

summary_table <- all_tests |>
  dplyr::group_by(Scope, Anchor, Domain) |>
  dplyr::summarise(
    Comparisons = dplyr::n(),
    RFRt_Better = sum(Preferred == "RF STA/LTA R(t)", na.rm = TRUE),
    RFmm_Better = sum(Preferred == "RF mm", na.rm = TRUE),
    Significant_RFRt = sum(Significant_Favored == "RF STA/LTA R(t)", na.rm = TRUE),
    Significant_RFmm = sum(Significant_Favored == "RF mm", na.rm = TRUE),
    .groups = "drop")
readr::write_csv(summary_table, file.path(OUT_DIR, "RF_HeadToHead_Summary.csv"), na = "")

# -----------------------------------------------------------------------------
# Anchor-focused A-C figures: same reporting domains as the current pipeline.
# No legend is needed because the x-axis explicitly names the two RF products.
# -----------------------------------------------------------------------------
plot_domain <- function(test_tbl, anchor, domain, scope) {
  dd <- test_tbl |>
    dplyr::filter(Anchor == anchor, Domain == domain) |>
    dplyr::select(Metric, Metric_Label, Mean_RF_mm, Mean_RF_Rt,
                  p_bonferroni, Significant_005) |>
    tidyr::pivot_longer(c(Mean_RF_mm, Mean_RF_Rt),
                        names_to = "Representation", values_to = "Value") |>
    dplyr::mutate(
      Representation = dplyr::recode(Representation,
                                      Mean_RF_mm = "RF mm",
                                      Mean_RF_Rt = "RF STA/LTA R(t)"),
      Representation = factor(Representation, levels = c("RF mm", "RF STA/LTA R(t)")))
  lab <- test_tbl |>
    dplyr::filter(Anchor == anchor, Domain == domain) |>
    dplyr::transmute(Metric_Label,
                     p_label = ifelse(is.finite(p_bonferroni),
                                      paste0("p(adj) = ", format.pval(p_bonferroni, digits = 2, eps = 0.001)),
                                      "p(adj) = NA"))
  ggplot2::ggplot(dd, ggplot2::aes(x = Representation, y = Value, group = Metric_Label)) +
    ggplot2::geom_line(linewidth = 0.55, colour = "grey45") +
    ggplot2::geom_point(size = 2.2) +
    ggplot2::facet_wrap(~ Metric_Label, scales = "free_y", nrow = 1) +
    ggplot2::geom_text(data = lab,
                       ggplot2::aes(x = 1.5, y = Inf, label = p_label),
                       inherit.aes = FALSE, vjust = 1.3, size = 2.5) +
    ggplot2::labs(x = NULL, y = ifelse(scope == "NCR", "Paired-year mean", "Paired-region mean"),
                  title = domain) +
    theme_pub() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
      legend.position = "none",
      plot.subtitle = ggplot2::element_blank(),
      plot.caption = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(10, 14, 12, 12))
}

build_three_domain <- function(test_tbl, anchor, scope) {
  p1 <- plot_domain(test_tbl, anchor, "Burden Accuracy", scope) +
    ggplot2::labs(title = "a    Burden Accuracy")
  p2 <- plot_domain(test_tbl, anchor, "Timeliness", scope) +
    ggplot2::labs(title = "b    Timeliness")
  p3 <- plot_domain(test_tbl, anchor, "False Alarms", scope) +
    ggplot2::labs(title = "c    False Alarms")
  p1 / p2 / p3
}

panel_specs <- list(
  list(domain = "Burden Accuracy", letter = "a", stem = "BurdenAccuracy"),
  list(domain = "Timeliness", letter = "b", stem = "Timeliness"),
  list(domain = "False Alarms", letter = "c", stem = "FalseAlarms")
)
for (cfg in PAIR_CONFIG) {
  for (sp in panel_specs) {
    pn <- plot_domain(ncr_tests, cfg$anchor, sp$domain, "NCR") +
      ggplot2::labs(title = paste0(sp$letter, "    ", sp$domain))
    save_pub(paste0("FigureH_NCR_", cfg$stem, "_panel_", sp$letter, "_", sp$stem),
             pn, NC_W_DOUBLE, ifelse(sp$domain == "False Alarms", 4.2, 4.8), OUT_DIR)
    pr <- plot_domain(regional_tests, cfg$anchor, sp$domain, "Regional") +
      ggplot2::labs(title = paste0(sp$letter, "    ", sp$domain))
    save_pub(paste0("FigureH_Regional_", cfg$stem, "_panel_", sp$letter, "_", sp$stem),
             pr, NC_W_DOUBLE, ifelse(sp$domain == "False Alarms", 4.2, 4.8), OUT_DIR)
  }
  save_pub(paste0("FigureH_NCR_RFmm_vs_RFRt_", cfg$stem),
           build_three_domain(ncr_tests, cfg$anchor, "NCR"), NC_W_DOUBLE, 9.0, OUT_DIR)
  save_pub(paste0("FigureH_Regional_RFmm_vs_RFRt_", cfg$stem),
           build_three_domain(regional_tests, cfg$anchor, "Regional"), NC_W_DOUBLE, 9.0, OUT_DIR)
}

writeLines(c(
  "Formal RF environmental-representation head-to-head comparison.",
  "P1 (RF mm) versus P3 (empirical RF STA/LTA R(t)) is evaluated on the Constant TA anchor.",
  "P2 (RF mm) versus P4 (empirical RF STA/LTA R(t)) is evaluated on the Outbreak Threshold anchor.",
  "NCR inference is paired by evaluable year and uses exact paired sign-flip tests.",
  "Regional inference treats region as the cluster: the 17 region-level paired summaries are tested with exact paired sign-flip tests.",
  "Paired 95% confidence intervals are percentile bootstrap intervals obtained by resampling the paired year or region units.",
  "Eight prespecified metrics are tested; p-values are Bonferroni-adjusted separately within each anchor and spatial scale.",
  "Higher is favorable for TAM, true alarms, PPV, sensitivity, mean lead time, warning persistence and ALY; lower is favorable for false alarms.",
  "Dominance probability bootstraps the paired unit after standardizing direction-adjusted metric differences.",
  "Lower/point/upper RF thresholds are compared using Stage C sensitivity/specificity as a stability diagnostic only; the point estimate remains the primary operational threshold for headline metrics.",
  "Overall preference is not forced. RF mm or RF STA/LTA R(t) is selected only when a prespecified majority, no-significant-deterioration, inferential/dominance, and bound-robustness rule is satisfied; otherwise the conclusion is No clear empirical superiority."
), file.path(OUT_DIR, "StageH_Methods_Note.txt"))

cat("Stage H complete: ", OUT_DIR, "\n", sep = "")
