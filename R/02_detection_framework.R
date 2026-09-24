# =============================================================================
# SHARED DETECTION-EVALUATION FRAMEWORK
# -----------------------------------------------------------------------------
# Single source of truth for the outbreak-detection evaluation used by BOTH
# Stage 3 (primary analysis) and Stage 6 (benchmarking against contemporary
# methods).
#
# The functions below were extracted VERBATIM from Stage3_QC_analysis.R. Stage 6
# sources this file so that every benchmarked method is scored with exactly the
# same anchor definitions, trigger classification and per-year metrics as the
# manuscript's primary detectors -- which is what makes the comparison fair.
#
# Anchors:
#   A1 = actionable early-warning window, peak-A1_LEAD_MAX .. peak-A1_LEAD_MIN
#   A2 = epidemic burden block, smallest contiguous run containing the peak
#        whose cumulative cases reach A2_BURDEN_FRAC of the annual total
#   A trigger is a True Alarm iff it falls in A1 OR A2.
#
# Required in the calling environment before sourcing:
#   A1_LEAD_MIN, A1_LEAD_MAX, A2_BURDEN_FRAC
#
# NOTE ON DUPLICATION: Stage 3 still carries its own copies of these functions.
# They are byte-identical today (this file was generated from them), but if
# Stage 3's definitions are ever edited, this file must be regenerated or the
# two will silently diverge. See the validation report.
# =============================================================================

compute_A1 <- function(df_in, year, lead_min = A1_LEAD_MIN, lead_max = A1_LEAD_MAX) {
  df_y <- df_in |> dplyr::filter(YR == year)
  if (nrow(df_y) == 0 || all(is.na(df_y$DC_QC)))
    return(list(peak_week = NA_integer_, A1_weeks = integer(0)))
  peak_idx <- which.max(df_y$DC_QC); peak_week <- df_y$WN[peak_idx]
  A1_weeks <- (peak_week - lead_max):(peak_week - lead_min)
  A1_weeks <- A1_weeks[A1_weeks >= 1L]
  list(peak_week = as.integer(peak_week), A1_weeks = as.integer(A1_weeks))
}

compute_A2 <- function(df_in, year, burden_frac = A2_BURDEN_FRAC) {
  df_y <- df_in |> dplyr::filter(YR == year) |> dplyr::arrange(WN)
  if (nrow(df_y) == 0 || all(is.na(df_y$DC_QC)))
    return(list(start_week = NA_integer_, end_week = NA_integer_, A2_weeks = integer(0)))
  cases <- ifelse(is.na(df_y$DC_QC), 0, df_y$DC_QC)
  weeks <- df_y$WN; n_w <- length(weeks); total <- sum(cases)
  if (total <= 0)
    return(list(start_week = NA_integer_, end_week = NA_integer_, A2_weeks = integer(0)))
  peak_idx <- which.max(cases); lo <- hi <- peak_idx; S <- cases[peak_idx]
  while (S / total < burden_frac && (lo > 1L || hi < n_w)) {
    L <- if (lo > 1L)  cases[lo - 1L] else -Inf
    R <- if (hi < n_w) cases[hi + 1L] else -Inf
    if (L >= R) { lo <- lo - 1L; S <- S + L } else { hi <- hi + 1L; S <- S + R }
  }
  list(start_week = as.integer(weeks[lo]), end_week = as.integer(weeks[hi]),
       A2_weeks = as.integer(weeks[lo:hi]))
}

compute_anchors_for_year <- function(df_in, year,
                                     lead_min = A1_LEAD_MIN, lead_max = A1_LEAD_MAX,
                                     burden_frac = A2_BURDEN_FRAC) {
  A1 <- compute_A1(df_in, year, lead_min, lead_max)
  A2 <- compute_A2(df_in, year, burden_frac)
  list(year = year, peak_week = A1$peak_week,
       A1_weeks = A1$A1_weeks, A2_weeks = A2$A2_weeks)
}

classify_trigger <- function(trigger_week, A1_weeks, A2_weeks) {
  in_A1 <- trigger_week %in% A1_weeks
  in_A2 <- trigger_week %in% A2_weeks
  is_true <- in_A1 || in_A2
  list(is_true = is_true, in_A1 = in_A1, in_A2 = in_A2)
}

evaluate_detector <- function(df_in, detector_col, evaluable_years,
                              lead_min = A1_LEAD_MIN, lead_max = A1_LEAD_MAX,
                              burden_frac = A2_BURDEN_FRAC) {
  per_trigger_rows <- list(); per_year_rows <- list()
  for (yr in evaluable_years) {
    anchors <- compute_anchors_for_year(df_in, yr,
                                        lead_min = lead_min, lead_max = lead_max,
                                        burden_frac = burden_frac)
    df_y <- df_in |> dplyr::filter(YR == yr) |> dplyr::arrange(WN)
    trig_idx <- which(df_y[[detector_col]] == 1L)
    year_true_count <- 0L; year_false_count <- 0L
    year_first_a1_true_week <- NA_integer_
    if (length(trig_idx) > 0L) {
      for (k in trig_idx) {
        wk <- df_y$WN[k]
        cls <- classify_trigger(wk, anchors$A1_weeks, anchors$A2_weeks)
        # Route attribute (route classifier preserved as a diagnostic field)
        route <- if (!cls$is_true) "FalseAlarm"
        else if (cls$in_A1 && cls$in_A2) "Both"
        else if (cls$in_A1) "Route1_EarlyWarning"
        else "Route2_WithinEpidemic"
        per_trigger_rows[[length(per_trigger_rows) + 1L]] <- data.frame(
          Year = yr, Week = as.integer(wk), Detector = detector_col,
          IsTrue = cls$is_true, Route = route,
          InA1 = cls$in_A1, InA2 = cls$in_A2,
          stringsAsFactors = FALSE)
        if (cls$is_true) {
          year_true_count <- year_true_count + 1L
          if (cls$in_A1 && (is.na(year_first_a1_true_week) || wk < year_first_a1_true_week))
            year_first_a1_true_week <- as.integer(wk)
        } else { year_false_count <- year_false_count + 1L }
      }
    }
    lead_time <- if (!is.na(year_first_a1_true_week) && !is.na(anchors$peak_week))
      as.integer(anchors$peak_week - year_first_a1_true_week) else NA_integer_
    per_year_rows[[length(per_year_rows) + 1L]] <- data.frame(
      Year = yr, Detector = detector_col, Peak_Week = anchors$peak_week,
      A1_start = if (length(anchors$A1_weeks)) min(anchors$A1_weeks) else NA_integer_,
      A1_end   = if (length(anchors$A1_weeks)) max(anchors$A1_weeks) else NA_integer_,
      A2_start = if (length(anchors$A2_weeks)) min(anchors$A2_weeks) else NA_integer_,
      A2_end   = if (length(anchors$A2_weeks)) max(anchors$A2_weeks) else NA_integer_,
      n_triggers = as.integer(length(trig_idx)),
      n_True = year_true_count, n_False = year_false_count,
      First_A1_True_Week = year_first_a1_true_week,
      Lead_Time_Weeks = lead_time, stringsAsFactors = FALSE)
  }
  per_trigger_df <- if (length(per_trigger_rows)) dplyr::bind_rows(per_trigger_rows) else
    data.frame(Year=integer(),Week=integer(),Detector=character(),IsTrue=logical(),
               Route=character(),InA1=logical(),InA2=logical(),
               stringsAsFactors=FALSE)
  per_year_df <- dplyr::bind_rows(per_year_rows)
  total_triggers <- sum(per_year_df$n_triggers, na.rm = TRUE)
  total_true     <- sum(per_year_df$n_True,     na.rm = TRUE)
  ppv <- if (total_triggers > 0) total_true / total_triggers else NA_real_
  evaluable_n  <- length(evaluable_years)
  # Sensitivity: A1-restricted (years with at least one True alarm in A1).
  years_with_A1 <- sum(!is.na(per_year_df$First_A1_True_Week))
  sensitivity   <- if (evaluable_n > 0) years_with_A1 / evaluable_n else NA_real_
  years_with_T  <- sum(per_year_df$n_True > 0L)
  # Same-denominator MLT: zero-coerce non-firing years.
  lead_zerocoerced <- ifelse(is.na(per_year_df$Lead_Time_Weeks), 0,
                             per_year_df$Lead_Time_Weeks)
  mean_lead             <- mean(lead_zerocoerced)
  mean_lead_conditional <- if (any(!is.na(per_year_df$Lead_Time_Weeks)))
    mean(per_year_df$Lead_Time_Weeks, na.rm = TRUE)
  else NA_real_
  if (is.nan(mean_lead))             mean_lead             <- NA_real_
  if (is.nan(mean_lead_conditional)) mean_lead_conditional <- NA_real_
  summary_row <- data.frame(
    Detector = detector_col, Total_Triggers = total_triggers,
    True_Alarms = total_true, False_Alarms = total_triggers - total_true,
    PPV = round(ppv, 3),
    Years_with_A1_True = years_with_A1,
    Years_with_True_Any = years_with_T,
    Years_Evaluable = evaluable_n,
    Sensitivity = round(sensitivity, 3),
    Mean_Lead_Time_wks = round(mean_lead, 2),
    Mean_Lead_Time_wks_conditional = round(mean_lead_conditional, 2),
    stringsAsFactors = FALSE)
  list(per_trigger = per_trigger_df, per_year = per_year_df, summary = summary_row)
}


# -----------------------------------------------------------------------------
# EIGHT-METRIC PER-YEAR COMPARTMENT MACHINERY
# -----------------------------------------------------------------------------
# Extracted verbatim from Stage3_QC_analysis.R (one documented parameterisation,
# noted below). These produce the eight study metrics used throughout the QC,
# Regional and Country analyses:
#   TAM, N_True_Alarms, PPV, Sensitivity, MLT, WP, ALY, N_False_Alarms
# Requires COMPARTMENT_ACTIONABLE_MIN / _MAX in the calling environment.

# Fallback bounds. A stage normally sets these to its own A1 lead window before
# sourcing this file; these defaults exist so classify_compartment() cannot fail
# with "object 'COMPARTMENT_ACTIONABLE_MIN' not found" if a stage forgets. They
# match the A1 window used throughout the manuscript (peak-8 .. peak-4).
if (!exists("COMPARTMENT_ACTIONABLE_MIN", inherits = TRUE)) {
  COMPARTMENT_ACTIONABLE_MIN <- 4L
}
if (!exists("COMPARTMENT_ACTIONABLE_MAX", inherits = TRUE)) {
  COMPARTMENT_ACTIONABLE_MAX <- 8L
}

classify_compartment <- function(lead_time) {
  if (is.na(lead_time))                                return(NA_character_)
  if (lead_time >= COMPARTMENT_ACTIONABLE_MIN &&
      lead_time <= COMPARTMENT_ACTIONABLE_MAX)         return("Actionable")
  return("Reactive")
}

augment_compartments <- function(trigger_df, per_year_df) {
  if (nrow(trigger_df) == 0) {
    trigger_df$Lead_Time   <- integer(0)
    trigger_df$Compartment <- character(0)
    return(trigger_df)
  }
  peak_lookup <- per_year_df |>
    dplyr::select(Year, Detector, Peak_Week) |> dplyr::distinct()
  trigger_df |>
    dplyr::left_join(peak_lookup, by = c("Year", "Detector")) |>
    dplyr::mutate(
      Lead_Time   = as.integer(Peak_Week - Week),
      Compartment = ifelse(
        !is.na(IsTrue) & IsTrue,
        vapply(Lead_Time, classify_compartment, character(1)),
        NA_character_
      )
    )
}

# NOTE (framework module): the original Stage 3 definition referenced the global
# `df` when joining weekly case counts for TAM. It is parameterised here as
# `case_df` (defaulting to `df`, so Stage 3 behaviour is byte-identical) purely
# so Stage 6 can pass its own equivalently-structured frame. No formula changed.
compute_per_year_compartments <- function(trigger_df_aug, per_year_df,
                                          case_df = df) {
  if (nrow(trigger_df_aug) == 0) {
    return(per_year_df |> dplyr::mutate(
      PPV_yr = NA_real_,
      n_Actionable_yr = 0L,
      n_Reactive_yr = 0L,
      n_TrueActionable_yr = 0L,
      Lead_Compartment = NA_character_,
      ALY_yr = 0,                ALY_yr_conditional = NA_real_,
      WP_yr_wks = 0,             WP_yr_wks_conditional = NA_real_,
      MLT_yr_wks = 0,            MLT_yr_wks_conditional = NA_real_,
      TAM_yr = NA_real_
    ))
  }
  per_year_compartments <- trigger_df_aug |>
    dplyr::group_by(Year, Detector) |>
    dplyr::summarise(
      n_Actionable_yr     = sum(Compartment == "Actionable", na.rm = TRUE),
      n_Reactive_yr       = sum(Compartment == "Reactive",   na.rm = TRUE),
      n_TrueActionable_yr = sum(IsTrue & Compartment == "Actionable", na.rm = TRUE),
      WP_yr_wks_conditional = ifelse(
        sum(IsTrue & Compartment == "Actionable", na.rm = TRUE) > 0,
        mean(Lead_Time[IsTrue & Compartment == "Actionable"], na.rm = TRUE),
        NA_real_
      ),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      WP_yr_wks_conditional = round(WP_yr_wks_conditional, 2),
      WP_yr_wks = ifelse(is.na(WP_yr_wks_conditional), 0, WP_yr_wks_conditional)
    )
  
  trigger_df_with_dc <- trigger_df_aug |>
    dplyr::left_join(
      case_df |> dplyr::transmute(Year = YR, Week = WN, DC = DC_QC),
      by = c("Year", "Week")
    )
  tam_per_year <- trigger_df_with_dc |>
    dplyr::filter(IsTrue) |>
    dplyr::group_by(Year, Detector) |>
    dplyr::summarise(TAM_yr = sum(DC, na.rm = TRUE), .groups = "drop")
  
  per_year_df |>
    dplyr::mutate(
      PPV_yr = ifelse(n_triggers > 0, round(n_True / n_triggers, 3), NA_real_),
      Lead_Compartment = vapply(Lead_Time_Weeks, classify_compartment, character(1))
    ) |>
    dplyr::left_join(per_year_compartments, by = c("Year", "Detector")) |>
    dplyr::left_join(tam_per_year,           by = c("Year", "Detector")) |>
    dplyr::mutate(
      n_Actionable_yr     = ifelse(is.na(n_Actionable_yr),     0L, n_Actionable_yr),
      n_Reactive_yr       = ifelse(is.na(n_Reactive_yr),       0L, n_Reactive_yr),
      n_TrueActionable_yr = ifelse(is.na(n_TrueActionable_yr), 0L, n_TrueActionable_yr),
      TAM_yr = ifelse(is.na(TAM_yr), 0, TAM_yr),
      WP_yr_wks             = ifelse(is.na(WP_yr_wks), 0, WP_yr_wks),
      WP_yr_wks_conditional = WP_yr_wks_conditional,
      ALY_yr_conditional = ifelse(n_True > 0,
                                  round(n_TrueActionable_yr / n_True, 3),
                                  NA_real_),
      ALY_yr             = ifelse(n_True > 0,
                                  round(n_TrueActionable_yr / n_True, 3),
                                  0),
      MLT_yr_wks_conditional = Lead_Time_Weeks,
      MLT_yr_wks             = ifelse(is.na(Lead_Time_Weeks), 0, Lead_Time_Weeks)
    )
}


# -----------------------------------------------------------------------------
# CONSTANT TA DETECTOR (Vaezi-style STA/LTA with hysteresis and frozen LTA)
# -----------------------------------------------------------------------------
# Extracted VERBATIM from Stage3_QC_analysis.R section B.8 and wrapped as a
# function so Stages 3-5 and Stage 6 all produce byte-identical triggers.
#
# The four features that make this detector what it is -- and that a naive
# STA/LTA reimplementation silently omits:
#   * LTA_WIN = 26 weeks, not 12.
#   * GUARD = 2: the LTA window ENDS two weeks before the STA window begins, so
#     the short and long windows never overlap.
#   * frozen_lta: once the detector latches ON, the long-term baseline is FROZEN
#     and stops updating. Without this the ongoing outbreak inflates its own
#     baseline, the ratio collapses, and the detector switches off early.
#   * MIN_OFF_RESET = 8: the baseline only unfreezes after 8 consecutive OFF
#     weeks, so brief dips inside an epidemic do not reset it.
# Plus a full state reset at each year boundary.
#
# @param dc  numeric weekly case counts, ordered by year then week
# @param yr  integer year for each element of dc
# @return list(trigger = integer 0/1, R = numeric STA/LTA ratio)
CONSTANT_TA_STA_WIN       <- 4L
CONSTANT_TA_LTA_WIN       <- 26L
CONSTANT_TA_GUARD         <- 2L
CONSTANT_TA_MIN_OFF_RESET <- 8L

build_constant_ta <- function(dc, yr, eta_on, eta_off,
                              STA_WIN = CONSTANT_TA_STA_WIN,
                              LTA_WIN = CONSTANT_TA_LTA_WIN,
                              GUARD = CONSTANT_TA_GUARD,
                              MIN_OFF_RESET = CONSTANT_TA_MIN_OFF_RESET) {
  MIN_T   <- STA_WIN + GUARD + LTA_WIN
  n_rows  <- length(dc)
  R_vaezi_v   <- rep(NA_real_, n_rows)
  triggered_v <- rep(FALSE, n_rows)
  is_on <- FALSE; frozen_lta <- NA_real_; consec_off <- 0L

  for (t in seq_len(n_rows)) {
    if (t > 1 && !is.na(yr[t]) && !is.na(yr[t - 1]) && yr[t] != yr[t - 1]) {
      is_on <- FALSE; frozen_lta <- NA_real_; consec_off <- 0L
    }
    if (!is_on) consec_off <- consec_off + 1L else consec_off <- 0L
    if (!is_on && consec_off >= MIN_OFF_RESET) frozen_lta <- NA_real_
    if (t < MIN_T) next
    sta_vals <- dc[(t - STA_WIN + 1):t]
    sta <- if (all(is.na(sta_vals))) NA_real_ else mean(sta_vals, na.rm = TRUE)
    if (!is_on || is.na(frozen_lta)) {
      lta_idx <- (t - STA_WIN - GUARD - LTA_WIN + 1):(t - STA_WIN - GUARD)
      if (length(lta_idx) == LTA_WIN && min(lta_idx) > 0) {
        lta_vals <- dc[lta_idx]
        frozen_lta <- if (all(is.na(lta_vals))) NA_real_ else
          mean(lta_vals, na.rm = TRUE)
      } else {
        frozen_lta <- NA_real_
      }
    }
    R_t <- if (!is.na(frozen_lta) && frozen_lta > 0 && !is.na(sta))
      sta / frozen_lta else NA_real_
    R_vaezi_v[t] <- R_t
    if (!is_on && !is.na(R_t) && R_t >= eta_on) { is_on <- TRUE; consec_off <- 0L }
    if (is_on && !is.na(R_t) && R_t < eta_off)  { is_on <- FALSE; frozen_lta <- NA_real_ }
    triggered_v[t] <- is_on
  }
  list(trigger = as.integer(triggered_v), R = R_vaezi_v)
}


# -----------------------------------------------------------------------------
# PAIRED WILCOXON SIGNED-RANK FOR HEAD-TO-HEAD COMPARISONS
# -----------------------------------------------------------------------------
# Consolidates the identical logic used in Stage 3 (pairing by YEAR) and
# Stages 4 and 5 (pairing by REGION / COUNTRY), so Stage 6 reports the same
# statistics computed the same way.
#
# Guards match the existing stages: the test runs only with at least 3 pairs and
# at least one non-zero difference, otherwise V and p are NA rather than a
# spurious value. `exact = FALSE` matches Stages 3-5 (normal approximation with
# continuity correction), which is what makes ties tractable.
#
# p_bonferroni is retained as an alias of the raw p-value: the study applies a
# per-pairwise alpha rather than a family-wise correction, and the column name is
# kept only for backward compatibility with the existing CSV schema.
HH_ALPHA_DEFAULT <- 0.05

hh_metric_category <- c(
  TAM            = "Epidemic_Burden_and_Alarm_Accuracy",
  N_True_Alarms  = "Epidemic_Burden_and_Alarm_Accuracy",
  PPV            = "Epidemic_Burden_and_Alarm_Accuracy",
  Sensitivity    = "Epidemic_Burden_and_Alarm_Accuracy",
  MLT            = "Timeliness",
  WP             = "Timeliness",
  ALY            = "Timeliness",
  N_False_Alarms = "False_Alarms"
)

hh_better_when <- function(metric) {
  if (identical(metric, "N_False_Alarms")) "lower" else "higher"
}

#' Cluster bootstrap CI of the median paired difference.
#' Mirrors hh_bootstrap_ci() in Stage 3, including its fixed seed.
hh_bootstrap_ci <- function(diff_vec, B = 1000L, alpha = 0.05,
                            seed = 20260101L) {
  diff_vec <- diff_vec[!is.na(diff_vec)]
  if (length(diff_vec) < 2L) return(c(lo = NA_real_, hi = NA_real_))
  set.seed(seed)
  reps <- vapply(seq_len(B), function(b) {
    idx <- sample(seq_along(diff_vec), size = length(diff_vec), replace = TRUE)
    stats::median(diff_vec[idx], na.rm = TRUE)
  }, numeric(1))
  c(lo = unname(stats::quantile(reps, alpha / 2,     na.rm = TRUE, type = 7)),
    hi = unname(stats::quantile(reps, 1 - alpha / 2, na.rm = TRUE, type = 7)))
}

#' Paired Wilcoxon signed-rank plus the accompanying summary statistics.
#'
#' @param a,b Aligned numeric vectors (already matched on the pairing key).
#' @return one-row data frame with V_statistic, p_value, p_pairwise,
#'   p_bonferroni, Significant_005, N_Pairs, Median_Diff_AminusB,
#'   Abs_Median_Diff, Bootstrap_CI_lo, Bootstrap_CI_hi.
hh_paired_wilcoxon <- function(a, b, boot_B = 1000L,
                               alpha = HH_ALPHA_DEFAULT) {
  keep <- is.finite(a) & is.finite(b)
  a <- a[keep]; b <- b[keep]
  n_pairs <- length(a)
  diffs <- a - b
  median_diff <- if (n_pairs > 0L) stats::median(diffs, na.rm = TRUE) else NA_real_

  if (n_pairs >= 3L && any(diffs != 0, na.rm = TRUE)) {
    wt <- tryCatch(suppressWarnings(stats::wilcox.test(
      a, b, paired = TRUE, alternative = "two.sided", exact = FALSE)),
      error = function(e) NULL)
    v_stat <- if (!is.null(wt)) as.numeric(wt$statistic) else NA_real_
    p_val  <- if (!is.null(wt)) as.numeric(wt$p.value)   else NA_real_
  } else {
    v_stat <- NA_real_; p_val <- NA_real_
  }

  ci <- hh_bootstrap_ci(diffs, B = boot_B, alpha = alpha)

  data.frame(
    N_Pairs             = n_pairs,
    Median_Diff_AminusB = median_diff,
    Abs_Median_Diff     = if (is.na(median_diff)) NA_real_ else abs(median_diff),
    Bootstrap_CI_lo     = unname(ci["lo"]),
    Bootstrap_CI_hi     = unname(ci["hi"]),
    V_statistic         = v_stat,
    p_value             = p_val,
    p_pairwise          = p_val,
    p_bonferroni        = p_val,
    Significant_005     = if (is.na(p_val)) NA else p_val < alpha,
    stringsAsFactors    = FALSE
  )
}


# -----------------------------------------------------------------------------
# COMPARATOR DETECTORS: FARRINGTON, EARS, EWARS
# -----------------------------------------------------------------------------
# Shared so Stage 3 (Figure 2, 14-detector dominance matrix) and Stage 6
# (benchmarking) build them identically. Definitions and provenance caveats are
# documented in scripts/Stage6_Benchmarking.R; in particular EWARS here is an
# EWARS-STYLE alarm-indicator model, not the WHO EWARS software.
#
# All are strictly prospective: the alarm for week w uses only data before w.

#' CDC EARS C2 statistic: (x_t - baseline mean) / baseline sd, where the
#' baseline is `win` weeks ending `lag_gap` weeks before t (guard band).
#' The 7-day CDC baseline is adapted to 7 weeks for weekly data.
ears_c2_stat <- function(x, lag_gap = 2L, win = 7L) {
  n <- length(x); s <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    hi <- i - 1L - lag_gap; lo <- hi - win + 1L
    if (lo < 1L) next
    b <- x[lo:hi]
    mu <- mean(b, na.rm = TRUE); sd <- stats::sd(b, na.rm = TRUE)
    if (!is.finite(sd) || sd <= 0) sd <- 1e-9
    s[i] <- (x[i] - mu) / sd
  }
  s
}

#' EARS C2 binary alarm at the CDC threshold of 3.
build_ears <- function(cases, threshold = 3) {
  sc <- ears_c2_stat(cases)
  as.integer(!is.na(sc) & sc > threshold)
}

# -----------------------------------------------------------------------------
# Faithful to the published algorithm:
#   1. Reference set: weeks within +/- w of the same seasonal week, for each of
#      the b preceding years, excluding the most recent weeks.
#   2. Quasi-Poisson GLM with a linear time trend, log link.
#   3. Anscombe-residual reweighting to downweight past outbreaks, then refit.
#   4. Trend retained only if significant and not extrapolating beyond the data.
#   5. Threshold via the 2/3 power transformation:
#         U = [ mu0^(2/3) + (2/3) * z * mu0^(-1/3) * sqrt(V) ]^(3/2)
#      with V = phi*mu0 + Var(mu_hat0), the variance of the prediction error.
#   6. limit54 rule: no alarm unless there are at least 5 cases in the last
#      4 weeks (suppresses alarms on very sparse counts).
farrington_native <- function(y, range_idx, b = 3L, w = 3L, freq = 52L,
                              alpha = 0.05, past_weeks_not_included = 3L,
                              limit54 = c(5, 4), use_trend = TRUE,
                              min_ref = 5L) {
  n <- length(y)
  z <- stats::qnorm(1 - alpha)
  alarm <- rep(NA_integer_, n)
  upper <- rep(NA_real_, n)
  nref  <- rep(NA_integer_, n)

  for (t0 in range_idx) {
    # --- 1. reference set (unavailable points are dropped, not fatal) --------
    ref <- unlist(lapply(seq_len(b), function(k) (t0 - k * freq) + (-w:w)))
    ref <- ref[ref >= 1L & ref <= (t0 - past_weeks_not_included)]
    ref <- sort(unique(ref))
    nref[t0] <- length(ref)
    if (length(ref) < min_ref) next

    yr <- as.numeric(y[ref])
    tr <- as.numeric(ref)
    if (all(!is.finite(yr))) next
    keep <- is.finite(yr); yr <- yr[keep]; tr <- tr[keep]
    if (length(yr) < min_ref) next

    # Degenerate reference set: no historical cases at all.
    if (sum(yr) <= 0) {
      upper[t0] <- 0
      alarm[t0] <- as.integer(y[t0] > 0)
      next
    }

    tc <- tr - t0          # centre time at the prediction point

    fit_glm <- function(with_trend, wts = NULL) {
      d <- data.frame(yy = yr, tt = tc)
      f <- if (with_trend) yy ~ tt else yy ~ 1
      suppressWarnings(try(
        stats::glm(f, family = stats::quasipoisson(link = "log"),
                   data = d, weights = wts),
        silent = TRUE))
    }

    m <- fit_glm(use_trend)
    if (inherits(m, "try-error") || !m$converged) m <- fit_glm(FALSE)
    if (inherits(m, "try-error") || !m$converged) next

    # --- 3. Anscombe-residual reweighting -----------------------------------
    mu <- stats::fitted(m)
    dfr <- max(1, stats::df.residual(m))
    phi <- max(1, sum(stats::residuals(m, type = "pearson")^2) / dfr)
    s <- (3 / 2) * (yr^(2/3) - mu^(2/3)) / (mu^(1/6) * sqrt(phi))
    wts <- ifelse(abs(s) > 1, 1 / (s^2), 1)
    wts[!is.finite(wts)] <- 1
    wts <- wts * length(wts) / sum(wts)          # normalise to mean 1

    m2 <- fit_glm(use_trend, wts = wts)
    if (inherits(m2, "try-error") || !m2$converged) m2 <- fit_glm(FALSE, wts = wts)
    if (inherits(m2, "try-error") || !m2$converged) m2 <- m

    # --- 4. keep the trend only if justified --------------------------------
    has_trend <- "tt" %in% names(stats::coef(m2))
    if (has_trend) {
      sm <- try(summary(m2), silent = TRUE)
      p_trend <- if (!inherits(sm, "try-error") &&
                     "tt" %in% rownames(sm$coefficients))
        sm$coefficients["tt", 4] else 1
      pred_try <- try(stats::predict(m2, newdata = data.frame(tt = 0),
                                     type = "response"), silent = TRUE)
      extrapolating <- inherits(pred_try, "try-error") ||
        !is.finite(pred_try) || pred_try > max(yr) * 2
      if (!is.finite(p_trend) || p_trend >= 0.05 || extrapolating || b < 3L) {
        m3 <- fit_glm(FALSE, wts = wts)
        if (!inherits(m3, "try-error") && m3$converged) m2 <- m3
      }
    }

    # --- 5. prediction and threshold ----------------------------------------
    # newdata always carries tt = 0 (the prediction point). For an
    # intercept-only fit the column is simply ignored, and nrow() still tells
    # predict() to return one value -- a zero-column data frame would not.
    pr <- try(stats::predict(m2, newdata = data.frame(tt = 0),
                             type = "response", se.fit = TRUE), silent = TRUE)
    if (inherits(pr, "try-error")) next
    mu0 <- as.numeric(pr$fit)
    se0 <- as.numeric(pr$se.fit)
    if (!is.finite(mu0) || mu0 <= 0) next

    dfr2 <- max(1, stats::df.residual(m2))
    phi2 <- max(1, sum(stats::residuals(m2, type = "pearson")^2) / dfr2)
    V <- phi2 * mu0 + se0^2                       # variance of prediction error
    if (!is.finite(V) || V < 0) V <- phi2 * mu0

    U <- (mu0^(2/3) + (2/3) * z * mu0^(-1/3) * sqrt(V))^(3/2)
    if (!is.finite(U)) next
    upper[t0] <- U

    # --- 6. limit54 sparse-count rule ---------------------------------------
    lo <- max(1L, t0 - limit54[2] + 1L)
    enough <- sum(y[lo:t0], na.rm = TRUE) >= limit54[1]
    alarm[t0] <- as.integer(enough && y[t0] > U)
  }

  list(alarm = alarm, upperbound = upper, n_reference = nref)
}

#' Farrington alarms over `eval_idx`, trying the surveillance package first and
#' falling back to the dependency-free farrington_native(). Returns a list with
#' the alarm vector and the execution path actually used.
build_farrington <- function(cases, yr, eval_idx, b = 3L, w = 3L, freq = 52L,
                             use_package = TRUE) {
  n <- length(cases)
  alarm <- rep(NA_integer_, n)
  path  <- NA_character_

  if (isTRUE(use_package) && requireNamespace("surveillance", quietly = TRUE) &&
      length(eval_idx) > 0L) {
    min_ok <- b * freq + w + 1L
    rng <- eval_idx[eval_idx >= min_ok]
    if (length(rng) > 0L) {
      att <- try({
        sts_obj <- surveillance::sts(
          observed = matrix(as.integer(round(cases)), ncol = 1),
          start = c(min(yr, na.rm = TRUE), 1), freq = freq)
        fr <- surveillance::farringtonFlexible(sts_obj, control = list(
          range = rng, b = b, w = w, weightsThreshold = 2.58,
          pastWeeksNotIncluded = w, pThresholdTrend = 0.05, trend = TRUE,
          thresholdMethod = "delta", alpha = 0.05))
        as.integer(surveillance::alarms(fr))
      }, silent = TRUE)
      if (!inherits(att, "try-error") && length(att) == length(rng)) {
        alarm[rng] <- att; path <- "farringtonFlexible"
      }
    }
  }

  unscored <- eval_idx[is.na(alarm[eval_idx])]
  if (length(unscored) > 0L) {
    att <- try(farrington_native(cases, unscored, b = b, w = w, freq = freq),
               silent = TRUE)
    if (!inherits(att, "try-error")) {
      got <- att$alarm[unscored]
      alarm[unscored] <- ifelse(is.na(got), 0L, got)
      path <- if (is.na(path)) "native" else paste0(path, "+native")
    }
  }
  alarm[is.na(alarm)] <- 0L
  list(alarm = alarm, path = if (is.na(path)) "none" else path)
}

#' EWARS-style alarm-indicator model. Logistic regression of an expanding-window
#' 75th-percentile outbreak label on lagged cases and lagged rainfall, refit for
#' each evaluation year on prior years only.
#'
#' The percentile label (rather than a per-week mean+2SD rule) is what allows
#' this to train on series that start only a few years before evaluation; see
#' STAGE6_SCALE_FIX.md.
build_ewars <- function(cases, rain, yr, eval_years, excluded_years = integer(0),
                        min_prior_weeks = 52L, label_q = 0.75,
                        pooled_train = NULL) {
  n <- length(cases)
  lagv <- function(x, k) c(rep(NA_real_, k), utils::head(x, -k))
  rn <- ifelse(is.na(rain), 0, rain)
  d <- data.frame(YR = yr,
                  lag_cases_4 = lagv(cases, 4L),
                  lag_rain_8  = lagv(rn, 8L),
                  lag_rain_12 = lagv(rn, 12L))
  tl <- rep(NA_integer_, n)
  for (i in seq_len(n)) {
    if (i <= min_prior_weeks) next
    thr <- stats::quantile(cases[seq_len(i - 1L)], label_q, na.rm = TRUE)
    if (!is.finite(thr)) next
    tl[i] <- as.integer(cases[i] > thr)
  }
  d$train_outbreak <- tl

  fit_one <- function(tr) {
    if (nrow(tr) < 40L) return(NULL)
    if (sum(tr$train_outbreak == 1L, na.rm = TRUE) < 5L ||
        sum(tr$train_outbreak == 0L, na.rm = TRUE) < 5L) return(NULL)
    f <- try(suppressWarnings(stats::glm(
      train_outbreak ~ log1p(lag_cases_4) + lag_rain_8 + lag_rain_12,
      data = tr, family = stats::binomial(),
      control = stats::glm.control(maxit = 100))), silent = TRUE)
    if (inherits(f, "try-error") || !isTRUE(f$converged)) return(NULL)
    fv <- stats::fitted(f)
    if (any(fv < 1e-8 | fv > 1 - 1e-8, na.rm = TRUE)) return(NULL)
    f
  }

  alarm <- rep(NA_integer_, n)
  srcs <- character(0)
  for (yy in eval_years) {
    te <- which(yr == yy)
    if (length(te) == 0L) next
    tr <- d[d$YR < yy & !d$YR %in% excluded_years, ]
    tr <- tr[stats::complete.cases(tr[, c("train_outbreak", "lag_cases_4",
                                          "lag_rain_8", "lag_rain_12")]), ]
    fit <- fit_one(tr); src <- "unit"
    if (is.null(fit) && !is.null(pooled_train)) {
      ptr <- pooled_train[pooled_train$YR < yy &
                            !pooled_train$YR %in% excluded_years, ]
      ptr <- ptr[stats::complete.cases(ptr[, c("train_outbreak", "lag_cases_4",
                                               "lag_rain_8", "lag_rain_12")]), ]
      fit <- fit_one(ptr); src <- "pooled"
    }
    if (is.null(fit)) { srcs <- c(srcs, "none"); next }
    pr <- try(suppressWarnings(stats::predict(fit, newdata = d[te, ],
                                              type = "response")), silent = TRUE)
    if (inherits(pr, "try-error")) { srcs <- c(srcs, "none"); next }
    pr <- as.numeric(pr)
    alarm[te] <- as.integer(!is.na(pr) & pr > 0.50)
    srcs <- c(srcs, src)
  }
  alarm[is.na(alarm)] <- 0L
  list(alarm = alarm, source = paste(unique(srcs), collapse = "/"))
}
