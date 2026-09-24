# =============================================================================
# STAGE C - FOUR OPERATIONAL RAINFALL THRESHOLD PRODUCTS
# =============================================================================
# P1: RF accumulation (mm) preceding Constant TA trigger
# P2: RF accumulation (mm) preceding Outbreak Threshold trigger
# P3: rainfall STA/LTA R(t) preceding Constant TA trigger
# P4: rainfall STA/LTA R(t) preceding Outbreak Threshold trigger
#
# OPERATIONAL CONCEPT
# -------------------
# Rainfall is treated as an upper-tail environmental activation signal:
# increasing rainfall represents increasing activation. Therefore all four
# products use RF >= threshold. Low rainfall can never become a "strong RF
# trigger" merely because a retrospective ROC routine reverses direction.
#
# Each RF threshold is evaluated on the antecedent environmental feature already
# available at the disease-trigger week. P1/P2 use cumulative RF_HDX during
# t-4..t-1; P3/P4 use the maximum rainfall STA/LTA R(t) during t-4..t-1. The
# disease-trigger week itself is excluded from the RF feature.
#
# A1/A2 remains the independent epidemiological truth standard. For each
# region-year and disease detector, the Stage C positive anchor is the FIRST
# detector-active trigger week that falls in A1 OR A2. Detector-active weeks
# outside both A1 and A2 are false-alarm comparators. Later active weeks inside
# A1/A2 are retained for audit but excluded from threshold fitting so one
# prolonged detector episode cannot contribute repeated positive observations.
#
# The point estimate is selected in the >= direction using A1/A2-classified
# disease-anchor performance. Candidate thresholds are first ranked by balanced
# accuracy; among near-optimal thresholds, lower A1/A2-defined false-alarm burden
# is preferred, followed by higher PPV, sensitivity, specificity, and a more
# conservative (higher) threshold.
#
# Bounds are stability bounds from re-running this same false-alarm-aware rule:
# leave-one-year-out (LOYO) for NCR and leave-one-region-out (LORO) for the
# pooled 17-region analysis. Operational interpretation is always:
#   RF < lower             : No RF activation
#   lower <= RF < point    : Precautionary RF trigger
#   point <= RF < upper    : Formal RF trigger
#   RF >= upper            : Strong RF trigger
# =============================================================================

REQUIRED_PACKAGES <- c(
  "readxl", "dplyr", "tidyr", "purrr", "tibble", "ISOweek",
  "pROC", "ggplot2", "patchwork", "readr"
)

.bootstrap <- function() {
  root <- Sys.getenv("TA_PROJECT_ROOT", unset = "")
  if (!nzchar(root) || !dir.exists(root)) {
    d <- normalizePath(getwd(), winslash = "/")
    for (i in seq_len(6)) {
      if (dir.exists(file.path(d, "R")) && dir.exists(file.path(d, "scripts"))) {
        root <- d
        break
      }
      p <- dirname(d)
      if (identical(p, d)) break
      d <- p
    }
  }
  if (!nzchar(root)) {
    stop("Set the working directory to the project root before sourcing.", call. = FALSE)
  }
  source(file.path(root, "R", "00_config.R"))
}
.bootstrap()
require_packages(REQUIRED_PACKAGES, "Stage C")
invisible(lapply(REQUIRED_PACKAGES, function(p) {
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}))
source(file.path(DIR_R, "01_publication_theme.R"), local = TRUE)
set.seed(GLOBAL_SEED)

OUT_DIR <- file.path(DIR_OUTPUT, "StageC_RF_Threshold_Derivation")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

.eta <- load_eta_thresholds()
CASE_ETA_ON <- .eta$ETA_ON
CASE_ETA_OFF <- .eta$ETA_OFF
.win <- load_sta_lta_windows()
STA_C <- .win$STA_WIN_CONSTANT
LTA_C <- .win$LTA_WIN_CONSTANT
CASE_GUARD <- 2L
MIN_OFF_RESET <- 8L

EVALUABLE_YEARS <- c(2018L, 2019L, 2022L, 2023L, 2024L, 2025L)
MIN_OBSERVED_CASE_WEEKS_PER_YEAR <- 40L
MIN_EVALUABLE_YEARS_PER_REGION <- 5L
CANONICAL_17 <- c(
  "BARMM", "CAR", "MIMAROPA", "NCR",
  "REGION I", "REGION II", "REGION III", "REGION IV-A",
  "REGION V", "REGION VI", "REGION VII", "REGION VIII",
  "REGION IX", "REGION X", "REGION XI", "REGION XII", "REGION XIII"
)

RF_WINDOW_CANDIDATES <- tibble::tribble(
  ~STA, ~LTA,
  2L, 13L,
  3L, 20L,
  4L, 26L
)
RF_GUARD <- 2L
RF_OPERATOR_FIXED <- ">="
OPERATIONAL_YOUDEN_TOLERANCE <- 0.01
RF_WARNING_HORIZON_WEEKS <- RF_PRECEDING_WEEKS

safe_quantile <- function(x, probs) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x)) return(rep(NA_real_, length(probs)))
  as.numeric(stats::quantile(x, probs = probs, na.rm = TRUE,
                             names = FALSE, type = 8))
}

# -----------------------------------------------------------------------------
# 1. DATA AND COVERAGE
# -----------------------------------------------------------------------------
df <- readxl::read_excel(DATA_FILE, sheet = SHEET_REGIONAL)
need <- c("REGION", "YR", "WN", "LC_DOH", "RF_HDX")
miss <- setdiff(need, names(df))
if (length(miss)) stop("Missing columns: ", paste(miss, collapse = ", "), call. = FALSE)

df <- df |>
  dplyr::transmute(
    REGION = as.character(REGION),
    YR = as.integer(YR),
    WN = as.integer(WN),
    Cases = suppressWarnings(as.numeric(LC_DOH)),
    RF = suppressWarnings(as.numeric(RF_HDX))
  ) |>
  dplyr::filter(!is.na(REGION), !is.na(YR), !is.na(WN),
                WN >= 1L, WN <= 53L, REGION %in% CANONICAL_17) |>
  dplyr::mutate(Date = ISOweek::ISOweek2date(sprintf("%d-W%02d-1", YR, WN))) |>
  dplyr::filter(!is.na(Date)) |>
  dplyr::arrange(REGION, Date)

coverage <- df |>
  dplyr::filter(YR %in% EVALUABLE_YEARS) |>
  dplyr::group_by(REGION, YR) |>
  dplyr::summarise(
    n_case_weeks = sum(is.finite(Cases)),
    n_rf_weeks = sum(is.finite(RF)),
    .groups = "drop"
  )

eligible_region_year <- coverage |>
  dplyr::filter(n_case_weeks >= MIN_OBSERVED_CASE_WEEKS_PER_YEAR,
                n_rf_weeks >= MIN_OBSERVED_CASE_WEEKS_PER_YEAR) |>
  dplyr::select(REGION, YR)

regions_in <- eligible_region_year |>
  dplyr::count(REGION, name = "n_eval") |>
  dplyr::filter(n_eval >= MIN_EVALUABLE_YEARS_PER_REGION) |>
  dplyr::pull(REGION) |>
  as.character()
if (!setequal(regions_in, CANONICAL_17) || length(unique(regions_in)) != 17L) {
  stop("Stage C coverage rule must retain all 17 regions. Missing: ",
       paste(setdiff(CANONICAL_17, regions_in), collapse = ", "), call. = FALSE)
}
readr::write_csv(coverage, file.path(OUT_DIR, "RF_Regional_Coverage_Audit.csv"), na = "")

# -----------------------------------------------------------------------------
# 2. DISEASE ANCHORS
# -----------------------------------------------------------------------------
add_disease_detectors <- function(x) {
  x <- x |> dplyr::arrange(Date)
  n <- nrow(x)
  dc <- x$Cases

  bl <- rep(NA_real_, n)
  years <- sort(unique(x$YR))
  for (y in years) {
    prior <- years[years < y]
    donors <- intersect(seq.int(y - 5L, y - 1L), prior)
    if (length(donors) < 3L) donors <- tail(prior, 5L)
    if (length(donors) < 1L) next
    idx <- which(x$YR == y)
    for (i in idx) {
      vals <- dc[x$YR %in% donors & x$WN == x$WN[i]]
      vals <- vals[is.finite(vals)]
      if (length(vals) == 1L) bl[i] <- vals[1]
      if (length(vals) >= 2L) bl[i] <- mean(vals) + 2 * stats::sd(vals)
    }
  }
  outbreak <- as.integer(is.finite(dc) & is.finite(bl) & dc > bl)

  ta <- rep(0L, n)
  ratio <- rep(NA_real_, n)
  is_on <- FALSE
  frozen <- NA_real_
  off <- 0L
  min_t <- STA_C + CASE_GUARD + LTA_C
  for (t in seq_len(n)) {
    if (t > 1L && x$YR[t] != x$YR[t - 1L]) {
      is_on <- FALSE
      frozen <- NA_real_
      off <- 0L
    }
    if (!is_on) off <- off + 1L else off <- 0L
    if (!is_on && off >= MIN_OFF_RESET) frozen <- NA_real_
    if (t < min_t) next
    fs <- t - STA_C - CASE_GUARD - LTA_C + 1L
    if (fs < 1L || !is_contiguous_week_span(x$Date, fs, t)) next
    sta_vals <- dc[(t - STA_C + 1L):t]
    if (all(!is.finite(sta_vals))) next
    sta <- mean(sta_vals, na.rm = TRUE)
    if (!is_on || !is.finite(frozen)) {
      li <- fs:(t - STA_C - CASE_GUARD)
      vv <- dc[li]
      frozen <- if (all(!is.finite(vv))) NA_real_ else mean(vv, na.rm = TRUE)
    }
    r <- if (is.finite(frozen) && frozen > 0 && is.finite(sta)) sta / frozen else NA_real_
    ratio[t] <- r
    if (!is_on && is.finite(r) && r >= CASE_ETA_ON) {
      is_on <- TRUE
      off <- 0L
    }
    if (is_on && is.finite(r) && r < CASE_ETA_OFF) {
      is_on <- FALSE
      frozen <- NA_real_
    }
    ta[t] <- as.integer(is_on)
  }

  x$Outbreak_Threshold <- bl
  x$Outbreak_Signal <- outbreak
  x$Constant_TA <- ta
  x$Case_Constant_TA_Rt <- ratio
  x$RF_Antecedent_4wk_mm <- preceding_contiguous_sum(x$RF, x$Date, RF_PRECEDING_WEEKS)
  x
}

det <- df |>
  dplyr::group_by(REGION) |>
  dplyr::group_modify(~ add_disease_detectors(.x)) |>
  dplyr::ungroup()

eval_det <- det |>
  dplyr::semi_join(eligible_region_year, by = c("REGION", "YR"))

# -----------------------------------------------------------------------------
# 3. UPPER-TAIL, FALSE-ALARM-AWARE THRESHOLD SELECTION
# -----------------------------------------------------------------------------
make_threshold_grid <- function(x, event) {
  x <- suppressWarnings(as.numeric(x))
  event <- suppressWarnings(as.integer(event))
  ok <- is.finite(x) & event %in% c(0L, 1L)
  x <- x[ok]
  event <- event[ok]
  if (!length(x)) return(numeric())
  q <- safe_quantile(x, seq(0, 1, by = 0.005))
  ev <- x[event == 1L]
  eps <- max(abs(max(x, na.rm = TRUE)) * 1e-9, 1e-9)
  sort(unique(c(q, ev, max(x, na.rm = TRUE) + eps)))
}

operational_threshold_search <- function(dat, youden_tolerance = OPERATIONAL_YOUDEN_TOLERANCE) {
  d <- dat |>
    dplyr::filter(is.finite(Predictor), Event %in% c(0L, 1L)) |>
    dplyr::mutate(Event = as.integer(Event))
  empty_sel <- data.frame(
    Threshold = NA_real_, Sensitivity = NA_real_, Specificity = NA_real_,
    Balanced_Accuracy = NA_real_, Youden_J = NA_real_, PPV = NA_real_, NPV = NA_real_,
    False_Alarm_Rate = NA_real_, False_Alarm_Weeks = NA_integer_,
    False_Alarm_Weeks_per_Evaluable_Unit = NA_real_, Trigger_Rate = NA_real_,
    TP = NA_integer_, FP = NA_integer_, FN = NA_integer_, TN = NA_integer_,
    stringsAsFactors = FALSE
  )
  empty <- list(
    auc = NA_real_, auc_ci_lower = NA_real_, auc_ci_upper = NA_real_,
    roc_youden = NA_real_, selected = empty_sel, grid = data.frame(), operator = ">="
  )
  if (nrow(d) < 3L || sum(d$Event == 1L) < 1L) return(empty)

  rr <- NULL
  if (length(unique(d$Event)) >= 2L) {
    rr <- tryCatch(
      pROC::roc(response = d$Event, predictor = d$Predictor,
                quiet = TRUE, direction = "<"),
      error = function(e) NULL
    )
  }
  auc <- if (is.null(rr)) NA_real_ else as.numeric(pROC::auc(rr))
  ci_auc <- if (is.null(rr)) c(NA_real_, NA_real_) else {
    z <- tryCatch(as.numeric(pROC::ci.auc(rr, conf.level = 0.95)),
                  error = function(e) c(NA_real_, NA_real_, NA_real_))
    if (length(z) >= 3L) c(z[1L], z[3L]) else c(NA_real_, NA_real_)
  }
  roc_youden <- if (is.null(rr)) NA_real_ else {
    cc <- tryCatch(
      pROC::coords(rr, "best", best.method = "youden", ret = "threshold",
                   transpose = FALSE),
      error = function(e) NULL
    )
    if (is.null(cc)) NA_real_ else as.numeric(unlist(cc, use.names = FALSE)[1L])
  }

  thresholds <- make_threshold_grid(d$Predictor, d$Event)
  if (!length(thresholds)) return(empty)
  n_units <- d |>
    dplyr::distinct(REGION, YR) |>
    nrow()
  if (!is.finite(n_units) || n_units < 1L) n_units <- 1L

  grid <- dplyr::bind_rows(lapply(thresholds, function(th) {
    pred <- is.finite(d$Predictor) & d$Predictor >= th
    tp <- sum(pred & d$Event == 1L)
    fp <- sum(pred & d$Event == 0L)
    fn <- sum(!pred & d$Event == 1L)
    tn <- sum(!pred & d$Event == 0L)
    sens <- if ((tp + fn) > 0L) tp / (tp + fn) else NA_real_
    spec <- if ((tn + fp) > 0L) tn / (tn + fp) else NA_real_
    ppv <- if ((tp + fp) > 0L) tp / (tp + fp) else NA_real_
    npv <- if ((tn + fn) > 0L) tn / (tn + fn) else NA_real_
    bal <- if (is.finite(sens) && is.finite(spec)) (sens + spec) / 2 else sens
    youden <- if (is.finite(sens) && is.finite(spec)) sens + spec - 1 else NA_real_
    far <- if ((fp + tn) > 0L) fp / (fp + tn) else 0
    data.frame(
      Threshold = th, Sensitivity = sens, Specificity = spec,
      Balanced_Accuracy = bal, Youden_J = youden, PPV = ppv, NPV = npv,
      False_Alarm_Rate = far, False_Alarm_Weeks = fp,
      False_Alarm_Weeks_per_Evaluable_Unit = fp / n_units,
      Trigger_Rate = mean(pred), TP = tp, FP = fp, FN = fn, TN = tn,
      stringsAsFactors = FALSE
    )
  }))

  eligible <- grid |> dplyr::filter(TP > 0L, is.finite(Sensitivity), Sensitivity > 0)
  if (!nrow(eligible)) return(empty)
  # If both A1/A2 truth classes are represented, use near-optimal balanced
  # accuracy then explicitly minimize A1/A2-defined false alarms. If only true
  # anchors are present in a small scope/fold, retain maximum sensitivity and
  # choose the highest threshold that preserves it; this avoids forcing a low
  # cutoff simply because no false-anchor class is available in that fold.
  if (length(unique(d$Event)) >= 2L && any(is.finite(eligible$Balanced_Accuracy))) {
    best_bal <- max(eligible$Balanced_Accuracy, na.rm = TRUE)
    near <- eligible |>
      dplyr::filter(is.finite(Balanced_Accuracy),
                    Balanced_Accuracy >= best_bal - youden_tolerance) |>
      dplyr::arrange(
        False_Alarm_Weeks_per_Evaluable_Unit,
        False_Alarm_Rate,
        dplyr::desc(PPV),
        dplyr::desc(Sensitivity),
        dplyr::desc(Specificity),
        dplyr::desc(Threshold)
      )
  } else {
    best_sens <- max(eligible$Sensitivity, na.rm = TRUE)
    near <- eligible |>
      dplyr::filter(Sensitivity >= best_sens - 1e-12) |>
      dplyr::arrange(dplyr::desc(Threshold))
  }
  selected <- if (nrow(near)) near[1L, , drop = FALSE] else empty_sel
  list(
    auc = auc, auc_ci_lower = ci_auc[1L], auc_ci_upper = ci_auc[2L],
    roc_youden = roc_youden, selected = selected, grid = grid,
    operator = ">="
  )
}

# A1/A2 helpers intentionally mirror Stage A/B. Stage C does not redefine truth.
A1_LEAD_MIN <- 4L
A1_LEAD_MAX <- 8L
compute_stagec_a1 <- function(x) {
  x <- x |> dplyr::arrange(WN)
  if (!nrow(x) || all(!is.finite(x$Cases))) return(list(peak = NA_integer_, weeks = integer()))
  cases <- ifelse(is.finite(x$Cases), x$Cases, -Inf)
  peak_i <- which.max(cases)
  peak <- as.integer(x$WN[peak_i])
  w <- (peak - A1_LEAD_MAX):(peak - A1_LEAD_MIN)
  list(peak = peak, weeks = as.integer(w[w >= 1L]))
}
compute_stagec_a2 <- function(x, burden_frac = A2_BURDEN_FRAC) {
  x <- x |> dplyr::arrange(WN)
  if (!nrow(x) || all(!is.finite(x$Cases))) return(list(start = NA_integer_, end = NA_integer_, weeks = integer()))
  cases <- ifelse(is.finite(x$Cases), x$Cases, 0)
  total <- sum(cases)
  if (!is.finite(total) || total <= 0) return(list(start = NA_integer_, end = NA_integer_, weeks = integer()))
  peak_i <- which.max(cases); lo <- hi <- peak_i; s <- cases[peak_i]
  while (s / total < burden_frac && (lo > 1L || hi < nrow(x))) {
    lv <- if (lo > 1L) cases[lo - 1L] else -Inf
    rv <- if (hi < nrow(x)) cases[hi + 1L] else -Inf
    if (lv >= rv) { lo <- lo - 1L; s <- s + lv } else { hi <- hi + 1L; s <- s + rv }
  }
  list(start = as.integer(x$WN[lo]), end = as.integer(x$WN[hi]),
       weeks = as.integer(x$WN[lo:hi]))
}

# Build every detector-active trigger week and classify it against A1/A2.
# Stage C then selects exactly one positive derivation anchor per region-year:
# the FIRST detector-active week that falls inside A1 OR A2. All detector-active
# weeks outside A1 and A2 are retained as false-alarm comparators. Later active
# weeks inside A1/A2 are preserved in the audit inventory but excluded from
# threshold fitting so a long detector episode cannot contribute repeated
# positive observations from the same region-year.
event_trigger_weeks_a1a2 <- function(x, signal, product_id, product_label) {
  x <- x |> dplyr::arrange(WN)
  idx <- which(as.integer(x[[signal]]) == 1L)
  if (!length(idx)) return(data.frame(
    Product_ID = integer(), Product = character(), WN = integer(),
    Date = as.Date(character()), Cases = numeric(), RF_Antecedent_4wk_mm = numeric(),
    Peak_Week = integer(), A1_Start = integer(), A1_End = integer(),
    A2_Start = integer(), A2_End = integer(), InA1 = logical(), InA2 = logical(),
    IsTrue = logical(), Alarm_Class = character(), stringsAsFactors = FALSE
  ))
  a1 <- compute_stagec_a1(x)
  a2 <- compute_stagec_a2(x)
  dplyr::bind_rows(lapply(idx, function(i) {
    wk <- as.integer(x$WN[i])
    ina1 <- wk %in% a1$weeks
    ina2 <- wk %in% a2$weeks
    data.frame(
      Product_ID = product_id, Product = product_label,
      WN = wk, Date = x$Date[i], Cases = x$Cases[i],
      RF_Antecedent_4wk_mm = x$RF_Antecedent_4wk_mm[i],
      Peak_Week = a1$peak,
      A1_Start = if (length(a1$weeks)) min(a1$weeks) else NA_integer_,
      A1_End = if (length(a1$weeks)) max(a1$weeks) else NA_integer_,
      A2_Start = a2$start, A2_End = a2$end,
      InA1 = ina1, InA2 = ina2, IsTrue = ina1 || ina2,
      Alarm_Class = if (ina1 || ina2) "True alarm" else "False alarm",
      stringsAsFactors = FALSE
    )
  }))
}

select_stagec_derivation_inventory <- function(trigger_inventory) {
  if (!nrow(trigger_inventory)) {
    trigger_inventory$Derivation_Role <- character()
    trigger_inventory$Use_For_Threshold_Derivation <- logical()
    return(trigger_inventory)
  }
  trigger_inventory |>
    dplyr::group_by(REGION, YR) |>
    dplyr::group_modify(~ {
      z <- .x |> dplyr::arrange(WN)
      true_rows <- z |> dplyr::filter(IsTrue)
      primary_true <- if (nrow(true_rows)) {
        true_rows |> dplyr::slice_min(WN, n = 1L, with_ties = FALSE)
      } else {
        true_rows
      }
      false_rows <- z |> dplyr::filter(!IsTrue)
      later_true <- if (nrow(true_rows) > 1L) {
        true_rows |> dplyr::filter(WN != min(true_rows$WN, na.rm = TRUE))
      } else {
        true_rows[0, , drop = FALSE]
      }
      dplyr::bind_rows(
        primary_true |> dplyr::mutate(
          Derivation_Role = "Primary true anchor: first trigger in A1 or A2",
          Use_For_Threshold_Derivation = TRUE),
        false_rows |> dplyr::mutate(
          Derivation_Role = "False-alarm comparator: trigger outside A1 and A2",
          Use_For_Threshold_Derivation = TRUE),
        later_true |> dplyr::mutate(
          Derivation_Role = "Later true trigger: audit only",
          Use_For_Threshold_Derivation = FALSE)
      ) |> dplyr::arrange(WN)
    }) |>
    dplyr::ungroup()
}

events_cta_all <- eval_det |>
  dplyr::group_by(REGION, YR) |>
  dplyr::group_modify(~ event_trigger_weeks_a1a2(
    .x, "Constant_TA", 1L, "RF preceding Constant TA trigger")) |>
  dplyr::ungroup()
events_ot_all <- eval_det |>
  dplyr::group_by(REGION, YR) |>
  dplyr::group_modify(~ event_trigger_weeks_a1a2(
    .x, "Outbreak_Signal", 2L, "RF preceding Outbreak Threshold trigger")) |>
  dplyr::ungroup()

events_cta_inventory <- select_stagec_derivation_inventory(events_cta_all)
events_ot_inventory <- select_stagec_derivation_inventory(events_ot_all)

events_cta <- events_cta_inventory |> dplyr::filter(Use_For_Threshold_Derivation)
events_ot <- events_ot_inventory |> dplyr::filter(Use_For_Threshold_Derivation)

events_cta_first_true <- events_cta |>
  dplyr::filter(IsTrue) |>
  dplyr::arrange(REGION, YR, WN)
events_ot_first_true <- events_ot |>
  dplyr::filter(IsTrue) |>
  dplyr::arrange(REGION, YR, WN)

# Explicit audit exports. The legacy first-event filenames now contain the
# primary A1/A2-true anchor (first qualifying trigger), not the first raw onset.
readr::write_csv(events_cta_first_true,
                 file.path(OUT_DIR, "Product1_First_Constant_TA_Events.csv"), na = "")
readr::write_csv(events_ot_first_true,
                 file.path(OUT_DIR, "Product2_First_Outbreak_Threshold_Events.csv"), na = "")
readr::write_csv(events_cta_first_true,
                 file.path(OUT_DIR, "Product1_First_A1A2_True_Constant_TA_Anchor.csv"), na = "")
readr::write_csv(events_ot_first_true,
                 file.path(OUT_DIR, "Product2_First_A1A2_True_Outbreak_Threshold_Anchor.csv"), na = "")
readr::write_csv(events_cta_all,
                 file.path(OUT_DIR, "Product1_All_Constant_TA_TriggerWeeks_A1A2.csv"), na = "")
readr::write_csv(events_ot_all,
                 file.path(OUT_DIR, "Product2_All_Outbreak_Threshold_TriggerWeeks_A1A2.csv"), na = "")
readr::write_csv(events_cta_inventory,
                 file.path(OUT_DIR, "Product1_StageC_Derivation_Inventory_A1A2.csv"), na = "")
readr::write_csv(events_ot_inventory,
                 file.path(OUT_DIR, "Product2_StageC_Derivation_Inventory_A1A2.csv"), na = "")
# Compatibility aliases retained for older downstream readers.
readr::write_csv(events_cta_inventory,
                 file.path(OUT_DIR, "Product1_All_Constant_TA_Onsets_A1A2.csv"), na = "")
readr::write_csv(events_ot_inventory,
                 file.path(OUT_DIR, "Product2_All_Outbreak_Threshold_Onsets_A1A2.csv"), na = "")

# Independent-support audit. P1/P3 share Constant-TA support; P2/P4 share
# Outbreak-Threshold support. The primary positive count is one first A1/A2
# qualifying trigger per region-year, not the number of active detector weeks.
stagec_support_one <- function(data_scope, all_events, inventory, anchor, scope_label) {
  units <- data_scope |> dplyr::distinct(REGION, YR)
  primary <- inventory |> dplyr::filter(Use_For_Threshold_Derivation, IsTrue)
  false_cmp <- inventory |> dplyr::filter(Use_For_Threshold_Derivation, !IsTrue)
  later_true <- inventory |> dplyr::filter(!Use_For_Threshold_Derivation, IsTrue)
  data.frame(
    Scope = scope_label,
    Anchor = anchor,
    N_Evaluable_Units = nrow(units),
    N_Evaluable_Years = dplyr::n_distinct(units$YR),
    N_Units_With_Any_Detector_Trigger = nrow(all_events |> dplyr::distinct(REGION, YR)),
    N_Units_With_A1A2_True_Trigger = nrow(all_events |> dplyr::filter(IsTrue) |> dplyr::distinct(REGION, YR)),
    N_Qualifying_Years = dplyr::n_distinct(primary$YR),
    N_Primary_True_Anchors = nrow(primary),
    N_False_Comparator_Triggers = nrow(false_cmp),
    N_Later_True_Triggers_AuditOnly = nrow(later_true),
    stringsAsFactors = FALSE
  )
}

support_reg_cta <- stagec_support_one(eval_det, events_cta_all, events_cta_inventory,
                                       "Constant TA", "Regional (pooled 17 regions)")
support_reg_ot <- stagec_support_one(eval_det, events_ot_all, events_ot_inventory,
                                     "Outbreak Threshold", "Regional (pooled 17 regions)")
support_ncr_cta <- stagec_support_one(
  eval_det |> dplyr::filter(REGION == "NCR"),
  events_cta_all |> dplyr::filter(REGION == "NCR"),
  events_cta_inventory |> dplyr::filter(REGION == "NCR"),
  "Constant TA", "NCR")
support_ncr_ot <- stagec_support_one(
  eval_det |> dplyr::filter(REGION == "NCR"),
  events_ot_all |> dplyr::filter(REGION == "NCR"),
  events_ot_inventory |> dplyr::filter(REGION == "NCR"),
  "Outbreak Threshold", "NCR")

stagec_support <- dplyr::bind_rows(
  support_reg_cta |> tidyr::uncount(2, .id = ".rep") |>
    dplyr::mutate(Product_ID = ifelse(.rep == 1L, 1L, 3L)) |> dplyr::select(-.rep),
  support_reg_ot |> tidyr::uncount(2, .id = ".rep") |>
    dplyr::mutate(Product_ID = ifelse(.rep == 1L, 2L, 4L)) |> dplyr::select(-.rep),
  support_ncr_cta |> tidyr::uncount(2, .id = ".rep") |>
    dplyr::mutate(Product_ID = ifelse(.rep == 1L, 1L, 3L)) |> dplyr::select(-.rep),
  support_ncr_ot |> tidyr::uncount(2, .id = ".rep") |>
    dplyr::mutate(Product_ID = ifelse(.rep == 1L, 2L, 4L)) |> dplyr::select(-.rep)
) |> dplyr::select(Scope, Product_ID, dplyr::everything())
readr::write_csv(stagec_support, file.path(OUT_DIR, "RF_Threshold_Derivation_Support.csv"), na = "")

key_of <- function(d) paste(d$REGION, d$YR, d$WN, sep = "|")
event_key_cta <- key_of(events_cta)
event_key_ot <- key_of(events_ot)
anchor_truth_lookup <- dplyr::bind_rows(
  events_cta |> dplyr::mutate(Anchor = "Constant TA"),
  events_ot |> dplyr::mutate(Anchor = "Outbreak Threshold")
) |>
  dplyr::mutate(Key = paste(REGION, YR, WN, sep = "|"))

rf_anchor_truth_product_id <- function(product_id) {
  pid <- suppressWarnings(as.integer(product_id)[1L])
  if (is.na(pid) || !pid %in% 1:4) {
    stop("RF product_id must be one of 1, 2, 3, or 4.", call. = FALSE)
  }
  # P1 and P3 share the Constant-TA truth inventory; P2 and P4 share the
  # Outbreak-Threshold truth inventory. The anchor table itself is stored under
  # Product_ID 1 (CTA) and Product_ID 2 (OT), so P3/P4 must be mapped back to
  # their disease-anchor truth IDs rather than searched as nonexistent IDs 3/4.
  if (pid %in% c(1L, 3L)) 1L else 2L
}

# Runtime regression for the four-product-to-two-anchor mapping. This executes
# before any RF STA/LTA window search, so a future mapping regression fails
# immediately instead of surfacing after the expensive derivation stages.
.rf_anchor_map_runtime <- vapply(1:4, rf_anchor_truth_product_id, integer(1))
if (!identical(unname(.rf_anchor_map_runtime), c(1L, 2L, 1L, 2L))) {
  stop("Internal RF anchor mapping regression: expected P1/P3 -> CTA and P2/P4 -> Outbreak Threshold.",
       call. = FALSE)
}
readr::write_csv(
  data.frame(
    RF_Product = 1:4,
    Anchor_Truth_Product_ID = unname(.rf_anchor_map_runtime),
    Disease_Anchor = c("Constant TA", "Outbreak Threshold", "Constant TA", "Outbreak Threshold")
  ),
  file.path(OUT_DIR, "RF_Product_Anchor_Truth_Mapping.csv")
)

label_future_anchor_window <- function(data_scope, event_keys, product_id, max_lead = RF_WARNING_HORIZON_WEEKS) {
  # Historical function name retained for compatibility. The returned rows are
  # now the actual disease-trigger weeks, labelled by A1/A2. Predictor values
  # are antecedent t-4..t-1 features already available at the anchor week.
  keys <- unique(as.character(event_keys))
  truth_pid <- rf_anchor_truth_product_id(product_id)
  truth <- anchor_truth_lookup |>
    dplyr::filter(Product_ID == truth_pid, Key %in% keys)
  if (!nrow(truth)) {
    stop("No A1/A2 anchor-truth rows were found for RF Product ", product_id,
         " (mapped disease-anchor truth Product_ID ", truth_pid, ").",
         call. = FALSE)
  }
  d <- data_scope |>
    dplyr::mutate(Key = paste(REGION, YR, WN, sep = "|")) |>
    dplyr::inner_join(
      truth |> dplyr::select(Key, Peak_Week, A1_Start, A1_End, A2_Start, A2_End,
                             InA1, InA2, IsTrue, Alarm_Class),
      by = "Key") |>
    dplyr::mutate(
      Anchor_Date = Date,
      Anchor_WN = WN,
      Signal_Date = Date - 7,
      Lead_to_Anchor_Weeks = 1L,
      Event = as.integer(IsTrue)
    )
  d
}

make_mm_input <- function(data_scope, event_keys, product_id) {
  data_scope |>
    dplyr::filter(YR %in% EVALUABLE_YEARS, is.finite(RF_Antecedent_4wk_mm)) |>
    label_future_anchor_window(event_keys, product_id, RF_WARNING_HORIZON_WEEKS) |>
    dplyr::transmute(
      Product_ID = product_id, REGION, YR, WN, Date, Signal_Date,
      Anchor_Date, Anchor_WN, Lead_to_Anchor_Weeks,
      Peak_Week, A1_Start, A1_End, A2_Start, A2_End,
      InA1, InA2, IsTrue, Alarm_Class,
      Predictor = RF_Antecedent_4wk_mm, Event
    )
}

roc_mm_cta <- make_mm_input(eval_det, event_key_cta, 1L)
roc_mm_ot <- make_mm_input(eval_det, event_key_ot, 2L)
roc_mm_cta_ncr <- roc_mm_cta |> dplyr::filter(REGION == "NCR")
roc_mm_ot_ncr <- roc_mm_ot |> dplyr::filter(REGION == "NCR")
readr::write_csv(dplyr::bind_rows(roc_mm_cta, roc_mm_ot),
                 file.path(OUT_DIR, "RF_mm_ROC_Input.csv"), na = "")

mm_summary_one <- function(rocdat, product_id, label, scope_label) {
  op <- operational_threshold_search(rocdat)
  s <- op$selected
  ev <- rocdat$Predictor[rocdat$Event == 1L & is.finite(rocdat$Predictor)]
  if (!nrow(s) || !is.finite(s$Threshold[1])) {
    stop("Could not derive an upper-tail A1/A2-consistent RF-mm threshold for Product ", product_id,
         " in scope ", scope_label, ".", call. = FALSE)
  }
  data.frame(
    Scope = scope_label,
    Product_ID = product_id, Product = label,
    Anchor = ifelse(product_id == 1L, "Constant TA", "Outbreak Threshold"),
    RF_Scale = "mm", Feature = "Cumulative RF_HDX during t-4..t-1; trigger week excluded",
    N_Events = nrow(rocdat),
    N_Positive_Warning_Weeks = sum(rocdat$Event == 1L),
    N_RF_Eligible_Qualifying_Units = dplyr::n_distinct(
      paste(rocdat$REGION[rocdat$Event == 1L], rocdat$YR[rocdat$Event == 1L], sep = "|")),
    N_RF_Eligible_Qualifying_Years = dplyr::n_distinct(rocdat$YR[rocdat$Event == 1L]),
    N_True_Anchor_Alarms = sum(rocdat$Event == 1L),
    N_False_Anchor_Alarms = sum(rocdat$Event == 0L),
    Truth_Definition = "Primary true anchor = first detector-active trigger week in A1 OR A2 per region-year; false comparators = detector-active trigger weeks outside A1 and A2",
    Event_Median = if (length(ev)) stats::median(ev) else NA_real_,
    ROC_AUC = op$auc, ROC_AUC_CI_Lower = op$auc_ci_lower,
    ROC_AUC_CI_Upper = op$auc_ci_upper, ROC_Youden = op$roc_youden,
    Sensitivity = s$Sensitivity[1], Specificity = s$Specificity[1],
    Balanced_Accuracy = s$Balanced_Accuracy[1], PPV = s$PPV[1],
    False_Alarm_Rate = s$False_Alarm_Rate[1],
    False_Alarm_Weeks = s$False_Alarm_Weeks[1],
    False_Alarm_Weeks_per_Evaluable_Unit = s$False_Alarm_Weeks_per_Evaluable_Unit[1],
    False_Alarm_Onsets = s$False_Alarm_Weeks[1],
    False_Alarms_per_Evaluable_Unit = s$False_Alarm_Weeks_per_Evaluable_Unit[1],
    Trigger_Rate = s$Trigger_Rate[1],
    Adopted_Threshold = s$Threshold[1], Units = "mm",
    Comparison_Operator = ">=", ROC_Direction = "<",
    RF_STA = NA_integer_, RF_LTA = NA_integer_, RF_Guard = NA_integer_,
    Adopted_Method = "Upper-tail A1/A2 anchor classification; near-optimal balanced accuracy with minimum false-alarm burden tie-break",
    Selection_Youden_Tolerance = OPERATIONAL_YOUDEN_TOLERANCE,
    Operational_Direction = "Higher rainfall activates (>=)",
    stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------------------------
# 4. RF STA/LTA PRODUCTS
# -----------------------------------------------------------------------------
add_rf_ratio_candidate <- function(data, sta, lta, guard) {
  data |>
    dplyr::group_by(REGION) |>
    dplyr::group_modify(~ {
      z <- .x |> dplyr::arrange(Date)
      z$RF_STA_LTA_Rt <- guarded_sta_lta_ratio(z$RF, z$Date, sta, lta, guard)
      z$RF_Antecedent_4wk_Rt_Max <- preceding_contiguous_max(
        z$RF_STA_LTA_Rt, z$Date, RF_PRECEDING_WEEKS)
      z
    }) |>
    dplyr::ungroup()
}

make_rt_input <- function(data_scope, sta, lta, event_keys, product_id) {
  add_rf_ratio_candidate(data_scope, sta, lta, RF_GUARD) |>
    dplyr::filter(YR %in% EVALUABLE_YEARS, is.finite(RF_Antecedent_4wk_Rt_Max)) |>
    label_future_anchor_window(event_keys, product_id, RF_WARNING_HORIZON_WEEKS) |>
    dplyr::transmute(
      Product_ID = product_id, REGION, YR, WN, Date, Signal_Date, Cases,
      Anchor_Date, Anchor_WN, Lead_to_Anchor_Weeks,
      Peak_Week, A1_Start, A1_End, A2_Start, A2_End,
      InA1, InA2, IsTrue, Alarm_Class,
      RF_STA_LTA_Rt, RF_Antecedent_4wk_Rt_Max,
      Predictor = RF_Antecedent_4wk_Rt_Max, Event
    )
}

score_rf_window <- function(data_scope, sta, lta, event_keys, product_id) {
  z <- make_rt_input(data_scope, sta, lta, event_keys, product_id)
  op <- operational_threshold_search(z)
  s <- op$selected
  if (!nrow(s)) s <- data.frame(
    Threshold = NA_real_, Sensitivity = NA_real_, Specificity = NA_real_,
    Balanced_Accuracy = NA_real_, PPV = NA_real_, False_Alarm_Rate = NA_real_,
    False_Alarm_Weeks_per_Evaluable_Unit = NA_real_, Trigger_Rate = NA_real_
  )
  data.frame(
    Product_ID = product_id, STA = sta, LTA = lta, Guard = RF_GUARD,
    ROC_AUC = op$auc, ROC_Youden = op$roc_youden,
    Selected_Threshold = s$Threshold[1], Sensitivity = s$Sensitivity[1],
    Specificity = s$Specificity[1], Balanced_Accuracy = s$Balanced_Accuracy[1],
    PPV = s$PPV[1], False_Alarm_Rate = s$False_Alarm_Rate[1],
    False_Alarm_Weeks_per_Evaluable_Unit = s$False_Alarm_Weeks_per_Evaluable_Unit[1],
    Trigger_Rate = s$Trigger_Rate[1], Comparison_Operator = ">=",
    N_Event = sum(z$Event == 1L), N_Control = sum(z$Event == 0L),
    stringsAsFactors = FALSE
  )
}

rank_rf_windows <- function(scored) {
  scored |>
    dplyr::arrange(
      dplyr::desc(Balanced_Accuracy),
      False_Alarm_Rate,
      False_Alarm_Weeks_per_Evaluable_Unit,
      dplyr::desc(PPV),
      dplyr::desc(Sensitivity),
      dplyr::desc(ROC_AUC),
      STA, LTA
    ) |>
    dplyr::mutate(Window_Rank = dplyr::row_number())
}

derive_rt_scope <- function(data_scope, event_keys, product_id, anchor_label, scope_label) {
  scored <- dplyr::bind_rows(lapply(seq_len(nrow(RF_WINDOW_CANDIDATES)), function(i) {
    score_rf_window(data_scope, RF_WINDOW_CANDIDATES$STA[i],
                    RF_WINDOW_CANDIDATES$LTA[i], event_keys, product_id)
  })) |>
    dplyr::mutate(Scope = scope_label, .before = 1)
  ranked <- rank_rf_windows(scored)
  if (!nrow(ranked) || !is.finite(ranked$Selected_Threshold[1])) {
    stop("Could not empirically select an upper-tail rainfall STA/LTA product for P",
         product_id, " in scope ", scope_label, ".", call. = FALSE)
  }
  sta <- as.integer(ranked$STA[1])
  lta <- as.integer(ranked$LTA[1])
  z <- make_rt_input(data_scope, sta, lta, event_keys, product_id)
  op <- operational_threshold_search(z)
  s <- op$selected
  event_vals <- z$Predictor[z$Event == 1L & is.finite(z$Predictor)]
  control_vals <- z$Predictor[z$Event == 0L & is.finite(z$Predictor)]

  methods <- data.frame(
    Scope = scope_label, Product_ID = product_id,
    Method = c(
      "Adopted upper-tail A1/A2 operating point",
      "Fixed-direction ROC Youden reference",
      "Non-event 90th percentile reference",
      "Event median reference"
    ),
    Threshold_Rt = c(
      s$Threshold[1], op$roc_youden,
      safe_quantile(control_vals, 0.90),
      if (length(event_vals)) stats::median(event_vals) else NA_real_
    ),
    Adopted = c(TRUE, FALSE, FALSE, FALSE),
    stringsAsFactors = FALSE
  )

  summary <- data.frame(
    Scope = scope_label, Product_ID = product_id,
    Product = paste0("RF STA/LTA R(t) preceding ", anchor_label, " trigger"),
    Anchor = anchor_label, RF_Scale = "R(t)",
    Feature = "Maximum rainfall STA/LTA R(t) during t-4..t-1; available by t-1",
    N_Events = dplyr::n_distinct(paste(z$REGION[z$Event == 1L], z$YR[z$Event == 1L], sep = "|")),
    N_Positive_Warning_Weeks = sum(z$Event == 1L),
    N_RF_Eligible_Qualifying_Units = dplyr::n_distinct(
      paste(z$REGION[z$Event == 1L], z$YR[z$Event == 1L], sep = "|")),
    N_RF_Eligible_Qualifying_Years = dplyr::n_distinct(z$YR[z$Event == 1L]),
    Event_Median = if (length(event_vals)) stats::median(event_vals) else NA_real_,
    ROC_AUC = op$auc, ROC_AUC_CI_Lower = op$auc_ci_lower,
    ROC_AUC_CI_Upper = op$auc_ci_upper, ROC_Youden = op$roc_youden,
    Sensitivity = s$Sensitivity[1], Specificity = s$Specificity[1],
    Balanced_Accuracy = s$Balanced_Accuracy[1], PPV = s$PPV[1],
    False_Alarm_Rate = s$False_Alarm_Rate[1],
    False_Alarm_Weeks = s$False_Alarm_Weeks[1],
    False_Alarm_Weeks_per_Evaluable_Unit = s$False_Alarm_Weeks_per_Evaluable_Unit[1],
    False_Alarm_Onsets = s$False_Alarm_Weeks[1],
    False_Alarms_per_Evaluable_Unit = s$False_Alarm_Weeks_per_Evaluable_Unit[1],
    Trigger_Rate = s$Trigger_Rate[1],
    Adopted_Threshold = s$Threshold[1], Units = "R(t)",
    Comparison_Operator = ">=", ROC_Direction = "<",
    RF_STA = sta, RF_LTA = lta, RF_Guard = RF_GUARD,
    Adopted_Method = "Upper-tail A1/A2 anchor classification; near-optimal balanced accuracy with minimum false-alarm burden tie-break",
    Selection_Youden_Tolerance = OPERATIONAL_YOUDEN_TOLERANCE,
    Operational_Direction = "Higher rainfall acceleration activates (>=)",
    stringsAsFactors = FALSE
  )

  events <- z |>
    dplyr::filter(Event == 1L) |>
    dplyr::mutate(
      Product_ID = product_id,
      Product = summary$Product[1],
      Comparison_Operator = ">=",
      RF_STA = sta, RF_LTA = lta, RF_Guard = RF_GUARD
    ) |>
    dplyr::select(
      Product_ID, Product, Comparison_Operator, REGION, YR, WN, Date, Cases,
      RF_STA, RF_LTA, RF_Guard,
      RF_Rt_at_Disease_Trigger = RF_STA_LTA_Rt,
      RF_Antecedent_4wk_Rt_Max
    )

  list(summary = summary, methods = methods, candidates = ranked,
       roc_input = z, events = events)
}

# Scope-specific point estimates ------------------------------------------------
reg_p1 <- mm_summary_one(roc_mm_cta, 1L,
                         "RF mm preceding Constant TA trigger",
                         "Regional (pooled 17 regions)")
reg_p2 <- mm_summary_one(roc_mm_ot, 2L,
                         "RF mm preceding Outbreak Threshold trigger",
                         "Regional (pooled 17 regions)")
reg_p3 <- derive_rt_scope(det, event_key_cta, 3L, "Constant TA",
                          "Regional (pooled 17 regions)")
reg_p4 <- derive_rt_scope(det, event_key_ot, 4L, "Outbreak Threshold",
                          "Regional (pooled 17 regions)")

ncr_det <- det |> dplyr::filter(REGION == "NCR")
ncr_event_key_cta <- event_key_cta[grepl("^NCR\\|", event_key_cta)]
ncr_event_key_ot <- event_key_ot[grepl("^NCR\\|", event_key_ot)]
ncr_p1 <- mm_summary_one(roc_mm_cta_ncr, 1L,
                         "RF mm preceding Constant TA trigger", "NCR")
ncr_p2 <- mm_summary_one(roc_mm_ot_ncr, 2L,
                         "RF mm preceding Outbreak Threshold trigger", "NCR")
ncr_p3 <- derive_rt_scope(ncr_det, ncr_event_key_cta, 3L, "Constant TA", "NCR")
ncr_p4 <- derive_rt_scope(ncr_det, ncr_event_key_ot, 4L, "Outbreak Threshold", "NCR")

readr::write_csv(reg_p3$events, file.path(OUT_DIR, "Product3_First_Constant_TA_RF_Rt_Events.csv"), na = "")
readr::write_csv(reg_p4$events, file.path(OUT_DIR, "Product4_First_Outbreak_Threshold_RF_Rt_Events.csv"), na = "")
readr::write_csv(dplyr::bind_rows(reg_p3$candidates, reg_p4$candidates,
                                  ncr_p3$candidates, ncr_p4$candidates),
                 file.path(OUT_DIR, "RF_STA_LTA_Window_Candidates.csv"), na = "")
readr::write_csv(dplyr::bind_rows(reg_p3$methods, reg_p4$methods,
                                  ncr_p3$methods, ncr_p4$methods),
                 file.path(OUT_DIR, "RF_Rt_Threshold_Empirical_Methods.csv"), na = "")
readr::write_csv(dplyr::bind_rows(
  reg_p3$roc_input |> dplyr::mutate(Scope = "Regional (pooled 17 regions)", .before = 1),
  reg_p4$roc_input |> dplyr::mutate(Scope = "Regional (pooled 17 regions)", .before = 1),
  ncr_p3$roc_input |> dplyr::mutate(Scope = "NCR", .before = 1),
  ncr_p4$roc_input |> dplyr::mutate(Scope = "NCR", .before = 1)
), file.path(OUT_DIR, "RF_Rt_ROC_Input.csv"), na = "")

# Full candidate-threshold audit. This export makes the A1/A2 false-alarm-aware
# point-estimate selection reproducible rather than leaving only the selected row.
threshold_grid_one <- function(rocdat, scope_label, product_id, anchor_label, rf_scale) {
  op <- operational_threshold_search(rocdat)
  g <- op$grid
  if (!nrow(g)) {
    stop("No operational threshold candidate grid for Product ", product_id,
         " in scope ", scope_label, ".", call. = FALSE)
  }
  chosen <- op$selected$Threshold[1]
  g |>
    dplyr::mutate(
      Scope = scope_label, Product_ID = product_id, Anchor = anchor_label,
      RF_Scale = rf_scale, Comparison_Operator = ">=",
      False_Alarm_Onsets = False_Alarm_Weeks,
      False_Alarms_per_Evaluable_Unit = False_Alarm_Weeks_per_Evaluable_Unit,
      Is_Selected_Point_Estimate = is.finite(Threshold) & is.finite(chosen) &
        abs(Threshold - chosen) <= .Machine$double.eps^0.5 * max(1, abs(chosen)),
      .before = 1
    )
}

threshold_candidate_grid <- dplyr::bind_rows(
  threshold_grid_one(roc_mm_cta, "Regional (pooled 17 regions)", 1L, "Constant TA", "mm"),
  threshold_grid_one(roc_mm_ot, "Regional (pooled 17 regions)", 2L, "Outbreak Threshold", "mm"),
  threshold_grid_one(reg_p3$roc_input, "Regional (pooled 17 regions)", 3L, "Constant TA", "R(t)"),
  threshold_grid_one(reg_p4$roc_input, "Regional (pooled 17 regions)", 4L, "Outbreak Threshold", "R(t)"),
  threshold_grid_one(roc_mm_cta_ncr, "NCR", 1L, "Constant TA", "mm"),
  threshold_grid_one(roc_mm_ot_ncr, "NCR", 2L, "Outbreak Threshold", "mm"),
  threshold_grid_one(ncr_p3$roc_input, "NCR", 3L, "Constant TA", "R(t)"),
  threshold_grid_one(ncr_p4$roc_input, "NCR", 4L, "Outbreak Threshold", "R(t)")
)
readr::write_csv(threshold_candidate_grid,
                 file.path(OUT_DIR, "RF_Operational_Threshold_Candidate_Grid.csv"), na = "")

# -----------------------------------------------------------------------------
# 5. THRESHOLD STABILITY BOUNDS
# -----------------------------------------------------------------------------
stability_mm <- function(rocdat, product_id, mode = c("region", "year")) {
  mode <- match.arg(mode)
  units <- if (mode == "region") CANONICAL_17 else sort(unique(rocdat$YR))
  dplyr::bind_rows(lapply(units, function(u) {
    train <- if (mode == "region") rocdat |> dplyr::filter(REGION != u) else rocdat |> dplyr::filter(YR != u)
    op <- operational_threshold_search(train)
    s <- op$selected
    data.frame(
      Product_ID = product_id,
      Held_Out_Region = if (mode == "region") as.character(u) else NA_character_,
      Held_Out_Year = if (mode == "year") as.integer(u) else NA_integer_,
      RF_STA = NA_integer_, RF_LTA = NA_integer_,
      Threshold = if (nrow(s)) s$Threshold[1] else NA_real_,
      AUC = op$auc,
      False_Alarm_Rate = if (nrow(s)) s$False_Alarm_Rate[1] else NA_real_,
      Comparison_Operator = ">=", stringsAsFactors = FALSE
    )
  }))
}

stability_rt <- function(data_scope, event_keys, product_id, mode = c("region", "year")) {
  # Re-run the full rainfall STA/LTA model-selection procedure inside every
  # stability fold. This lets the uncertainty bounds reflect both threshold
  # uncertainty and uncertainty in the empirically selected RF STA/LTA window.
  mode <- match.arg(mode)
  units <- if (mode == "region") CANONICAL_17 else sort(unique(data_scope$YR))
  dplyr::bind_rows(lapply(units, function(u) {
    train_scope <- if (mode == "region") {
      data_scope |> dplyr::filter(REGION != u)
    } else {
      data_scope |> dplyr::filter(YR != u)
    }
    scored <- dplyr::bind_rows(lapply(seq_len(nrow(RF_WINDOW_CANDIDATES)), function(i) {
      score_rf_window(
        train_scope,
        RF_WINDOW_CANDIDATES$STA[i], RF_WINDOW_CANDIDATES$LTA[i],
        event_keys, product_id
      )
    }))
    ranked <- rank_rf_windows(scored)
    if (!nrow(ranked) || !is.finite(ranked$Selected_Threshold[1])) {
      return(data.frame(
        Product_ID = product_id,
        Held_Out_Region = if (mode == "region") as.character(u) else NA_character_,
        Held_Out_Year = if (mode == "year") as.integer(u) else NA_integer_,
        RF_STA = NA_integer_, RF_LTA = NA_integer_,
        Threshold = NA_real_, AUC = NA_real_, False_Alarm_Rate = NA_real_,
        Comparison_Operator = ">=", stringsAsFactors = FALSE
      ))
    }
    best <- ranked[1L, , drop = FALSE]
    data.frame(
      Product_ID = product_id,
      Held_Out_Region = if (mode == "region") as.character(u) else NA_character_,
      Held_Out_Year = if (mode == "year") as.integer(u) else NA_integer_,
      RF_STA = as.integer(best$STA[1]), RF_LTA = as.integer(best$LTA[1]),
      Threshold = as.numeric(best$Selected_Threshold[1]),
      AUC = as.numeric(best$ROC_AUC[1]),
      False_Alarm_Rate = as.numeric(best$False_Alarm_Rate[1]),
      Comparison_Operator = ">=", stringsAsFactors = FALSE
    )
  }))
}

loro_all <- dplyr::bind_rows(
  stability_mm(roc_mm_cta, 1L, "region"),
  stability_mm(roc_mm_ot, 2L, "region"),
  stability_rt(eval_det, event_key_cta, 3L, "region"),
  stability_rt(eval_det, event_key_ot, 4L, "region")
)
loyo_ncr <- dplyr::bind_rows(
  stability_mm(roc_mm_cta_ncr, 1L, "year"),
  stability_mm(roc_mm_ot_ncr, 2L, "year"),
  stability_rt(eval_det |> dplyr::filter(REGION == "NCR"), event_key_cta, 3L, "year"),
  stability_rt(eval_det |> dplyr::filter(REGION == "NCR"), event_key_ot, 4L, "year")
)
readr::write_csv(loro_all, file.path(OUT_DIR, "RF_Threshold_LORO_Stability.csv"), na = "")
readr::write_csv(loyo_ncr, file.path(OUT_DIR, "RF_Threshold_LOYO_Stability_NCR.csv"), na = "")

add_bounds <- function(summary_tbl, stability_tbl, method_label) {
  out <- summary_tbl
  out$Threshold_Lower <- NA_real_
  out$Threshold_Upper <- NA_real_
  out$Precautionary_Activation_Bound <- NA_real_
  out$Strong_Activation_Bound <- NA_real_
  out$Bounds_Method <- method_label
  out$N_LORO_Thresholds <- 0L
  out$N_Stability_Thresholds <- 0L
  for (i in seq_len(nrow(out))) {
    pid <- out$Product_ID[i]
    pt <- out$Adopted_Threshold[i]
    vals <- stability_tbl$Threshold[stability_tbl$Product_ID == pid]
    vals <- suppressWarnings(as.numeric(vals))
    vals <- vals[is.finite(vals)]
    if (length(vals) >= 2L) {
      b <- safe_quantile(vals, c(0.025, 0.975))
      lo <- min(b[1], pt)
      hi <- max(b[2], pt)
    } else if (length(vals) == 1L) {
      lo <- min(vals[1], pt)
      hi <- max(vals[1], pt)
    } else {
      lo <- pt
      hi <- pt
    }
    out$Threshold_Lower[i] <- lo
    out$Threshold_Upper[i] <- hi
    out$Precautionary_Activation_Bound[i] <- lo
    out$Strong_Activation_Bound[i] <- hi
    out$N_LORO_Thresholds[i] <- length(vals)
    out$N_Stability_Thresholds[i] <- length(vals)
  }
  out
}

summary_regional <- dplyr::bind_rows(reg_p1, reg_p2, reg_p3$summary, reg_p4$summary) |>
  dplyr::arrange(Product_ID) |>
  add_bounds(loro_all,
             "2.5th-97.5th percentiles of A1/A2 false-alarm-aware LORO upper-tail thresholds with RF STA/LTA window re-selection") |>
  dplyr::left_join(
    stagec_support |> dplyr::filter(Scope == "Regional (pooled 17 regions)"),
    by = c("Scope", "Product_ID", "Anchor"))
summary_ncr <- dplyr::bind_rows(ncr_p1, ncr_p2, ncr_p3$summary, ncr_p4$summary) |>
  dplyr::arrange(Product_ID) |>
  add_bounds(loyo_ncr,
             "2.5th-97.5th percentiles of A1/A2 false-alarm-aware LOYO upper-tail thresholds with RF STA/LTA window re-selection") |>
  dplyr::left_join(
    stagec_support |> dplyr::filter(Scope == "NCR"),
    by = c("Scope", "Product_ID", "Anchor"))

for (zname in c("summary_regional", "summary_ncr")) {
  z <- get(zname)
  if (!identical(as.integer(z$Product_ID), 1:4) ||
      any(!is.finite(z$Adopted_Threshold)) ||
      any(z$Comparison_Operator != ">=") ||
      any(!is.finite(z$Threshold_Lower)) ||
      any(!is.finite(z$Threshold_Upper)) ||
      any(z$Threshold_Lower > z$Adopted_Threshold) ||
      any(z$Threshold_Upper < z$Adopted_Threshold)) {
    stop("Stage C failed to derive four finite upper-tail RF products for ", zname, ".",
         call. = FALSE)
  }
}

summary_regional <- summary_regional |>
  dplyr::mutate(
    Primary_Operational_Rule = paste0("RF >= point estimate ",
                                      format(round(Adopted_Threshold, 6), trim = TRUE), " ", Units),
    Precautionary_Activation_Rule = paste0("RF >= lower bound ",
                                            format(round(Threshold_Lower, 6), trim = TRUE), " ", Units),
    Strong_Activation_Rule = paste0("RF >= upper bound ",
                                     format(round(Threshold_Upper, 6), trim = TRUE), " ", Units)
  )
summary_ncr <- summary_ncr |>
  dplyr::mutate(
    Primary_Operational_Rule = paste0("RF >= point estimate ",
                                      format(round(Adopted_Threshold, 6), trim = TRUE), " ", Units),
    Precautionary_Activation_Rule = paste0("RF >= lower bound ",
                                            format(round(Threshold_Lower, 6), trim = TRUE), " ", Units),
    Strong_Activation_Rule = paste0("RF >= upper bound ",
                                     format(round(Threshold_Upper, 6), trim = TRUE), " ", Units)
  )


# -----------------------------------------------------------------------------
# 5b. INDEPENDENT PER-REGION RF THRESHOLDS
# -----------------------------------------------------------------------------
# These thresholds are distinct from the pooled Regional thresholds. Each
# region is calibrated independently using the same first-A1/A2-qualifying
# trigger definition as NCR. P1/P3 share Constant-TA truth; P2/P4 share
# Outbreak-Threshold truth. Bounds are region-specific LOYO stability bounds.
derive_independent_region_bundle <- function(region_name) {
  scope_label <- paste0("Independent region: ", region_name)
  dreg <- det |> dplyr::filter(REGION == region_name)
  ereg <- eval_det |> dplyr::filter(REGION == region_name)
  mm1 <- roc_mm_cta |> dplyr::filter(REGION == region_name)
  mm2 <- roc_mm_ot  |> dplyr::filter(REGION == region_name)
  keys1 <- event_key_cta[startsWith(event_key_cta, paste0(region_name, "|"))]
  keys2 <- event_key_ot[startsWith(event_key_ot, paste0(region_name, "|"))]

  p1 <- mm_summary_one(mm1, 1L, "RF mm preceding Constant TA trigger", scope_label)
  p2 <- mm_summary_one(mm2, 2L, "RF mm preceding Outbreak Threshold trigger", scope_label)
  p3 <- derive_rt_scope(dreg, keys1, 3L, "Constant TA", scope_label)
  p4 <- derive_rt_scope(dreg, keys2, 4L, "Outbreak Threshold", scope_label)

  stab <- dplyr::bind_rows(
    stability_mm(mm1, 1L, "year"),
    stability_mm(mm2, 2L, "year"),
    stability_rt(dreg, keys1, 3L, "year"),
    stability_rt(dreg, keys2, 4L, "year")
  ) |> dplyr::mutate(REGION = region_name, Scope = scope_label, .before = 1)

  support_cta <- stagec_support_one(
    ereg,
    events_cta_all |> dplyr::filter(REGION == region_name),
    events_cta_inventory |> dplyr::filter(REGION == region_name),
    "Constant TA", scope_label)
  support_ot <- stagec_support_one(
    ereg,
    events_ot_all |> dplyr::filter(REGION == region_name),
    events_ot_inventory |> dplyr::filter(REGION == region_name),
    "Outbreak Threshold", scope_label)
  support <- dplyr::bind_rows(
    support_cta |> tidyr::uncount(2, .id = ".rep") |>
      dplyr::mutate(Product_ID = ifelse(.rep == 1L, 1L, 3L)) |> dplyr::select(-.rep),
    support_ot |> tidyr::uncount(2, .id = ".rep") |>
      dplyr::mutate(Product_ID = ifelse(.rep == 1L, 2L, 4L)) |> dplyr::select(-.rep)
  ) |> dplyr::select(Scope, Product_ID, dplyr::everything())

  sm <- dplyr::bind_rows(p1, p2, p3$summary, p4$summary) |>
    dplyr::arrange(Product_ID) |>
    add_bounds(stab,
      "2.5th-97.5th percentiles of region-specific LOYO upper-tail thresholds with RF STA/LTA window re-selection") |>
    dplyr::left_join(support, by = c("Scope", "Product_ID", "Anchor")) |>
    dplyr::mutate(
      REGION = region_name,
      Threshold_Derivation_Scope = "Independent region",
      Primary_Operational_Rule = paste0("RF >= point estimate ",
        format(round(Adopted_Threshold, 6), trim = TRUE), " ", Units),
      Precautionary_Activation_Rule = paste0("RF >= lower bound ",
        format(round(Threshold_Lower, 6), trim = TRUE), " ", Units),
      Strong_Activation_Rule = paste0("RF >= upper bound ",
        format(round(Threshold_Upper, 6), trim = TRUE), " ", Units),
      .before = 1)

  if (!identical(as.integer(sm$Product_ID), 1:4) ||
      any(!is.finite(sm$Adopted_Threshold)) ||
      any(!is.finite(sm$Threshold_Lower)) || any(!is.finite(sm$Threshold_Upper)) ||
      any(sm$Comparison_Operator != ">=") ||
      any(sm$Threshold_Lower > sm$Adopted_Threshold) ||
      any(sm$Threshold_Upper < sm$Adopted_Threshold)) {
    stop("Independent per-region RF derivation failed for ", region_name,
         ": Products 1-4 require finite ordered upper-tail thresholds.", call. = FALSE)
  }
  list(summary = sm, stability = stab,
       candidates = dplyr::bind_rows(p3$candidates, p4$candidates) |>
         dplyr::mutate(REGION = region_name, .before = 1),
       support = support |> dplyr::mutate(REGION = region_name, .before = 1))
}

per_region_bundles <- lapply(CANONICAL_17, derive_independent_region_bundle)
summary_per_region <- dplyr::bind_rows(lapply(per_region_bundles, `[[`, "summary"))
stability_per_region <- dplyr::bind_rows(lapply(per_region_bundles, `[[`, "stability"))
candidates_per_region <- dplyr::bind_rows(lapply(per_region_bundles, `[[`, "candidates"))
support_per_region <- dplyr::bind_rows(lapply(per_region_bundles, `[[`, "support"))

if (nrow(summary_per_region) != 17L * 4L ||
    dplyr::n_distinct(summary_per_region$REGION) != 17L ||
    any(table(summary_per_region$REGION) != 4L)) {
  stop("Independent per-region RF derivation must produce exactly four RF products for all 17 regions.",
       call. = FALSE)
}
readr::write_csv(summary_per_region,
                 file.path(OUT_DIR, "RF_Threshold_Derivation_Summary_PerRegion.csv"), na = "")
readr::write_csv(stability_per_region,
                 file.path(OUT_DIR, "RF_Threshold_LOYO_Stability_PerRegion.csv"), na = "")
readr::write_csv(candidates_per_region,
                 file.path(OUT_DIR, "RF_STA_LTA_Window_Candidates_PerRegion.csv"), na = "")
readr::write_csv(support_per_region,
                 file.path(OUT_DIR, "RF_Threshold_Derivation_Support_PerRegion.csv"), na = "")

readr::write_csv(summary_regional, file.path(OUT_DIR, "RF_Threshold_Derivation_Summary.csv"), na = "")
readr::write_csv(summary_ncr, file.path(OUT_DIR, "RF_Threshold_Derivation_Summary_NCR.csv"), na = "")
readr::write_csv(dplyr::bind_rows(summary_ncr, summary_regional),
                 file.path(OUT_DIR, "RF_Threshold_Derivation_Summary_AllScopes.csv"), na = "")
readr::write_csv(summary_regional |>
                   dplyr::select(Scope, Product_ID, Product, Anchor, RF_Scale, Units,
                                 Threshold_Lower, Adopted_Threshold, Threshold_Upper,
                                 Precautionary_Activation_Bound, Strong_Activation_Bound,
                                 Comparison_Operator, Primary_Operational_Rule,
                                 Precautionary_Activation_Rule, Strong_Activation_Rule,
                                 Bounds_Method, N_Stability_Thresholds),
                 file.path(OUT_DIR, "RF_Threshold_Operational_Bounds.csv"), na = "")
readr::write_csv(summary_ncr |>
                   dplyr::select(Scope, Product_ID, Product, Anchor, RF_Scale, Units,
                                 Threshold_Lower, Adopted_Threshold, Threshold_Upper,
                                 Precautionary_Activation_Bound, Strong_Activation_Bound,
                                 Comparison_Operator, Primary_Operational_Rule,
                                 Precautionary_Activation_Rule, Strong_Activation_Rule,
                                 Bounds_Method, N_Stability_Thresholds),
                 file.path(OUT_DIR, "RF_Threshold_Operational_Bounds_NCR.csv"), na = "")

# -----------------------------------------------------------------------------
# 6. BOUND SENSITIVITY INCLUDING FALSE-ALARM BURDEN
# -----------------------------------------------------------------------------
rf_perf_at_threshold <- function(dat, threshold) {
  d <- dat |>
    dplyr::filter(is.finite(Predictor), Event %in% c(0L, 1L))
  if (!nrow(d)) {
    return(data.frame(
      N = 0L, Sensitivity = NA_real_, Specificity = NA_real_,
      Balanced_Accuracy = NA_real_, PPV = NA_real_, False_Alarm_Rate = NA_real_,
      False_Alarm_Weeks = NA_integer_, False_Alarm_Weeks_per_Evaluable_Unit = NA_real_,
      Trigger_Rate = NA_real_, stringsAsFactors = FALSE
    ))
  }
  pred <- d$Predictor >= threshold
  tp <- sum(pred & d$Event == 1L)
  fp <- sum(pred & d$Event == 0L)
  fn <- sum(!pred & d$Event == 1L)
  tn <- sum(!pred & d$Event == 0L)
  sens <- if ((tp + fn) > 0L) tp / (tp + fn) else NA_real_
  spec <- if ((tn + fp) > 0L) tn / (tn + fp) else NA_real_
  ppv <- if ((tp + fp) > 0L) tp / (tp + fp) else NA_real_
  n_units <- d |> dplyr::distinct(REGION, YR) |> nrow()
  if (n_units < 1L) n_units <- 1L
  data.frame(
    N = nrow(d), Sensitivity = sens, Specificity = spec,
    Balanced_Accuracy = if (is.finite(sens) && is.finite(spec)) (sens + spec) / 2 else sens,
    PPV = ppv,
    False_Alarm_Rate = if ((fp + tn) > 0L) fp / (fp + tn) else 0,
    False_Alarm_Weeks = fp,
    False_Alarm_Weeks_per_Evaluable_Unit = fp / n_units,
    Trigger_Rate = mean(pred), stringsAsFactors = FALSE
  )
}

roc_inputs_regional <- list(
  `1` = roc_mm_cta, `2` = roc_mm_ot,
  `3` = reg_p3$roc_input, `4` = reg_p4$roc_input
)
roc_inputs_ncr <- list(
  `1` = roc_mm_cta_ncr, `2` = roc_mm_ot_ncr,
  `3` = ncr_p3$roc_input, `4` = ncr_p4$roc_input
)

bound_perf_for_scope <- function(summary_tbl, inputs, scope_label) {
  dplyr::bind_rows(lapply(1:4, function(pid) {
    ss <- summary_tbl |> dplyr::filter(Product_ID == pid) |> dplyr::slice(1)
    levels_tbl <- data.frame(
      Bound = c("Lower", "Point", "Upper"),
      Operational_Role = c("Precautionary RF trigger", "Formal RF trigger", "Strong RF trigger"),
      Threshold = c(ss$Threshold_Lower, ss$Adopted_Threshold, ss$Threshold_Upper),
      stringsAsFactors = FALSE
    )
    dplyr::bind_rows(lapply(seq_len(nrow(levels_tbl)), function(j) {
      pp <- rf_perf_at_threshold(inputs[[as.character(pid)]], levels_tbl$Threshold[j])
      cbind(data.frame(
        Scope = scope_label, Product_ID = pid, RF_Scale = ss$RF_Scale,
        Comparison_Operator = ">=", Bound = levels_tbl$Bound[j],
        Operational_Role = levels_tbl$Operational_Role[j],
        Threshold = levels_tbl$Threshold[j], stringsAsFactors = FALSE
      ), pp)
    }))
  }))
}

rf_bound_perf_regional <- bound_perf_for_scope(summary_regional, roc_inputs_regional,
                                                "Regional (pooled 17 regions)")
rf_bound_perf_ncr <- bound_perf_for_scope(summary_ncr, roc_inputs_ncr, "NCR")
readr::write_csv(rf_bound_perf_regional,
                 file.path(OUT_DIR, "RF_Threshold_Bound_Sensitivity.csv"), na = "")
readr::write_csv(rf_bound_perf_ncr,
                 file.path(OUT_DIR, "RF_Threshold_Bound_Sensitivity_NCR.csv"), na = "")
readr::write_csv(dplyr::bind_rows(rf_bound_perf_ncr, rf_bound_perf_regional),
                 file.path(OUT_DIR, "RF_Threshold_Bound_Sensitivity_AllScopes.csv"), na = "")

# -----------------------------------------------------------------------------
# 7. DIRECT STAGE-C HEAD-TO-HEAD: RF mm VS RF STA/LTA R(t)
# -----------------------------------------------------------------------------
stagec_h2h_metric_names <- c(
  "AUC", "Sensitivity", "Specificity", "Balanced_Accuracy", "PPV", "False_Alarm_Rate"
)

stagec_h2h_metrics <- function(event, predictor, threshold) {
  event <- suppressWarnings(as.integer(event))
  predictor <- suppressWarnings(as.numeric(predictor))
  keep <- is.finite(predictor) & event %in% c(0L, 1L)
  event <- event[keep]
  predictor <- predictor[keep]
  if (!length(event)) {
    return(data.frame(
      N = 0L, N_Events = 0L, AUC = NA_real_, AUC_CI_Lower = NA_real_,
      AUC_CI_Upper = NA_real_, Sensitivity = NA_real_, Specificity = NA_real_,
      Balanced_Accuracy = NA_real_, PPV = NA_real_, False_Alarm_Rate = NA_real_,
      Trigger_Rate = NA_real_, stringsAsFactors = FALSE
    ))
  }
  rr <- if (length(unique(event)) >= 2L) tryCatch(
    pROC::roc(event, predictor, direction = "<", quiet = TRUE),
    error = function(e) NULL) else NULL
  auc <- if (is.null(rr)) NA_real_ else as.numeric(pROC::auc(rr))
  ci_auc <- if (is.null(rr)) c(NA_real_, NA_real_) else {
    z <- tryCatch(as.numeric(pROC::ci.auc(rr, conf.level = 0.95)),
                  error = function(e) c(NA_real_, NA_real_, NA_real_))
    if (length(z) >= 3L) c(z[1L], z[3L]) else c(NA_real_, NA_real_)
  }
  pred <- predictor >= threshold
  tp <- sum(pred & event == 1L)
  fp <- sum(pred & event == 0L)
  fn <- sum(!pred & event == 1L)
  tn <- sum(!pred & event == 0L)
  sens <- if ((tp + fn) > 0L) tp / (tp + fn) else NA_real_
  spec <- if ((tn + fp) > 0L) tn / (tn + fp) else NA_real_
  ppv <- if ((tp + fp) > 0L) tp / (tp + fp) else NA_real_
  data.frame(
    N = length(event), N_Events = sum(event == 1L),
    AUC = auc, AUC_CI_Lower = ci_auc[1L], AUC_CI_Upper = ci_auc[2L],
    Sensitivity = sens, Specificity = spec,
    Balanced_Accuracy = if (is.finite(sens) && is.finite(spec)) (sens + spec) / 2 else sens,
    PPV = ppv,
    False_Alarm_Rate = if ((fp + tn) > 0L) fp / (fp + tn) else 0,
    Trigger_Rate = mean(pred), stringsAsFactors = FALSE
  )
}

stagec_h2h_exact_signflip <- function(d, exact_max_n = 20L, B = 100000L, seed = 1L) {
  d <- suppressWarnings(as.numeric(d))
  d <- d[is.finite(d)]
  n <- length(d)
  if (!n) return(NA_real_)
  if (all(abs(d) <= sqrt(.Machine$double.eps))) return(1)
  obs <- abs(mean(d))
  tol <- sqrt(.Machine$double.eps)
  if (n <= exact_max_n) {
    idx <- 0:(2^n - 1)
    perm_sum <- numeric(length(idx))
    for (j in seq_len(n)) {
      sgn <- ifelse(bitwAnd(idx, bitwShiftL(1L, j - 1L)) == 0L, -1, 1)
      perm_sum <- perm_sum + sgn * d[j]
    }
    return(mean(abs(perm_sum / n) >= (obs - tol)))
  }
  set.seed(seed)
  perm <- replicate(B, abs(mean(d * sample(c(-1, 1), n, replace = TRUE))))
  (1 + sum(perm >= (obs - tol))) / (B + 1)
}

stagec_h2h_boot_ci <- function(d, B = 5000L, conf = 0.95, seed = 1L) {
  d <- suppressWarnings(as.numeric(d))
  d <- d[is.finite(d)]
  n <- length(d)
  if (!n) return(c(NA_real_, NA_real_))
  if (n == 1L) return(c(d, d))
  set.seed(seed)
  z <- replicate(B, mean(d[sample.int(n, n, replace = TRUE)]))
  alpha <- (1 - conf) / 2
  as.numeric(stats::quantile(z, probs = c(alpha, 1 - alpha),
                             na.rm = TRUE, names = FALSE, type = 8))
}

make_common <- function(mm, rt) {
  mm2 <- mm |>
    dplyr::transmute(REGION, YR, WN, Event,
                     RF_mm = suppressWarnings(as.numeric(Predictor)))
  rt2 <- rt |>
    dplyr::transmute(REGION, YR, WN, Event,
                     RF_Rt = suppressWarnings(as.numeric(Predictor)))
  dplyr::inner_join(mm2, rt2, by = c("REGION", "YR", "WN", "Event")) |>
    dplyr::filter(is.finite(RF_mm), is.finite(RF_Rt), Event %in% c(0L, 1L)) |>
    dplyr::arrange(REGION, YR, WN)
}

scope_h2h_row <- function(common, scope_name, anchor, mm_product, rt_product,
                          summary_tbl) {
  smm <- summary_tbl |> dplyr::filter(Product_ID == mm_product) |> dplyr::slice(1)
  srt <- summary_tbl |> dplyr::filter(Product_ID == rt_product) |> dplyr::slice(1)
  mm_thr <- smm$Adopted_Threshold[1]
  rt_thr <- srt$Adopted_Threshold[1]
  mm_m <- stagec_h2h_metrics(common$Event, common$RF_mm, mm_thr)
  rt_m <- stagec_h2h_metrics(common$Event, common$RF_Rt, rt_thr)
  roc_mm <- tryCatch(pROC::roc(common$Event, common$RF_mm, direction = "<", quiet = TRUE),
                     error = function(e) NULL)
  roc_rt <- tryCatch(pROC::roc(common$Event, common$RF_Rt, direction = "<", quiet = TRUE),
                     error = function(e) NULL)
  p_delong <- if (!is.null(roc_mm) && !is.null(roc_rt)) {
    tryCatch(as.numeric(pROC::roc.test(roc_mm, roc_rt, paired = TRUE,
                                       method = "delong")$p.value),
             error = function(e) NA_real_)
  } else NA_real_
  data.frame(
    Scope = scope_name, Anchor = anchor,
    RF_mm_Product = mm_product, RF_Rt_Product = rt_product,
    RF_mm_Adopted_Threshold = mm_thr, RF_Rt_Adopted_Threshold = rt_thr,
    RF_mm_Operator = ">=", RF_Rt_Operator = ">=",
    N_Common_Weeks = nrow(common),
    N_Events = dplyr::n_distinct(paste(common$REGION[common$Event == 1L], common$YR[common$Event == 1L], sep = "|")),
    N_Positive_Warning_Weeks = sum(common$Event == 1L),
    AUC_RF_mm = mm_m$AUC, AUC_RF_mm_CI_Lower = mm_m$AUC_CI_Lower,
    AUC_RF_mm_CI_Upper = mm_m$AUC_CI_Upper,
    AUC_RF_Rt = rt_m$AUC, AUC_RF_Rt_CI_Lower = rt_m$AUC_CI_Lower,
    AUC_RF_Rt_CI_Upper = rt_m$AUC_CI_Upper,
    Delta_AUC_Rt_minus_mm = rt_m$AUC - mm_m$AUC,
    p_Paired_DeLong = p_delong,
    Sensitivity_RF_mm = mm_m$Sensitivity, Sensitivity_RF_Rt = rt_m$Sensitivity,
    Specificity_RF_mm = mm_m$Specificity, Specificity_RF_Rt = rt_m$Specificity,
    Balanced_Accuracy_RF_mm = mm_m$Balanced_Accuracy,
    Balanced_Accuracy_RF_Rt = rt_m$Balanced_Accuracy,
    PPV_RF_mm = mm_m$PPV, PPV_RF_Rt = rt_m$PPV,
    False_Alarm_Rate_RF_mm = mm_m$False_Alarm_Rate,
    False_Alarm_Rate_RF_Rt = rt_m$False_Alarm_Rate,
    Trigger_Rate_RF_mm = mm_m$Trigger_Rate, Trigger_Rate_RF_Rt = rt_m$Trigger_Rate,
    stringsAsFactors = FALSE
  )
}

h2h_pair <- function(anchor, mm_product, rt_product) {
  reg_common <- make_common(roc_inputs_regional[[as.character(mm_product)]],
                            roc_inputs_regional[[as.character(rt_product)]])
  ncr_common <- make_common(roc_inputs_ncr[[as.character(mm_product)]],
                            roc_inputs_ncr[[as.character(rt_product)]])
  if (!nrow(reg_common) || !nrow(ncr_common)) {
    stop("Stage C RF head-to-head has insufficient paired rows for ", anchor, ".",
         call. = FALSE)
  }
  overall <- dplyr::bind_rows(
    scope_h2h_row(reg_common, "Pooled 17 regions", anchor, mm_product, rt_product,
                  summary_regional),
    scope_h2h_row(ncr_common, "NCR", anchor, mm_product, rt_product,
                  summary_ncr)
  )

  smm <- summary_regional |> dplyr::filter(Product_ID == mm_product) |> dplyr::slice(1)
  srt <- summary_regional |> dplyr::filter(Product_ID == rt_product) |> dplyr::slice(1)
  by_region <- dplyr::bind_rows(lapply(CANONICAL_17, function(r) {
    z <- reg_common |> dplyr::filter(REGION == r)
    mm_m <- stagec_h2h_metrics(z$Event, z$RF_mm, smm$Adopted_Threshold[1])
    rt_m <- stagec_h2h_metrics(z$Event, z$RF_Rt, srt$Adopted_Threshold[1])
    data.frame(
      Anchor = anchor, REGION = r,
      N_Common_Weeks = nrow(z),
      N_Events = dplyr::n_distinct(z$YR[z$Event == 1L]),
      N_Positive_Warning_Weeks = sum(z$Event == 1L),
      AUC_RF_mm = mm_m$AUC, AUC_RF_Rt = rt_m$AUC,
      Sensitivity_RF_mm = mm_m$Sensitivity, Sensitivity_RF_Rt = rt_m$Sensitivity,
      Specificity_RF_mm = mm_m$Specificity, Specificity_RF_Rt = rt_m$Specificity,
      Balanced_Accuracy_RF_mm = mm_m$Balanced_Accuracy,
      Balanced_Accuracy_RF_Rt = rt_m$Balanced_Accuracy,
      PPV_RF_mm = mm_m$PPV, PPV_RF_Rt = rt_m$PPV,
      False_Alarm_Rate_RF_mm = mm_m$False_Alarm_Rate,
      False_Alarm_Rate_RF_Rt = rt_m$False_Alarm_Rate,
      stringsAsFactors = FALSE
    )
  }))
  list(reg_common = reg_common, ncr_common = ncr_common,
       overall = overall, by_region = by_region)
}

h2h_cta <- h2h_pair("Constant TA", 1L, 3L)
h2h_ot <- h2h_pair("Outbreak Threshold", 2L, 4L)
stagec_h2h_common <- dplyr::bind_rows(
  h2h_cta$reg_common |> dplyr::mutate(Scope = "Pooled 17 regions", Anchor = "Constant TA", .before = 1),
  h2h_cta$ncr_common |> dplyr::mutate(Scope = "NCR", Anchor = "Constant TA", .before = 1),
  h2h_ot$reg_common |> dplyr::mutate(Scope = "Pooled 17 regions", Anchor = "Outbreak Threshold", .before = 1),
  h2h_ot$ncr_common |> dplyr::mutate(Scope = "NCR", Anchor = "Outbreak Threshold", .before = 1)
)
stagec_h2h_overall <- dplyr::bind_rows(h2h_cta$overall, h2h_ot$overall) |>
  dplyr::group_by(Scope) |>
  dplyr::mutate(
    p_Paired_DeLong_Bonferroni = pmin(1, p_Paired_DeLong * sum(is.finite(p_Paired_DeLong))),
    AUC_Significant_005 = is.finite(p_Paired_DeLong_Bonferroni) & p_Paired_DeLong_Bonferroni < 0.05
  ) |>
  dplyr::ungroup()
stagec_h2h_region <- dplyr::bind_rows(h2h_cta$by_region, h2h_ot$by_region)

stagec_h2h_region_long <- stagec_h2h_region |>
  tidyr::pivot_longer(
    cols = -c(Anchor, REGION, N_Common_Weeks, N_Events),
    names_to = c("Metric", "Representation"),
    names_pattern = "(.*)_RF_(mm|Rt)$",
    values_to = "Value"
  ) |>
  dplyr::filter(Metric %in% stagec_h2h_metric_names) |>
  tidyr::pivot_wider(names_from = Representation, values_from = Value,
                     names_prefix = "RF_")

stagec_h2h_tests <- stagec_h2h_region_long |>
  dplyr::group_by(Anchor, Metric) |>
  dplyr::group_modify(~ {
    z <- .x |> dplyr::filter(is.finite(RF_mm), is.finite(RF_Rt))
    d <- z$RF_Rt - z$RF_mm
    seed <- as.integer((sum(utf8ToInt(paste(.y$Anchor, .y$Metric))) %% 1000000L) + 3001L)
    ci <- stagec_h2h_boot_ci(d, B = 5000L, seed = seed)
    lower_better <- identical(as.character(.y$Metric[1]), "False_Alarm_Rate")
    mean_diff <- if (nrow(z)) mean(d) else NA_real_
    preferred <- if (!is.finite(mean_diff) || mean_diff == 0) {
      "Tie"
    } else if (lower_better) {
      if (mean_diff < 0) "RF STA/LTA R(t)" else "RF mm"
    } else {
      if (mean_diff > 0) "RF STA/LTA R(t)" else "RF mm"
    }
    data.frame(
      N_Paired_Regions = nrow(z),
      Mean_RF_mm = if (nrow(z)) mean(z$RF_mm) else NA_real_,
      Mean_RF_Rt = if (nrow(z)) mean(z$RF_Rt) else NA_real_,
      Mean_Difference_Rt_minus_mm = mean_diff,
      Difference_CI_Lower = ci[1L], Difference_CI_Upper = ci[2L],
      p_exact_signflip = stagec_h2h_exact_signflip(d, seed = seed + 19L),
      Preferred = preferred, stringsAsFactors = FALSE
    )
  }) |>
  dplyr::ungroup() |>
  dplyr::group_by(Anchor) |>
  dplyr::mutate(
    p_bonferroni = pmin(1, p_exact_signflip * sum(is.finite(p_exact_signflip))),
    Significant_005 = is.finite(p_bonferroni) & p_bonferroni < 0.05
  ) |>
  dplyr::ungroup()

stagec_h2h_loro <- dplyr::bind_rows(
  loro_all |> dplyr::filter(Product_ID == 1L) |>
    dplyr::transmute(Anchor = "Constant TA", Held_Out_Region,
                     RF_mm_Threshold = Threshold, RF_mm_AUC = AUC),
  loro_all |> dplyr::filter(Product_ID == 2L) |>
    dplyr::transmute(Anchor = "Outbreak Threshold", Held_Out_Region,
                     RF_mm_Threshold = Threshold, RF_mm_AUC = AUC)
) |>
  dplyr::left_join(
    dplyr::bind_rows(
      loro_all |> dplyr::filter(Product_ID == 3L) |>
        dplyr::transmute(Anchor = "Constant TA", Held_Out_Region,
                         RF_Rt_Threshold = Threshold, RF_Rt_AUC = AUC),
      loro_all |> dplyr::filter(Product_ID == 4L) |>
        dplyr::transmute(Anchor = "Outbreak Threshold", Held_Out_Region,
                         RF_Rt_Threshold = Threshold, RF_Rt_AUC = AUC)
    ), by = c("Anchor", "Held_Out_Region")) |>
  dplyr::mutate(AUC_Difference_Rt_minus_mm = RF_Rt_AUC - RF_mm_AUC)

stagec_h2h_loro_tests <- stagec_h2h_loro |>
  dplyr::group_by(Anchor) |>
  dplyr::group_modify(~ {
    d <- .x$AUC_Difference_Rt_minus_mm
    d <- d[is.finite(d)]
    seed <- as.integer((sum(utf8ToInt(.y$Anchor)) %% 1000000L) + 5001L)
    ci <- stagec_h2h_boot_ci(d, B = 5000L, seed = seed)
    data.frame(
      N_Paired_LORO = length(d),
      Mean_AUC_Difference_Rt_minus_mm = if (length(d)) mean(d) else NA_real_,
      Difference_CI_Lower = ci[1L], Difference_CI_Upper = ci[2L],
      p_exact_signflip = stagec_h2h_exact_signflip(d, seed = seed + 31L),
      Preferred = if (!length(d)) NA_character_ else if (mean(d) > 0) "RF STA/LTA R(t)" else if (mean(d) < 0) "RF mm" else "Tie",
      stringsAsFactors = FALSE
    )
  }) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    p_bonferroni = pmin(1, p_exact_signflip * sum(is.finite(p_exact_signflip))),
    Significant_005 = is.finite(p_bonferroni) & p_bonferroni < 0.05
  )

stagec_h2h_decision <- dplyr::bind_rows(lapply(c("Constant TA", "Outbreak Threshold"), function(a) {
  pooled <- stagec_h2h_overall |> dplyr::filter(Anchor == a, Scope == "Pooled 17 regions") |> dplyr::slice(1)
  ncr <- stagec_h2h_overall |> dplyr::filter(Anchor == a, Scope == "NCR") |> dplyr::slice(1)
  far_test <- stagec_h2h_tests |> dplyr::filter(Anchor == a, Metric == "False_Alarm_Rate") |> dplyr::slice(1)
  evidence <- character()
  if (nrow(pooled) && isTRUE(pooled$AUC_Significant_005[1])) evidence <- c(evidence, "pooled paired DeLong AUC")
  if (nrow(ncr) && isTRUE(ncr$AUC_Significant_005[1])) evidence <- c(evidence, "NCR paired DeLong AUC")
  if (nrow(far_test) && isTRUE(far_test$Significant_005[1])) evidence <- c(evidence, "17-region false-alarm rate")
  data.frame(
    Anchor = a,
    Derivation_Stage_Conclusion = if (!length(evidence)) "No clear empirical superiority" else "See significant evidence; final adoption requires Stage H operational performance",
    Significant_Evidence = if (length(evidence)) paste(evidence, collapse = "; ") else "None after specified tests",
    Note = "Stage C compares upper-tail threshold derivation/discrimination including false-alarm rate; Stage H evaluates downstream operational performance.",
    stringsAsFactors = FALSE
  )
}))

readr::write_csv(stagec_h2h_common, file.path(OUT_DIR, "RF_Threshold_HeadToHead_CommonRows.csv"), na = "")
readr::write_csv(stagec_h2h_overall, file.path(OUT_DIR, "RF_Threshold_HeadToHead_Overall.csv"), na = "")
readr::write_csv(stagec_h2h_region, file.path(OUT_DIR, "RF_Threshold_HeadToHead_ByRegion.csv"), na = "")
readr::write_csv(stagec_h2h_tests, file.path(OUT_DIR, "RF_Threshold_HeadToHead_ExactPairedTests.csv"), na = "")
readr::write_csv(stagec_h2h_loro, file.path(OUT_DIR, "RF_Threshold_HeadToHead_LORO.csv"), na = "")
readr::write_csv(stagec_h2h_loro_tests, file.path(OUT_DIR, "RF_Threshold_HeadToHead_LORO_Tests.csv"), na = "")
readr::write_csv(stagec_h2h_decision, file.path(OUT_DIR, "RF_Threshold_HeadToHead_Decision.csv"), na = "")

# -----------------------------------------------------------------------------
# 8. SCOPE-SPECIFIC HANDOFF
# -----------------------------------------------------------------------------
value_for <- function(tbl, pid, col) tbl[[col]][tbl$Product_ID == pid][1]
writeLines(c(
  "# Auto-generated by StageC_RF_Threshold_Derivation.R -- do not edit.",
  "# All RF products are upper-tail operational activations (RF >= threshold).",
  "# Regional variables retain legacy names; NCR variables carry _NCR suffix.",
  sprintf("RF_THRESHOLD_PRODUCT1_MM <- %.12f", value_for(summary_regional, 1L, "Adopted_Threshold")),
  sprintf("RF_THRESHOLD_PRODUCT2_MM <- %.12f", value_for(summary_regional, 2L, "Adopted_Threshold")),
  sprintf("RF_THRESHOLD_PRODUCT3_RT <- %.12f", value_for(summary_regional, 3L, "Adopted_Threshold")),
  sprintf("RF_THRESHOLD_PRODUCT4_RT <- %.12f", value_for(summary_regional, 4L, "Adopted_Threshold")),
  sprintf("RF_THRESHOLD_LOWER_PRODUCT1_MM <- %.12f", value_for(summary_regional, 1L, "Threshold_Lower")),
  sprintf("RF_THRESHOLD_UPPER_PRODUCT1_MM <- %.12f", value_for(summary_regional, 1L, "Threshold_Upper")),
  sprintf("RF_THRESHOLD_LOWER_PRODUCT2_MM <- %.12f", value_for(summary_regional, 2L, "Threshold_Lower")),
  sprintf("RF_THRESHOLD_UPPER_PRODUCT2_MM <- %.12f", value_for(summary_regional, 2L, "Threshold_Upper")),
  sprintf("RF_THRESHOLD_LOWER_PRODUCT3_RT <- %.12f", value_for(summary_regional, 3L, "Threshold_Lower")),
  sprintf("RF_THRESHOLD_UPPER_PRODUCT3_RT <- %.12f", value_for(summary_regional, 3L, "Threshold_Upper")),
  sprintf("RF_THRESHOLD_LOWER_PRODUCT4_RT <- %.12f", value_for(summary_regional, 4L, "Threshold_Lower")),
  sprintf("RF_THRESHOLD_UPPER_PRODUCT4_RT <- %.12f", value_for(summary_regional, 4L, "Threshold_Upper")),
  sprintf("RF_STA_PRODUCT3 <- %dL", value_for(summary_regional, 3L, "RF_STA")),
  sprintf("RF_LTA_PRODUCT3 <- %dL", value_for(summary_regional, 3L, "RF_LTA")),
  sprintf("RF_GUARD_PRODUCT3 <- %dL", value_for(summary_regional, 3L, "RF_Guard")),
  sprintf("RF_STA_PRODUCT4 <- %dL", value_for(summary_regional, 4L, "RF_STA")),
  sprintf("RF_LTA_PRODUCT4 <- %dL", value_for(summary_regional, 4L, "RF_LTA")),
  sprintf("RF_GUARD_PRODUCT4 <- %dL", value_for(summary_regional, 4L, "RF_Guard")),
  'RF_OPERATOR_PRODUCT1 <- ">="', 'RF_OPERATOR_PRODUCT2 <- ">="',
  'RF_OPERATOR_PRODUCT3 <- ">="', 'RF_OPERATOR_PRODUCT4 <- ">="',
  sprintf("RF_THRESHOLD_PRODUCT1_MM_NCR <- %.12f", value_for(summary_ncr, 1L, "Adopted_Threshold")),
  sprintf("RF_THRESHOLD_PRODUCT2_MM_NCR <- %.12f", value_for(summary_ncr, 2L, "Adopted_Threshold")),
  sprintf("RF_THRESHOLD_PRODUCT3_RT_NCR <- %.12f", value_for(summary_ncr, 3L, "Adopted_Threshold")),
  sprintf("RF_THRESHOLD_PRODUCT4_RT_NCR <- %.12f", value_for(summary_ncr, 4L, "Adopted_Threshold")),
  sprintf("RF_THRESHOLD_LOWER_PRODUCT1_MM_NCR <- %.12f", value_for(summary_ncr, 1L, "Threshold_Lower")),
  sprintf("RF_THRESHOLD_UPPER_PRODUCT1_MM_NCR <- %.12f", value_for(summary_ncr, 1L, "Threshold_Upper")),
  sprintf("RF_THRESHOLD_LOWER_PRODUCT2_MM_NCR <- %.12f", value_for(summary_ncr, 2L, "Threshold_Lower")),
  sprintf("RF_THRESHOLD_UPPER_PRODUCT2_MM_NCR <- %.12f", value_for(summary_ncr, 2L, "Threshold_Upper")),
  sprintf("RF_THRESHOLD_LOWER_PRODUCT3_RT_NCR <- %.12f", value_for(summary_ncr, 3L, "Threshold_Lower")),
  sprintf("RF_THRESHOLD_UPPER_PRODUCT3_RT_NCR <- %.12f", value_for(summary_ncr, 3L, "Threshold_Upper")),
  sprintf("RF_THRESHOLD_LOWER_PRODUCT4_RT_NCR <- %.12f", value_for(summary_ncr, 4L, "Threshold_Lower")),
  sprintf("RF_THRESHOLD_UPPER_PRODUCT4_RT_NCR <- %.12f", value_for(summary_ncr, 4L, "Threshold_Upper")),
  sprintf("RF_STA_PRODUCT3_NCR <- %dL", value_for(summary_ncr, 3L, "RF_STA")),
  sprintf("RF_LTA_PRODUCT3_NCR <- %dL", value_for(summary_ncr, 3L, "RF_LTA")),
  sprintf("RF_GUARD_PRODUCT3_NCR <- %dL", value_for(summary_ncr, 3L, "RF_Guard")),
  sprintf("RF_STA_PRODUCT4_NCR <- %dL", value_for(summary_ncr, 4L, "RF_STA")),
  sprintf("RF_LTA_PRODUCT4_NCR <- %dL", value_for(summary_ncr, 4L, "RF_LTA")),
  sprintf("RF_GUARD_PRODUCT4_NCR <- %dL", value_for(summary_ncr, 4L, "RF_Guard")),
  'RF_OPERATOR_PRODUCT1_NCR <- ">="', 'RF_OPERATOR_PRODUCT2_NCR <- ">="',
  'RF_OPERATOR_PRODUCT3_NCR <- ">="', 'RF_OPERATOR_PRODUCT4_NCR <- ">="'
), file.path(OUT_DIR, "rf_thresholds_derived.R"))

# -----------------------------------------------------------------------------
# 9. STAGE C FIGURES: SCOPE-SEPARATED, HORIZONTAL X-AXES
# -----------------------------------------------------------------------------
# NCR and Regional figures are deliberately exported as separate files. No
# mixed-scope Stage C figure is produced, so each scope retains adequate axis,
# label and legend space.
product_labels_mm <- c(
  "P1\nRF mm preceding\nConstant TA",
  "P2\nRF mm preceding\nOutbreak Threshold"
)
product_labels_rt <- c(
  "P3\nRF STA/LTA R(t) preceding\nConstant TA",
  "P4\nRF STA/LTA R(t) preceding\nOutbreak Threshold"
)

make_event_distribution_panel <- function(scope_name, mm_events, rt3_events, rt4_events) {
  mm_plot <- dplyr::bind_rows(
    mm_events[[1]] |> dplyr::mutate(Product = product_labels_mm[1]),
    mm_events[[2]] |> dplyr::mutate(Product = product_labels_mm[2])
  ) |> dplyr::filter(is.finite(RF_Antecedent_4wk_mm))
  mm_plot$Product <- factor(mm_plot$Product, levels = product_labels_mm)
  p_mm <- ggplot2::ggplot(mm_plot, ggplot2::aes(x = Product, y = RF_Antecedent_4wk_mm)) +
    ggplot2::geom_boxplot(outlier.shape = NA, width = 0.55) +
    ggplot2::geom_jitter(width = 0.10, height = 0, alpha = 0.50, size = 1.6) +
    ggplot2::labs(x = NULL, y = "Antecedent rainfall (mm)",
                  title = paste0(scope_name, ": RF mm preceding first A1/A2-qualifying trigger")) +
    theme_pub() +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, vjust = 1),
      plot.margin = ggplot2::margin(8, 12, 12, 12))

  rt_plot <- dplyr::bind_rows(
    rt3_events |> dplyr::mutate(Product = product_labels_rt[1]),
    rt4_events |> dplyr::mutate(Product = product_labels_rt[2])
  ) |> dplyr::filter(is.finite(RF_Antecedent_4wk_Rt_Max))
  rt_plot$Product <- factor(rt_plot$Product, levels = product_labels_rt)
  p_rt <- ggplot2::ggplot(rt_plot, ggplot2::aes(x = Product, y = RF_Antecedent_4wk_Rt_Max)) +
    ggplot2::geom_boxplot(outlier.shape = NA, width = 0.55) +
    ggplot2::geom_jitter(width = 0.10, height = 0, alpha = 0.50, size = 1.6) +
    ggplot2::labs(x = NULL, y = "Antecedent rainfall STA/LTA R(t)",
                  title = paste0(scope_name, ": RF STA/LTA preceding first A1/A2-qualifying trigger")) +
    theme_pub() +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, vjust = 1),
      plot.margin = ggplot2::margin(8, 12, 12, 12))
  list(mm = p_mm, rt = p_rt)
}

reg_mm_cta_pre <- roc_mm_cta |> dplyr::filter(Event == 1L) |>
  dplyr::mutate(RF_Antecedent_4wk_mm = Predictor)
reg_mm_ot_pre <- roc_mm_ot |> dplyr::filter(Event == 1L) |>
  dplyr::mutate(RF_Antecedent_4wk_mm = Predictor)
ncr_mm_cta_pre <- roc_mm_cta_ncr |> dplyr::filter(Event == 1L) |>
  dplyr::mutate(RF_Antecedent_4wk_mm = Predictor)
ncr_mm_ot_pre <- roc_mm_ot_ncr |> dplyr::filter(Event == 1L) |>
  dplyr::mutate(RF_Antecedent_4wk_mm = Predictor)

fig_products_reg <- make_event_distribution_panel(
  "Regional", list(reg_mm_cta_pre, reg_mm_ot_pre),
  reg_p3$roc_input |> dplyr::filter(Event == 1L),
  reg_p4$roc_input |> dplyr::filter(Event == 1L))
fig_products_ncr <- make_event_distribution_panel(
  "NCR", list(ncr_mm_cta_pre, ncr_mm_ot_pre),
  ncr_p3$roc_input |> dplyr::filter(Event == 1L),
  ncr_p4$roc_input |> dplyr::filter(Event == 1L))

fig_products_ncr_full <- patchwork::wrap_plots(
  list(fig_products_ncr$mm, fig_products_ncr$rt), ncol = 1, nrow = 2)
fig_products_reg_full <- patchwork::wrap_plots(
  list(fig_products_reg$mm, fig_products_reg$rt), ncol = 1, nrow = 2)
save_pub("FigureC_NCR_RF_Threshold_Products", fig_products_ncr_full,
         NC_W_DOUBLE, 8.0, OUT_DIR)
save_pub("FigureC_Regional_RF_Threshold_Products", fig_products_reg_full,
         NC_W_DOUBLE, 8.0, OUT_DIR)

window_plot <- dplyr::bind_rows(reg_p3$candidates, reg_p4$candidates,
                                ncr_p3$candidates, ncr_p4$candidates) |>
  dplyr::mutate(
    Product = factor(Product_ID, levels = c(3, 4),
                     labels = c("P3: RF STA/LTA preceding\nConstant TA",
                                "P4: RF STA/LTA preceding\nOutbreak Threshold")),
    Window = paste0(STA, "/", LTA)
  )

make_window_plot <- function(scope_name, title_scope) {
  z <- window_plot |> dplyr::filter(Scope == scope_name)
  ggplot2::ggplot(z,
                  ggplot2::aes(x = Window, y = Balanced_Accuracy,
                               group = Product, shape = Product)) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::geom_line(linewidth = 0.60) +
    ggplot2::labs(x = "Rainfall STA/LTA window", y = "Balanced accuracy",
                  title = paste0(title_scope, ": rainfall STA/LTA windows preceding disease triggers")) +
    ggplot2::guides(shape = ggplot2::guide_legend(
      nrow = 2, byrow = TRUE, direction = "horizontal")) +
    theme_pub() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
      legend.position = "bottom", legend.direction = "horizontal",
      legend.box = "horizontal", legend.justification = "center",
      plot.margin = ggplot2::margin(8, 12, 14, 12))
}
p_win_ncr <- make_window_plot("NCR", "NCR")
p_win_reg <- make_window_plot("Regional (pooled 17 regions)", "Regional")
save_pub("FigureC_NCR_RF_STA_LTA_Window_Selection", p_win_ncr,
         NC_W_DOUBLE, 5.0, OUT_DIR)
save_pub("FigureC_Regional_RF_STA_LTA_Window_Selection", p_win_reg,
         NC_W_DOUBLE, 5.0, OUT_DIR)

threshold_plot_data <- dplyr::bind_rows(summary_ncr, summary_regional) |>
  dplyr::mutate(
    Product_Label = factor(
      dplyr::case_when(
        Product_ID == 1L ~ "P1\nRF mm preceding\nConstant TA",
        Product_ID == 2L ~ "P2\nRF mm preceding\nOutbreak Threshold",
        Product_ID == 3L ~ "P3\nRF STA/LTA preceding\nConstant TA",
        TRUE ~ "P4\nRF STA/LTA preceding\nOutbreak Threshold"
      ),
      levels = c("P1\nRF mm preceding\nConstant TA",
                 "P2\nRF mm preceding\nOutbreak Threshold",
                 "P3\nRF STA/LTA preceding\nConstant TA",
                 "P4\nRF STA/LTA preceding\nOutbreak Threshold")
    ),
    Scale_Label = ifelse(RF_Scale == "mm", "RF accumulation (mm)",
                         "Rainfall STA/LTA R(t)")
  )

make_bounds_plot <- function(scope_name, title_scope) {
  z <- threshold_plot_data |> dplyr::filter(Scope == scope_name)
  ggplot2::ggplot(z, ggplot2::aes(x = Product_Label, y = Adopted_Threshold)) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = Threshold_Lower, ymax = Threshold_Upper),
      width = 0.15, linewidth = 0.65, na.rm = TRUE) +
    ggplot2::geom_point(size = 2.5, na.rm = TRUE) +
    ggplot2::facet_wrap(~ Scale_Label, ncol = 1, scales = "free_y") +
    ggplot2::labs(x = NULL, y = "Fixed operational RF threshold",
                  title = paste0(title_scope, ": RF thresholds preceding disease triggers")) +
    theme_pub() +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, vjust = 1),
      plot.margin = ggplot2::margin(8, 12, 14, 12))
}
p_bounds_ncr <- make_bounds_plot("NCR", "NCR")
p_bounds_reg <- make_bounds_plot("Regional (pooled 17 regions)", "Regional")
save_pub("FigureC_NCR_RF_Operational_Threshold_Bounds", p_bounds_ncr,
         NC_W_DOUBLE, 7.0, OUT_DIR)
save_pub("FigureC_Regional_RF_Operational_Threshold_Bounds", p_bounds_reg,
         NC_W_DOUBLE, 7.0, OUT_DIR)

stagec_h2h_make_figure <- function(anchor, scope_name, title_scope) {
  ov <- stagec_h2h_overall |> dplyr::filter(Anchor == anchor, Scope == scope_name)
  if (nrow(ov) != 1L) {
    stop("Stage C figure assembly expected one head-to-head row for ", title_scope,
         " / ", anchor, "; found ", nrow(ov), ".", call. = FALSE)
  }
  rep_levels <- c("RF accumulation\n(mm)", "RF STA/LTA\nR(t)")
  auc_long <- dplyr::bind_rows(
    ov |> dplyr::transmute(Representation = rep_levels[1], AUC = AUC_RF_mm,
                            Lower = AUC_RF_mm_CI_Lower, Upper = AUC_RF_mm_CI_Upper),
    ov |> dplyr::transmute(Representation = rep_levels[2], AUC = AUC_RF_Rt,
                            Lower = AUC_RF_Rt_CI_Lower, Upper = AUC_RF_Rt_CI_Upper)
  )
  auc_long$Representation <- factor(auc_long$Representation, levels = rep_levels)
  p_a <- ggplot2::ggplot(auc_long, ggplot2::aes(x = Representation, y = AUC)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = Lower, ymax = Upper),
                           width = 0.12, na.rm = TRUE) +
    ggplot2::geom_point(size = 2.5, na.rm = TRUE) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(x = NULL, y = "ROC AUC", title = "a    Discrimination") +
    theme_pub() +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5))

  op_long <- dplyr::bind_rows(
    ov |> dplyr::transmute(Representation = rep_levels[1],
                            Sensitivity = Sensitivity_RF_mm, PPV = PPV_RF_mm,
                            Specificity = Specificity_RF_mm),
    ov |> dplyr::transmute(Representation = rep_levels[2],
                            Sensitivity = Sensitivity_RF_Rt, PPV = PPV_RF_Rt,
                            Specificity = Specificity_RF_Rt)
  ) |>
    tidyr::pivot_longer(c(Sensitivity, PPV, Specificity),
                        names_to = "Metric", values_to = "Value") |>
    dplyr::mutate(
      Metric = factor(Metric, levels = c("Sensitivity", "PPV", "Specificity")),
      Representation = factor(Representation, levels = rep_levels)
    )
  p_b <- ggplot2::ggplot(
    op_long, ggplot2::aes(x = Metric, y = Value, group = Representation,
                          shape = Representation, linetype = Representation)) +
    ggplot2::geom_line(linewidth = 0.55, na.rm = TRUE) +
    ggplot2::geom_point(size = 2.1, na.rm = TRUE) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(x = NULL, y = "Performance", title = "b    Point-threshold performance") +
    ggplot2::guides(
      shape = ggplot2::guide_legend(nrow = 1, byrow = TRUE, direction = "horizontal"),
      linetype = "none") +
    theme_pub() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
      legend.position = "bottom", legend.direction = "horizontal",
      legend.box = "horizontal", legend.justification = "center")

  fa_long <- dplyr::bind_rows(
    ov |> dplyr::transmute(Representation = rep_levels[1],
                            False_Alarm_Rate = False_Alarm_Rate_RF_mm),
    ov |> dplyr::transmute(Representation = rep_levels[2],
                            False_Alarm_Rate = False_Alarm_Rate_RF_Rt)
  )
  fa_long$Representation <- factor(fa_long$Representation, levels = rep_levels)
  p_c <- ggplot2::ggplot(fa_long,
                          ggplot2::aes(x = Representation, y = False_Alarm_Rate)) +
    ggplot2::geom_col(width = 0.60) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(x = NULL, y = "False-alarm rate", title = "c    False-alarm burden") +
    theme_pub() +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5))

  patchwork::wrap_plots(list(p_a, p_b, p_c), ncol = 1, nrow = 3,
                        heights = c(1.0, 1.15, 0.9)) +
    patchwork::plot_annotation(
      title = paste0(title_scope, ": RF representation preceding ", anchor))
}

fig_stagec_ncr_cta <- stagec_h2h_make_figure("Constant TA", "NCR", "NCR")
fig_stagec_ncr_ot <- stagec_h2h_make_figure("Outbreak Threshold", "NCR", "NCR")
fig_stagec_reg_cta <- stagec_h2h_make_figure(
  "Constant TA", "Pooled 17 regions", "Regional")
fig_stagec_reg_ot <- stagec_h2h_make_figure(
  "Outbreak Threshold", "Pooled 17 regions", "Regional")
save_pub("FigureC_NCR_RF_HeadToHead_ConstantTA", fig_stagec_ncr_cta,
         NC_W_DOUBLE, 9.2, OUT_DIR)
save_pub("FigureC_NCR_RF_HeadToHead_OutbreakThreshold", fig_stagec_ncr_ot,
         NC_W_DOUBLE, 9.2, OUT_DIR)
save_pub("FigureC_Regional_RF_HeadToHead_ConstantTA", fig_stagec_reg_cta,
         NC_W_DOUBLE, 9.2, OUT_DIR)
save_pub("FigureC_Regional_RF_HeadToHead_OutbreakThreshold", fig_stagec_reg_ot,
         NC_W_DOUBLE, 9.2, OUT_DIR)

# -----------------------------------------------------------------------------
# 10. METHODS NOTE
# -----------------------------------------------------------------------------
writeLines(c(
  "All four RF products are operational upper-tail activation thresholds: RF >= threshold.",
  "Low rainfall is never interpreted as a stronger RF trigger.",
  "Products 1 and 2 use cumulative RF_HDX during t-4..t-1; Products 3 and 4 use the maximum rainfall STA/LTA R(t) during t-4..t-1.",
  "The disease-trigger week is excluded from every environmental predictor, so the RF feature is available by t-1 and precedes the Constant TA or Outbreak Threshold trigger week at t.",
  "For each region-year and disease anchor, the primary true derivation anchor is the first detector-active trigger week that falls in A1 OR A2. Detector-active weeks outside both A1 and A2 are false-alarm comparators; later active weeks inside A1/A2 are retained for audit but excluded from threshold fitting to avoid repeated positive observations from one detector episode.",
  "The point estimate is chosen in the >= direction from this A1/A2 derivation inventory. Thresholds within 0.01 of the best balanced accuracy are ranked by lower false-alarm burden, lower false-alarm rate, higher PPV, higher sensitivity and specificity, then by a more conservative higher threshold. The RF feature itself is calculated only from t-4..t-1, excluding the disease-trigger week.",
  "For rainfall STA/LTA products, candidate windows 2/13, 3/20, and 4/26 with a 2-week guard are compared using balanced accuracy, false-alarm burden, PPV, sensitivity, and ROC AUC. Rainfall windows are independent of the leptospirosis case STA/LTA windows, and the RF window is re-selected inside every LOYO/LORO stability fold.",
  "NCR thresholds are derived only from NCR history and are fixed across NCR evaluation years. NCR stability bounds use leave-one-year-out re-derivation.",
  "Regional thresholds are derived from the pooled 17-region history and are fixed across all 17 regional reruns. Regional stability bounds use leave-one-region-out re-derivation.",
  "Lower and upper bounds are empirical stability bounds, not separate annually re-estimated thresholds.",
  "Operational interpretation is fixed: RF below the lower bound = no activation; lower to below point = precautionary activation; point to below upper = formal RF trigger; at or above upper = strong RF trigger.",
  "Stage C head-to-head comparisons pair RF mm versus RF STA/LTA on the same primary-true/false-comparator derivation trigger inventory and explicitly report false-alarm rate in addition to discrimination and point-threshold performance.",
  "Stage H remains the downstream operational head-to-head after each scope-specific RF product is applied to leptospirosis detection."
), file.path(OUT_DIR, "StageC_Methods_Note.txt"))

cat("Stage C complete: upper-tail, A1/A2-consistent RF thresholds derived for NCR and Regional scopes.\n")
cat("Outputs: ", OUT_DIR, "\n", sep = "")
