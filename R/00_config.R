# =============================================================================
# PROJECT CONFIGURATION
# -----------------------------------------------------------------------------
# Portable paths and shared constants. Sourced first by every stage script.
#
# Replaces the hard-coded "C:/Users/User/Desktop/..." paths that appeared in all
# five original scripts. The project root is discovered automatically, so the
# project runs unmodified on Windows, macOS and Linux from any location.
# =============================================================================


# -----------------------------------------------------------------------------
# 1. PROJECT ROOT DISCOVERY
# -----------------------------------------------------------------------------
# Resolution order:
#   1. TA_PROJECT_ROOT environment variable (explicit override)
#   2. here::here() if the 'here' package is installed
#   3. Walk up from the working directory looking for the project marker
.find_project_root <- function() {
  env <- Sys.getenv("TA_PROJECT_ROOT", unset = "")
  if (nzchar(env) && dir.exists(env)) return(normalizePath(env, winslash = "/"))

  if (requireNamespace("here", quietly = TRUE)) {
    r <- try(here::here(), silent = TRUE)
    if (!inherits(r, "try-error") && dir.exists(file.path(r, "R"))) {
      return(normalizePath(r, winslash = "/"))
    }
  }

  # Walk upward for a directory containing both R/ and scripts/
  d <- normalizePath(getwd(), winslash = "/")
  for (i in seq_len(6)) {
    if (dir.exists(file.path(d, "R")) && dir.exists(file.path(d, "scripts"))) {
      return(d)
    }
    parent <- dirname(d)
    if (identical(parent, d)) break
    d <- parent
  }

  stop("Could not locate the project root.\n",
       "Set the working directory to the project folder (the one containing ",
       "R/ and scripts/), or set TA_PROJECT_ROOT.", call. = FALSE)
}

PROJECT_ROOT <- .find_project_root()


# -----------------------------------------------------------------------------
# 2. STANDARD DIRECTORIES
# -----------------------------------------------------------------------------
DIR_R       <- file.path(PROJECT_ROOT, "R")
DIR_SCRIPTS <- file.path(PROJECT_ROOT, "scripts")
DIR_DATA    <- file.path(PROJECT_ROOT, "data")
DIR_OUTPUT  <- file.path(PROJECT_ROOT, "outputs")

for (d in c(DIR_DATA, DIR_OUTPUT)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}


# -----------------------------------------------------------------------------
# 3. INPUT DATA
# -----------------------------------------------------------------------------
DATA_FILE <- file.path(DIR_DATA, "Dengue_Environmental_Vars_Dataset.xlsx")

# This adaptation evaluates leptospirosis, which in the source workbook is
# only available on the Regional Data sheet (LC_DOH), so only that sheet is
# read. QC Data and Multi-Setting Data are dengue-only in this workbook and
# are not touched by any script in this pipeline.
SHEET_REGIONAL <- "Regional Data"  # REGION, YR, WN, DC_DOH, LC_DOH, RF_HDX, flags

if (!file.exists(DATA_FILE)) {
  stop("Input dataset not found:\n  ", DATA_FILE,
       "\nPlace 'Dengue_Environmental_Vars_Dataset.xlsx' in the data/ directory.",
       call. = FALSE)
}


# -----------------------------------------------------------------------------
# 4. OPTIONAL GEOMETRY FOR THE STAGEB MAP
# -----------------------------------------------------------------------------
# StageB draws a Philippine regional choropleth. Geometry is resolved in this
# order (see scripts/StageB_Regional_Lepto_Analysis.R):
#   1. A shapefile / GeoPackage placed in data/geometry/  (fully offline)
#   2. geodata::gadm()          (downloads, cached under data/geometry/)
#   3. rnaturalearth::ne_states() (downloads)
#
# Supplying local geometry is strongly recommended for a reproducible archive:
# it removes the only remaining network dependency and pins the boundary
# vintage, which matters because GADM revises Philippine regions periodically.
DIR_GEOMETRY <- file.path(DIR_DATA, "geometry")
if (!dir.exists(DIR_GEOMETRY)) {
  dir.create(DIR_GEOMETRY, recursive = TRUE, showWarnings = FALSE)
}

# Set to an explicit file path to force a specific boundary file.
PH_SHAPEFILE <- {
  cands <- list.files(DIR_GEOMETRY, pattern = "\\.(shp|gpkg|geojson)$",
                      full.names = TRUE, recursive = TRUE)
  if (length(cands) > 0) cands[1] else NULL
}


# -----------------------------------------------------------------------------
# 5. ANALYSIS CONSTANTS SHARED ACROSS STAGES
# -----------------------------------------------------------------------------
# Years excluded from threshold derivation and evaluation.
#   2020, 2021 - COVID-19 surveillance disruption; the source Regional Data
#                sheet carries no rows for these two years at all, for either
#                disease, so they are excluded by absence as well as by name.
# Unlike the dengue pipeline this adaptation does NOT drop 2025. The dengue
# exclusion of 2025 reflected a truncated, still-updating QC extract; the
# Regional Data sheet has no such truncation, LC_DOH carries a full 52-week
# panel for 2025 in every region, so 2025 is kept as an evaluable year here.
# Each stage script still carries its own local copy of this constant (the
# scripts are sourced into separate environments), updated to match.
EXCLUDED_YEARS <- c(2020L, 2021L)

# Canonical Philippine regions used by pooled and independent-region branches.
CANONICAL_17 <- c(
  "BARMM", "CAR", "MIMAROPA", "NCR",
  paste("REGION", c("I", "II", "III", "IV-A", "V", "VI", "VII", "VIII",
                    "IX", "X", "XI", "XII", "XIII"))
)

# Reproducibility
GLOBAL_SEED <- 12345L
set.seed(GLOBAL_SEED)
options(scipen = 999, stringsAsFactors = FALSE)


# -----------------------------------------------------------------------------
# 6. ETA THRESHOLDS (Stage 1 -> StageA / StageB / StageC handoff)
# -----------------------------------------------------------------------------
# Stage 1 derives eta_ON / eta_OFF empirically and writes them to
# outputs/Stage1_eta_thresholds_NCR/eta_thresholds_derived.R.
#
# In the original code the downstream stages hard-coded ETA_ON <- 1.33 /
# ETA_OFF <- 0.73 and never read Stage 1's output, so re-running Stage 1
# could not change them -- a silent reproducibility break. load_eta_
# thresholds() now prefers the derived values and falls back to the
# published constants with a warning.
ETA_ON_PUBLISHED  <- 1.33
ETA_OFF_PUBLISHED <- 0.73

load_eta_thresholds <- function(verbose = TRUE) {
  derived <- file.path(DIR_OUTPUT, "Stage1_eta_thresholds_NCR",
                       "eta_thresholds_derived.R")

  r <- read_value_file(derived,
                       required = c("ETA_ON_ADOPTED", "ETA_OFF_ADOPTED"),
                       label = "eta_thresholds_derived.R")

  if (r$ok) {
    on  <- r$values[["ETA_ON_ADOPTED"]]
    off <- r$values[["ETA_OFF_ADOPTED"]]
    if (is.numeric(on) && is.numeric(off) && is.finite(on) && is.finite(off)) {
      if (verbose) {
        message(sprintf(
          "[eta] Using Stage 1 derived thresholds: eta_ON = %.2f, eta_OFF = %.2f",
          on, off))
      }
      return(list(ETA_ON = on, ETA_OFF = off, source = "stage1"))
    }
  }

  # The fallback previously fired SILENTLY because of the emptyenv bug above,
  # so Stage 1's derived thresholds never actually reached Stage 3. The reason
  # is now surfaced in the warning so a masked failure cannot recur unnoticed.
  if (verbose) {
    warning("Falling back to published eta thresholds ",
            sprintf("(eta_ON = %.2f, eta_OFF = %.2f).\n  Reason: ",
                    ETA_ON_PUBLISHED, ETA_OFF_PUBLISHED),
            if (nzchar(r$message)) r$message else "unknown",
            "\n  Run Stage 1 for a fully reproducible chain.", call. = FALSE)
  }
  list(ETA_ON = ETA_ON_PUBLISHED, ETA_OFF = ETA_OFF_PUBLISHED,
       source = "published")
}

# -----------------------------------------------------------------------------
# 6a. READING AUTO-GENERATED VALUE FILES
# -----------------------------------------------------------------------------
# Shared reader for the small sourceable file Stage 1 emits.
#
# ROOT CAUSE OF A PAST FAILURE IN THE ORIGINAL PIPELINE (fixed here, carried
# over unchanged):
#   These files were previously sourced into new.env(parent = emptyenv()).
#   In R an assignment such as `X <- 0.74` is a call to the function `<-`,
#   resolved lexically through the environment's parent chain. With emptyenv()
#   as the parent that chain is empty, so `<-` itself cannot be found and
#   sourcing failed with "could not find function \"<-\"". The generated files
#   were always correct; the environment they were read into was not.
#   Parent is now baseenv(), which exposes base functions while still isolating
#   the file from the caller's workspace.
read_value_file <- function(path, required, label = basename(path)) {
  if (!file.exists(path)) {
    return(list(ok = FALSE, reason = "missing", values = NULL,
                message = paste0("File not found: ", path)))
  }
  e <- new.env(parent = baseenv())
  res <- tryCatch({ sys.source(path, envir = e); NULL },
                  error = function(err) conditionMessage(err))
  if (!is.null(res)) {
    return(list(ok = FALSE, reason = "source_error", values = NULL,
                message = paste0("Could not source ", label, ": ", res)))
  }
  missing <- required[!vapply(required, exists, logical(1),
                              envir = e, inherits = FALSE)]
  if (length(missing) > 0) {
    return(list(ok = FALSE, reason = "missing_object", values = NULL,
                message = paste0(label, " does not define: ",
                                 paste(missing, collapse = ", "),
                                 "\n  It defines: ",
                                 paste(ls(e), collapse = ", "))))
  }
  vals <- mget(required, envir = e)
  list(ok = TRUE, reason = "ok", values = vals, message = "")
}

#' Validate a value is a single finite proportion strictly inside (0, 1).
is_valid_proportion <- function(x) {
  is.numeric(x) && length(x) == 1L && is.finite(x) && x > 0 && x < 1
}

# -----------------------------------------------------------------------------
# 6a-i. STA/LTA WINDOWS (Stage 1 -> StageA / StageB / StageC handoff)
# -----------------------------------------------------------------------------
# NEW for this leptospirosis adaptation (not present in the original dengue
# pipeline, which fixed these windows a priori). Stage 1 now runs an
# empirical window-optimisation step and writes the winning (STA, LTA)
# pair for each of the two transmission-acceleration detectors to
# outputs/Stage1_eta_thresholds_NCR/sta_lta_windows_derived.R. If that file
# is absent, this falls back to the published windows (Constant TA 4/26,
# Continuous TA 3/12) with a warning, the same pattern load_eta_thresholds()
# uses.
STA_WIN_CONSTANT_PUBLISHED   <- 4L
LTA_WIN_CONSTANT_PUBLISHED   <- 26L
STA_WIN_CONTINUOUS_PUBLISHED <- 3L
LTA_WIN_CONTINUOUS_PUBLISHED <- 12L

load_sta_lta_windows <- function(verbose = TRUE) {
  derived <- file.path(DIR_OUTPUT, "Stage1_eta_thresholds_NCR",
                       "sta_lta_windows_derived.R")

  r <- read_value_file(derived,
                       required = c("STA_WIN_CONSTANT_ADOPTED", "LTA_WIN_CONSTANT_ADOPTED",
                                    "STA_WIN_CONTINUOUS_ADOPTED", "LTA_WIN_CONTINUOUS_ADOPTED"),
                       label = "sta_lta_windows_derived.R")

  if (r$ok) {
    sc <- r$values[["STA_WIN_CONSTANT_ADOPTED"]]
    lc <- r$values[["LTA_WIN_CONSTANT_ADOPTED"]]
    sn <- r$values[["STA_WIN_CONTINUOUS_ADOPTED"]]
    ln <- r$values[["LTA_WIN_CONTINUOUS_ADOPTED"]]
    ok_numeric <- all(vapply(list(sc, lc, sn, ln),
                             function(x) is.numeric(x) && length(x) == 1L && is.finite(x),
                             logical(1)))
    if (ok_numeric) {
      if (verbose) {
        message(sprintf(
          "[sta_lta] Using Stage 1 derived windows: Constant TA = %d/%d, Continuous TA = %d/%d",
          as.integer(sc), as.integer(lc), as.integer(sn), as.integer(ln)))
      }
      return(list(STA_WIN_CONSTANT = as.integer(sc), LTA_WIN_CONSTANT = as.integer(lc),
                  STA_WIN_CONTINUOUS = as.integer(sn), LTA_WIN_CONTINUOUS = as.integer(ln),
                  source = "stage1"))
    }
  }

  if (verbose) {
    warning("Falling back to published STA/LTA windows ",
            sprintf("(Constant TA = %d/%d, Continuous TA = %d/%d).\n  Reason: ",
                    STA_WIN_CONSTANT_PUBLISHED, LTA_WIN_CONSTANT_PUBLISHED,
                    STA_WIN_CONTINUOUS_PUBLISHED, LTA_WIN_CONTINUOUS_PUBLISHED),
            if (nzchar(r$message)) r$message else "unknown",
            "\n  Run Stage 1 for a fully reproducible chain.", call. = FALSE)
  }
  list(STA_WIN_CONSTANT = STA_WIN_CONSTANT_PUBLISHED, LTA_WIN_CONSTANT = LTA_WIN_CONSTANT_PUBLISHED,
      STA_WIN_CONTINUOUS = STA_WIN_CONTINUOUS_PUBLISHED, LTA_WIN_CONTINUOUS = LTA_WIN_CONTINUOUS_PUBLISHED,
      source = "published")
}

# -----------------------------------------------------------------------------
# 6b. EPIDEMIC BURDEN FRACTION  (fixed at 70%)
# -----------------------------------------------------------------------------
# The A2 anchor ("Epidemic Burden" block), and hence the Activation Threshold,
# uses a FIXED burden fraction of 0.70, matching the published manuscript and
# retained unchanged for this leptospirosis adaptation.
#
# This value is defined once here and read by StageA, StageB and StageC. It
# is a plain constant: no file is read, nothing is sourced, and no stage
# depends on the output of any other stage for it.
A2_BURDEN_FRAC <- 0.70
A2_BURDEN_PCT  <- sprintf("%.0f%%", 100 * A2_BURDEN_FRAC)   # "70%" for captions

# -----------------------------------------------------------------------------
# 7. DEPENDENCY CHECK
# -----------------------------------------------------------------------------
# The original scripts called install.packages() at run time. That mutates the
# user's library without consent, fails on offline / HPC nodes, and breaks
# non-interactive execution. Dependencies are now verified, never installed;
# R/install_dependencies.R performs installation as a deliberate one-off step.
require_packages <- function(pkgs, purpose = NULL) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Missing required package(s)",
         if (!is.null(purpose)) paste0(" for ", purpose) else "", ":\n  ",
         paste(missing, collapse = ", "),
         "\n\nInstall them with:\n  source(\"R/install_dependencies.R\")",
         call. = FALSE)
  }
  invisible(TRUE)
}

message("[config] Project root: ", PROJECT_ROOT)

# -----------------------------------------------------------------------------
# 8. CALENDAR-CONTIGUOUS ROLLING HELPERS AND RAINFALL-GATE HANDOFF
# -----------------------------------------------------------------------------
# The Regional Data sheet has no rows for 2020-2021. Row-adjacent rolling
# windows can therefore bridge directly from 2019 to 2022 unless calendar
# continuity is enforced explicitly. All revised TA and RF computations use
# these helpers.
is_contiguous_week_span <- function(dates, start_idx, end_idx) {
  if (length(dates) == 0L || start_idx < 1L || end_idx > length(dates) ||
      start_idx > end_idx) return(FALSE)
  dd <- dates[start_idx:end_idx]
  if (any(is.na(dd))) return(FALSE)
  if (length(dd) <= 1L) return(TRUE)
  all(as.numeric(diff(dd)) == 7)
}

rolling_contiguous_mean <- function(values, dates, k) {
  n <- length(values)
  out <- rep(NA_real_, n)
  if (n == 0L || k < 1L) return(out)
  for (i in seq_len(n)) {
    s <- i - k + 1L
    if (s < 1L || !is_contiguous_week_span(dates, s, i)) next
    vv <- values[s:i]
    if (all(is.na(vv))) next
    out[i] <- mean(vv, na.rm = TRUE)
  }
  out
}

# Sum rainfall in the k COMPLETE weeks immediately BEFORE each index: t-k..t-1.
# The trigger week itself is deliberately excluded. A value is returned only
# when those k weeks and the trigger week form one uninterrupted 7-day sequence.
preceding_contiguous_sum <- function(values, dates, k = 4L) {
  n <- length(values)
  out <- rep(NA_real_, n)
  if (n == 0L || k < 1L) return(out)
  for (i in seq_len(n)) {
    s <- i - k
    e <- i - 1L
    if (s < 1L || e < s) next
    if (!is_contiguous_week_span(dates, s, i)) next
    vv <- values[s:e]
    if (any(!is.finite(vv))) next
    out[i] <- sum(vv)
  }
  out
}

# Maximum of a weekly feature in the k COMPLETE weeks before each index.
# This is used by the RF STA/LTA products so the environmental ratio always
# precedes the disease trigger and never uses trigger-week rainfall.
preceding_contiguous_max <- function(values, dates, k = 4L) {
  n <- length(values)
  out <- rep(NA_real_, n)
  if (n == 0L || k < 1L) return(out)
  for (i in seq_len(n)) {
    s <- i - k
    e <- i - 1L
    if (s < 1L || e < s) next
    if (!is_contiguous_week_span(dates, s, i)) next
    vv <- values[s:e]
    if (any(!is.finite(vv))) next
    out[i] <- max(vv)
  }
  out
}

# Guarded STA/LTA ratio for rainfall. Unlike the disease Constant-TA state
# machine, this function returns the raw weekly environmental ratio only. The
# short-term mean ends at week t; the long-term mean ends before the guard band.
# Product 3/4 then use the maximum of this ratio during t-4..t-1, so no rainfall
# observed in the disease-trigger week contributes to the environmental gate.
guarded_sta_lta_ratio <- function(values, dates, sta_win, lta_win, guard = 2L) {
  sta_win <- as.integer(sta_win); lta_win <- as.integer(lta_win); guard <- as.integer(guard)
  if (sta_win < 1L || lta_win < 1L || guard < 0L)
    stop("Invalid STA/LTA window specification.", call. = FALSE)
  n <- length(values)
  out <- rep(NA_real_, n)
  min_t <- sta_win + guard + lta_win
  if (n == 0L) return(out)
  for (t in seq_len(n)) {
    if (t < min_t) next
    full_start <- t - sta_win - guard - lta_win + 1L
    lta_end <- t - sta_win - guard
    sta_start <- t - sta_win + 1L
    if (full_start < 1L || lta_end < full_start ||
        !is_contiguous_week_span(dates, full_start, t)) next
    sta_vals <- values[sta_start:t]
    lta_vals <- values[full_start:lta_end]
    if (any(!is.finite(sta_vals)) || any(!is.finite(lta_vals))) next
    sta <- mean(sta_vals)
    lta <- mean(lta_vals)
    if (is.finite(sta) && is.finite(lta) && lta > 0) out[t] <- sta / lta
  }
  out
}

RF_PRECEDING_WEEKS <- 4L

# Normalize the operational RF comparison operator. The current framework is
# deliberately upper-tail only: increasing rainfall represents increasing
# environmental activation, so Products 1-4 must all use RF >= threshold.
normalize_rf_operator <- function(operator) {
  op <- trimws(unname(as.character(operator))[1L])
  if (length(op) != 1L || is.na(op) || !nzchar(op) || op != ">=") {
    stop("Unsupported operational RF comparison operator: ",
         paste(as.character(operator), collapse = ", "),
         ". Products 1-4 require >=.", call. = FALSE)
  }
  op
}

apply_rf_gate <- function(feature, threshold, operator) {
  normalize_rf_operator(operator)
  thr <- suppressWarnings(as.numeric(unname(threshold)[1L]))
  if (length(thr) != 1L || !is.finite(thr)) {
    stop("RF threshold must be one finite numeric value.", call. = FALSE)
  }
  is.finite(feature) & feature >= thr
}

# Convert weekly RF threshold exceedances into an antecedent qualification for
# a disease-trigger week. A threshold signal is valid only when it occurs 1 to
# max_lead complete weeks before the disease week. The earliest qualifying RF
# signal is retained so figures can display the actual first RF threshold
# crossing rather than forcing a marker at t-1. `signal_ok[i]` corresponds to
# the RF feature available at dates[i] - 7 days.
rf_precedence_qualification <- function(signal_ok, dates, max_lead = RF_PRECEDING_WEEKS) {
  signal_ok <- as.logical(signal_ok)
  dates <- as.Date(dates)
  n <- length(signal_ok)
  if (length(dates) != n) stop("RF precedence inputs must have equal length.", call. = FALSE)
  qualified <- rep(FALSE, n)
  source_index <- rep(NA_integer_, n)
  lead_weeks <- rep(NA_integer_, n)
  for (t in seq_len(n)) {
    if (is.na(dates[t])) next
    lo <- max(1L, t - as.integer(max_lead) + 1L)
    cand <- lo:t
    signal_dates <- dates[cand] - 7
    lead <- as.integer(round(as.numeric(dates[t] - signal_dates) / 7))
    valid <- !is.na(signal_ok[cand]) & signal_ok[cand] &
      !is.na(signal_dates) & lead >= 1L & lead <= as.integer(max_lead)
    if (any(valid)) {
      valid_idx <- cand[valid]
      valid_lead <- lead[valid]
      # First/earliest qualifying RF signal in the allowed pre-anchor window.
      k <- which.max(valid_lead)
      qualified[t] <- TRUE
      source_index[t] <- valid_idx[k]
      lead_weeks[t] <- valid_lead[k]
    }
  }
  list(qualified = qualified, source_index = source_index, lead_weeks = lead_weeks)
}


# Classify RF observations against one fixed lower/point/upper threshold triplet.
# The current operational definition is strictly upper-tail: below lower = no
# activation; lower-to-point = precautionary; point-to-upper = formal; >= upper
# = strong. The point estimate is the primary headline RF trigger.
classify_rf_activation <- function(feature, point, lower, upper, operator) {
  op <- normalize_rf_operator(operator)
  if (op != ">=") {
    stop("Operational RF activation requires the upper-tail >= rule.", call. = FALSE)
  }
  x <- suppressWarnings(as.numeric(feature))
  pt <- suppressWarnings(as.numeric(unname(point)[1L]))
  lo <- suppressWarnings(as.numeric(unname(lower)[1L]))
  hi <- suppressWarnings(as.numeric(unname(upper)[1L]))
  if (any(!is.finite(c(pt, lo, hi))) || lo > pt || hi < pt) {
    stop("RF activation bounds must be finite and ordered lower <= point <= upper.", call. = FALSE)
  }
  out <- rep(NA_character_, length(x))
  ok <- is.finite(x)
  out[ok] <- "No RF activation"
  out[ok & x >= lo] <- "Precautionary RF trigger"
  out[ok & x >= pt] <- "Formal RF trigger"
  out[ok & x >= hi] <- "Strong RF trigger"
  out
}

# Four environmental products are derived in Stage C:
#   P1 RF mm preceding first Constant-TA trigger
#   P2 RF mm preceding first Outbreak-Threshold trigger
#   P3 empirically derived RF STA/LTA R(t) preceding first Constant-TA trigger
#   P4 empirically derived RF STA/LTA R(t) preceding first Outbreak-Threshold trigger
# P3/P4 have rainfall-specific STA/LTA windows and activation cutoffs. They do
# not borrow either the case windows or eta_ON from the leptospirosis detector.
load_rf_thresholds <- function(scope = c("Regional", "NCR", "IndependentRegion"), region = NULL, verbose = TRUE) {
  scope <- match.arg(scope)
  if (scope == "IndependentRegion") {
    reg <- trimws(as.character(region)[1L])
    if (!nzchar(reg)) stop("IndependentRegion RF thresholds require an explicit region name.", call. = FALSE)
    csv_path <- file.path(DIR_OUTPUT, "StageC_RF_Threshold_Derivation",
                          "RF_Threshold_Derivation_Summary_PerRegion.csv")
    if (!file.exists(csv_path)) {
      stop("Independent-region RF thresholds are unavailable: ", csv_path,
           "\nRun Stage C first.", call. = FALSE)
    }
    z <- utils::read.csv(csv_path, stringsAsFactors = FALSE, check.names = FALSE)
    need <- c("REGION", "Product_ID", "RF_Scale", "Adopted_Threshold",
              "Threshold_Lower", "Threshold_Upper", "Comparison_Operator",
              "RF_STA", "RF_LTA", "RF_Guard")
    miss <- setdiff(need, names(z))
    if (length(miss)) stop("Independent-region RF threshold table is missing: ",
                           paste(miss, collapse=", "), call. = FALSE)
    z <- z[z$REGION == reg, , drop=FALSE]
    if (nrow(z) != 4L || !identical(sort(as.integer(z$Product_ID)), 1:4)) {
      stop("Independent-region RF threshold table must contain Products 1-4 for ", reg,
           ". Found ", nrow(z), " rows.", call. = FALSE)
    }
    z <- z[match(1:4, as.integer(z$Product_ID)), , drop=FALSE]
    mk <- function(i) {
      sc <- tolower(as.character(z$RF_Scale[i]))
      list(scale = if (sc == "mm") "mm" else "rt",
           threshold = as.numeric(z$Adopted_Threshold[i]),
           lower = as.numeric(z$Threshold_Lower[i]),
           upper = as.numeric(z$Threshold_Upper[i]),
           operator = normalize_rf_operator(z$Comparison_Operator[i]),
           sta = if (sc == "mm") NA_integer_ else as.integer(z$RF_STA[i]),
           lta = if (sc == "mm") NA_integer_ else as.integer(z$RF_LTA[i]),
           guard = if (sc == "mm") NA_integer_ else as.integer(z$RF_Guard[i]))
    }
    out <- list(product1=mk(1), product2=mk(2), product3=mk(3), product4=mk(4),
                source=csv_path, scope=scope, region=reg)
    for (i in 1:4) {
      cfg <- out[[paste0("product", i)]]
      if (any(!is.finite(c(cfg$lower, cfg$threshold, cfg$upper))) ||
          cfg$lower > cfg$threshold || cfg$upper < cfg$threshold) {
        stop("Invalid independent-region RF threshold triplet for ", reg,
             " Product ", i, ".", call. = FALSE)
      }
      out[[paste0("product", i)]]$activation_bound <- cfg$lower
      out[[paste0("product", i)]]$strong_bound <- cfg$upper
    }
    if (verbose) message("[rf:IndependentRegion] Using region-specific P1-P4 thresholds for ", reg)
    return(out)
  }
  path <- file.path(DIR_OUTPUT, "StageC_RF_Threshold_Derivation",
                    "rf_thresholds_derived.R")
  base_required <- c(
    "RF_THRESHOLD_PRODUCT1_MM", "RF_THRESHOLD_PRODUCT2_MM",
    "RF_THRESHOLD_PRODUCT3_RT", "RF_THRESHOLD_PRODUCT4_RT",
    "RF_THRESHOLD_LOWER_PRODUCT1_MM", "RF_THRESHOLD_UPPER_PRODUCT1_MM",
    "RF_THRESHOLD_LOWER_PRODUCT2_MM", "RF_THRESHOLD_UPPER_PRODUCT2_MM",
    "RF_THRESHOLD_LOWER_PRODUCT3_RT", "RF_THRESHOLD_UPPER_PRODUCT3_RT",
    "RF_THRESHOLD_LOWER_PRODUCT4_RT", "RF_THRESHOLD_UPPER_PRODUCT4_RT",
    "RF_STA_PRODUCT3", "RF_LTA_PRODUCT3", "RF_GUARD_PRODUCT3",
    "RF_STA_PRODUCT4", "RF_LTA_PRODUCT4", "RF_GUARD_PRODUCT4",
    "RF_OPERATOR_PRODUCT1", "RF_OPERATOR_PRODUCT2",
    "RF_OPERATOR_PRODUCT3", "RF_OPERATOR_PRODUCT4"
  )
  ncr_required <- paste0(base_required, "_NCR")
  r <- read_value_file(path, required = unique(c(base_required, ncr_required)),
                       label = "rf_thresholds_derived.R")
  if (!r$ok) {
    stop("RF-gated analysis requested but Stage C threshold handoff is unavailable.\n",
         r$message, "\nRun scripts/StageC_RF_Threshold_Derivation.R first.",
         call. = FALSE)
  }
  suffix <- if (scope == "NCR") "_NCR" else ""
  getv <- function(base) r$values[[paste0(base, suffix)]]
  num_base <- c(
    "RF_THRESHOLD_PRODUCT1_MM", "RF_THRESHOLD_PRODUCT2_MM",
    "RF_THRESHOLD_PRODUCT3_RT", "RF_THRESHOLD_PRODUCT4_RT",
    "RF_THRESHOLD_LOWER_PRODUCT1_MM", "RF_THRESHOLD_UPPER_PRODUCT1_MM",
    "RF_THRESHOLD_LOWER_PRODUCT2_MM", "RF_THRESHOLD_UPPER_PRODUCT2_MM",
    "RF_THRESHOLD_LOWER_PRODUCT3_RT", "RF_THRESHOLD_UPPER_PRODUCT3_RT",
    "RF_THRESHOLD_LOWER_PRODUCT4_RT", "RF_THRESHOLD_UPPER_PRODUCT4_RT",
    "RF_STA_PRODUCT3", "RF_LTA_PRODUCT3", "RF_GUARD_PRODUCT3",
    "RF_STA_PRODUCT4", "RF_LTA_PRODUCT4", "RF_GUARD_PRODUCT4"
  )
  v <- lapply(num_base, function(nm) suppressWarnings(as.numeric(getv(nm))))
  names(v) <- num_base
  scalar_finite <- vapply(v, function(x) length(x) == 1L && is.finite(x), logical(1))
  if (!all(scalar_finite)) {
    stop("Stage C RF threshold handoff contains missing/non-finite numeric values for ",
         scope, ": ", paste(names(v)[!scalar_finite], collapse = ", "), call. = FALSE)
  }
  ops <- unname(vapply(1:4, function(i) {
    normalize_rf_operator(getv(paste0("RF_OPERATOR_PRODUCT", i)))
  }, character(1)))
  if (any(ops != ">=")) {
    stop("Operational RF thresholds must use the upper-tail >= rule for ", scope,
         ". Re-run Stage C with the current pipeline.", call. = FALSE)
  }
  out <- list(
    product1 = list(scale = "mm", threshold = v$RF_THRESHOLD_PRODUCT1_MM,
                    lower = v$RF_THRESHOLD_LOWER_PRODUCT1_MM,
                    upper = v$RF_THRESHOLD_UPPER_PRODUCT1_MM,
                    operator = ">=", sta = NA_integer_, lta = NA_integer_, guard = NA_integer_),
    product2 = list(scale = "mm", threshold = v$RF_THRESHOLD_PRODUCT2_MM,
                    lower = v$RF_THRESHOLD_LOWER_PRODUCT2_MM,
                    upper = v$RF_THRESHOLD_UPPER_PRODUCT2_MM,
                    operator = ">=", sta = NA_integer_, lta = NA_integer_, guard = NA_integer_),
    product3 = list(scale = "rt", threshold = v$RF_THRESHOLD_PRODUCT3_RT,
                    lower = v$RF_THRESHOLD_LOWER_PRODUCT3_RT,
                    upper = v$RF_THRESHOLD_UPPER_PRODUCT3_RT,
                    operator = ">=", sta = as.integer(v$RF_STA_PRODUCT3),
                    lta = as.integer(v$RF_LTA_PRODUCT3), guard = as.integer(v$RF_GUARD_PRODUCT3)),
    product4 = list(scale = "rt", threshold = v$RF_THRESHOLD_PRODUCT4_RT,
                    lower = v$RF_THRESHOLD_LOWER_PRODUCT4_RT,
                    upper = v$RF_THRESHOLD_UPPER_PRODUCT4_RT,
                    operator = ">=", sta = as.integer(v$RF_STA_PRODUCT4),
                    lta = as.integer(v$RF_LTA_PRODUCT4), guard = as.integer(v$RF_GUARD_PRODUCT4)),
    source = path, scope = scope
  )
  for (i in 1:4) {
    cfg <- out[[paste0("product", i)]]
    if (!is.finite(cfg$lower) || !is.finite(cfg$upper) ||
        cfg$lower > cfg$threshold || cfg$upper < cfg$threshold) {
      stop("Stage C RF threshold bounds are invalid for ", scope, " Product ", i,
           ": expected lower <= point <= upper.", call. = FALSE)
    }
    out[[paste0("product", i)]]$activation_bound <- cfg$lower
    out[[paste0("product", i)]]$strong_bound <- cfg$upper
  }
  if (verbose) {
    message(sprintf(
      "[rf:%s] P1 >= %.3f mm [%.3f-%.3f]; P2 >= %.3f mm [%.3f-%.3f]; P3 >= %.3f R(t) [%.3f-%.3f; %d/%d,g=%d]; P4 >= %.3f R(t) [%.3f-%.3f; %d/%d,g=%d]",
      scope,
      out$product1$threshold, out$product1$lower, out$product1$upper,
      out$product2$threshold, out$product2$lower, out$product2$upper,
      out$product3$threshold, out$product3$lower, out$product3$upper,
      out$product3$sta, out$product3$lta, out$product3$guard,
      out$product4$threshold, out$product4$lower, out$product4$upper,
      out$product4$sta, out$product4$lta, out$product4$guard))
  }
  out
}

get_rf_gate <- function(scope = c("Regional", "NCR", "IndependentRegion"), region = NULL, verbose = TRUE) {
  scope <- match.arg(scope)
  product <- suppressWarnings(as.integer(Sys.getenv("TA_RF_GATE_PRODUCT", unset = "0")))
  if (is.na(product) || !product %in% 0:4) {
    stop("TA_RF_GATE_PRODUCT must be 0, 1, 2, 3, or 4.", call. = FALSE)
  }
  suffix <- Sys.getenv("TA_OUTPUT_SUFFIX", unset = "")
  if (product == 0L) {
    return(list(active = FALSE, product = 0L, scale = "none", threshold = NA_real_,
                lower = NA_real_, upper = NA_real_, activation_bound = NA_real_,
                strong_bound = NA_real_, operator = NA_character_,
                threshold_mm = NA_real_, threshold_rt = NA_real_, sta = NA_integer_,
                lta = NA_integer_, guard = NA_integer_, suffix = suffix,
                scope = scope, label = "Disease-only"))
  }
  th <- load_rf_thresholds(scope = scope, region = region, verbose = verbose)
  cfg <- th[[paste0("product", product)]]
  labels <- c(
    "RF Product 1: mm threshold preceding Constant TA trigger",
    "RF Product 2: mm threshold preceding Outbreak Threshold crossing",
    "RF Product 3: empirical RF STA/LTA R(t) preceding Constant TA trigger",
    "RF Product 4: empirical RF STA/LTA R(t) preceding Outbreak Threshold crossing"
  )
  list(active = TRUE, product = product, scale = cfg$scale,
       threshold = cfg$threshold, lower = cfg$lower, upper = cfg$upper,
       activation_bound = cfg$lower, strong_bound = cfg$upper,
       operator = ">=", threshold_mm = if (cfg$scale == "mm") cfg$threshold else NA_real_,
       threshold_rt = if (cfg$scale == "rt") cfg$threshold else NA_real_,
       sta = cfg$sta, lta = cfg$lta, guard = cfg$guard,
       suffix = suffix, scope = scope, label = labels[product])
}

# Regional year-cluster bootstrap replicates used by Stage B and final validation.
REGIONAL_BOOT_N_CI <- 1000L
