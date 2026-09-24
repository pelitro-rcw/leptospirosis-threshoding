# ==============================================================================
# STAGE B - REGIONAL LEPTOSPIROSIS ANALYSIS
# 17 Philippine regions; FigureB4 and FigureB5, plus tables
# ------------------------------------------------------------------------------
# Adapted from the published Stage 4 (Regional Analysis) script, which
# evaluated dengue outbreak detection. This version evaluates a
# four-method framework (Continuous TA, Constant TA, Outbreak Threshold,
# Alarm Threshold; see Item 6 of the Delphi review) against weekly
# leptospirosis case counts
# (LC_DOH) across the 17 administrative regions of the Philippines
# (2018-2025, excluding the 2020-2021 pandemic disruption). Everything
# below this header carries the original script's statistical machinery
# unchanged; only the input column, the year window, and the regional
# case-volume inclusion threshold (recalibrated for leptospirosis' smaller
# case counts; see docs/Delphi_Panel_Methodological_Review.md, Item 4) were
# adjusted. The peak used to anchor the A1 and A2 anchors below is the
# leptospirosis peak (max weekly LC_DOH), not the dengue peak.
#
# Unlike the dengue Stage 4, 2025 is NOT excluded here: the Regional Data
# sheet's 2025 rows are a complete 52-week panel for every disease, not a
# truncated in-progress extract, so there is no truncation concern to guard
# against for either year window.
#
# OUTPUTS
# -------
# Figure 4 (six multipanel figures, one PDF/PNG per file):
#   Per-metric regional dominance for the five reported operational metrics
#   plus a cross-metric Method Summary figure:
#
#     FigureB4_TAM
#     FigureB4_N_True_Alarms
#     FigureB4_Sensitivity
#     FigureB4_Mean_Lead_Time
#     FigureB4_Warning_Persistence
#     FigureB4_Method_Summary
#
#   Each per-metric figure has three sub-panels:
#     a. Regional dominance matrix  - 6 methods (rows) by 17 regions
#                                    (columns); cell shading by within-region
#                                    min-max-normalised dominance score; right
#                                    column reports the per-method count of
#                                    regions swept at score >= 0.75.
#     b. Per-detector dot plot      - one dot per region per detector;
#                                    bootstrap 95% CI bars; group-mean
#                                    crossbar.
#     c. Wilcoxon detector-paired   - two pairwise contrasts, both against
#                                    Outbreak Threshold: Constant TA vs
#                                    Outbreak Threshold, and Continuous TA
#                                    vs Outbreak Threshold; paired by region
#                                    (n = 17). Alarm Threshold is reported in
#                                    panels a and b but is not part of
#                                    either head-to-head, per instruction.
#
# Figure 5 (single composite figure):
#   a. Choropleth detector map of the Philippines (consensus winner per
#      region, four-tier classification scheme).
#   b. Per-region metric table with embedded per-metric significance.
#   c. Per-detector dot plot of regional dominance probabilities.
#
# Tables (CSV):
#   StageB_Regional_Framework_Metrics.csv
#   StageB_Regional_Framework_Metrics_with_CIs.csv
#   StageB_Regional_8Metric_Summary.csv
#   StageB_Regional_8Metric_Summary_with_CIs.csv
#   StageB_Regional_Dominance_Matrix.csv
#   StageB_Regional_Dominance_Probabilities.csv
#   StageB_Regional_Wilcoxon_PerMetric.csv
#   StageB_Regional_Wilcoxon_ConstantTA_vs_Comparators.csv  # vs Outbreak Threshold
#   StageB_Method_Summary_Long.csv
#   StageB_Method_Summary_Aggregate.csv
#   StageB_Regional_Bootstrap_Replicates.csv
#   StageB_Detector_Map_RegionTable.csv
#   StageB_Detector_Map_PerMetricSignificance.csv
#   StageB_Detector_Map_CrossRegion_Wilcoxon.csv
#   StageB_Detector_Map_PanelC_DotPlot.csv
#   StageB_Detector_Map_PanelD_HeatmapData.csv
#   StageB_Detector_Map_JoinAudit.csv
#   StageB_FigureB4_Legend.txt
#   StageB_FigureB5_Legend.txt
#
# FRAMEWORK
# ---------
# Two-anchor true-alarm rule:
#   T = A1 union A2
#     A1 (Actionable Window):  4-8 weeks pre-peak.
#     A2 (Epidemic Burden):    contiguous weeks centred on the seasonal
#                              peak that cumulate to 70% of seasonal cases.
#
# Two-bucket compartment scheme on True alarms:
#   Actionable: 4 <= lead <= 8 (within A1).
#   Reactive:   any other lead value on a True alarm.
#
# Same-denominator early warning timeliness:
#   Mean_Lead_Time and WP per (region, detector) divide by the count of
#   evaluable years; years that do not produce a qualifying trigger
#   contribute zero rather than NA.
#
# A1-restricted Sensitivity:
#   years_with_A1_true / years_evaluable.
#
# Year-cluster bootstrap (Cameron, Gelbach & Miller 2008):
#   B = 1000 replicates per region; year is the cluster unit because
#   within-year weekly observations are not independent.
#
# Regional inclusion criteria:
#   >= 40 observed weekly case values per evaluable region-year;
#   >= 5 evaluable years per region.
#   This is deliberately coverage-based rather than incidence-based, so all
#   17 canonical regions with adequate surveillance remain represented even
#   when weekly leptospirosis counts are low.
#
# REQUIREMENTS
# ------------
# R packages (install explicitly with R/install_dependencies.R if missing):
#   readxl, dplyr, tidyr, purrr, ggplot2, zoo, ISOweek, scales, tibble,
#   grid, patchwork, ggrepel, rlang, sf, readr, stringr.
# Optional geometry sources (one of):
#   geodata, rnaturalearth (with rnaturalearthdata).
#
# Input: regional weekly leptospirosis dataset specified by `INPUT_PATH`,
# sheet `Regional Data`, columns REGION, YR, WN, LC_DOH.
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. SETUP
# ------------------------------------------------------------------------------
SCRIPT_TITLE   <- "STAGE B - Regional leptospirosis analysis (17 Philippine regions; FigureB4 and FigureB5, plus tables)"
SCRIPT_VERSION <- "1.0 (final)"

cat("\n=============================================================\n")
cat(SCRIPT_TITLE, "\n", sep = "")
cat("Version: ", SCRIPT_VERSION, "\n", sep = "")
cat("=============================================================\n\n")

REQUIRED_PACKAGES <- c(
  "readxl", "dplyr", "tidyr", "purrr", "ggplot2", "zoo", "ISOweek",
  "scales", "tibble", "grid", "patchwork", "ggrepel", "rlang",
  "sf", "readr", "stringr", "ggspatial", "prettymapr", "viridisLite"
)
# Required for unpacking the bundled GADM PackedSpatVector RDS, but deliberately
# not attached because terra masks several functions from other plotting/data
# packages. Namespace-qualified terra:: calls are used in the geometry loader.
GEOMETRY_REQUIRED_PACKAGES <- c("terra")
# Geometry sources for the Figure 5 map. Only needed when no local boundary
# file is present in data/geometry/.
OPTIONAL_GEOM_PACKAGES <- c("geodata", "rnaturalearth", "rnaturalearthdata")

# --- Project bootstrap -------------------------------------------------------
# Portable paths + shared publication theme. Replaces the previous inline
# install.packages() loops and the hard-coded Desktop path.
.bootstrap <- function() {
  root <- Sys.getenv("TA_PROJECT_ROOT", unset = "")
  if (!nzchar(root) || !dir.exists(root)) {
    d <- normalizePath(getwd(), winslash = "/")
    for (i in seq_len(6)) {
      if (dir.exists(file.path(d, "R")) && dir.exists(file.path(d, "scripts"))) {
        root <- d; break
      }
      pp <- dirname(d); if (identical(pp, d)) break; d <- pp
    }
  }
  if (!nzchar(root)) {
    stop("Set the working directory to the project root before sourcing.",
         call. = FALSE)
  }
  source(file.path(root, "R", "00_config.R"))
}
.bootstrap()

require_packages(c(REQUIRED_PACKAGES, GEOMETRY_REQUIRED_PACKAGES), purpose = "Stage B")
invisible(lapply(REQUIRED_PACKAGES, function(p)
  suppressWarnings(suppressPackageStartupMessages(library(p, character.only = TRUE)))))

source(file.path(DIR_R, "01_publication_theme.R"), local = TRUE)

# Shared detection framework: build_ears(), build_farrington(), build_ewars(),
# the anchor/compartment helpers and hh_paired_wilcoxon().
# local = TRUE is REQUIRED -- run_all.R runs each stage with
# sys.source(envir = env); without it these functions land in globalenv and
# cannot see this stage's own constants.
source(file.path(DIR_R, "02_detection_framework.R"), local = TRUE)

# The map labels use ggrepel's text halo (bg.colour / bg.r), added in 0.9.0.
# Checked up front so an outdated install fails immediately with a clear
# message, rather than partway through building Figure 5.
if (utils::packageVersion("ggrepel") < "0.9.0") {
  stop("ggrepel >= 0.9.0 is required (map label halos use bg.colour). ",
       "Installed: ", as.character(utils::packageVersion("ggrepel")),
       "\nRun: source(\"R/install_dependencies.R\")", call. = FALSE)
}

set.seed(GLOBAL_SEED)
options(scipen = 999)

# ------------------------------------------------------------------------------
# 1. USER PARAMETERS
# ------------------------------------------------------------------------------
# Path to the input regional dengue dataset. Edit as appropriate, or use
# `file.choose()` interactively.
INPUT_PATH <- DATA_FILE
SHEET_NAME <- SHEET_REGIONAL

RF_GATE <- get_rf_gate(scope = "Regional")
OUTPUT_DIR <- file.path(DIR_OUTPUT, paste0("StageB_Regional_Lepto_Analysis", RF_GATE$suffix))
REGION_SCOPE_LABEL <- "17 Philippine regions (leptospirosis)"
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# PH_SHAPEFILE is resolved in R/00_config.R: any .shp/.gpkg/.geojson placed in
# data/geometry/ is picked up automatically, which removes the last network
# dependency and pins the administrative boundary vintage.

base_family_global <- PUB_FAMILY

# Delegates to save_pub(): vector PDF + 600 dpi PNG, dimensions clamped to the
# journal print area.
save_plot_file <- function(stem, plot_obj, width, height, dpi = 600,
                           max_height = NC_H_SUPP) {
  if (isTRUE(RF_GATE$active)) {
    rf_legend_theme <- ggplot2::theme(
      legend.position = "bottom", legend.direction = "horizontal",
      legend.box = "horizontal", legend.justification = "center",
      legend.box.just = "center",
      legend.text = ggplot2::element_text(size = 6.1, lineheight = 0.92,
                                          margin = ggplot2::margin(l = 2, r = 6)),
      legend.key.width = grid::unit(11, "pt"),
      legend.key.height = grid::unit(8.5, "pt"),
      legend.spacing.x = grid::unit(5, "pt"),
      legend.spacing.y = grid::unit(2, "pt"),
      legend.margin = ggplot2::margin(3, 6, 3, 6),
      legend.box.spacing = grid::unit(4, "pt"),
      plot.margin = ggplot2::margin(14, 24, 18, 14)
    )
    if (inherits(plot_obj, "patchwork")) plot_obj <- plot_obj & rf_legend_theme
    else if (inherits(plot_obj, "ggplot")) plot_obj <- plot_obj + rf_legend_theme
  }
  save_pub(stem = sub("\\.(pdf|png)$", "", stem), plot = plot_obj,
           width = width, height = height, dir = OUTPUT_DIR, dpi = dpi,
           max_height = max_height)
}

# ------------------------------------------------------------------------------
# 2. FRAMEWORK CONSTANTS
# ------------------------------------------------------------------------------
# Constant TA (Vaezi-style hysteresis)
# eta_ON / eta_OFF are read from Stage 1's derivation rather than hard-coded,
# so re-running Stage 1 propagates here. If Stage 1 output is absent,
# load_eta_thresholds() falls back to the published constants (1.33 / 0.73)
# with a warning.
.eta    <- load_eta_thresholds()
ETA_ON  <- .eta$ETA_ON
ETA_OFF <- .eta$ETA_OFF

# STA/LTA windows for both TA detectors are also read from Stage 1, which
# now runs an empirical window-optimisation step (not present in the
# original dengue pipeline, which fixed these at 4/26 and 3/12) and
# selects the best-performing pair from a candidate grid on the NCR anchor
# series. Falls back to the published windows with a warning if Stage 1's
# window-optimisation output is absent. This is the single source of the
# STA/LTA windows for every stage in this pipeline: StageA, StageB and
# StageC all call the same loader, so the same optimised windows are used
# throughout. See Item 7 of the Delphi review.
.win <- load_sta_lta_windows()
STA_WIN_CONSTANT   <- .win$STA_WIN_CONSTANT
LTA_WIN_CONSTANT   <- .win$LTA_WIN_CONSTANT
STA_WIN_CONTINUOUS <- .win$STA_WIN_CONTINUOUS
LTA_WIN_CONTINUOUS <- .win$LTA_WIN_CONTINUOUS
GUARD              <- 2L  # guard band, Constant TA only; not empirically varied

# Continuous TA uses the same activation threshold as Constant TA's eta_ON,
# per the original design (a single ratio cut, no hysteresis).
ETA_ON_CLASSIC <- ETA_ON

# Vaezi STA/LTA OFF-state reset (consecutive-off weeks forcing LTA recompute).
MIN_OFF_RESET  <- 8L

# Anchor parameters (two-anchor framework: A1 and A2; no A3).
A1_LEAD_MIN    <- 4L
A1_LEAD_MAX    <- 8L
# A2_BURDEN_FRAC (0.70) and A2_BURDEN_PCT ("70%") are defined in R/00_config.R.
# Fixed by design: this stage reads no file and depends on no other stage for it.

# Lead-time compartment thresholds.
COMPARTMENT_ACTIONABLE_MIN <- 4L  # lower bound of A1
COMPARTMENT_ACTIONABLE_MAX <- 8L  # upper bound of A1

# Year inclusion. LC_DOH has no coverage before 2018 (blank in the source
# sheet), and 2020-2021 are absent from the sheet entirely for every
# disease. 2025 is kept (see header note): its LC_DOH panel is complete,
# not truncated, so there is no reason to drop it here as the dengue script
# did.
EXCLUDED_YEARS  <- c(2020L, 2021L)
TARGET_YEARS    <- 2018L:2025L
EVALUABLE_YEARS <- setdiff(TARGET_YEARS, EXCLUDED_YEARS)  # n = 6

# Regional inclusion is coverage-based, not incidence-based. All 17 regions
# have complete leptospirosis surveillance coverage in the evaluable years;
# low-incidence regions must not disappear from the national map merely because
# their annual peak is small. A region-year is evaluable when at least 40 weekly
# case observations are present, and a region must contribute at least 5 such
# years. With the supplied workbook this retains all 17 canonical regions.
MIN_OBSERVED_CASE_WEEKS_PER_YEAR <- 40L
MIN_EVALUABLE_YEARS_PER_REGION   <- 5L

# Bootstrap configuration
BOOT_N_CI <- REGIONAL_BOOT_N_CI
# N_BOOTS is an alias used throughout the pairwise-test and dominance-summary
# helpers below; defined here (immediately after BOOT_N_CI) rather than near
# its heaviest use further down, since it is referenced as an eager default
# argument value (B = N_BOOTS) well before that point in the script.
N_BOOTS   <- BOOT_N_CI

# Canonical list of the 17 Philippine administrative regions. Defined here
# (moved up from its original position in the "REGION CANONICALISER" section
# further down) because Supplementary Table 13's construction, earlier in
# this script, references it directly.
CANONICAL_17 <- c("BARMM", "CAR", "MIMAROPA", "NCR",
                  "REGION I", "REGION II", "REGION III", "REGION IV-A",
                  "REGION V", "REGION VI", "REGION VII", "REGION VIII",
                  "REGION IX", "REGION X", "REGION XI", "REGION XII",
                  "REGION XIII")

# Significance and Bonferroni configuration
SIG_LEVEL          <- 0.05
HH_ALPHA           <- 0.05
BONF_LEVEL1_FACTOR <- 5L  # 5 metrics within (region x ordered pair)
BONF_LEVEL2_FACTOR <- 2L  # 2 pairs within metric (Constant TA vs OT, Continuous TA vs OT)

# Within-region dominance score threshold
DOMINANCE_THRESHOLD <- 0.75
CONFIDENCE_HI       <- 0.90

# ------------------------------------------------------------------------------
# 3. COLOR SYSTEM AND THEMES
# ------------------------------------------------------------------------------
TYPE_COLORS <- c(
  "Surveillance-Guideline Percentile Thresholds" = "#F8766D",
  "Retrospective Thresholds" = "#00BA38",
  "Acceleration Measures" = "#619CFF"
)
type_order <- c(
  "Surveillance-Guideline Percentile Thresholds",
  "Retrospective Thresholds",
  "Acceleration Measures"
)

# Built on the shared publication theme so typography matches every other
# stage; only the dashboard legend placement is overridden.
theme_dashboard <- function(base_size = PUB_BASE, base_family = PUB_FAMILY) {
  theme_pub(base_size, base_family) %+replace%
    ggplot2::theme(
      legend.position  = "bottom",
      legend.direction = "horizontal"
    )
}

# ------------------------------------------------------------------------------
# 4. DATA LOADING
# ------------------------------------------------------------------------------
if (!file.exists(INPUT_PATH)) {
  stop(paste0("Input file not found at:\n", INPUT_PATH))
}

df_raw <- readxl::read_excel(INPUT_PATH, sheet = SHEET_NAME)

required_cols <- c("REGION", "YR", "WN", "LC_DOH", "RF_HDX")
missing_cols  <- setdiff(required_cols, names(df_raw))
if (length(missing_cols) > 0) {
  stop(paste0("Missing required columns in sheet '", SHEET_NAME, "': ",
              paste(missing_cols, collapse = ", ")))
}

df_all <- df_raw |>
  dplyr::mutate(
    REGION = as.character(REGION),
    YR     = suppressWarnings(as.integer(YR)),
    WN     = suppressWarnings(as.integer(WN)),
    LC_DOH = suppressWarnings(as.numeric(LC_DOH)),
    RF_HDX = suppressWarnings(as.numeric(RF_HDX))
  ) |>
  dplyr::filter(!is.na(REGION), !is.na(YR), !is.na(WN),
                WN >= 1, WN <= 53) |>
  dplyr::mutate(
    ISOweek = sprintf("%d-W%02d", YR, WN),
    Date    = ISOweek::ISOweek2date(paste0(ISOweek, "-1"))
  ) |>
  dplyr::filter(!is.na(Date)) |>
  dplyr::arrange(REGION, Date) |>
  dplyr::filter(YR %in% TARGET_YEARS)

if (nrow(df_all) == 0) {
  stop("No rows remain after filtering to target years.")
}

# Anchor functions reference a column named DC_QC. Add an alias so the
# regional pipeline can use them verbatim.
df_all <- df_all |> dplyr::mutate(DC_QC = LC_DOH)

# ------------------------------------------------------------------------------
# 5. SURGE GENERATION (PER REGION)
# ------------------------------------------------------------------------------
# This leptospirosis adaptation computes the six requested comparison
# detectors: Constant TA, Continuous TA, Outbreak Threshold, Alarm Threshold,
# and the 75th/90th historical percentile comparators. Additional detector
# families from the dengue reference are not part of the active benchmark.
# #
# Both transmission-acceleration measures now use an empirically DERIVED
# STA/LTA window rather than a fixed one: STA_WIN_CONTINUOUS / LTA_WIN_
# CONTINUOUS and STA_WIN_CONSTANT / LTA_WIN_CONSTANT are read from Stage 1's
# window-optimisation output (falling back to the published windows, 3/12
# and 4/26, with a warning if Stage 1 has not been run). See Item 7 of the
# Delphi review and the "STA/LTA WINDOW OPTIMISATION" section of Stage 1.
#
# Calendar continuity is enforced explicitly for STA/LTA and RF rolling windows,
# so no detector lookback can bridge the missing 2020-2021 surveillance gap.
#
# Per-region detector pipeline, retained detectors only:
#   - Week-specific baseline (mean, sd, Alarm Threshold = mean+1sd,
#     Outbreak Threshold = mean+2sd) computed from up to 5 prior donor years.
#   - Continuous TA: STA_WIN_CONTINUOUS / LTA_WIN_CONTINUOUS moving-average
#     ratio, no hysteresis.
#   - Constant TA: Vaezi STA/LTA hysteresis at STA_WIN_CONSTANT /
#     LTA_WIN_CONSTANT, with year-boundary reset, MIN_OFF_RESET
#     consecutive-off weeks forcing fresh LTA recompute, and >= for ON
#     activation.
# ------------------------------------------------------------------------------

safe_quantile <- function(x, probs, na.rm = TRUE) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(rep(NA_real_, length(probs)))
  suppressWarnings(as.numeric(stats::quantile(
    x, probs = probs, na.rm = na.rm, names = FALSE
  )))
}

create_surges <- function(df_region) {
  df_region <- df_region |> dplyr::arrange(Date)
  n <- nrow(df_region)

  # 5.1 Per-region rolling donor map -----------------------------------------
  region_target_years <- sort(unique(df_region$YR))
  donor_pool_region   <- region_target_years
  rolling_map <- list()
  for (y in region_target_years) {
    donors <- intersect(seq(y - 5L, y - 1L), donor_pool_region)
    if (length(donors) < 3L) {
      donors <- tail(donor_pool_region[donor_pool_region < y], 5L)
    }
    rolling_map[[as.character(y)]] <- sort(unique(donors))
  }

  # 5.2 Rolling weekly baseline (week-specific, donor-year based) ------------
  # Compute the same donor-year baseline family used in the NCR comparison:
  # Alarm = mean + 1 SD, Outbreak = mean + 2 SD, plus the requested 75th and
  # 90th percentile comparators. These names describe the analytical
  # comparators in this pipeline and do not assert external validation.
  bl_rows <- list()
  for (y in names(rolling_map)) {
    donors <- rolling_map[[y]]; y_int <- as.integer(y)
    if (length(donors) == 0L) next
    weeks_y <- df_region |> dplyr::filter(YR == y_int) |>
      dplyr::pull(WN) |> unique() |> sort()
    if (length(weeks_y) == 0L) next
    for (w in weeks_y) {
      vals <- df_region |> dplyr::filter(YR %in% donors, WN == w) |>
        dplyr::pull(LC_DOH)
      vals <- vals[is.finite(vals)]
      n_vals <- length(vals)
      if (n_vals == 0L) {
        mean_val <- NA_real_; sd_val <- NA_real_; p75 <- NA_real_; p90 <- NA_real_
      } else if (n_vals == 1L) {
        mean_val <- vals[1]; sd_val <- NA_real_; p75 <- vals[1]; p90 <- vals[1]
      } else {
        mean_val <- mean(vals, na.rm = TRUE)
        sd_val   <- stats::sd(vals, na.rm = TRUE)
        p75 <- safe_quantile(vals, 0.75); p90 <- safe_quantile(vals, 0.90)
      }
      bl_rows[[length(bl_rows) + 1L]] <- data.frame(
        YR = y_int, WN = as.integer(w),
        bl_mean = mean_val, bl_sd = sd_val, bl_p75 = p75, bl_p90 = p90,
        bl_alarm    = if (is.na(sd_val)) mean_val else mean_val + 1 * sd_val,
        bl_outbreak = if (is.na(sd_val)) mean_val else mean_val + 2 * sd_val,
        bl_n = n_vals,
        stringsAsFactors = FALSE
      )
    }
  }
  thresholds_region <- if (length(bl_rows) > 0L) {
    dplyr::bind_rows(bl_rows) |> dplyr::arrange(YR, WN)
  } else {
    data.frame(
      YR = integer(), WN = integer(),
      bl_mean = numeric(), bl_sd = numeric(), bl_p75 = numeric(), bl_p90 = numeric(),
      bl_alarm = numeric(), bl_outbreak = numeric(), bl_n = integer(),
      stringsAsFactors = FALSE
    )
  }
  df_region <- df_region |>
    dplyr::left_join(thresholds_region, by = c("YR", "WN"))

  # 5.3 Continuous TA: STA_WIN_CONTINUOUS / LTA_WIN_CONTINUOUS ratio ---------
  dc_ma_sta <- rolling_contiguous_mean(df_region$LC_DOH, df_region$Date, STA_WIN_CONTINUOUS)
  dc_ma_lta <- rolling_contiguous_mean(df_region$LC_DOH, df_region$Date, LTA_WIN_CONTINUOUS)

  # 5.4 Primary detector signals ---------------------------------------------
  df_region <- df_region |>
    dplyr::mutate(
      surge_mean_2sd = as.integer(!is.na(LC_DOH) & !is.na(bl_outbreak) &
                                    LC_DOH > bl_outbreak),
      surge_mean_1sd = as.integer(!is.na(LC_DOH) & !is.na(bl_alarm) &
                                    LC_DOH > bl_alarm),
      surge_who_75 = as.integer(!is.na(LC_DOH) & !is.na(bl_p75) & LC_DOH > bl_p75),
      surge_who_90 = as.integer(!is.na(LC_DOH) & !is.na(bl_p90) & LC_DOH > bl_p90),
      surge_sta_lta  = as.integer(!is.na(dc_ma_sta) & !is.na(dc_ma_lta) &
                                    dc_ma_lta > 0 &
                                    (dc_ma_sta / dc_ma_lta) > ETA_ON_CLASSIC)
    )

  # 5.5 Vaezi-style STA/LTA hysteresis (Constant TA) -------------------------
  MIN_T <- STA_WIN_CONSTANT + GUARD + LTA_WIN_CONSTANT
  dc <- df_region$LC_DOH
  triggered_v <- rep(FALSE, n)
  R_vaezi_v   <- rep(NA_real_, n)
  is_on <- FALSE; frozen_lta <- NA_real_; consec_off <- 0L
  for (t in seq_len(n)) {
    if (t > 1L && !is.na(df_region$YR[t]) &&
        !is.na(df_region$YR[t - 1L]) &&
        df_region$YR[t] != df_region$YR[t - 1L]) {
      is_on <- FALSE; frozen_lta <- NA_real_; consec_off <- 0L
    }
    if (!is_on) consec_off <- consec_off + 1L else consec_off <- 0L
    if (!is_on && consec_off >= MIN_OFF_RESET) frozen_lta <- NA_real_
    if (t < MIN_T) next

    sta_start <- t - STA_WIN_CONSTANT + 1L
    full_start <- t - STA_WIN_CONSTANT - GUARD - LTA_WIN_CONSTANT + 1L
    if (full_start < 1L || !is_contiguous_week_span(df_region$Date, full_start, t)) next
    sta_vals <- dc[sta_start:t]
    sta <- if (all(is.na(sta_vals))) NA_real_ else mean(sta_vals, na.rm = TRUE)

    if (!is_on || is.na(frozen_lta)) {
      lta_idx <- full_start:(t - STA_WIN_CONSTANT - GUARD)
      if (length(lta_idx) == LTA_WIN_CONSTANT) {
        lta_vals <- dc[lta_idx]
        frozen_lta <- if (all(is.na(lta_vals))) NA_real_ else mean(lta_vals, na.rm = TRUE)
      } else {
        frozen_lta <- NA_real_
      }
    }
    R_t <- if (!is.na(frozen_lta) && frozen_lta > 0 && !is.na(sta))
      sta / frozen_lta else NA_real_
    R_vaezi_v[t] <- R_t

    if (!is_on && !is.na(R_t) && R_t >= ETA_ON) { is_on <- TRUE; consec_off <- 0L }
    if ( is_on && !is.na(R_t) && R_t < ETA_OFF) { is_on <- FALSE; frozen_lta <- NA_real_ }
    triggered_v[t] <- is_on
  }
  df_region$surge_sta_lta_vaezi <- as.integer(triggered_v)
  df_region$R_vaezi             <- R_vaezi_v

  df_region$ungated_surge_sta_lta_vaezi <- df_region$surge_sta_lta_vaezi
  df_region$ungated_surge_mean_2sd <- df_region$surge_mean_2sd

  df_region$RF_Antecedent_4wk_mm <- preceding_contiguous_sum(
    df_region$RF_HDX, df_region$Date, RF_PRECEDING_WEEKS)
  df_region$RF_STA_LTA_Rt <- NA_real_
  df_region$RF_Antecedent_4wk_Rt_Max <- NA_real_
  df_region$RF_Gate_Feature <- NA_real_
  df_region$RF_signal_formal_ok <- FALSE
  df_region$RF_signal_precautionary_ok <- FALSE
  df_region$RF_signal_strong_ok <- FALSE
  df_region$RF_gate_ok <- FALSE
  df_region$RF_gate_precautionary_ok <- FALSE
  df_region$RF_gate_strong_ok <- FALSE
  df_region$RF_gate_source_index <- NA_integer_
  df_region$RF_gate_lead_weeks <- NA_integer_
  df_region$RF_gate_source_WN <- NA_integer_
  df_region$RF_gate_source_window_start_WN <- NA_integer_
  df_region$RF_gate_source_window_end_WN <- NA_integer_
  df_region$RF_gate_source_feature <- NA_real_
  df_region$RF_gate_source_status <- NA_character_
  df_region$RF_Activation_Status <- NA_character_
  if (RF_GATE$active) {
    RF_GATE$operator <- normalize_rf_operator(RF_GATE$operator)
    if (identical(RF_GATE$scale, "mm")) {
      gate_feature <- df_region$RF_Antecedent_4wk_mm
    } else if (identical(RF_GATE$scale, "rt")) {
      df_region$RF_STA_LTA_Rt <- guarded_sta_lta_ratio(
        df_region$RF_HDX, df_region$Date,
        RF_GATE$sta, RF_GATE$lta, RF_GATE$guard)
      df_region$RF_Antecedent_4wk_Rt_Max <- preceding_contiguous_max(
        df_region$RF_STA_LTA_Rt, df_region$Date, RF_PRECEDING_WEEKS)
      gate_feature <- df_region$RF_Antecedent_4wk_Rt_Max
    } else {
      stop("Unsupported RF gate scale: ", RF_GATE$scale, call. = FALSE)
    }
    signal_formal <- apply_rf_gate(gate_feature, RF_GATE$threshold, RF_GATE$operator)
    signal_precautionary <- apply_rf_gate(gate_feature, RF_GATE$activation_bound, RF_GATE$operator)
    signal_strong <- apply_rf_gate(gate_feature, RF_GATE$strong_bound, RF_GATE$operator)
    # The RF feature at disease week t already summarizes t-4..t-1. No extra
    # temporal shift is applied. This exactly matches Stage C derivation.
    .idx <- seq_len(nrow(df_region))
    make_direct_q <- function(sig) list(
      qualified = as.logical(sig),
      source_index = ifelse(as.logical(sig), .idx, NA_integer_),
      lead_weeks = ifelse(as.logical(sig), 1L, NA_integer_)
    )
    q_formal <- make_direct_q(signal_formal)
    q_precautionary <- make_direct_q(signal_precautionary)
    q_strong <- make_direct_q(signal_strong)
    df_region$RF_Gate_Feature <- suppressWarnings(as.numeric(gate_feature))
    df_region$RF_signal_formal_ok <- as.logical(signal_formal)
    df_region$RF_signal_precautionary_ok <- as.logical(signal_precautionary)
    df_region$RF_signal_strong_ok <- as.logical(signal_strong)
    df_region$RF_gate_ok <- as.logical(q_formal$qualified)
    df_region$RF_gate_precautionary_ok <- as.logical(q_precautionary$qualified)
    df_region$RF_gate_strong_ok <- as.logical(q_strong$qualified)
    df_region$RF_gate_source_index <- q_formal$source_index
    df_region$RF_gate_lead_weeks <- q_formal$lead_weeks
    df_region$RF_Activation_Status <- classify_rf_activation(
      gate_feature, RF_GATE$threshold, RF_GATE$lower, RF_GATE$upper, RF_GATE$operator)
    .src_ok <- which(!is.na(q_formal$source_index))
    if (length(.src_ok)) {
      .src <- q_formal$source_index[.src_ok]
      .sig_date <- df_region$Date[.src] - 7
      .sig_wn <- as.integer(df_region$WN[.src] - 1L)
      .same_year <- as.integer(format(.sig_date, "%Y")) == df_region$YR[.src_ok]
      .sig_wn[!.same_year | .sig_wn < 1L] <- NA_integer_
      df_region$RF_gate_source_WN[.src_ok] <- .sig_wn
      df_region$RF_gate_source_window_start_WN[.src_ok] <- pmax(1L, as.integer(df_region$WN[.src] - RF_PRECEDING_WEEKS))
      df_region$RF_gate_source_window_end_WN[.src_ok] <- pmax(1L, as.integer(df_region$WN[.src] - 1L))
      df_region$RF_gate_source_feature[.src_ok] <- suppressWarnings(as.numeric(gate_feature[.src]))
      df_region$RF_gate_source_status[.src_ok] <- classify_rf_activation(
        gate_feature[.src], RF_GATE$threshold, RF_GATE$lower, RF_GATE$upper, RF_GATE$operator)
    }
    gate_cols <- c("surge_sta_lta_vaezi","surge_sta_lta","surge_mean_2sd",
                   "surge_mean_1sd","surge_who_75","surge_who_90")
    for (nm in gate_cols) df_region[[nm]] <- as.integer(df_region[[nm]] == 1L & df_region$RF_gate_ok)
  }
  df_region
}

cat("Generating per-region surge signals for the six requested detectors...\n")
df_all <- df_all |>
  dplyr::group_by(REGION) |>
  dplyr::group_modify(~ create_surges(.x)) |>
  dplyr::ungroup()
if (RF_GATE$active) {
  utils::write.csv(data.frame(
    RF_Product=RF_GATE$product, Label=RF_GATE$label, Threshold_Scope=RF_GATE$scope, RF_Scale=RF_GATE$scale,
    Preceding_Weeks=RF_PRECEDING_WEEKS, Threshold=RF_GATE$threshold,
    Lower_Bound=RF_GATE$lower, Upper_Bound=RF_GATE$upper,
    Precautionary_Activation_Bound=RF_GATE$activation_bound,
    Strong_Activation_Bound=RF_GATE$strong_bound,
    Operator=RF_GATE$operator, Primary_Metrics_Threshold="Point estimate",
    Units=ifelse(RF_GATE$scale == "mm", "mm", "R(t)"),
    RF_STA=RF_GATE$sta, RF_LTA=RF_GATE$lta, RF_Guard=RF_GATE$guard,
    Trigger_Week_Excluded=TRUE,
    Overlay_Definition=paste0(
      "The RF feature is computed from t-4..t-1 and evaluated at the disease-trigger week; trigger-week rainfall is excluded.")),
    file.path(OUTPUT_DIR,"RF_Gate_Metadata.csv"), row.names=FALSE)

  rf_operational_weekly_regional <- df_all |>
    dplyr::transmute(
      REGION, YR, WN, Date,
      RF_Product = RF_GATE$product,
      RF_Scale = RF_GATE$scale,
      RF_Feature = RF_Gate_Feature,
      Threshold_Lower = RF_GATE$lower,
      Threshold_Point = RF_GATE$threshold,
      Threshold_Upper = RF_GATE$upper,
      Comparison_Operator = RF_GATE$operator,
      Precautionary_RF_Signal = RF_signal_precautionary_ok,
      Formal_RF_Signal = RF_signal_formal_ok,
      Strong_RF_Signal = RF_signal_strong_ok,
      Antecedent_t4_to_t1_Meets_Point_Threshold = RF_gate_ok,
      RF_Activation_Status
    )
  readr::write_csv(rf_operational_weekly_regional,
                   file.path(OUTPUT_DIR, "RF_Operational_Weekly_Regional.csv"), na = "")
  writeLines(c(
    paste0("Primary RF product: ", RF_GATE$label),
    "Primary RF-specific outputs are RF_to_Anchor_Trigger_Detail.csv, RF_to_Anchor_Regional_Summary.csv, FigureB6/B6A RF overlays, and FigureB7 RF-to-anchor map.",
    "FigureB5 and the six-detector tables are retained as diagnostic case-detector context and should not be interpreted as the primary RF product visualization.",
    "P1/P3 map RF preceding Constant TA; P2/P4 map RF preceding Outbreak Threshold. A1/A2 remains the true/false alarm standard."
  ), file.path(OUTPUT_DIR, "RF_Product_Output_Guide.txt"))

  rf_anchor_col <- if (RF_GATE$product %in% c(1L, 3L))
    "ungated_surge_sta_lta_vaezi" else "ungated_surge_mean_2sd"
  rf_anchor_label <- if (RF_GATE$product %in% c(1L, 3L))
    "Constant TA" else "Outbreak Threshold"
  rf_scale_label <- if (identical(RF_GATE$scale, "mm")) "RF mm" else "RF STA/LTA R(t)"

  regional_rf_overlay <- df_all |>
    dplyr::group_by(REGION, YR) |>
    dplyr::group_modify(~ {
      z <- .x |> dplyr::arrange(WN)
      av <- as.integer(z[[rf_anchor_col]] == 1L)
      starts <- which(av == 1L & c(TRUE, head(av, -1L) == 0L))
      if (!length(starts)) {
        return(data.frame(
          RF_Product = RF_GATE$product, RF_Scale = rf_scale_label,
          Anchor = rf_anchor_label, RF_Trigger_WN = NA_integer_,
          Anchor_Trigger_WN = NA_integer_,
          RF_Warning_Horizon_Start_WN = NA_integer_, RF_Warning_Horizon_End_WN = NA_integer_,
          RF_Window_Start_WN = NA_integer_, RF_Window_End_WN = NA_integer_,
          RF_Feature_At_RF_Trigger = NA_real_, RF_Max_Preceding_Feature = NA_real_, RF_Feature_At_Anchor = NA_real_,
          RF_Threshold = RF_GATE$threshold,
          RF_Threshold_Lower = RF_GATE$lower, RF_Threshold_Upper = RF_GATE$upper,
          RF_Precautionary_Bound = RF_GATE$activation_bound,
          RF_Strong_Bound = RF_GATE$strong_bound,
          RF_Activation_Status = NA_character_,
          RF_Operator = RF_GATE$operator,
          RF_Lead_Weeks = NA_integer_, stringsAsFactors = FALSE))
      }
      i <- starts[1L]
      anchor_wk <- as.integer(z$WN[i])
      rf_wk <- suppressWarnings(as.integer(z$RF_gate_source_WN[i]))
      src_feature <- suppressWarnings(as.numeric(z$RF_gate_source_feature[i]))
      signal_week <- as.integer(z$WN - 1L)
      warning_idx <- which(signal_week >= max(1L, anchor_wk - RF_PRECEDING_WEEKS) &
                           signal_week <= anchor_wk - 1L & is.finite(z$RF_Gate_Feature))
      max_pre <- if (length(warning_idx)) max(z$RF_Gate_Feature[warning_idx], na.rm = TRUE) else NA_real_
      anchor_feature <- suppressWarnings(as.numeric(z$RF_Gate_Feature[i]))
      display_status <- if (is.finite(anchor_feature))
        classify_rf_activation(anchor_feature, RF_GATE$threshold, RF_GATE$lower, RF_GATE$upper, RF_GATE$operator) else NA_character_
      data.frame(
        RF_Product = RF_GATE$product, RF_Scale = rf_scale_label,
        Anchor = rf_anchor_label, RF_Trigger_WN = rf_wk,
        Anchor_Trigger_WN = anchor_wk,
        RF_Warning_Horizon_Start_WN = as.integer(max(1L, anchor_wk - RF_PRECEDING_WEEKS)),
        RF_Warning_Horizon_End_WN = as.integer(max(1L, anchor_wk - 1L)),
        RF_Window_Start_WN = suppressWarnings(as.integer(z$RF_gate_source_window_start_WN[i])),
        RF_Window_End_WN = suppressWarnings(as.integer(z$RF_gate_source_window_end_WN[i])),
        RF_Feature_At_RF_Trigger = src_feature,
        RF_Max_Preceding_Feature = max_pre,
        RF_Feature_At_Anchor = anchor_feature,
        RF_Threshold = RF_GATE$threshold,
        RF_Threshold_Lower = RF_GATE$lower, RF_Threshold_Upper = RF_GATE$upper,
        RF_Precautionary_Bound = RF_GATE$activation_bound,
        RF_Strong_Bound = RF_GATE$strong_bound,
        RF_Activation_Status = display_status,
        RF_Operator = RF_GATE$operator,
        RF_Lead_Weeks = suppressWarnings(as.integer(z$RF_gate_lead_weeks[i])),
        stringsAsFactors = FALSE)
    }) |>
    dplyr::ungroup()
  utils::write.csv(regional_rf_overlay,
                   file.path(OUTPUT_DIR, "RF_Trigger_Overlay_Regional.csv"),
                   row.names = FALSE, na = "")

  # RF-focused regional figures. Weekly leptospirosis cases are light-blue
  # bars. The right axis is the RF gate's own environmental scale: mm for
  # Products 1/2 and empirically derived rainfall STA/LTA R(t) for Products 3/4.
  # No leptospirosis Constant/Continuous TA ratio curve is drawn.
  rf_is_mm <- identical(RF_GATE$scale, "mm")
  rf_axis_name <- if (rf_is_mm)
    "Antecedent rainfall, t-4 to t-1 (mm)" else
    "Antecedent rainfall STA/LTA R(t)"
  rf_series_label <- if (rf_is_mm)
    "Antecedent rainfall (mm)" else "Antecedent rainfall STA/LTA R(t)"
  rf_threshold_label <- if (rf_is_mm)
    "Adopted RF threshold (point estimate)" else "Adopted RF R(t) threshold (point estimate)"
  rf_bounds_label <- if (rf_is_mm)
    "RF threshold lower/upper bounds" else "RF R(t) lower/upper bounds"
  anchor_legend_label <- if (RF_GATE$product %in% c(1L, 3L))
    "Constant TA trigger"
  else "Outbreak threshold trigger"
  anchor_colour <- if (RF_GATE$product %in% c(1L, 3L)) "#D55E00" else "black"

  build_regional_rf_plot <- function(plot_y, overlay_y, title_text, facet_formula) {
    rf_feature_col <- if (rf_is_mm) "RF_Antecedent_4wk_mm" else "RF_Antecedent_4wk_Rt_Max"
    plot_y$.rf_feature_display <- suppressWarnings(as.numeric(plot_y[[rf_feature_col]]))

    case_upper <- suppressWarnings(max(plot_y$LC_DOH,
      if (RF_GATE$product %in% c(2L, 4L)) plot_y$bl_outbreak else NA_real_, na.rm = TRUE))
    if (!is.finite(case_upper) || case_upper <= 0) case_upper <- 1
    case_upper <- case_upper * 1.15
    rf_upper <- suppressWarnings(max(plot_y$.rf_feature_display, RF_GATE$threshold,
                                      RF_GATE$lower, RF_GATE$upper, na.rm = TRUE))
    if (!is.finite(rf_upper) || rf_upper <= 0) {
      rf_upper <- max(abs(c(RF_GATE$threshold, RF_GATE$lower, RF_GATE$upper)), 1, na.rm = TRUE)
    }
    rf_upper <- rf_upper * 1.15
    rf_scale <- case_upper / rf_upper

    overlay_y <- overlay_y |>
      dplyr::mutate(
        .op = vapply(as.character(RF_Operator), normalize_rf_operator, character(1)),
        .op_symbol = ifelse(.op == ">=", "\u2265", "\u2264"),
        .units = ifelse(rf_is_mm, "mm", "R(t)"),
        .digits = ifelse(rf_is_mm, 1L, 2L),
        .display_rf = dplyr::if_else(is.finite(RF_Feature_At_RF_Trigger),
                                      RF_Feature_At_RF_Trigger, RF_Max_Preceding_Feature),
        .gate_label = dplyr::if_else(
          is.finite(.display_rf),
          paste0("RF = ",
                 format(round(.display_rf, ifelse(rf_is_mm, 1L, 2L)), trim = TRUE),
                 " ", .units, "\nthreshold = ",
                 format(round(RF_Threshold, ifelse(rf_is_mm, 1L, 2L)), trim = TRUE),
                 " ", .units, " [",
                 format(round(RF_Threshold_Lower, ifelse(rf_is_mm, 1L, 2L)), trim = TRUE),
                 "-",
                 format(round(RF_Threshold_Upper, ifelse(rf_is_mm, 1L, 2L)), trim = TRUE),
                 " ", .units, "]; ", RF_Activation_Status),
          ""),
        .label_x = pmax(3, RF_Warning_Horizon_Start_WN),
        .label_y = case_upper * 0.94
      )

    p <- ggplot2::ggplot(plot_y, ggplot2::aes(x = WN)) +
      ggplot2::geom_col(
        ggplot2::aes(y = LC_DOH, fill = "Leptospirosis cases"),
        width = 0.84, alpha = 0.90, na.rm = TRUE) +
      ggplot2::geom_line(
        ggplot2::aes(y = .rf_feature_display * rf_scale,
                     colour = rf_series_label, linetype = rf_series_label),
        linewidth = 0.72, na.rm = TRUE) +
      ggplot2::geom_hline(
        data = data.frame(.ythr = RF_GATE$threshold * rf_scale, .lab = rf_threshold_label),
        ggplot2::aes(yintercept = .ythr, colour = .lab, linetype = .lab),
        linewidth = 0.66, inherit.aes = FALSE, na.rm = TRUE) +
      ggplot2::geom_hline(
        data = data.frame(.ythr = c(RF_GATE$lower, RF_GATE$upper) * rf_scale,
                          .lab = rf_bounds_label),
        ggplot2::aes(yintercept = .ythr, colour = .lab, linetype = .lab),
        linewidth = 0.48, inherit.aes = FALSE, na.rm = TRUE) +
      ggplot2::geom_rect(
        data = overlay_y,
        ggplot2::aes(xmin = RF_Warning_Horizon_Start_WN - 0.5,
                     xmax = RF_Warning_Horizon_End_WN + 0.5, ymin = 0, ymax = Inf),
        inherit.aes = FALSE, fill = "#0072B2", alpha = 0.04, na.rm = TRUE) +
      ggplot2::geom_vline(
        data = overlay_y,
        ggplot2::aes(xintercept = RF_Trigger_WN,
                     colour = "RF trigger (precedes anchor)", linetype = "RF trigger (precedes anchor)"),
        linewidth = 0.62, na.rm = TRUE) +
      ggplot2::geom_vline(
        data = overlay_y,
        ggplot2::aes(xintercept = Anchor_Trigger_WN,
                     colour = anchor_legend_label, linetype = anchor_legend_label),
        linewidth = 0.68, na.rm = TRUE)

    if (RF_GATE$product %in% c(2L, 4L) && "bl_outbreak" %in% names(plot_y)) {
      p <- p + ggplot2::geom_line(
        ggplot2::aes(y = bl_outbreak,
                     colour = "Outbreak threshold (mean + 2SD)",
                     linetype = "Outbreak threshold (mean + 2SD)"),
        linewidth = 0.62, na.rm = TRUE)
    }

    p <- p + ggplot2::geom_label(
      data = overlay_y |> dplyr::filter(nzchar(.gate_label)),
      ggplot2::aes(x = .label_x, y = .label_y, label = .gate_label),
      inherit.aes = FALSE, hjust = 0, vjust = 1,
      size = pub_text_size(PUB_ANNOT - 0.5), family = base_family_global,
      fontface = "bold", linewidth = 0.20, fill = "white", colour = "black",
      na.rm = TRUE)

    breaks <- c(rf_series_label, rf_threshold_label, rf_bounds_label,
                "RF trigger (precedes anchor)", anchor_legend_label)
    vals <- c("#0072B2", "#2166AC", "#67A9CF", "#009E73", anchor_colour)
    ltys <- c("solid", "dashed", "dotted", "dotdash", "longdash")
    if (RF_GATE$product %in% c(2L, 4L)) {
      breaks <- c(breaks, "Outbreak threshold (mean + 2SD)")
      vals <- c(vals, "black")
      ltys <- c(ltys, "dashed")
    }

    p +
      facet_formula +
      ggplot2::scale_fill_manual(
        name = NULL, values = c("Leptospirosis cases" = "#A6CEE3"),
        breaks = "Leptospirosis cases",
        guide = ggplot2::guide_legend(order = 1, nrow = 1, direction = "horizontal")) +
      ggplot2::scale_colour_manual(
        name = NULL, values = stats::setNames(vals, breaks), limits = breaks, breaks = breaks,
        drop = FALSE, guide = ggplot2::guide_legend(order = 2, nrow = 2, byrow = TRUE, direction = "horizontal",
                                                       keywidth = grid::unit(11, "pt"), keyheight = grid::unit(8.5, "pt"))) +
      ggplot2::scale_linetype_manual(
        name = NULL, values = stats::setNames(ltys, breaks), breaks = breaks,
        drop = FALSE, guide = "none") +
      ggplot2::scale_x_continuous(
        name = "Week number", breaks = c(1, 13, 26, 39, 52),
        labels = paste0("W", c(1, 13, 26, 39, 52))) +
      ggplot2::scale_y_continuous(
        name = "Weekly leptospirosis cases", breaks = scales::pretty_breaks(n = 5),
        sec.axis = ggplot2::sec_axis(~ . / rf_scale, name = rf_axis_name,
                                     breaks = scales::pretty_breaks(n = 5))) +
      ggplot2::coord_cartesian(ylim = c(0, case_upper), clip = "off") +
      ggplot2::labs(title = title_text, subtitle = NULL, caption = NULL) +
      theme_pub(base_size = PUB_BASE - 0.5, base_family = base_family_global) +
      ggplot2::theme(
        legend.position = "bottom", legend.direction = "horizontal",
        # Environmental keys wrap horizontally; guide blocks are stacked so
        # the legend uses the available width instead of producing long gaps.
        legend.box = "horizontal", legend.justification = "center",
        legend.box.just = "center",
        legend.text = ggplot2::element_text(size = 6.1, lineheight = 0.92,
                                            margin = ggplot2::margin(l = 2, r = 6)),
        legend.key.width = grid::unit(11, "pt"),
        legend.key.height = grid::unit(8.5, "pt"),
        legend.spacing.x = grid::unit(5, "pt"),
        legend.spacing.y = grid::unit(2, "pt"),
        legend.margin = ggplot2::margin(3, 6, 3, 6),
        legend.box.spacing = grid::unit(4, "pt"),
        strip.text = ggplot2::element_text(size = 6.5, face = "bold"),
        axis.text.x = ggplot2::element_text(size = 5.8),
        axis.text.y = ggplot2::element_text(size = 5.8),
        axis.title.y.right = ggplot2::element_text(margin = ggplot2::margin(l = 7)),
        plot.subtitle = ggplot2::element_blank(), plot.caption = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(12, 24, 18, 14))
  }

  # Year-wise 17-region comparison figures: fixed paired axes so the right-axis
  # environmental units remain valid across every facet.
  for (yy in EVALUABLE_YEARS) {
    plot_y <- df_all |>
      dplyr::filter(YR == yy, REGION %in% CANONICAL_17) |>
      dplyr::mutate(REGION = factor(REGION, levels = CANONICAL_17))
    overlay_y <- regional_rf_overlay |>
      dplyr::filter(YR == yy) |>
      dplyr::mutate(REGION = factor(REGION, levels = CANONICAL_17))
    if (nrow(plot_y) == 0L) next

    p_rf_overlay <- build_regional_rf_plot(
      plot_y, overlay_y,
      paste0(rf_scale_label, " preceding ", rf_anchor_label, ": regional overlay, ", yy),
      ggplot2::facet_wrap(~ REGION, ncol = 4, scales = "fixed", drop = FALSE))
    save_plot_file(paste0("FigureB6_RF_Trigger_Overlay_Regional_", yy),
                   p_rf_overlay, NC_W_DOUBLE, 9.2, max_height = NC_H_SUPP)
  }

  # Region-wise multi-year figures preserve a valid region-specific secondary
  # axis and make low-incidence regions readable without losing any of the 17.
  for (rr in CANONICAL_17) {
    plot_r <- df_all |>
      dplyr::filter(REGION == rr, YR %in% EVALUABLE_YEARS) |>
      dplyr::mutate(YR = factor(YR, levels = EVALUABLE_YEARS))
    overlay_r <- regional_rf_overlay |>
      dplyr::filter(REGION == rr, YR %in% EVALUABLE_YEARS) |>
      dplyr::mutate(YR = factor(YR, levels = EVALUABLE_YEARS))
    if (nrow(plot_r) == 0L) next
    safe_rr <- gsub("[^A-Za-z0-9]+", "_", rr)
    p_rf_region <- build_regional_rf_plot(
      plot_r, overlay_r,
      paste0(rf_scale_label, " preceding ", rf_anchor_label, ": ", rr),
      ggplot2::facet_wrap(~ YR, ncol = 2, scales = "fixed", drop = FALSE))
    save_plot_file(paste0("FigureB6A_RF_Overlay_", safe_rr),
                   p_rf_region, NC_W_DOUBLE, 8.8, max_height = NC_H_SUPP)
  }
}

# ------------------------------------------------------------------------------
# 6. METHOD MAPPING
# ------------------------------------------------------------------------------
# The regional benchmark uses the same six requested detectors as NCR:
# two transmission-acceleration detectors, two mean/SD historical thresholds,
# and two historical percentile comparators.
surge_defs <- c(
  "Continuous Transmission Acceleration" = "surge_sta_lta",
  "Constant Transmission Acceleration" = "surge_sta_lta_vaezi",
  "Outbreak Threshold" = "surge_mean_2sd",
  "Alarm Threshold" = "surge_mean_1sd",
  "WHO 75th Percentile Threshold" = "surge_who_75",
  "WHO 90th Percentile Threshold" = "surge_who_90"
)

# NAME-KEYED, not positional: a positional vector here is what produced the
# "Size 14 / Size 11" tibble error when detectors were added in the original
# Stage 3.
PARADIGM_LOOKUP <- c(
  "Continuous Transmission Acceleration" = "Acceleration Measures",
  "Constant Transmission Acceleration" = "Acceleration Measures",
  "Outbreak Threshold" = "Retrospective Thresholds",
  "Alarm Threshold" = "Retrospective Thresholds",
  "WHO 75th Percentile Threshold" = "Surveillance-Guideline Percentile Thresholds",
  "WHO 90th Percentile Threshold" = "Surveillance-Guideline Percentile Thresholds"
)
.untyped <- setdiff(names(surge_defs), names(PARADIGM_LOOKUP))
if (length(.untyped) > 0) {
  stop("Detector(s) with no paradigm: ", paste(.untyped, collapse = ", "),
       call. = FALSE)
}
method_type_map <- tibble::tibble(
  Method   = names(surge_defs),
  Paradigm = unname(PARADIGM_LOOKUP[names(surge_defs)])
)

# Canonical method display order (paradigm-grouped):
#   Acceleration Measures (2), Retrospective Thresholds (2).
method_order <- c(
  "Continuous Transmission Acceleration",
  "Constant Transmission Acceleration",
  "Outbreak Threshold",
  "Alarm Threshold",
  "WHO 75th Percentile Threshold",
  "WHO 90th Percentile Threshold"
)
stopifnot(setequal(method_order, names(surge_defs)))
stopifnot(length(method_order) == length(surge_defs))

# Multi-line wrapped detector tick labels.
method_two_line <- c(
  "Continuous Transmission Acceleration" = "Continuous\nTransmission\nAcceleration",
  "Constant Transmission Acceleration" = "Constant\nTransmission\nAcceleration",
  "Outbreak Threshold" = "Outbreak\nThreshold",
  "Alarm Threshold" = "Alarm\nThreshold",
  "WHO 75th Percentile Threshold" = "WHO 75th\nPercentile Threshold",
  "WHO 90th Percentile Threshold" = "WHO 90th\nPercentile Threshold"
)
.missing_two_line_keys <- setdiff(method_order, names(method_two_line))
if (length(.missing_two_line_keys) > 0L) {
  stop("method_two_line is missing keys for: ",
       paste(.missing_two_line_keys, collapse = ", "),
       ". Update the named vector to include all method_order entries.")
}

# ------------------------------------------------------------------------------
# 7. ANCHORS AND TRIGGER CLASSIFICATION
# ------------------------------------------------------------------------------
peak_index_whichmax <- function(dc_vec) {
  if (length(dc_vec) == 0) return(NA_integer_)
  valid_dc <- ifelse(is.na(dc_vec), -Inf, dc_vec)
  if (all(!is.finite(valid_dc))) return(NA_integer_)
  pk <- which.max(valid_dc)
  if (length(pk) == 0 || !is.finite(valid_dc[pk])) return(NA_integer_)
  as.integer(pk)
}

compute_A1 <- function(df_in, year,
                       lead_min = A1_LEAD_MIN, lead_max = A1_LEAD_MAX) {
  df_y <- df_in |> dplyr::filter(YR == year)
  if (nrow(df_y) == 0 || all(is.na(df_y$DC_QC)))
    return(list(peak_week = NA_integer_, A1_weeks = integer(0)))
  peak_idx  <- which.max(df_y$DC_QC)
  peak_week <- df_y$WN[peak_idx]
  A1_weeks  <- (peak_week - lead_max):(peak_week - lead_min)
  A1_weeks  <- A1_weeks[A1_weeks >= 1L]
  list(peak_week = as.integer(peak_week),
       A1_weeks  = as.integer(A1_weeks))
}

compute_A2 <- function(df_in, year, burden_frac = A2_BURDEN_FRAC) {
  df_y <- df_in |> dplyr::filter(YR == year) |> dplyr::arrange(WN)
  if (nrow(df_y) == 0 || all(is.na(df_y$DC_QC)))
    return(list(start_week = NA_integer_, end_week = NA_integer_,
                A2_weeks = integer(0)))
  cases <- ifelse(is.na(df_y$DC_QC), 0, df_y$DC_QC)
  weeks <- df_y$WN; n_w <- length(weeks); total <- sum(cases)
  if (total <= 0)
    return(list(start_week = NA_integer_, end_week = NA_integer_,
                A2_weeks = integer(0)))
  peak_idx <- which.max(cases); lo <- hi <- peak_idx; S <- cases[peak_idx]
  while (S / total < burden_frac && (lo > 1L || hi < n_w)) {
    L <- if (lo > 1L)  cases[lo - 1L] else -Inf
    R <- if (hi < n_w) cases[hi + 1L] else -Inf
    if (L >= R) { lo <- lo - 1L; S <- S + L }
    else        { hi <- hi + 1L; S <- S + R }
  }
  list(start_week = as.integer(weeks[lo]),
       end_week   = as.integer(weeks[hi]),
       A2_weeks   = as.integer(weeks[lo:hi]))
}

compute_anchors_for_year <- function(df_in, year,
                                     lead_min = A1_LEAD_MIN,
                                     lead_max = A1_LEAD_MAX,
                                     burden_frac = A2_BURDEN_FRAC) {
  A1 <- compute_A1(df_in, year, lead_min, lead_max)
  A2 <- compute_A2(df_in, year, burden_frac)
  list(year = year, peak_week = A1$peak_week,
       A1_weeks = A1$A1_weeks, A2_weeks = A2$A2_weeks)
}

# True alarm iff t in A1 OR t in A2 (T = A1 union A2).
classify_trigger <- function(trigger_week, A1_weeks, A2_weeks) {
  in_A1 <- trigger_week %in% A1_weeks
  in_A2 <- trigger_week %in% A2_weeks
  is_true <- in_A1 || in_A2
  list(is_true = is_true, in_A1 = in_A1, in_A2 = in_A2)
}

# Two-bucket compartment classifier on True alarms.
classify_compartment <- function(lead_time) {
  if (is.na(lead_time)) return(NA_character_)
  if (lead_time >= COMPARTMENT_ACTIONABLE_MIN &&
      lead_time <= COMPARTMENT_ACTIONABLE_MAX) return("Actionable")
  return("Reactive")
}

# ------------------------------------------------------------------------------
# 8. REGIONAL INCLUSION CRITERIA
# ------------------------------------------------------------------------------
coverage_per_region_year <- df_all |>
  dplyr::filter(YR %in% EVALUABLE_YEARS) |>
  dplyr::group_by(REGION, YR) |>
  dplyr::summarise(
    n_case_weeks = sum(is.finite(LC_DOH)),
    annual_peak = ifelse(any(is.finite(LC_DOH)), max(LC_DOH, na.rm = TRUE), NA_real_),
    .groups = "drop"
  )

evaluable_region_year <- coverage_per_region_year |>
  dplyr::filter(n_case_weeks >= MIN_OBSERVED_CASE_WEEKS_PER_YEAR) |>
  dplyr::select(REGION, YR)

region_year_count <- evaluable_region_year |>
  dplyr::count(REGION, name = "n_evaluable_years")

regions_in <- region_year_count |>
  dplyr::filter(n_evaluable_years >= MIN_EVALUABLE_YEARS_PER_REGION) |>
  dplyr::pull(REGION) |>
  as.character()

regions_excluded <- setdiff(CANONICAL_17, regions_in)
if (length(regions_excluded) > 0L) {
  stop("Regional coverage validation failed: the following canonical region(s) lack the required surveillance coverage: ",
       paste(regions_excluded, collapse = ", "), call. = FALSE)
}
if (!setequal(regions_in, CANONICAL_17) || length(unique(regions_in)) != 17L) {
  stop("Regional analysis must retain exactly the 17 canonical Philippine regions.", call. = FALSE)
}
regions_in <- CANONICAL_17

cat("Included regions: 17/17; coverage rule = >= ",
    MIN_OBSERVED_CASE_WEEKS_PER_YEAR, " observed case weeks/year in >= ",
    MIN_EVALUABLE_YEARS_PER_REGION, " evaluable years.\n", sep = "")

# ------------------------------------------------------------------------------
# 9. PER-REGION FRAMEWORK METRICS
# ------------------------------------------------------------------------------
build_trigger_detail_for_region <- function(df_region, evaluable_years_for_region) {
  rows <- list()
  for (yr in evaluable_years_for_region) {
    anchors <- compute_anchors_for_year(df_region, yr)
    df_y <- df_region |> dplyr::filter(YR == yr) |> dplyr::arrange(WN)
    if (nrow(df_y) == 0) next

    for (method_name in names(surge_defs)) {
      col_name <- surge_defs[[method_name]]
      trig_idx <- which(df_y[[col_name]] == 1L)
      if (length(trig_idx) == 0L) next
      for (k in trig_idx) {
        wk  <- df_y$WN[k]
        cls <- classify_trigger(wk, anchors$A1_weeks, anchors$A2_weeks)
        lead_time <- if (!is.na(anchors$peak_week))
          as.integer(anchors$peak_week - wk) else NA_integer_
        # Compartment is meaningful only for True alarms.
        comp <- if (isTRUE(cls$is_true)) classify_compartment(lead_time)
                else NA_character_
        case_count <- df_y$LC_DOH[k]
        rows[[length(rows) + 1L]] <- data.frame(
          Year = yr, Method = method_name, Week = as.integer(wk),
          Peak_Week = anchors$peak_week, Lead_Time = lead_time,
          Compartment = comp,
          IsTrue = cls$is_true,
          InA1 = cls$in_A1, InA2 = cls$in_A2,
          DC = if (is.na(case_count)) 0 else case_count,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (length(rows) == 0L) {
    return(data.frame(
      Year = integer(), Method = character(), Week = integer(),
      Peak_Week = integer(), Lead_Time = integer(),
      Compartment = character(), IsTrue = logical(),
      InA1 = logical(), InA2 = logical(),
      DC = numeric(), stringsAsFactors = FALSE
    ))
  }
  dplyr::bind_rows(rows)
}

compute_yearly_lead_for_region <- function(trig_aug, evaluable_years_for_region) {
  if (nrow(trig_aug) == 0)
    return(data.frame(
      Method = character(), Year = integer(),
      First_A1_True_Week = integer(), Lead_Time_Yr = numeric(),
      stringsAsFactors = FALSE
    ))
  trig_aug |>
    dplyr::filter(Year %in% evaluable_years_for_region, InA1, IsTrue) |>
    dplyr::group_by(Method, Year) |>
    dplyr::slice_min(Week, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::transmute(
      Method,
      Year = as.integer(Year),
      First_A1_True_Week = as.integer(Week),
      Lead_Time_Yr = as.numeric(Lead_Time)
    )
}

compute_method_metrics_region <- function(trig_aug, yearly_lead,
                                          evaluable_years_for_region) {
  rows <- list()
  n_eval <- length(evaluable_years_for_region)
  for (m in names(surge_defs)) {
    trig_sub <- trig_aug |>
      dplyr::filter(Method == m, Year %in% evaluable_years_for_region)
    lead_sub <- yearly_lead |>
      dplyr::filter(Method == m, Year %in% evaluable_years_for_region)

    total      <- nrow(trig_sub)
    true_n     <- sum(trig_sub$IsTrue, na.rm = TRUE)
    false_n    <- total - true_n
    reactive_n <- sum(trig_sub$Compartment == "Reactive", na.rm = TRUE)
    truact_n   <- sum(trig_sub$IsTrue & trig_sub$Compartment == "Actionable",
                      na.rm = TRUE)
    truact_lt  <- trig_sub$Lead_Time[
      trig_sub$IsTrue & trig_sub$Compartment == "Actionable"
    ]

    ppv <- if (total > 0)  true_n / total else NA_real_

    # A1-restricted Sensitivity.
    years_with_a1_true <- length(unique(
      trig_sub$Year[trig_sub$IsTrue & trig_sub$Compartment == "Actionable"]
    ))
    sens <- if (n_eval > 0) years_with_a1_true / n_eval else NA_real_

    # Mean Lead Time: same-denominator headline + conditional variant.
    mean_lead_conditional <- if (nrow(lead_sub) > 0)
      mean(lead_sub$Lead_Time_Yr, na.rm = TRUE) else NA_real_
    if (is.nan(mean_lead_conditional)) mean_lead_conditional <- NA_real_
    mean_lead <- if (n_eval > 0)
      sum(lead_sub$Lead_Time_Yr, na.rm = TRUE) / n_eval else NA_real_
    if (is.nan(mean_lead)) mean_lead <- NA_real_

    # Warning Persistence: same-denominator headline + conditional variant.
    wp_conditional <- if (length(truact_lt) > 0)
      mean(truact_lt, na.rm = TRUE) else NA_real_
    if (is.nan(wp_conditional)) wp_conditional <- NA_real_
    if (n_eval > 0L) {
      per_year_wp <- vapply(evaluable_years_for_region, function(y) {
        lt <- trig_sub$Lead_Time[
          trig_sub$Year == y &
            trig_sub$IsTrue &
            trig_sub$Compartment == "Actionable"
        ]
        if (length(lt) > 0L) mean(lt, na.rm = TRUE) else 0
      }, numeric(1))
      wp <- mean(per_year_wp, na.rm = TRUE)
      if (is.nan(wp)) wp <- NA_real_
    } else {
      wp <- NA_real_
    }

    aly       <- if (true_n > 0) truact_n / true_n else NA_real_
    n_true_yr <- if (n_eval > 0) true_n  / n_eval else NA_real_
    n_false_yr <- if (n_eval > 0) false_n / n_eval else NA_real_

    # TAM (True-Alarm Magnitude): per-year T-restricted alarm-on case sum,
    # averaged across evaluable years.
    tam_per_year <- vapply(evaluable_years_for_region, function(y) {
      ix <- which(trig_sub$Year == y & trig_sub$IsTrue)
      if (length(ix) == 0) 0 else sum(trig_sub$DC[ix], na.rm = TRUE)
    }, numeric(1))
    tam <- mean(tam_per_year, na.rm = TRUE)

    rows[[m]] <- data.frame(
      Method = m,
      # Reported metrics
      TAM            = tam,
      N_True_Alarms  = n_true_yr,
      Sensitivity    = sens,
      Mean_Lead_Time = mean_lead,
      WP             = wp,
      # Diagnostic / wide-CSV-only columns
      Mean_Lead_Time_conditional = mean_lead_conditional,
      WP_conditional             = wp_conditional,
      PPV            = ppv,
      ALY            = aly,
      N_False_Alarms = n_false_yr,
      Total_Triggers = total,
      True_Alarms    = true_n,
      False_Alarms   = false_n,
      n_Reactive            = reactive_n,
      n_TrueActionable      = truact_n,
      N_Years_with_A1_True  = years_with_a1_true,
      N_Years_Evaluable     = n_eval,
      stringsAsFactors = FALSE
    )
  }
  dplyr::bind_rows(rows)
}

# Dedicated RF-to-anchor inventory for RF product reruns. Unlike the six-detector
# diagnostic outputs, this table evaluates only the disease anchor that the RF
# product is designed to precede: P1/P3 -> Constant TA; P2/P4 -> Outbreak
# Threshold. A1/A2 remains the independent true/false-alarm standard.
RF_TO_ANCHOR_DETAIL <- data.frame()
RF_TO_ANCHOR_SUMMARY <- data.frame()
if (isTRUE(RF_GATE$active)) {
  rf_anchor_col_detail <- if (RF_GATE$product %in% c(1L, 3L))
    "ungated_surge_sta_lta_vaezi" else "ungated_surge_mean_2sd"
  rf_anchor_name_detail <- if (RF_GATE$product %in% c(1L, 3L))
    "Constant TA" else "Outbreak Threshold"
  rf_representation_detail <- if (RF_GATE$scale == "mm") "RF mm" else "RF STA/LTA R(t)"

  detail_rows <- list()
  for (reg in regions_in) {
    yrs_reg <- evaluable_region_year |> dplyr::filter(REGION == reg) |> dplyr::pull(YR) |> sort()
    dr <- df_all |> dplyr::filter(REGION == reg)
    for (yr in yrs_reg) {
      dy <- dr |> dplyr::filter(YR == yr) |> dplyr::arrange(WN)
      if (!nrow(dy)) next
      anchors <- compute_anchors_for_year(dr, yr)
      idx <- which(as.integer(dy[[rf_anchor_col_detail]]) == 1L)
      if (!length(idx)) next
      for (ii in idx) {
        wk <- as.integer(dy$WN[ii])
        cls <- classify_trigger(wk, anchors$A1_weeks, anchors$A2_weeks)
        detail_rows[[length(detail_rows) + 1L]] <- data.frame(
          REGION = reg, YR = yr, WN = wk, Date = dy$Date[ii],
          RF_Product = RF_GATE$product,
          RF_Representation = rf_representation_detail,
          Disease_Anchor = rf_anchor_name_detail,
          RF_Feature = suppressWarnings(as.numeric(dy$RF_Gate_Feature[ii])),
          RF_Threshold_Lower = RF_GATE$lower,
          RF_Threshold_Point = RF_GATE$threshold,
          RF_Threshold_Upper = RF_GATE$upper,
          RF_Activation_Status = as.character(dy$RF_Activation_Status[ii]),
          RF_Supported = isTRUE(dy$RF_gate_ok[ii]),
          InA1 = isTRUE(cls$in_A1), InA2 = isTRUE(cls$in_A2),
          IsTrue_A1A2 = isTRUE(cls$is_true),
          Alarm_Class = if (isTRUE(cls$is_true)) "True alarm" else "False alarm",
          Peak_Week = anchors$peak_week,
          Lead_Time = if (is.finite(anchors$peak_week)) anchors$peak_week - wk else NA_integer_,
          stringsAsFactors = FALSE)
      }
    }
  }
  RF_TO_ANCHOR_DETAIL <- dplyr::bind_rows(detail_rows)
  if (!nrow(RF_TO_ANCHOR_DETAIL)) {
    stop("RF product has no disease-anchor trigger weeks to map: Product ", RF_GATE$product,
         " -> ", rf_anchor_name_detail, call. = FALSE)
  }
  RF_TO_ANCHOR_SUMMARY <- RF_TO_ANCHOR_DETAIL |>
    dplyr::group_by(REGION) |>
    dplyr::summarise(
      RF_Product = dplyr::first(RF_Product),
      RF_Representation = dplyr::first(RF_Representation),
      Disease_Anchor = dplyr::first(Disease_Anchor),
      N_Anchor_Trigger_Weeks = dplyr::n(),
      N_True_A1A2_Anchor_Weeks = sum(IsTrue_A1A2, na.rm = TRUE),
      N_False_A1A2_Anchor_Weeks = sum(!IsTrue_A1A2, na.rm = TRUE),
      N_RF_Supported_Anchor_Weeks = sum(RF_Supported, na.rm = TRUE),
      N_RF_Supported_True_Alarms = sum(RF_Supported & IsTrue_A1A2, na.rm = TRUE),
      N_RF_Supported_False_Alarms = sum(RF_Supported & !IsTrue_A1A2, na.rm = TRUE),
      RF_Supported_PPV = ifelse(N_RF_Supported_Anchor_Weeks > 0,
                                N_RF_Supported_True_Alarms / N_RF_Supported_Anchor_Weeks, NA_real_),
      RF_Qualification_Rate = N_RF_Supported_Anchor_Weeks / N_Anchor_Trigger_Weeks,
      N_Years_With_RF_Supported_True_Alarm = dplyr::n_distinct(YR[RF_Supported & IsTrue_A1A2]),
      .groups = "drop")
  .rf_anchor_required_cols <- c(
    "REGION", "RF_Product", "RF_Representation", "Disease_Anchor",
    "N_Anchor_Trigger_Weeks", "N_RF_Supported_Anchor_Weeks",
    "N_RF_Supported_True_Alarms", "N_RF_Supported_False_Alarms",
    "RF_Supported_PPV", "RF_Qualification_Rate"
  )
  .rf_anchor_missing_cols <- setdiff(.rf_anchor_required_cols, names(RF_TO_ANCHOR_SUMMARY))
  if (length(.rf_anchor_missing_cols)) {
    stop("RF-to-anchor summary schema is incomplete. Missing: ",
         paste(.rf_anchor_missing_cols, collapse = ", "), call. = FALSE)
  }
  .expected_rf_anchor <- if (RF_GATE$product %in% c(1L, 3L)) "Constant TA" else "Outbreak Threshold"
  if (any(as.integer(RF_TO_ANCHOR_SUMMARY$RF_Product) != RF_GATE$product) ||
      any(as.character(RF_TO_ANCHOR_SUMMARY$Disease_Anchor) != .expected_rf_anchor)) {
    stop("RF-to-anchor product mapping failed: P1/P3 must map to Constant TA and P2/P4 to Outbreak Threshold.",
         call. = FALSE)
  }
  .ppv_ok <- is.na(RF_TO_ANCHOR_SUMMARY$RF_Supported_PPV) |
    (is.finite(RF_TO_ANCHOR_SUMMARY$RF_Supported_PPV) &
       RF_TO_ANCHOR_SUMMARY$RF_Supported_PPV >= 0 & RF_TO_ANCHOR_SUMMARY$RF_Supported_PPV <= 1)
  if (any(!.ppv_ok)) {
    stop("RF-to-anchor summary contains an invalid RF-supported PPV outside [0,1].", call. = FALSE)
  }

  if (!setequal(as.character(RF_TO_ANCHOR_SUMMARY$REGION), CANONICAL_17)) {
    missing <- setdiff(CANONICAL_17, as.character(RF_TO_ANCHOR_SUMMARY$REGION))
    # Regions with an anchor detector but no active weeks are explicitly added
    # rather than silently disappearing from the RF map.
    if (length(missing)) {
      RF_TO_ANCHOR_SUMMARY <- dplyr::bind_rows(
        RF_TO_ANCHOR_SUMMARY,
        data.frame(REGION = missing, RF_Product = RF_GATE$product,
                   RF_Representation = rf_representation_detail,
                   Disease_Anchor = rf_anchor_name_detail,
                   N_Anchor_Trigger_Weeks = 0L, N_True_A1A2_Anchor_Weeks = 0L,
                   N_False_A1A2_Anchor_Weeks = 0L, N_RF_Supported_Anchor_Weeks = 0L,
                   N_RF_Supported_True_Alarms = 0L, N_RF_Supported_False_Alarms = 0L,
                   RF_Supported_PPV = NA_real_, RF_Qualification_Rate = 0,
                   N_Years_With_RF_Supported_True_Alarm = 0L,
                   stringsAsFactors = FALSE))
    }
  }
  RF_TO_ANCHOR_SUMMARY <- RF_TO_ANCHOR_SUMMARY |>
    dplyr::mutate(REGION = factor(REGION, levels = CANONICAL_17)) |>
    dplyr::arrange(REGION) |>
    dplyr::mutate(REGION = as.character(REGION))
  readr::write_csv(RF_TO_ANCHOR_DETAIL,
                   file.path(OUTPUT_DIR, "RF_to_Anchor_Trigger_Detail.csv"), na = "")
  readr::write_csv(RF_TO_ANCHOR_SUMMARY,
                   file.path(OUTPUT_DIR, "RF_to_Anchor_Regional_Summary.csv"), na = "")
}

# Per-region pipeline driver
regional_metrics_list <- list()
for (reg in regions_in) {
  yrs_reg <- evaluable_region_year |>
    dplyr::filter(REGION == reg) |>
    dplyr::pull(YR) |> sort()
  if (length(yrs_reg) < MIN_EVALUABLE_YEARS_PER_REGION) next
  df_reg     <- df_all |> dplyr::filter(REGION == reg)
  trig_aug   <- build_trigger_detail_for_region(df_reg, yrs_reg)
  yearly_ld  <- compute_yearly_lead_for_region(trig_aug, yrs_reg)
  metrics    <- compute_method_metrics_region(trig_aug, yearly_ld, yrs_reg)
  metrics$REGION <- reg
  regional_metrics_list[[reg]] <- metrics
}

regional_metrics <- dplyr::bind_rows(regional_metrics_list) |>
  dplyr::left_join(method_type_map, by = "Method") |>
  dplyr::mutate(
    Method   = factor(Method,   levels = method_order),
    Paradigm = factor(Paradigm, levels = type_order)
  )

if (nrow(regional_metrics) == 0)
  stop("No regional metrics were computed. Check inclusion criteria and data.")

# ------------------------------------------------------------------------------
# 10. YEAR-CLUSTER BOOTSTRAP (Cameron-Gelbach-Miller, B = 1000)
# ------------------------------------------------------------------------------
# For each region, resample evaluable years with replacement (year is the
# cluster unit because within-year weekly observations are not independent).
# On each replicate, recompute all framework metrics from the cached trigger
# detail subset. Trigger columns and anchors are deterministic given case
# data and are cached outside the bootstrap loop for speed.
# ------------------------------------------------------------------------------
cat("\n=== Year-cluster bootstrap (B =", BOOT_N_CI, "per region) ===\n")
cat("This computation may take a few minutes for ",
    length(regions_in), " regions x ", BOOT_N_CI, " replicates.\n", sep = "")

# Anchor cache (one entry per region/year).
anchor_cache <- list()
for (reg in regions_in) {
  yrs_reg <- evaluable_region_year |>
    dplyr::filter(REGION == reg) |>
    dplyr::pull(YR) |> sort()
  if (length(yrs_reg) < MIN_EVALUABLE_YEARS_PER_REGION) next
  df_reg <- df_all |> dplyr::filter(REGION == reg)
  anchor_cache[[reg]] <- list()
  for (yr in yrs_reg) {
    anchor_cache[[reg]][[as.character(yr)]] <- compute_anchors_for_year(df_reg, yr)
  }
}

# Bootstrap-specific trigger detail builder using the anchor cache.
build_trigger_detail_cached <- function(df_region, evaluable_years_for_region,
                                        region_name) {
  rows <- list()
  for (yr in evaluable_years_for_region) {
    anchors <- anchor_cache[[region_name]][[as.character(yr)]]
    if (is.null(anchors)) next
    df_y <- df_region |> dplyr::filter(YR == yr) |> dplyr::arrange(WN)
    if (nrow(df_y) == 0) next
    for (method_name in names(surge_defs)) {
      col_name <- surge_defs[[method_name]]
      trig_idx <- which(df_y[[col_name]] == 1L)
      if (length(trig_idx) == 0L) next
      for (k in trig_idx) {
        wk  <- df_y$WN[k]
        cls <- classify_trigger(wk, anchors$A1_weeks, anchors$A2_weeks)
        lead_time <- if (!is.na(anchors$peak_week))
          as.integer(anchors$peak_week - wk) else NA_integer_
        comp <- if (isTRUE(cls$is_true)) classify_compartment(lead_time)
                else NA_character_
        case_count <- df_y$LC_DOH[k]
        rows[[length(rows) + 1L]] <- data.frame(
          Year = yr, Method = method_name, Week = as.integer(wk),
          Peak_Week = anchors$peak_week, Lead_Time = lead_time,
          Compartment = comp,
          IsTrue = cls$is_true,
          InA1 = cls$in_A1, InA2 = cls$in_A2,
          DC = if (is.na(case_count)) 0 else case_count,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (length(rows) == 0L) {
    return(data.frame(
      Year = integer(), Method = character(), Week = integer(),
      Peak_Week = integer(), Lead_Time = integer(),
      Compartment = character(), IsTrue = logical(),
      InA1 = logical(), InA2 = logical(),
      DC = numeric(), stringsAsFactors = FALSE
    ))
  }
  dplyr::bind_rows(rows)
}

# Pre-compute the FULL (unsampled) trigger detail per region.
full_trig_detail <- list()
for (reg in names(anchor_cache)) {
  df_reg  <- df_all |> dplyr::filter(REGION == reg)
  yrs_reg <- as.integer(names(anchor_cache[[reg]]))
  full_trig_detail[[reg]] <- build_trigger_detail_cached(df_reg, yrs_reg, reg)
}

# Single-replicate metric computation.
bootstrap_metrics_one_replicate <- function(reg, yrs_resampled, trig_full_override = NULL) {
  trig_full <- if (is.null(trig_full_override)) full_trig_detail[[reg]] else trig_full_override
  # RF-gated products can legitimately have zero triggers in an entire region.
  # A zero-trigger bootstrap replicate is still an evaluable replicate and must
  # return the full method schema; returning NULL drops Method and causes the
  # downstream group_by(Method) failure.
  if (is.null(trig_full)) {
    trig_full <- data.frame(
      Year = integer(), Method = character(), Week = integer(),
      Peak_Week = integer(), Lead_Time = integer(),
      Compartment = character(), IsTrue = logical(),
      InA1 = logical(), InA2 = logical(), DC = numeric(),
      stringsAsFactors = FALSE
    )
  }

  trig_sub <- dplyr::bind_rows(
    lapply(yrs_resampled, function(y) trig_full |> dplyr::filter(Year == y))
  )

  if (nrow(trig_sub) == 0L) {
    yearly_ld <- data.frame(
      Method = character(), Year = integer(),
      Lead_Time_Yr = numeric(), stringsAsFactors = FALSE
    )
  } else {
    yearly_ld <- trig_sub |>
      dplyr::filter(InA1, IsTrue) |>
      dplyr::group_by(Method, Year) |>
      dplyr::slice_min(Week, n = 1, with_ties = FALSE) |>
      dplyr::ungroup() |>
      dplyr::transmute(Method, Year = as.integer(Year),
                       Lead_Time_Yr = as.numeric(Lead_Time))
    yr_mult <- table(yrs_resampled)
    yearly_ld <- yearly_ld |>
      dplyr::mutate(.mult = as.integer(unname(yr_mult[as.character(Year)]))) |>
      dplyr::filter(!is.na(.mult), .mult > 0L) |>
      tidyr::uncount(.mult)
  }

  rows <- list()
  n_eval <- length(yrs_resampled)
  for (m in names(surge_defs)) {
    trig_m <- trig_sub |> dplyr::filter(Method == m)
    lead_m <- yearly_ld |> dplyr::filter(Method == m)

    total      <- nrow(trig_m)
    true_n     <- sum(trig_m$IsTrue, na.rm = TRUE)
    false_n    <- total - true_n
    reactive_n <- sum(trig_m$Compartment == "Reactive", na.rm = TRUE)
    truact_n   <- sum(trig_m$IsTrue & trig_m$Compartment == "Actionable",
                      na.rm = TRUE)
    truact_lt  <- trig_m$Lead_Time[
      trig_m$IsTrue & trig_m$Compartment == "Actionable"
    ]

    ppv <- if (total > 0)  true_n / total else NA_real_
    years_with_a1_true <- if (n_eval > 0L) {
      sum(vapply(yrs_resampled, function(y) {
        any(trig_m$Year == y & trig_m$IsTrue &
              trig_m$Compartment == "Actionable", na.rm = TRUE)
      }, logical(1)))
    } else 0L
    sens <- if (n_eval > 0) years_with_a1_true / n_eval else NA_real_

    mean_lead <- if (n_eval > 0)
      sum(lead_m$Lead_Time_Yr, na.rm = TRUE) / n_eval else NA_real_
    if (is.nan(mean_lead)) mean_lead <- NA_real_

    if (n_eval > 0L) {
      per_year_wp <- vapply(yrs_resampled, function(y) {
        lt <- trig_m$Lead_Time[
          trig_m$Year == y &
            trig_m$IsTrue &
            trig_m$Compartment == "Actionable"
        ]
        if (length(lt) > 0L) mean(lt, na.rm = TRUE) else 0
      }, numeric(1))
      wp <- mean(per_year_wp, na.rm = TRUE)
      if (is.nan(wp)) wp <- NA_real_
    } else {
      wp <- NA_real_
    }

    aly       <- if (true_n > 0) truact_n / true_n else NA_real_
    n_true_yr <- if (n_eval > 0) true_n  / n_eval else NA_real_
    n_false_yr <- if (n_eval > 0) false_n / n_eval else NA_real_

    # Cluster-bootstrap TAM: per-year sums computed once on the unduplicated
    # trigger detail, then averaged across the resampled year vector.
    trig_m_unique <- trig_full |> dplyr::filter(Method == m)
    yr_to_sum <- if (sum(trig_m_unique$IsTrue, na.rm = TRUE) > 0L) {
      tapply(
        trig_m_unique$DC[trig_m_unique$IsTrue],
        trig_m_unique$Year[trig_m_unique$IsTrue],
        sum, na.rm = TRUE
      )
    } else {
      stats::setNames(numeric(0), character(0))
    }
    tam_per_year <- vapply(yrs_resampled, function(y) {
      v <- yr_to_sum[as.character(y)]
      if (is.null(v) || length(v) == 0L || is.na(v)) 0 else as.numeric(v)
    }, numeric(1))
    tam <- mean(tam_per_year, na.rm = TRUE)

    rows[[m]] <- data.frame(
      Method = m,
      TAM = tam, N_True_Alarms = n_true_yr,
      Sensitivity = sens,
      Mean_Lead_Time = mean_lead, WP = wp,
      PPV = ppv, ALY = aly, N_False_Alarms = n_false_yr,
      stringsAsFactors = FALSE
    )
  }
  dplyr::bind_rows(rows)
}

# Zero-trigger bootstrap regression test ------------------------------------------------
# RF-gated reruns can legitimately leave a region with no qualifying trigger.
# Such a replicate must still return one row per detector with the full schema.
.bootstrap_required_cols <- c(
  "Method", "TAM", "N_True_Alarms", "Sensitivity", "Mean_Lead_Time",
  "WP", "PPV", "ALY", "N_False_Alarms"
)
.zero_trigger_detail <- data.frame(
  Year = integer(), Method = character(), Week = integer(),
  Peak_Week = integer(), Lead_Time = integer(),
  Compartment = character(), IsTrue = logical(),
  InA1 = logical(), InA2 = logical(), DC = numeric(),
  stringsAsFactors = FALSE
)
.test_years <- if (length(names(anchor_cache))) {
  as.integer(names(anchor_cache[[names(anchor_cache)[1L]]]))
} else {
  EVALUABLE_YEARS
}
if (!length(.test_years)) .test_years <- 2019L
.zero_boot <- bootstrap_metrics_one_replicate(
  reg = if (length(names(anchor_cache))) names(anchor_cache)[1L] else "NCR",
  yrs_resampled = .test_years,
  trig_full_override = .zero_trigger_detail
)
if (!all(.bootstrap_required_cols %in% names(.zero_boot)) ||
    nrow(.zero_boot) != length(surge_defs) ||
    !setequal(as.character(.zero_boot$Method), names(surge_defs))) {
  stop(
    "Regional bootstrap zero-trigger regression failed: expected a complete ",
    length(surge_defs), "-method schema including Method. This protects RF Products 1-4 ",
    "from empty-replicate group_by(Method) failures.",
    call. = FALSE
  )
}
rm(.zero_trigger_detail, .zero_boot, .test_years)
cat("Regional bootstrap zero-trigger regression: PASS (complete six-method schema).\n")

# Bootstrap loop
boot_results <- list()
set.seed(20260101L)
for (reg in names(anchor_cache)) {
  yrs_reg <- as.integer(names(anchor_cache[[reg]]))
  if (length(yrs_reg) < 2L) {
    cat("  Skipping region '", reg, "' (only ",
        length(yrs_reg), " evaluable year(s)).\n", sep = "")
    next
  }
  reps <- vector("list", BOOT_N_CI)
  for (b in seq_len(BOOT_N_CI)) {
    yrs_resampled <- sample(yrs_reg, size = length(yrs_reg), replace = TRUE)
    reps[[b]] <- bootstrap_metrics_one_replicate(reg, yrs_resampled)
  }
  reps_df <- dplyr::bind_rows(reps, .id = ".rep")
  missing_boot_cols <- setdiff(c(".rep", .bootstrap_required_cols), names(reps_df))
  expected_boot_rows <- BOOT_N_CI * length(surge_defs)
  if (length(missing_boot_cols) || nrow(reps_df) != expected_boot_rows) {
    stop(
      "Regional bootstrap schema failure for ", reg, ": missing column(s) [",
      paste(missing_boot_cols, collapse = ", "), "]; expected ", expected_boot_rows,
      " replicate-method rows but found ", nrow(reps_df), ".",
      call. = FALSE
    )
  }
  .rep_counts <- table(reps_df$.rep)
  if (length(.rep_counts) != BOOT_N_CI || any(as.integer(.rep_counts) != length(surge_defs))) {
    stop("Regional bootstrap replicate completeness failure for ", reg,
         ": every replicate must contain exactly one row for each of the ",
         length(surge_defs), " methods.", call. = FALSE)
  }

  ci_summary <- reps_df |>
    dplyr::group_by(Method) |>
    dplyr::summarise(
      TAM_lo            = safe_quantile(TAM, 0.025),
      TAM_hi            = safe_quantile(TAM, 0.975),
      N_True_Alarms_lo  = safe_quantile(N_True_Alarms, 0.025),
      N_True_Alarms_hi  = safe_quantile(N_True_Alarms, 0.975),
      PPV_lo            = safe_quantile(PPV, 0.025),
      PPV_hi            = safe_quantile(PPV, 0.975),
      Sens_lo           = safe_quantile(Sensitivity, 0.025),
      Sens_hi           = safe_quantile(Sensitivity, 0.975),
      MLT_lo            = safe_quantile(Mean_Lead_Time, 0.025),
      MLT_hi            = safe_quantile(Mean_Lead_Time, 0.975),
      WP_lo             = safe_quantile(WP, 0.025),
      WP_hi             = safe_quantile(WP, 0.975),
      ALY_lo            = safe_quantile(ALY, 0.025),
      ALY_hi            = safe_quantile(ALY, 0.975),
      N_False_Alarms_lo = safe_quantile(N_False_Alarms, 0.025),
      N_False_Alarms_hi = safe_quantile(N_False_Alarms, 0.975),
      .groups = "drop"
    ) |>
    dplyr::mutate(REGION = reg)

  boot_results[[reg]] <- list(replicates = reps_df, ci = ci_summary)
  cat("  Region '", reg, "' bootstrap complete (",
      BOOT_N_CI, " replicates).\n", sep = "")
}

boot_ci <- dplyr::bind_rows(lapply(boot_results, function(x) x$ci))

# Merge CIs into the metrics table
regional_metrics <- regional_metrics |>
  dplyr::mutate(REGION = as.character(REGION)) |>
  dplyr::left_join(
    boot_ci |> dplyr::mutate(Method = as.character(Method)),
    by = c("REGION", "Method")
  ) |>
  dplyr::mutate(REGION = factor(REGION, levels = unique(REGION)))

# ------------------------------------------------------------------------------
# 11. DOMINANCE PROBABILITY (composite mean-rank across the 5 reported metrics)
# ------------------------------------------------------------------------------
# For each replicate, rank the six requested target detectors on each of the five
# reported metrics (rank 1 = best; ties get average rank; all metrics are
# higher-is-better), then compute the mean rank across the five metrics.
# The replicate winner is the detector with the lowest mean rank. The
# Dominance_Probability for a region is the proportion of replicates in
# which the observed point-estimate winner wins on the same composite.
#
# Properties:
#   - Equal weight to each of the five reported metrics (Borda-style).
#   - Insensitive to absolute scale across metrics.
#   - Ties handled by ties.method = "average".
#   - PPV is computed as a diagnostic but is NOT used in this composite.
# ------------------------------------------------------------------------------
COMPOSITE_DOMINANCE_METRICS <- c(
  "TAM", "N_True_Alarms", "Sensitivity", "Mean_Lead_Time", "WP"
)

.composite_winner_per_rep <- function(reps_target,
                                      composite_metrics = COMPOSITE_DOMINANCE_METRICS) {
  ok <- Reduce(`&`, lapply(composite_metrics,
                            function(m) is.finite(reps_target[[m]])))
  reps_target <- reps_target[ok, , drop = FALSE]
  if (nrow(reps_target) == 0L) return(reps_target[0, , drop = FALSE])
  reps_target <- reps_target |>
    dplyr::group_by(.rep) |>
    dplyr::mutate(
      r_TAM  = rank(-TAM,            ties.method = "average"),
      r_NTA  = rank(-N_True_Alarms,  ties.method = "average"),
      r_Sens = rank(-Sensitivity,    ties.method = "average"),
      r_MLT  = rank(-Mean_Lead_Time, ties.method = "average"),
      r_WP   = rank(-WP,             ties.method = "average"),
      mean_rank = (r_TAM + r_NTA + r_Sens + r_MLT + r_WP) / 5
    ) |>
    dplyr::slice_min(mean_rank, n = 1, with_ties = FALSE) |>
    dplyr::ungroup()
  reps_target
}

dominance_results <- list()
# All six requested detectors participate in regional dominance probability.
TARGET_DETECTORS <- c(
  "Constant Transmission Acceleration",
  "Continuous Transmission Acceleration",
  "Outbreak Threshold",
  "Alarm Threshold",
  "WHO 75th Percentile Threshold",
  "WHO 90th Percentile Threshold"
)
TARGET_P_COLS <- c(
  "Constant Transmission Acceleration" = "P_ConstantTA",
  "Continuous Transmission Acceleration" = "P_ContinuousTA",
  "Outbreak Threshold" = "P_OutbreakThreshold",
  "Alarm Threshold" = "P_AlarmThreshold",
  "WHO 75th Percentile Threshold" = "P_WHO75",
  "WHO 90th Percentile Threshold" = "P_WHO90"
)
for (reg in names(boot_results)) {
  reps_df     <- boot_results[[reg]]$replicates
  reps_target <- reps_df |> dplyr::filter(Method %in% TARGET_DETECTORS)

  rep_winners <- .composite_winner_per_rep(reps_target,
                                           COMPOSITE_DOMINANCE_METRICS)

  obs_pt <- regional_metrics |>
    dplyr::filter(REGION == reg, Method %in% TARGET_DETECTORS) |>
    dplyr::select(Method, dplyr::all_of(COMPOSITE_DOMINANCE_METRICS))
  if (nrow(obs_pt) >= 1L &&
      all(vapply(COMPOSITE_DOMINANCE_METRICS,
                 function(m) is.finite(obs_pt[[m]][1]),
                 logical(1)))) {
    obs_pt <- obs_pt |>
      dplyr::mutate(
        r_TAM  = rank(-TAM,            ties.method = "average"),
        r_NTA  = rank(-N_True_Alarms,  ties.method = "average"),
        r_Sens = rank(-Sensitivity,    ties.method = "average"),
        r_MLT  = rank(-Mean_Lead_Time, ties.method = "average"),
        r_WP   = rank(-WP,             ties.method = "average"),
        mean_rank = (r_TAM + r_NTA + r_Sens + r_MLT + r_WP) / 5
      ) |>
      dplyr::arrange(mean_rank)
    obs_winner <- as.character(obs_pt$Method[1])
  } else {
    obs_winner <- character(0)
  }

  win_tab    <- table(rep_winners$Method)
  total_reps <- sum(win_tab)

  dom_prob <- if (length(obs_winner) > 0L && obs_winner %in% names(win_tab)) {
    as.numeric(win_tab[obs_winner]) / total_reps
  } else {
    NA_real_
  }

  dominance_results[[reg]] <- data.frame(
    REGION = reg,
    Observed_Winner = if (length(obs_winner) > 0L) obs_winner else NA_character_,
    Dominance_Probability = dom_prob,
    Bootstrap_N    = total_reps,
    Composite_Rule = "mean_rank_5metrics",
    stringsAsFactors = FALSE
  )
  # One P_* column per target detector, added programmatically.
  for (dn in TARGET_DETECTORS) {
    dominance_results[[reg]][[unname(TARGET_P_COLS[dn])]] <-
      if (dn %in% names(win_tab)) as.numeric(win_tab[dn]) / total_reps else 0
  }
}
dominance_df <- dplyr::bind_rows(dominance_results) |>
  dplyr::mutate(dplyr::across(dplyr::starts_with("P_"),
                              ~ ifelse(is.na(.), 0, .)))

cat("\n=== Bootstrap dominance probability (composite mean-rank, 5 metrics) ===\n")
print(as.data.frame(dominance_df), row.names = FALSE)
cat("\n")


# ------------------------------------------------------------------------------
# HEAD-TO-HEAD CONSENSUS: EACH TA DETECTOR VS EPIDEMIC THRESHOLD
# ------------------------------------------------------------------------------
# REVISION: this block previously ran a single anchor (Constant TA) against
# four comparators (Outbreak Threshold, Farrington, EWARS, EARS). With those
# three contemporary comparators removed and Outbreak Threshold the only
# remaining comparator, the design is now two independent single-comparator
# head-to-head tests, exactly as instructed:
#   Constant Transmission Acceleration    vs  Outbreak Threshold
#   Continuous Transmission Acceleration  vs  Outbreak Threshold
# Each is its own scientific question (consistent with the "we do NOT
# multiply p-values across ... pairs" principle stated elsewhere in this
# script for the DETECTOR_PAIRS Wilcoxon panel), so the two families are
# Bonferroni-corrected separately rather than pooled into one k. See Item 8
# of the Delphi review.
#
# This block is CSV-only. It does not alter the FigureB4/FigureB5 objects
# built from DETECTOR_PAIRS above (which report the same two comparisons at
# the display level, via a different, region-paired Wilcoxon test). This
# block instead uses the within-region bootstrap win-count test (2x2
# contingency, one-sided) that the original Supplementary Table 13 used, so
# both a region-paired and a within-region-bootstrap view of the same two
# comparisons are available.
#
# Bonferroni:
#   strict cross-region, per family: k = n_included_regions (computed below
#                                    from the regions that actually pass the
#                                    inclusion filter, not hardcoded;
#                                    17 in this coverage-complete leptospirosis dataset
#   within-region:                   k = 1 (a single comparator per anchor,
#                                    so this equals the raw p-value; kept as
#                                    a column for interface consistency with
#                                    the original four-comparator table)
# ------------------------------------------------------------------------------

ST13_K_STRICT <- length(unique(as.character(dominance_df$REGION)))
ST13_K_WITHIN <- 1L

.st13_pair_test <- function(p_anchor, p_comp, anchor_label, B = N_BOOTS) {
  if (!is.finite(p_anchor) || !is.finite(p_comp)) {
    return(list(p_raw = NA_real_, winner = NA_character_))
  }
  n_anchor <- round(p_anchor * B)
  n_comp   <- round(p_comp   * B)
  n_pair   <- n_anchor + n_comp
  if (n_pair <= 0L || n_anchor == n_comp) {
    return(list(p_raw = if (n_pair <= 0L) NA_real_ else 1.0,
                winner = NA_character_))
  }
  if (n_anchor > n_comp) {
    p_raw <- stats::binom.test(n_anchor, n_pair, p = 0.5,
                               alternative = "greater")$p.value
    winner <- anchor_label
  } else {
    p_raw <- stats::binom.test(n_comp, n_pair, p = 0.5,
                               alternative = "greater")$p.value
    winner <- "Outbreak Threshold"
  }
  list(p_raw = p_raw, winner = winner)
}

ST13_ANCHORS <- c("Constant Transmission Acceleration",
                  "Continuous Transmission Acceleration")

st13_rows <- list()
for (anchor in ST13_ANCHORS) {
  p_col <- unname(TARGET_P_COLS[anchor])
  for (i in seq_len(nrow(dominance_df))) {
    rr  <- dominance_df[i, , drop = FALSE]
    reg <- as.character(rr$REGION)
    p_anchor <- as.numeric(rr[[p_col]])
    p_comp   <- as.numeric(rr[["P_OutbreakThreshold"]])
    tt <- .st13_pair_test(p_anchor, p_comp, anchor, B = N_BOOTS)
    st13_rows[[length(st13_rows) + 1L]] <- data.frame(
      REGION = reg,
      Anchor = anchor,
      Comparator = "Outbreak Threshold",
      P_Anchor = p_anchor,
      P_Comparator = p_comp,
      Pair_Winner = tt$winner,
      p_raw = tt$p_raw,
      p_bonf_strict = ifelse(is.na(tt$p_raw), NA_real_, min(1, tt$p_raw * ST13_K_STRICT)),
      p_bonf_within = ifelse(is.na(tt$p_raw), NA_real_, min(1, tt$p_raw * ST13_K_WITHIN)),
      stringsAsFactors = FALSE
    )
  }
}

st13_long <- dplyr::bind_rows(st13_rows) |>
  dplyr::mutate(
    Sig_Strict = !is.na(p_bonf_strict) & p_bonf_strict < HH_ALPHA,
    Sig_Within = !is.na(p_bonf_within) & p_bonf_within < HH_ALPHA,
    # Three-tier result per anchor-vs-Outbreak Threshold pair (simplified
    # from the original four-tier scheme, which needed multiple comparators
    # per anchor to have "strong" vs "partial" tiers be meaningful):
    #   significant_win : anchor wins, Bonferroni-significant
    #   lead_only       : anchor's point estimate favours it, not significant
    #   no_advantage    : Outbreak Threshold wins, or the pair is untestable
    Result_Strict = dplyr::case_when(
      is.na(p_bonf_strict) ~ "not_testable",
      Sig_Strict & Pair_Winner == Anchor ~ "significant_win",
      !is.na(Pair_Winner) & Pair_Winner == Anchor ~ "lead_only",
      TRUE ~ "no_advantage"
    ),
    Result_Within = dplyr::case_when(
      is.na(p_bonf_within) ~ "not_testable",
      Sig_Within & Pair_Winner == Anchor ~ "significant_win",
      !is.na(Pair_Winner) & Pair_Winner == Anchor ~ "lead_only",
      TRUE ~ "no_advantage"
    )
  )

st13_wide <- st13_long |>
  dplyr::select(REGION, Anchor, P_Anchor, P_Comparator, Pair_Winner,
                p_raw, p_bonf_strict, p_bonf_within,
                Result_Strict, Result_Within) |>
  tidyr::pivot_wider(
    id_cols = REGION,
    names_from = Anchor,
    values_from = c(P_Anchor, p_raw, p_bonf_strict, p_bonf_within,
                    Result_Strict, Result_Within),
    names_glue = "{.value}__{Anchor}"
  ) |>
  dplyr::mutate(REGION = factor(REGION, levels = CANONICAL_17)) |>
  dplyr::arrange(REGION) |>
  dplyr::mutate(REGION = as.character(REGION))

utils::write.csv(
  st13_long,
  file.path(OUTPUT_DIR, "StageB_Regional_TA_vs_OutbreakThreshold_Consensus_Long.csv"),
  row.names = FALSE
)
utils::write.csv(
  st13_wide,
  file.path(OUTPUT_DIR, "StageB_Regional_TA_vs_OutbreakThreshold_Consensus.csv"),
  row.names = FALSE
)

cat("Saved Constant TA / Continuous TA vs Outbreak Threshold consensus outputs:\n")
cat("  StageB_Regional_TA_vs_OutbreakThreshold_Consensus_Long.csv\n")
cat("  StageB_Regional_TA_vs_OutbreakThreshold_Consensus.csv\n")
cat("  Strict Bonferroni k (per family) = ", ST13_K_STRICT,
    "; within-family Bonferroni k = ", ST13_K_WITHIN, "\n", sep = "")

# Join dominance probability into regional_metrics (region-level property,
# inherited by every method row in that region).
regional_metrics <- regional_metrics |>
  dplyr::mutate(REGION = as.character(REGION)) |>
  dplyr::left_join(
    dominance_df |>
      dplyr::select(REGION, Observed_Winner, Dominance_Probability,
                    dplyr::all_of(unname(TARGET_P_COLS)), Bootstrap_N),
    by = "REGION"
  ) |>
  dplyr::mutate(REGION = factor(REGION, levels = unique(REGION)))

# ------------------------------------------------------------------------------
# 12. ORDER REGIONS AND METHODS; ESTABLISH METRIC METADATA
# ------------------------------------------------------------------------------
region_ppv_order <- regional_metrics |>
  dplyr::group_by(REGION) |>
  dplyr::summarise(mean_ppv = mean(PPV, na.rm = TRUE), .groups = "drop") |>
  dplyr::arrange(dplyr::desc(mean_ppv)) |>
  dplyr::pull(REGION)

regional_metrics <- regional_metrics |>
  dplyr::mutate(REGION = factor(REGION, levels = rev(region_ppv_order))) |>
  dplyr::mutate(
    Method   = factor(as.character(Method),   levels = method_order),
    Paradigm = factor(as.character(Paradigm), levels = type_order)
  )

# Five reported head-to-head metrics (Figure 4 figure set).
HH_METRICS_REGIONAL <- c(
  "TAM", "N_True_Alarms", "Sensitivity",
  "Mean_Lead_Time", "WP"
)

hh_metric_full_name <- c(
  TAM            = "True-Alarm Magnitude",
  N_True_Alarms  = "Number of True Alarms",
  Sensitivity    = "Sensitivity",
  Mean_Lead_Time = "Mean Lead Time",
  WP             = "Warning Persistence"
)

hh_metric_y_label <- c(
  TAM            = "True-Alarm Magnitude (cases per year)",
  N_True_Alarms  = "Number of True Alarms (per year)",
  Sensitivity    = "Sensitivity",
  Mean_Lead_Time = "Mean Lead Time (weeks before peak)",
  WP             = "Warning Persistence (weeks before peak)"
)

hh_metric_ci_cols <- list(
  TAM            = c("TAM_lo",            "TAM_hi"),
  N_True_Alarms  = c("N_True_Alarms_lo",  "N_True_Alarms_hi"),
  Sensitivity    = c("Sens_lo",           "Sens_hi"),
  Mean_Lead_Time = c("MLT_lo",            "MLT_hi"),
  WP             = c("WP_lo",             "WP_hi")
)

# All five reported metrics are higher-is-better.
hh_metric_higher_better <- c(
  TAM = TRUE,  N_True_Alarms = TRUE,  Sensitivity = TRUE,
  Mean_Lead_Time = TRUE, WP = TRUE
)

# Operational reference threshold (NA = no fixed threshold).
hh_metric_threshold <- c(
  TAM            = NA_real_,
  N_True_Alarms  = NA_real_,
  Sensitivity    = 1.00,
  Mean_Lead_Time = 4,
  WP             = 4
)

hh_metric_is_proportion <- c(
  TAM = FALSE, N_True_Alarms = FALSE, Sensitivity = TRUE,
  Mean_Lead_Time = FALSE, WP = FALSE
)

# Sub-letter numeric subscript used in per-metric panel headers.
hh_metric_subscript <- c(
  TAM            = "1",
  N_True_Alarms  = "2",
  Sensitivity    = "3",
  Mean_Lead_Time = "4",
  WP             = "5"
)

# Filename slug per metric (Figure 4 outputs).
hh_metric_slug <- c(
  TAM            = "TAM",
  N_True_Alarms  = "N_True_Alarms",
  Sensitivity    = "Sensitivity",
  Mean_Lead_Time = "Mean_Lead_Time",
  WP             = "Warning_Persistence"
)

# Detector-pair definitions for the Wilcoxon side panel.
# REVISION: per instruction, the head-to-head comparisons are now exactly
# Constant TA vs Outbreak Threshold and Continuous TA vs Outbreak Threshold.
# The original third display-level pair (Constant TA vs Continuous TA) is
# dropped, since it is not one of the two comparisons requested. See Item 8
# of the Delphi review.
DETECTOR_PAIRS <- list(
  list(A = "Constant Transmission Acceleration",
       B = "Outbreak Threshold",
       short = "Constant TA  vs  OT"),
  list(A = "Continuous Transmission Acceleration",
       B = "Outbreak Threshold",
       short = "Continuous TA  vs  OT")
)
HH_ALPHA_BONF <- HH_ALPHA / length(DETECTOR_PAIRS)  # reference value

# ADDITIONAL_DETECTOR_PAIRS previously extended a single Constant-TA-anchored
# comparison out to four comparators (Outbreak Threshold, Farrington, EWARS,
# EARS). With only Outbreak Threshold remaining as a comparator, and with
# Continuous TA now also getting its own head-to-head against it, this list
# is identical in content to DETECTOR_PAIRS above; it is kept as a separate
# name only so the Section 20 consensus-classification code below (which
# refers to ADDITIONAL_DETECTOR_PAIRS by name) does not need restructuring.
ADDITIONAL_DETECTOR_PAIRS <- DETECTOR_PAIRS

# ------------------------------------------------------------------------------
# 13. WILCOXON DETECTOR-PAIRED TESTS (PER METRIC)
# ------------------------------------------------------------------------------
# For each metric, two pairwise Wilcoxon signed-rank tests paired by
# REGION (n = 17), both against Outbreak Threshold: Constant TA vs Outbreak
# Threshold, and Continuous TA vs Outbreak Threshold. Significance is
# evaluated PER PAIRWISE COMPARISON at alpha = 0.05; we do NOT multiply
# p-values across the two pairs (each pair answers a distinct scientific
# question).
# ------------------------------------------------------------------------------
compute_detector_paired_wilcoxon <- function(metric_id, regional_metrics_df, detector_pairs = DETECTOR_PAIRS) {
  metric_long <- regional_metrics_df |>
    dplyr::transmute(
      REGION = as.character(REGION),
      Method = as.character(Method),
      Score  = .data[[metric_id]]
    )
  metric_wide <- metric_long |>
    tidyr::pivot_wider(names_from = Method, values_from = Score) |>
    as.data.frame()

  rows <- list()
  for (pp in detector_pairs) {
    a_vals <- if (pp$A %in% names(metric_wide)) metric_wide[[pp$A]]
              else rep(NA_real_, nrow(metric_wide))
    b_vals <- if (pp$B %in% names(metric_wide)) metric_wide[[pp$B]]
              else rep(NA_real_, nrow(metric_wide))
    paired_keep <- !is.na(a_vals) & !is.na(b_vals)
    a_paired <- a_vals[paired_keep]
    b_paired <- b_vals[paired_keep]
    n_pairs  <- length(a_paired)
    diffs    <- a_paired - b_paired
    median_diff <- if (n_pairs > 0L) stats::median(diffs, na.rm = TRUE)
                   else NA_real_
    abs_median_diff <- if (!is.na(median_diff)) abs(median_diff) else NA_real_

    if (n_pairs >= 3L && any(diffs != 0, na.rm = TRUE)) {
      wt <- suppressWarnings(stats::wilcox.test(
        a_paired, b_paired,
        paired = TRUE,
        alternative = "two.sided",
        exact = FALSE
      ))
      v_stat <- as.numeric(wt$statistic)
      p_val  <- as.numeric(wt$p.value)
    } else {
      v_stat <- NA_real_
      p_val  <- NA_real_
    }
    p_pairwise <- p_val
    p_bonf     <- p_val
    sig        <- if (!is.na(p_pairwise)) p_pairwise < HH_ALPHA else NA

    rows[[length(rows) + 1L]] <- data.frame(
      Metric = metric_id,
      Comparison = pp$short,
      Detector_A = pp$A,
      Detector_B = pp$B,
      N_Regions_Paired = n_pairs,
      Median_Diff_AminusB = median_diff,
      Abs_Median_Diff = abs_median_diff,
      V_statistic = v_stat,
      p_value = p_val,
      p_pairwise = p_pairwise,
      p_bonferroni = p_bonf,
      Significant_005 = sig,
      stringsAsFactors = FALSE
    )
  }
  dplyr::bind_rows(rows)
}

all_wilcoxon_results <- dplyr::bind_rows(
  lapply(HH_METRICS_REGIONAL, function(m) {
    compute_detector_paired_wilcoxon(m, regional_metrics)
  })
)

cat("\n=== Wilcoxon detector-paired (5 metrics x 2 pairs = 10 rows) ===\n")
print(all_wilcoxon_results, row.names = FALSE, digits = 3)
cat("\n")

# ADDITIONAL_DETECTOR_PAIRS is now identical to DETECTOR_PAIRS (see Section
# 13's definition above): with the contemporary comparators removed, there
# is nothing left to add beyond the two comparisons already computed. This
# block is kept only so the CSV filename referenced in the header docstring
# still exists; its content now duplicates all_wilcoxon_results exactly.
additional_wilcoxon_results <- dplyr::bind_rows(
  lapply(HH_METRICS_REGIONAL, function(m) {
    compute_detector_paired_wilcoxon(
      metric_id = m,
      regional_metrics_df = regional_metrics,
      detector_pairs = ADDITIONAL_DETECTOR_PAIRS
    )
  })
)

cat("\n=== Additional Wilcoxon head-to-head (duplicates the pairs above; 5 metrics x 2 pairs = 10 rows) ===\n")
print(additional_wilcoxon_results, row.names = FALSE, digits = 3)
cat("\n")

# ------------------------------------------------------------------------------
# 14. REGIONAL DOMINANCE MATRIX (LONG)
# ------------------------------------------------------------------------------
# For each metric, compute per-(region, detector) normalized dominance score
# in [0, 1] using min-max normalization within (region, metric) across the
# 11 detectors. All five reported metrics are higher-is-better, so no
# directional flip is needed.
# ------------------------------------------------------------------------------
rescale_01_safe_local <- function(v) {
  vmin <- min(v, na.rm = TRUE)
  vmax <- max(v, na.rm = TRUE)
  if (!is.finite(vmin) || !is.finite(vmax) || vmax == vmin)
    return(rep(NA_real_, length(v)))
  (v - vmin) / (vmax - vmin)
}

compute_regional_dominance_long <- function() {
  rows <- list()
  for (m in HH_METRICS_REGIONAL) {
    higher_better <- hh_metric_higher_better[[m]]
    sub <- regional_metrics |>
      dplyr::transmute(
        REGION = as.character(REGION),
        Method = as.character(Method),
        Paradigm = as.character(Paradigm),
        Value = .data[[m]]
      )
    sub <- sub |>
      dplyr::group_by(REGION) |>
      dplyr::mutate(
        raw_norm = rescale_01_safe_local(Value),
        Dominance_Score = if (higher_better) raw_norm else 1 - raw_norm
      ) |>
      dplyr::ungroup() |>
      dplyr::mutate(
        Metric = m,
        Is_Dominant = !is.na(Dominance_Score) &
          Dominance_Score >= DOMINANCE_THRESHOLD
      )
    rows[[m]] <- sub
  }
  dplyr::bind_rows(rows)
}
regional_dominance_long <- compute_regional_dominance_long()

# Single-color (blue) cell fill scaled by dominance score.
DOMINANCE_BLUE_LIGHT <- "#FFFFFF"
DOMINANCE_BLUE_DARK  <- "#0B2447"
mix_blue_intensity <- function(score) {
  if (is.na(score)) return("#F2F2F2")
  rgb_lo <- grDevices::col2rgb(DOMINANCE_BLUE_LIGHT) / 255
  rgb_hi <- grDevices::col2rgb(DOMINANCE_BLUE_DARK)  / 255
  r <- rgb_lo[1] * (1 - score) + rgb_hi[1] * score
  g <- rgb_lo[2] * (1 - score) + rgb_hi[2] * score
  b <- rgb_lo[3] * (1 - score) + rgb_hi[3] * score
  grDevices::rgb(r, g, b, maxColorValue = 1)
}
regional_dominance_long$cell_fill <- vapply(
  regional_dominance_long$Dominance_Score, mix_blue_intensity, character(1)
)

# ------------------------------------------------------------------------------
# 15. METHOD-ACROSS-METRICS AGGREGATION
# ------------------------------------------------------------------------------
compute_method_summary_long <- function() {
  regional_dominance_long |>
    dplyr::filter(Metric %in% HH_METRICS_REGIONAL) |>
    dplyr::group_by(Method, Metric) |>
    dplyr::summarise(
      Sweep_Count = sum(Is_Dominant, na.rm = TRUE),
      N_Regions   = sum(!is.na(Dominance_Score)),
      .groups     = "drop"
    ) |>
    dplyr::mutate(
      Method = factor(as.character(Method), levels = method_order),
      Metric = factor(as.character(Metric), levels = HH_METRICS_REGIONAL)
    )
}
method_summary_long <- compute_method_summary_long()

method_summary_agg <- method_summary_long |>
  dplyr::group_by(Method) |>
  dplyr::summarise(
    Mean_Sweep = mean(Sweep_Count, na.rm = TRUE),
    SD_Sweep   = stats::sd(Sweep_Count, na.rm = TRUE),
    .groups    = "drop"
  ) |>
  dplyr::mutate(
    grand_mean = mean(Mean_Sweep, na.rm = TRUE),
    grand_sd   = stats::sd(Mean_Sweep, na.rm = TRUE),
    Z_Sweep    = ifelse(is.finite(grand_sd) & grand_sd > 0,
                        (Mean_Sweep - grand_mean) / grand_sd,
                        NA_real_),
    Sweep_Class = dplyr::case_when(
      is.na(Z_Sweep)  ~ "Insufficient",
      Z_Sweep >=  1.0 ~ "Sweeper",
      Z_Sweep <= -1.0 ~ "Non-sweeper",
      TRUE            ~ "Average"
    )
  ) |>
  dplyr::arrange(dplyr::desc(Mean_Sweep))

method_summary_order <- as.character(method_summary_agg$Method)

MAX_REGIONS_POSSIBLE <- {
  .mrp <- suppressWarnings(max(method_summary_long$N_Regions, na.rm = TRUE))
  if (!is.finite(.mrp) || .mrp <= 0) 17L else as.integer(.mrp)
}

mix_blue_count <- function(count, max_count = MAX_REGIONS_POSSIBLE) {
  if (is.na(count) || !is.finite(count)) return("#FFFFFF")
  frac <- pmax(0, pmin(1, count / max(1, max_count)))
  mix_blue_intensity(frac)
}
method_summary_long$cell_fill <- vapply(
  method_summary_long$Sweep_Count, mix_blue_count, character(1)
)

# ------------------------------------------------------------------------------
# 16. FIGURE 4 BUILDERS - PER-METRIC PANELS AND METHOD SUMMARY
# ------------------------------------------------------------------------------

# 16a. Regional dominance matrix (per metric) ---------------------------------
build_regional_dominance_matrix <- function(metric_id) {
  metric_label <- hh_metric_full_name[[metric_id]]
  this_metric_df <- regional_dominance_long |>
    dplyr::filter(Metric == metric_id)

  # Per-method count of regions swept at score >= threshold.
  per_method_count <- this_metric_df |>
    dplyr::group_by(Method) |>
    dplyr::summarise(
      N_Regions_Swept = sum(Is_Dominant, na.rm = TRUE),
      N_Regions_Total = sum(!is.na(Dominance_Score)),
      Mean_Score      = mean(Dominance_Score, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(Sweep_Label = paste0(N_Regions_Swept, "/", N_Regions_Total))

  method_sweep_order <- per_method_count |>
    dplyr::arrange(dplyr::desc(N_Regions_Swept),
                   dplyr::desc(Mean_Score)) |>
    dplyr::pull(Method) |> as.character()

  region_score_order <- this_metric_df |>
    dplyr::group_by(REGION) |>
    dplyr::summarise(mean_score = mean(Dominance_Score, na.rm = TRUE),
                     .groups = "drop") |>
    dplyr::arrange(dplyr::desc(mean_score)) |>
    dplyr::pull(REGION) |> as.character()

  this_metric_df <- this_metric_df |>
    dplyr::mutate(
      REGION = factor(as.character(REGION), levels = region_score_order),
      Method = factor(as.character(Method), levels = rev(method_sweep_order))
    )
  per_method_count <- per_method_count |>
    dplyr::mutate(Method = factor(as.character(Method),
                                  levels = rev(method_sweep_order)))

  N_REGIONS     <- length(region_score_order)
  SWEEP_COUNT_X <- N_REGIONS + 1L

  ggplot2::ggplot(this_metric_df) +
    ggplot2::geom_tile(
      ggplot2::aes(x = REGION, y = Method),
      fill = this_metric_df$cell_fill,
      colour = "white", linewidth = 0.5, show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = per_method_count,
      ggplot2::aes(x = SWEEP_COUNT_X, y = Method, label = Sweep_Label),
      size = 2.6, fontface = "bold", family = base_family_global,
      colour = "#0B2447", inherit.aes = FALSE, show.legend = FALSE
    ) +
    ggplot2::geom_vline(xintercept = SWEEP_COUNT_X - 0.5,
                        linetype = "solid",
                        linewidth = 0.50, colour = "grey25") +
    ggplot2::annotate(
      "text", x = SWEEP_COUNT_X, y = length(method_sweep_order) + 1.20,
      label = paste0("Regions Swept\n(score \u2265 ",
                     sprintf("%.2f", DOMINANCE_THRESHOLD), ")"),
      size = 2.2, family = base_family_global, fontface = "bold",
      colour = "#0B2447", lineheight = 0.92, vjust = 0
    ) +
    ggplot2::scale_x_discrete(
      limits = c(region_score_order, "__SWEEP_COUNT__"),
      labels = {
        raw <- as.character(c(region_score_order, " "))
        ifelse(is.na(raw) | raw == "", " ", raw)
      },
      position = "top",
      expand = ggplot2::expansion(add = c(0.04, 0.50))
    ) +
    ggplot2::scale_y_discrete(
      limits = rev(method_sweep_order),
      labels = {
        raw <- as.character(unname(method_two_line[rev(method_sweep_order)]))
        ifelse(is.na(raw) | raw == "", " ", raw)
      },
      expand = ggplot2::expansion(add = c(0.30, 1.45))
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(
      title = paste0("a", hh_metric_subscript[[metric_id]],
                     "    Regional dominance matrix on ", metric_label,
                     " (rows = methods, columns = regions)"),
      x = NULL, y = NULL
    ) +
    theme_dashboard(base_size = PUB_BASE, base_family = base_family_global) +
    theme_bold_axes() +
    ggplot2::theme(
      axis.text.x.top = ggplot2::element_text(color = "black", size = 6.8,
                                              hjust = 0.5, vjust = 0,
                                              margin = ggplot2::margin(b = 4)),
      axis.text.y     = ggplot2::element_text(color = "black", size = 6.8,
                                              lineheight = 0.85,
                                              hjust = 1, vjust = 0.5,
                                              margin = ggplot2::margin(r = 4)),
      panel.grid    = ggplot2::element_blank(),
      panel.border  = ggplot2::element_blank(),
      legend.position = "none",
      plot.margin   = ggplot2::margin(40, 14, 8, 8)
    )
}

# 16b. Dominance score legend strip -------------------------------------------
build_dominance_legend <- function() {
  legend_score_seq <- seq(0, 1, length.out = 100L)
  legend_df <- data.frame(
    x = legend_score_seq,
    y = 1L,
    fill_col = vapply(legend_score_seq, mix_blue_intensity, character(1))
  )
  ggplot2::ggplot(legend_df) +
    ggplot2::geom_tile(ggplot2::aes(x = x, y = y),
                       fill = legend_df$fill_col,
                       width = 1 / nrow(legend_df),
                       height = 1) +
    ggplot2::scale_x_continuous(
      breaks = c(0, 0.25, 0.50, DOMINANCE_THRESHOLD, 1.00),
      labels = c("0.00\nweak", "0.25", "0.50\nmoderate",
                 sprintf("%.2f\ndominant", DOMINANCE_THRESHOLD),
                 "1.00\nsweep"),
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0))) +
    ggplot2::labs(
      title = paste0("Dominance score (min-max normalized within region; ",
                     "darker = stronger dominance)"),
      x = NULL, y = NULL
    ) +
    ggplot2::theme_void(base_family = base_family_global) +
    ggplot2::theme(
      plot.title  = ggplot2::element_text(face = "bold", size = 7.5,
                                          hjust = 0,
                                          margin = ggplot2::margin(b = 4)),
      axis.text.x = ggplot2::element_text(size = PUB_AXIS_TXT, lineheight = 0.95,
                                          margin = ggplot2::margin(t = 2)),
      plot.margin = ggplot2::margin(10, 12, 8, 8)
    )
}

# 16c. Per-detector dot plot --------------------------------------------------
build_regional_dot_panel <- function(metric_id) {
  metric_label  <- hh_metric_full_name[[metric_id]]
  ylab          <- hh_metric_y_label[[metric_id]]
  higher_better <- hh_metric_higher_better[[metric_id]]
  threshold     <- hh_metric_threshold[[metric_id]]
  is_proportion <- hh_metric_is_proportion[[metric_id]]
  ci_cols       <- hh_metric_ci_cols[[metric_id]]

  plot_df <- regional_metrics |>
    dplyr::transmute(
      Method, REGION,
      Value = .data[[metric_id]]
    ) |>
    dplyr::mutate(REGION_chr = as.character(REGION),
                  Method_chr = as.character(Method))

  if (all(ci_cols %in% names(boot_ci))) {
    plot_df <- plot_df |>
      dplyr::left_join(
        boot_ci |>
          dplyr::transmute(
            REGION_chr = as.character(REGION),
            Method_chr = as.character(Method),
            ci_lo = .data[[ci_cols[1]]],
            ci_hi = .data[[ci_cols[2]]]
          ),
        by = c("REGION_chr", "Method_chr")
      )
  } else {
    plot_df$ci_lo <- NA_real_
    plot_df$ci_hi <- NA_real_
  }

  fav_layer <- NULL
  ref_line  <- NULL
  if (!is.na(threshold)) {
    if (higher_better) {
      fav_layer <- ggplot2::annotate(
        "rect", xmin = -Inf, xmax = Inf,
        ymin = threshold, ymax = Inf,
        fill = "#cfe2f3", alpha = 0.30
      )
    } else {
      fav_layer <- ggplot2::annotate(
        "rect", xmin = -Inf, xmax = Inf,
        ymin = -Inf, ymax = threshold,
        fill = "#cfe2f3", alpha = 0.30
      )
    }
    ref_line <- ggplot2::geom_hline(
      yintercept = threshold, linetype = "dashed",
      linewidth = 0.35, colour = "grey45"
    )
  }

  pj <- ggplot2::position_jitter(width = 0.18, height = 0, seed = 12345)
  p  <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Method, y = Value))
  if (!is.null(fav_layer)) p <- p + fav_layer
  if (!is.null(ref_line))  p <- p + ref_line

  p <- p +
    ggplot2::geom_linerange(
      ggplot2::aes(ymin = ci_lo, ymax = ci_hi),
      position = pj,
      colour = "grey55", alpha = 0.55, linewidth = 0.30, na.rm = TRUE
    ) +
    ggplot2::geom_point(
      ggplot2::aes(fill = Method),
      shape = 21, position = pj,
      colour = "grey20", stroke = 0.30, size = 2.2, alpha = 0.92,
      show.legend = FALSE, na.rm = TRUE
    ) +
    ggplot2::stat_summary(
      fun = mean, geom = "crossbar",
      width = 0.50, linewidth = 0.40,
      colour = "black", fatten = 0,
      fill = NA, na.rm = TRUE
    ) +
    ggplot2::scale_fill_manual(
      values = stats::setNames(rep("#FFFFFF", length(method_order)),
                               as.character(method_order)),
      drop = FALSE
    ) +
    ggplot2::scale_x_discrete(
      limits = method_order,
      labels = {
        raw <- as.character(unname(method_two_line[as.character(method_order)]))
        ifelse(is.na(raw) | raw == "", " ", raw)
      }
    ) +
    ggplot2::labs(
      title = paste0("b", hh_metric_subscript[[metric_id]],
                     "    Per-detector regional values on ", metric_label),
      x = NULL, y = ylab
    ) +
    theme_dashboard(base_size = PUB_BASE, base_family = base_family_global) +
    theme_bold_axes() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = PUB_AXIS_TXT, lineheight = 0.85,
                                          hjust = 0.5, vjust = 1,
                                          margin = ggplot2::margin(t = 2)),
      axis.text.y = ggplot2::element_text(size = PUB_AXIS_TXT),
      axis.title.y = ggplot2::element_text(size = PUB_AXIS_TIT,
                                           margin = ggplot2::margin(r = 4)),
      legend.position = "none",
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(10, 12, 22, 8)
    )

  if (is_proportion) {
    p <- p + ggplot2::scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      limits = c(0, 1.05),
      expand = ggplot2::expansion(mult = c(0, 0))
    )
  }
  p
}

# 16d. Wilcoxon detector-paired sig bars --------------------------------------
build_regional_sig_panel <- function(metric_id, all_wilcoxon_results) {
  res <- all_wilcoxon_results |>
    dplyr::filter(Metric == metric_id) |>
    dplyr::mutate(Comparison = factor(Comparison,
                                      levels = vapply(DETECTOR_PAIRS,
                                                      function(p) p$short,
                                                      character(1))))
  res <- res |>
    dplyr::mutate(
      p_label = ifelse(
        is.na(p_bonferroni),
        "n/a",
        ifelse(p_bonferroni < 0.001,
               "p < 0.001",
               paste0("p = ", sprintf("%.3f", p_bonferroni)))
      ),
      sig_star = ifelse(!is.na(Significant_005) & Significant_005, "*", ""),
      bar_fill = ifelse(!is.na(Significant_005) & Significant_005,
                        "#1F2D5C", "#C7CCD9")
    )

  x_max <- max(res$Abs_Median_Diff, na.rm = TRUE)
  if (!is.finite(x_max) || x_max <= 0) x_max <- 1
  x_lim_hi <- x_max * 1.45

  ggplot2::ggplot(res, ggplot2::aes(y = Comparison, x = Abs_Median_Diff)) +
    ggplot2::geom_col(
      fill = res$bar_fill,
      width = 0.55, na.rm = TRUE
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = paste0(p_label, " ", sig_star),
        x = ifelse(is.na(Abs_Median_Diff),
                   x_lim_hi * 0.05,
                   pmin(Abs_Median_Diff + x_lim_hi * 0.04, x_lim_hi * 0.95))
      ),
      hjust = 0,
      family = base_family_global, size = pub_text_size(PUB_ANNOT),
      fontface = "bold",
      colour = "black", na.rm = TRUE   # p-values in bold black
    ) +
    # CROPPING FIX. The p-value labels are drawn hjust = 0 starting as far right
    # as 0.95 * x_lim_hi, so with zero expansion their text ran past the panel
    # edge and was cut. The right-hand expansion reserves room for them.
    ggplot2::scale_x_continuous(
      limits = c(0, x_lim_hi),
      expand = ggplot2::expansion(mult = c(0, 0.30))
    ) +
    ggplot2::scale_y_discrete(
      limits = rev(vapply(DETECTOR_PAIRS, function(p) p$short, character(1)))
    ) +
    ggplot2::labs(
      title = paste0("c", hh_metric_subscript[[metric_id]],
                     "    Wilcoxon detector-paired (n = ",
                     {
                       n_max <- suppressWarnings(max(res$N_Regions_Paired,
                                                     na.rm = TRUE))
                       if (!is.finite(n_max)) "?" else as.character(n_max)
                     },
                     " regions)"),
      x = "|Median diff|",
      y = NULL
    ) +
    theme_dashboard(base_size = PUB_BASE - 0.5, base_family = base_family_global) +
    theme_bold_axes() +
    ggplot2::theme(
      axis.text.y  = ggplot2::element_text(size = PUB_AXIS_TXT, face = "bold",
                                           lineheight = 0.95),
      axis.text.x  = ggplot2::element_text(size = PUB_AXIS_TXT),
      axis.title.x = ggplot2::element_text(size = PUB_AXIS_TIT),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),
      legend.position    = "none"
    )
}

# 16e. Compose one per-metric multipanel figure -------------------------------
build_regional_multipanel <- function(metric_id) {
  top          <- build_regional_dominance_matrix(metric_id)
  legend_strip <- build_dominance_legend()
  bot_l        <- build_regional_dot_panel(metric_id)
  bot_r        <- build_regional_sig_panel(metric_id, all_wilcoxon_results)
  bottom       <- (bot_l | bot_r) + patchwork::plot_layout(widths = c(7, 3))

  (top / legend_strip / bottom) +
    patchwork::plot_layout(heights = c(1.55, 0.18, 1.00)) +
    patchwork::plot_annotation(
      title = paste0("Regional generalisability: ",
                     hh_metric_full_name[[metric_id]]),
      # Panel letters are INLINE in each sub-panel's own title
      # ("a<sub>metric</sub>    Regional dominance matrix on ...").
      # No patchwork tag is used: tag_levels would also letter the legend
      # strip, and a separate tag would print a second letter above the title.
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = PUB_TITLE + 1,
                                           family = base_family_global,
                                           margin = ggplot2::margin(b = 6)),
        plot.tag = ggplot2::element_text(size = PUB_TAG, face = "bold",
                                         colour = "black",
                                         family = base_family_global,
                                       hjust = 0, vjust = 1),
      # Same left edge as the title: tags anchor to the plot, and
      # plot.title.position = "plot" moves the title off the panel edge to match.
      plot.tag.position     = "topleft",
      plot.title.position   = "plot",
        plot.margin = ggplot2::margin(10, 12, 8, 8)
      )
    ) &
    # Applied to EVERY sub-panel with `&`. plot.title.position = "plot" keeps
    # each inline-lettered title flush with the plot's left edge rather than
    # indented past the y-axis labels.
    ggplot2::theme(
      plot.margin           = ggplot2::margin(10, 12, 8, 8),
      plot.tag              = ggplot2::element_text(size = PUB_TAG,
                                                    face = "bold",
                                                    colour = "black",
                                                    family = base_family_global,
                                                    hjust = 0, vjust = 1),
      plot.tag.position     = "topleft",
      plot.title.position   = "plot"
    )
}

# 16f. Method-across-metrics summary panels -----------------------------------
build_method_summary_matrix <- function() {
  df <- method_summary_long |>
    dplyr::mutate(
      Method = factor(as.character(Method), levels = rev(method_summary_order)),
      Metric = factor(as.character(Metric), levels = HH_METRICS_REGIONAL)
    )
  agg <- method_summary_agg |>
    dplyr::mutate(
      Method = factor(as.character(Method), levels = rev(method_summary_order)),
      Summary_Label = sprintf("%.1f (z=%+.2f)", Mean_Sweep, Z_Sweep)
    )

  N_METRICS <- length(HH_METRICS_REGIONAL)
  SUMMARY_X <- N_METRICS + 1L
  N_METHODS <- length(method_summary_order)
  HEADER_Y  <- N_METHODS + 1L

  metric_full_names <- vapply(HH_METRICS_REGIONAL,
                              function(m) hh_metric_full_name[[m]],
                              character(1))
  metric_cell_labels <- c(
    "TAM"            = "True-Alarm\nMagnitude",
    "N_True_Alarms"  = "Number of\nTrue Alarms",
    "Sensitivity"    = "Sensitivity",
    "Mean_Lead_Time" = "Mean\nLead Time",
    "WP"             = "Warning\nPersistence"
  )
  for (m in HH_METRICS_REGIONAL) {
    if (is.na(metric_cell_labels[m]) || is.null(metric_cell_labels[[m]])) {
      metric_cell_labels[m] <- metric_full_names[[m]]
    }
  }

  header_df <- data.frame(
    x_pos = c(seq_len(N_METRICS), SUMMARY_X),
    y_pos = HEADER_Y,
    label = c(unname(metric_cell_labels[HH_METRICS_REGIONAL]),
              "Mean Sweep\n(z-score)"),
    stringsAsFactors = FALSE
  )

  ggplot2::ggplot(df) +
    ggplot2::geom_tile(
      ggplot2::aes(x = as.integer(Metric), y = as.integer(Method)),
      fill = df$cell_fill, colour = "white", linewidth = 0.5
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = as.integer(Metric), y = as.integer(Method),
                   label = sprintf("%d", Sweep_Count),
                   colour = Sweep_Count >= round(MAX_REGIONS_POSSIBLE * 0.55)),
      size = 2.6, family = base_family_global, fontface = "bold",
      show.legend = FALSE
    ) +
    ggplot2::scale_colour_manual(
      values = c("FALSE" = "grey15", "TRUE" = "white"),
      guide  = "none"
    ) +
    ggplot2::geom_text(
      data = agg,
      ggplot2::aes(x = SUMMARY_X, y = as.integer(Method), label = Summary_Label),
      size = 2.5, fontface = "bold", family = base_family_global,
      colour = "#0B2447", inherit.aes = FALSE
    ) +
    # Header background tiles removed (no background shading).
    ggplot2::geom_text(
      data = header_df,
      ggplot2::aes(x = x_pos, y = y_pos, label = label),
      size = 2.4, family = base_family_global, fontface = "bold",
      colour = "#0B2447", lineheight = 0.85, inherit.aes = FALSE
    ) +
    ggplot2::geom_vline(xintercept = SUMMARY_X - 0.5,
                        linetype = "solid",
                        linewidth = 0.50, colour = "grey25") +
    ggplot2::scale_x_continuous(
      breaks = NULL,
      limits = c(0.5, SUMMARY_X + 0.5),
      expand = ggplot2::expansion(add = c(0, 0))
    ) +
    ggplot2::scale_y_continuous(
      breaks = c(seq_len(N_METHODS), HEADER_Y),
      labels = c(unname(method_two_line[rev(method_summary_order)]), ""),
      limits = c(0.5, HEADER_Y + 0.5),
      expand = ggplot2::expansion(add = c(0.05, 0.30))
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(
      title = "Regional dominance matrix on 5 operational metrics",
      x = NULL, y = NULL
    ) +
    theme_dashboard(base_size = PUB_BASE, base_family = base_family_global) +
    theme_bold_axes() +
    ggplot2::theme(
      axis.text.x.top = ggplot2::element_blank(),
      axis.text.x     = ggplot2::element_blank(),
      axis.ticks.x    = ggplot2::element_blank(),
      axis.text.y     = ggplot2::element_text(color = "black", size = 6.8,
                                              lineheight = 0.85,
                                              hjust = 1, vjust = 0.5,
                                              margin = ggplot2::margin(r = 4)),
      panel.grid    = ggplot2::element_blank(),
      panel.border  = ggplot2::element_blank(),
      legend.position = "none",
      plot.margin   = ggplot2::margin(20, 14, 8, 8)
    )
}

build_method_summary_dot <- function() {
  df_dots <- method_summary_long |>
    dplyr::mutate(
      Method = factor(as.character(Method), levels = rev(method_summary_order))
    )
  df_means <- method_summary_agg |>
    dplyr::mutate(
      Method = factor(as.character(Method), levels = rev(method_summary_order))
    )

  ggplot2::ggplot() +
    ggplot2::geom_jitter(
      data = df_dots,
      ggplot2::aes(x = Sweep_Count, y = Method),
      width = 0, height = 0.18,
      shape = 21, fill = "#3673B6", colour = "white",
      stroke = 0.4, size = 2.4, alpha = 0.85
    ) +
    ggplot2::geom_segment(
      data = df_means,
      ggplot2::aes(x = Mean_Sweep, xend = Mean_Sweep,
                   y    = as.numeric(Method) - 0.30,
                   yend = as.numeric(Method) + 0.30),
      colour = "#0B2447", linewidth = 1.0
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0, MAX_REGIONS_POSSIBLE + 0.5),
      breaks = pretty(c(0, MAX_REGIONS_POSSIBLE), n = 5),
      expand = ggplot2::expansion(mult = c(0.02, 0.05))
    ) +
    ggplot2::scale_y_discrete(
      limits = rev(method_summary_order),
      labels = unname(method_two_line[rev(method_summary_order)])
    ) +
    ggplot2::labs(
      title = "Per-method sweep counts across the 5 metrics",
      x = "Regions swept (count)", y = NULL
    ) +
    theme_dashboard(base_size = PUB_BASE, base_family = base_family_global) +
    theme_bold_axes() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(color = "black", size = 6.8,
                                          lineheight = 0.85,
                                          hjust = 1, vjust = 0.5,
                                          margin = ggplot2::margin(r = 4)),
            panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      legend.position = "none",
      plot.margin = ggplot2::margin(10, 12, 8, 8)
    )
}

# Composite Wilcoxon panel for the method summary figure (DUAL-STATISTIC).
#
#   p-value      : paired Wilcoxon signed-rank on within-region Dominance_Score
#                  paired by (REGION x METRIC). n_pairs ~ 17 x 5 = 85.
#                  Properly powered; p < 0.05 is reachable.
#   Bar length   : |median(diff in regions-swept count)| across the 5 reported
#                  metrics. Interpretable in operational units of "regions".
#
# The dual specification answers two distinct questions in the same panel:
# (i) is the difference statistically real (n=85 powered test); and (ii) how
# big is the practical difference in operational regions (n=5 sweep-count
# median). Both numbers are also written to the console for auditability.
build_method_summary_sig <- function() {
  # REVISION: trimmed from three pairs to the two instructed, both against
  # Outbreak Threshold. Constant TA vs Continuous TA is no longer computed
  # or plotted here. See Item 8 of the Delphi review.
  TARGETS <- c("Constant Transmission Acceleration",
               "Continuous Transmission Acceleration",
               "Outbreak Threshold")
  PAIRS_LIST <- list(
    list(a = "Constant Transmission Acceleration",
         b = "Outbreak Threshold",
         label = "Constant TA vs Outbreak Threshold"),
    list(a = "Continuous Transmission Acceleration",
         b = "Outbreak Threshold",
         label = "Continuous TA vs Outbreak Threshold")
  )

  pivot_score <- regional_dominance_long |>
    dplyr::filter(as.character(Method) %in% TARGETS,
                  as.character(Metric) %in% HH_METRICS_REGIONAL) |>
    dplyr::select(REGION, Metric, Method, Dominance_Score) |>
    dplyr::mutate(Method = as.character(Method)) |>
    tidyr::pivot_wider(names_from = Method, values_from = Dominance_Score)

  pivot_count <- method_summary_long |>
    dplyr::filter(as.character(Method) %in% TARGETS,
                  as.character(Metric) %in% HH_METRICS_REGIONAL) |>
    dplyr::select(Method, Metric, Sweep_Count) |>
    dplyr::mutate(Method = as.character(Method)) |>
    tidyr::pivot_wider(names_from = Method, values_from = Sweep_Count)

  rows <- list()
  for (pr in PAIRS_LIST) {
    a_score <- pivot_score[[pr$a]]; b_score <- pivot_score[[pr$b]]
    keep_s  <- !is.na(a_score) & !is.na(b_score)
    a_score <- a_score[keep_s]; b_score <- b_score[keep_s]
    if (length(a_score) >= 2L && any(a_score != b_score)) {
      tt <- suppressWarnings(stats::wilcox.test(a_score, b_score,
                                                paired = TRUE,
                                                exact = FALSE))
      pval           <- tt$p.value
      med_score_diff <- stats::median(a_score - b_score, na.rm = TRUE)
    } else {
      pval           <- NA_real_
      med_score_diff <- NA_real_
    }
    n_score_pairs <- length(a_score)

    a_count <- pivot_count[[pr$a]]; b_count <- pivot_count[[pr$b]]
    keep_c  <- !is.na(a_count) & !is.na(b_count)
    a_count <- a_count[keep_c]; b_count <- b_count[keep_c]
    med_count_diff <- if (length(a_count) >= 1L) {
      stats::median(a_count - b_count, na.rm = TRUE)
    } else {
      NA_real_
    }

    rows[[length(rows) + 1L]] <- data.frame(
      Pair               = pr$label,
      a_method           = pr$a,
      b_method           = pr$b,
      Median_Score_Diff  = med_score_diff,
      Median_Count_Diff  = med_count_diff,
      P_Value            = pval,
      N_Score_Pairs      = n_score_pairs,
      N_Metric_Pairs     = length(a_count),
      Sig                = !is.na(pval) & pval < 0.05,
      stringsAsFactors   = FALSE
    )
  }
  res <- dplyr::bind_rows(rows)
  res$Pair_short <- factor(res$Pair, levels = res$Pair)

  res$Bar_Mag_raw <- ifelse(is.na(res$Median_Count_Diff), 0,
                            abs(res$Median_Count_Diff))
  res$Sig_plot    <- ifelse(is.na(res$Sig), FALSE, res$Sig)
  MIN_VISIBLE_BAR <- 0.5
  res$Bar_Mag <- ifelse(res$Sig_plot & res$Bar_Mag_raw < MIN_VISIBLE_BAR,
                        MIN_VISIBLE_BAR, res$Bar_Mag_raw)

  cat("\n--- Method Summary Panel c: dual-statistic Wilcoxon ---\n")
  cat("    Test: paired Wilcoxon on Dominance_Score (17 REGIONS x 5 METRICS, n ~ 85)\n")
  cat("    Bar:  |median(diff in sweep counts)| across the 5 metrics, units 'regions'\n")
  print(res[, c("Pair", "Median_Score_Diff", "Median_Count_Diff",
                "P_Value", "N_Score_Pairs", "N_Metric_Pairs", "Sig")],
        row.names = FALSE)

  n_typical <- {
    .nm <- suppressWarnings(median(res$N_Score_Pairs, na.rm = TRUE))
    if (!is.finite(.nm)) "?" else as.character(round(.nm))
  }
  res$Sig_star <- ifelse(res$Sig_plot, "*", "")

  ggplot2::ggplot(res) +
    ggplot2::geom_col(
      ggplot2::aes(x = Bar_Mag, y = Pair_short, fill = Sig_plot),
      width = 0.65
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = Bar_Mag, y = Pair_short,
                   label = paste0(
                     ifelse(is.na(P_Value), "n/a",
                            ifelse(P_Value < 0.001,
                                   "p < 0.001",
                                   sprintf("p = %.3g", P_Value))),
                     " ", Sig_star
                   )),
      hjust = -0.10, size = 2.5, family = base_family_global,
      fontface = "bold", colour = "grey15"
    ) +
    ggplot2::scale_fill_manual(
      values = c("FALSE" = "#BFC6CE", "TRUE" = "#3673B6"),
      guide  = "none"
    ) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.30))
    ) +
    ggplot2::labs(
      title = paste0("c    Wilcoxon detector-paired (n \u2248 ",
                     n_typical, ", 17 regions \u00D7 5 metrics)"),
      x = "Median diff in regions swept (across 5 metrics)",
      y = NULL
    ) +
    theme_dashboard(base_size = PUB_BASE, base_family = base_family_global) +
    theme_bold_axes() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(color = "black", size = 7.0,
                                          hjust = 1, vjust = 0.5,
                                          margin = ggplot2::margin(r = 4)),
            panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(10, 14, 8, 8)
    )
}

build_method_summary_multipanel <- function() {
  top          <- build_method_summary_matrix()
  legend_strip <- build_dominance_legend()
  bot_l        <- build_method_summary_dot()
  bot_r        <- build_method_summary_sig()
  bottom       <- (bot_l | bot_r) + patchwork::plot_layout(widths = c(7, 3))

  (top / legend_strip / bottom) +
    patchwork::plot_layout(heights = c(1.50, 0.18, 1.00)) +
    patchwork::plot_annotation(
      title = paste0("Figure 4 | Method-across-metrics summary ",
                     "(6 methods x 5 reported metrics)"),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = PUB_TITLE + 1,
                                           family = base_family_global,
                                           margin = ggplot2::margin(b = 6)),
        plot.margin = ggplot2::margin(10, 12, 8, 8)
      )
    ) &
    # Applied to EVERY sub-panel with `&`. plot.title.position = "plot" keeps
    # each inline-lettered title flush with the plot's left edge rather than
    # indented past the y-axis labels.
    ggplot2::theme(
      plot.margin           = ggplot2::margin(10, 12, 8, 8),
      plot.tag              = ggplot2::element_text(size = PUB_TAG,
                                                    face = "bold",
                                                    colour = "black",
                                                    family = base_family_global,
                                                    hjust = 0, vjust = 1),
      plot.tag.position     = "topleft",
      plot.title.position   = "plot"
    )
}

# ------------------------------------------------------------------------------
# 17. SAVE FIGURE 4 (PER-METRIC + METHOD SUMMARY)
# ------------------------------------------------------------------------------
# Figure 4 set, authored at final print size. Previously 16 x 13.5 in, which
# production reduced by ~0.44x, printing axis text near 3.5 pt.
fig4_w <- NC_W_DOUBLE
# Composite stacks a dominance matrix, a legend strip and two sub-panels;
# the taller supplementary canvas prevents the rows from compressing.
fig4_h <- NC_H_SUPP

cat("\n=== Building 5 per-metric multipanel figures (Figure 4 set) ===\n")
for (m in HH_METRICS_REGIONAL) {
  fig <- build_regional_multipanel(m)
  slug <- hh_metric_slug[[m]]
  save_plot_file(paste0("FigureB4_", slug), fig, fig4_w, fig4_h)
}

cat("\n=== Building Method-Across-Metrics Summary figure (Figure 4 set) ===\n")
{
  fig_summary <- build_method_summary_multipanel()
  save_plot_file("FigureB4_Method_Summary", fig_summary, fig4_w, fig4_h)
}

# ------------------------------------------------------------------------------
# 18. WRITE FIGURE 4 CSV TABLES
# ------------------------------------------------------------------------------
# Wide regional metrics table (carries the 5 reported metrics plus diagnostic
# columns). RFR / RFR_RX columns dropped (not used by this pipeline).
utils::write.csv(
  regional_metrics |> dplyr::mutate(REGION = as.character(REGION),
                                     Method = as.character(Method)),
  file.path(OUTPUT_DIR, "StageB_Regional_Framework_Metrics.csv"),
  row.names = FALSE
)

# 8-metric harmonized regional CSV.
regional_8metric_csv <- regional_metrics |>
  dplyr::mutate(REGION = as.character(REGION),
                Method = as.character(Method),
                Paradigm = as.character(Paradigm)) |>
  dplyr::transmute(
    REGION, Method, Paradigm,
    TAM, N_True_Alarms, PPV, Sensitivity,
    Mean_Lead_Time, WP, ALY,
    Mean_Lead_Time_conditional, WP_conditional,
    N_False_Alarms,
    Total_Triggers, True_Alarms, False_Alarms,
    n_TrueActionable, n_Reactive,
    N_Years_with_A1_True, N_Years_Evaluable
  )
utils::write.csv(
  regional_8metric_csv,
  file.path(OUTPUT_DIR, "StageB_Regional_8Metric_Summary.csv"),
  row.names = FALSE
)

# 8-metric CSV with bootstrap 95% CIs joined per (REGION, Method).
regional_8metric_with_ci <- regional_8metric_csv |>
  dplyr::left_join(
    boot_ci |> dplyr::mutate(REGION = as.character(REGION),
                              Method = as.character(Method)) |>
      dplyr::select(REGION, Method,
                    TAM_lo, TAM_hi,
                    N_True_Alarms_lo, N_True_Alarms_hi,
                    PPV_lo, PPV_hi,
                    Sens_lo, Sens_hi,
                    MLT_lo, MLT_hi,
                    WP_lo, WP_hi,
                    ALY_lo, ALY_hi,
                    N_False_Alarms_lo, N_False_Alarms_hi),
    by = c("REGION", "Method")
  )
utils::write.csv(
  regional_8metric_with_ci,
  file.path(OUTPUT_DIR, "StageB_Regional_8Metric_Summary_with_CIs.csv"),
  row.names = FALSE
)

# Wide regional metrics table with bootstrap CIs.
metrics_with_ci <- regional_metrics |>
  dplyr::mutate(REGION = as.character(REGION),
                Method = as.character(Method)) |>
  dplyr::left_join(
    boot_ci |>
      dplyr::mutate(REGION = as.character(REGION),
                    Method = as.character(Method)),
    by = c("REGION", "Method")
  )
utils::write.csv(
  metrics_with_ci,
  file.path(OUTPUT_DIR, "StageB_Regional_Framework_Metrics_with_CIs.csv"),
  row.names = FALSE
)

# Dominance probability (region-level).
utils::write.csv(
  dominance_df,
  file.path(OUTPUT_DIR, "StageB_Regional_Dominance_Probabilities.csv"),
  row.names = FALSE
)

# Regional Dominance Matrix (long).
utils::write.csv(
  regional_dominance_long |>
    dplyr::transmute(
      REGION = as.character(REGION),
      Method = as.character(Method),
      Paradigm = as.character(Paradigm),
      Metric = Metric,
      Value = Value,
      Dominance_Score = Dominance_Score,
      Is_Dominant_at_threshold_075 = Is_Dominant
    ),
  file.path(OUTPUT_DIR, "StageB_Regional_Dominance_Matrix.csv"),
  row.names = FALSE
)

# Wilcoxon per-metric detector-paired results.
utils::write.csv(
  all_wilcoxon_results,
  file.path(OUTPUT_DIR, "StageB_Regional_Wilcoxon_PerMetric.csv"),
  row.names = FALSE
)

# Additional cross-region head-to-head comparisons requested for CSV output.
# Kept separate so the existing Figure 4 and Figure 5 inference objects remain
# unchanged.
utils::write.csv(
  additional_wilcoxon_results,
  file.path(OUTPUT_DIR, "StageB_Regional_Wilcoxon_ConstantTA_vs_Comparators.csv"),
  row.names = FALSE
)

# Method-across-metrics summary CSVs.
utils::write.csv(
  method_summary_long |>
    dplyr::mutate(
      Method = as.character(Method),
      Metric = as.character(Metric)
    ) |>
    dplyr::select(Method, Metric, Sweep_Count, N_Regions),
  file.path(OUTPUT_DIR, "StageB_Method_Summary_Long.csv"),
  row.names = FALSE
)
utils::write.csv(
  method_summary_agg |>
    dplyr::mutate(Method = as.character(Method)) |>
    dplyr::select(Method, Mean_Sweep, SD_Sweep, Z_Sweep, Sweep_Class),
  file.path(OUTPUT_DIR, "StageB_Method_Summary_Aggregate.csv"),
  row.names = FALSE
)

# Bootstrap replicates (long; one row per REGION x .rep x Method).
all_reps_df <- dplyr::bind_rows(
  lapply(names(boot_results), function(reg) {
    reps <- boot_results[[reg]]$replicates
    reps$REGION <- reg
    reps
  })
) |>
  dplyr::select(REGION, .rep, Method,
                TAM, N_True_Alarms, Sensitivity,
                Mean_Lead_Time, WP, PPV, ALY, N_False_Alarms)
utils::write.csv(
  all_reps_df,
  file.path(OUTPUT_DIR, "StageB_Regional_Bootstrap_Replicates.csv"),
  row.names = FALSE
)
cat("Saved bootstrap replicates: ",
    file.path(OUTPUT_DIR, "StageB_Regional_Bootstrap_Replicates.csv"),
    "\n", sep = "")

# ==============================================================================
# FIGURE 5 - DETECTOR MAP, METRIC TABLE, AND DOT PLOT
# ==============================================================================
# Composes a single composite figure summarising regional dominance:
#   a. Choropleth detector map of the 17 administrative regions, encoding the
#      consensus winner under a four-tier classification (strong, partial,
#      lead_only, contested) with Bonferroni-adjusted head-to-head testing.
#   b. Per-region metric table with embedded per-metric significance below
#      each cell value.
#   c. Per-detector dot plot showing the bootstrap dominance probabilities
#      for each region and detector.
#
# Per-region per-metric significance is also computed for six ORDERED
# detector pairs (each ordered direction is its own one-sided test of the
# explicit signed difference X - Y). The full 6 x 5 x 17 table is saved as
# StageB_Detector_Map_PanelD_HeatmapData.csv for downstream reporting; it
# is NOT rendered as a figure panel.
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# 19. REGION CANONICALISER
# ------------------------------------------------------------------------------
.normalize_key <- function(x) {
  s <- toupper(trimws(as.character(x)))
  s <- gsub("[().,]", " ", s)
  s <- gsub("\\s+", " ", s)
  trimws(s)
}

region_recoder <- c(
  "NATIONAL CAPITAL REGION" = "NCR", "METROPOLITAN MANILA" = "NCR",
  "METRO MANILA" = "NCR", "MANILA" = "NCR", "NCR" = "NCR",
  "CORDILLERA ADMINISTRATIVE REGION" = "CAR",
  "CORDILLERA ADMINISTRATIVE REGION CAR" = "CAR",
  "CORDILLERA" = "CAR", "CAR" = "CAR",
  "ILOCOS REGION" = "REGION I", "ILOCOS" = "REGION I",
  "REGION 1" = "REGION I", "REGION I" = "REGION I",
  "CAGAYAN VALLEY" = "REGION II", "REGION 2" = "REGION II",
  "REGION II" = "REGION II",
  "CENTRAL LUZON" = "REGION III", "REGION 3" = "REGION III",
  "REGION III" = "REGION III",
  "CALABARZON" = "REGION IV-A", "REGION 4-A" = "REGION IV-A",
  "REGION 4A" = "REGION IV-A", "REGION IV-A" = "REGION IV-A",
  "REGION IVA" = "REGION IV-A", "REGION IV A" = "REGION IV-A",
  "MIMAROPA REGION" = "MIMAROPA", "MIMAROPA" = "MIMAROPA",
  "REGION 4-B" = "MIMAROPA", "REGION 4B" = "MIMAROPA",
  "REGION IV-B" = "MIMAROPA", "REGION IV B" = "MIMAROPA",
  "SOUTHWESTERN TAGALOG REGION" = "MIMAROPA",
  "BICOL REGION" = "REGION V", "BICOL" = "REGION V",
  "REGION 5" = "REGION V", "REGION V" = "REGION V",
  "WESTERN VISAYAS" = "REGION VI", "REGION 6" = "REGION VI",
  "REGION VI" = "REGION VI",
  "CENTRAL VISAYAS" = "REGION VII", "REGION 7" = "REGION VII",
  "REGION VII" = "REGION VII",
  "EASTERN VISAYAS" = "REGION VIII", "REGION 8" = "REGION VIII",
  "REGION VIII" = "REGION VIII",
  "ZAMBOANGA PENINSULA" = "REGION IX", "REGION 9" = "REGION IX",
  "REGION IX" = "REGION IX",
  "NORTHERN MINDANAO" = "REGION X", "REGION 10" = "REGION X",
  "REGION X" = "REGION X",
  "DAVAO REGION" = "REGION XI", "DAVAO" = "REGION XI",
  "REGION 11" = "REGION XI", "REGION XI" = "REGION XI",
  "SOCCSKSARGEN" = "REGION XII", "REGION 12" = "REGION XII",
  "REGION XII" = "REGION XII",
  "CARAGA" = "REGION XIII", "CARAGA REGION" = "REGION XIII",
  "REGION 13" = "REGION XIII", "REGION XIII" = "REGION XIII",
  "AUTONOMOUS REGION IN MUSLIM MINDANAO" = "BARMM",
  "AUTONOMOUS REGION OF MUSLIM MINDANAO" = "BARMM",
  "ARMM" = "BARMM", "BANGSAMORO" = "BARMM",
  "BANGSAMORO AUTONOMOUS REGION IN MUSLIM MINDANAO" = "BARMM",
  "BARMM" = "BARMM"
)

# Province-to-region lookup for GADM/Natural Earth dissolves
province_to_region <- c(
  # CAR
  "ABRA" = "CAR", "APAYAO" = "CAR", "BENGUET" = "CAR", "IFUGAO" = "CAR",
  "KALINGA" = "CAR", "MOUNTAIN PROVINCE" = "CAR", "BAGUIO" = "CAR",
  "BAGUIO CITY" = "CAR", "CITY OF BAGUIO" = "CAR",
  # Region I
  "ILOCOS NORTE" = "REGION I", "ILOCOS SUR" = "REGION I",
  "LA UNION" = "REGION I", "PANGASINAN" = "REGION I",
  # Region II
  "BATANES" = "REGION II", "CAGAYAN" = "REGION II", "ISABELA" = "REGION II",
  "NUEVA VIZCAYA" = "REGION II", "QUIRINO" = "REGION II",
  # Region III
  "AURORA" = "REGION III", "BATAAN" = "REGION III", "BULACAN" = "REGION III",
  "NUEVA ECIJA" = "REGION III", "PAMPANGA" = "REGION III",
  "TARLAC" = "REGION III", "ZAMBALES" = "REGION III",
  # Region IV-A
  "BATANGAS" = "REGION IV-A", "CAVITE" = "REGION IV-A",
  "LAGUNA" = "REGION IV-A", "QUEZON" = "REGION IV-A",
  "RIZAL" = "REGION IV-A",
  # MIMAROPA
  "MARINDUQUE" = "MIMAROPA", "OCCIDENTAL MINDORO" = "MIMAROPA",
  "ORIENTAL MINDORO" = "MIMAROPA", "MINDORO OCCIDENTAL" = "MIMAROPA",
  "MINDORO ORIENTAL" = "MIMAROPA", "PALAWAN" = "MIMAROPA",
  "ROMBLON" = "MIMAROPA",
  # Region V
  "ALBAY" = "REGION V", "CAMARINES NORTE" = "REGION V",
  "CAMARINES SUR" = "REGION V", "CATANDUANES" = "REGION V",
  "MASBATE" = "REGION V", "SORSOGON" = "REGION V",
  # Region VI
  "AKLAN" = "REGION VI", "ANTIQUE" = "REGION VI", "CAPIZ" = "REGION VI",
  "GUIMARAS" = "REGION VI", "ILOILO" = "REGION VI",
  "NEGROS OCCIDENTAL" = "REGION VI",
  # Region VII
  "BOHOL" = "REGION VII", "CEBU" = "REGION VII",
  "NEGROS ORIENTAL" = "REGION VII", "SIQUIJOR" = "REGION VII",
  # Region VIII
  "BILIRAN" = "REGION VIII", "EASTERN SAMAR" = "REGION VIII",
  "LEYTE" = "REGION VIII", "NORTHERN SAMAR" = "REGION VIII",
  "SAMAR" = "REGION VIII", "WESTERN SAMAR" = "REGION VIII",
  "SAMAR WESTERN" = "REGION VIII", "SOUTHERN LEYTE" = "REGION VIII",
  # Region IX
  "ZAMBOANGA DEL NORTE" = "REGION IX", "ZAMBOANGA DEL SUR" = "REGION IX",
  "ZAMBOANGA SIBUGAY" = "REGION IX", "CITY OF ISABELA" = "REGION IX",
  "ISABELA CITY" = "REGION IX",
  # Region X
  "BUKIDNON" = "REGION X", "CAMIGUIN" = "REGION X",
  "LANAO DEL NORTE" = "REGION X", "MISAMIS OCCIDENTAL" = "REGION X",
  "MISAMIS ORIENTAL" = "REGION X",
  # Region XI
  "COMPOSTELA VALLEY" = "REGION XI", "DAVAO DE ORO" = "REGION XI",
  "DAVAO DEL NORTE" = "REGION XI", "DAVAO DEL SUR" = "REGION XI",
  "DAVAO OCCIDENTAL" = "REGION XI", "DAVAO ORIENTAL" = "REGION XI",
  # Region XII
  "COTABATO" = "REGION XII", "NORTH COTABATO" = "REGION XII",
  "COTABATO NORTH" = "REGION XII", "SARANGANI" = "REGION XII",
  "SOUTH COTABATO" = "REGION XII", "COTABATO SOUTH" = "REGION XII",
  "SULTAN KUDARAT" = "REGION XII",
  # Region XIII
  "AGUSAN DEL NORTE" = "REGION XIII", "AGUSAN DEL SUR" = "REGION XIII",
  "DINAGAT ISLANDS" = "REGION XIII", "ISLANDS OF DINAGAT" = "REGION XIII",
  "SURIGAO DEL NORTE" = "REGION XIII", "SURIGAO DEL SUR" = "REGION XIII",
  # BARMM
  "BASILAN" = "BARMM", "LANAO DEL SUR" = "BARMM", "MAGUINDANAO" = "BARMM",
  "MAGUINDANAO DEL NORTE" = "BARMM", "MAGUINDANAO DEL SUR" = "BARMM",
  "SULU" = "BARMM", "TAWI-TAWI" = "BARMM", "TAWI TAWI" = "BARMM",
  # NCR cities
  "QUEZON CITY" = "NCR", "CALOOCAN" = "NCR",
  "LAS PINAS" = "NCR", "LAS PINAS CITY" = "NCR", "MAKATI" = "NCR",
  "MAKATI CITY" = "NCR", "MALABON" = "NCR", "MANDALUYONG" = "NCR",
  "MARIKINA" = "NCR", "MUNTINLUPA" = "NCR", "NAVOTAS" = "NCR",
  "PARANAQUE" = "NCR", "PASAY" = "NCR", "PASIG" = "NCR",
  "SAN JUAN" = "NCR", "TAGUIG" = "NCR", "VALENZUELA" = "NCR",
  "PATEROS" = "NCR"
)

canonical_region <- function(x) {
  s   <- .normalize_key(x)
  out <- unname(region_recoder[s])
  is_unmatched <- is.na(out) | !(out %in% CANONICAL_17)
  if (any(is_unmatched)) {
    prov_lookup <- unname(province_to_region[s[is_unmatched]])
    out[is_unmatched] <- prov_lookup
  }
  out[is.na(out)] <- s[is.na(out)]
  out
}

# Apply canonicaliser to dominance and metrics tables before joining map.
dominance_df_canonical <- dominance_df |>
  dplyr::mutate(REGION = canonical_region(REGION))
regional_metrics_canonical <- regional_metrics |>
  dplyr::mutate(REGION = canonical_region(REGION))
boot_ci_canonical <- boot_ci |>
  dplyr::mutate(REGION = canonical_region(REGION))

# ------------------------------------------------------------------------------
# 20. WINNER-ROW TABLE AND HEAD-TO-HEAD CONSENSUS TIER CLASSIFICATION
# ------------------------------------------------------------------------------
# REVISION: generalised from the original three-pair round robin (every
# detector tested against both others) to the two comparisons instructed:
#   Constant Transmission Acceleration    vs  Outbreak Threshold
#   Continuous Transmission Acceleration  vs  Outbreak Threshold
# Constant TA and Continuous TA are no longer compared directly against
# each other. This makes the number of "available pairs" asymmetric across
# detectors: each TA detector has exactly one defined pair (vs Outbreak
# Threshold), while Outbreak Threshold has two (vs each TA). The tier
# definitions below generalise directly from the original scheme to that
# asymmetry, rather than assuming every detector has the same pair count:
#   strong    : Observed_Winner has a significant win in EVERY pair it is
#               defined in, and no significant loss. For a TA detector
#               (one defined pair) this means winning that one pair; for
#               Outbreak Threshold (two defined pairs) this means winning
#               both.
#   partial   : Observed_Winner has at least one significant win but not
#               all of its defined pairs (only reachable for Outbreak
#               Threshold, which is the only detector with more than one
#               defined pair).
#   lead_only : Observed_Winner has zero significant wins and zero
#               significant losses, and its dominance probability is
#               >= 0.50.
#   contested : a rival is significant against the Observed_Winner in any
#               of its defined pairs, or (zero wins/losses and dominance
#               probability < 0.50).
# See Item 8 of the Delphi review.
# ------------------------------------------------------------------------------

winner_rows <- dominance_df_canonical |>
  dplyr::select(REGION, Observed_Winner, Dominance_Probability,
                dplyr::all_of(unname(TARGET_P_COLS))) |>
  dplyr::left_join(
    regional_metrics_canonical |>
      dplyr::select(REGION, Method, TAM, N_True_Alarms,
                    Sensitivity, Mean_Lead_Time, WP) |>
      dplyr::mutate(Method = as.character(Method)),
    by = c("REGION" = "REGION", "Observed_Winner" = "Method")
  )

n_missing_metric <- sum(is.na(winner_rows$TAM))
if (n_missing_metric > 0) {
  warning("Could not match metric values for ", n_missing_metric, " region(s).")
}

# Pair test (chi-squared with Fisher fallback); unchanged from the original.
.pair_test <- function(name_X, p_X, name_Y, p_Y, N) {
  if (any(is.na(c(p_X, p_Y, N))) || N <= 0L)
    return(list(p_one_sided = NA_real_, pair_winner = NA_character_))
  n_X <- as.integer(round(p_X * N))
  n_Y <- as.integer(round(p_Y * N))
  if (n_X == n_Y) return(list(p_one_sided = 1.0, pair_winner = NA_character_))
  if (n_X > n_Y) {
    n_higher <- n_X; n_lower <- n_Y; pair_winner <- name_X
  } else {
    n_higher <- n_Y; n_lower <- n_X; pair_winner <- name_Y
  }
  tab <- matrix(c(n_higher, N - n_higher, n_lower, N - n_lower),
                nrow = 2L, byrow = TRUE)
  expected_cells <- as.numeric(rowSums(tab) %o% colSums(tab) / sum(tab))
  use_fisher <- any(expected_cells < 5)
  if (use_fisher) {
    p_two <- tryCatch(stats::fisher.test(tab)$p.value,
                      error = function(e) NA_real_)
  } else {
    p_two <- tryCatch(
      suppressWarnings(stats::chisq.test(tab, correct = FALSE))$p.value,
      error = function(e) NA_real_)
  }
  if (is.na(p_two)) return(list(p_one_sided = NA_real_,
                                pair_winner = NA_character_))
  p_one <- pmin(pmax(p_two / 2, .Machine$double.eps), 1)
  list(p_one_sided = p_one, pair_winner = pair_winner)
}

N_PAIRS  <- length(DETECTOR_PAIRS)  # 2: ConstantTA-vs-ET, ContinuousTA-vs-ET
n_regions <- nrow(winner_rows)
total_corr_strict <- N_PAIRS * n_regions
total_corr_within <- N_PAIRS

p_const_vs_ot_raw <- numeric(n_regions)
p_cont_vs_ot_raw  <- numeric(n_regions)
w_const_vs_ot     <- character(n_regions)
w_cont_vs_ot      <- character(n_regions)
for (i in seq_len(n_regions)) {
  r  <- winner_rows[i, , drop = FALSE]
  t1 <- .pair_test("ConstantTA",   r$P_ConstantTA,
                   "OutbreakThr",  r$P_OutbreakThreshold, N_BOOTS)
  t2 <- .pair_test("ContinuousTA", r$P_ContinuousTA,
                   "OutbreakThr",  r$P_OutbreakThreshold, N_BOOTS)
  p_const_vs_ot_raw[i] <- t1$p_one_sided
  p_cont_vs_ot_raw[i]  <- t2$p_one_sided
  w_const_vs_ot[i]     <- t1$pair_winner
  w_cont_vs_ot[i]      <- t2$pair_winner
}

# Helper: long detector names -> short codes.
.detector_short <- function(d) {
  dplyr::case_when(
    d == "Constant Transmission Acceleration"   ~ "ConstantTA",
    d == "Continuous Transmission Acceleration" ~ "ContinuousTA",
    d == "Outbreak Threshold"                   ~ "OutbreakThr",
    d == "Alarm Threshold"                      ~ "AlarmThr",
    d == "WHO 75th Percentile Threshold"          ~ "WHO75",
    d == "WHO 90th Percentile Threshold"          ~ "WHO90",
    TRUE                                        ~ NA_character_
  )
}

winner_rows <- winner_rows |>
  dplyr::mutate(
    p_const_vs_ot_bonf        = pmin(p_const_vs_ot_raw * total_corr_strict, 1),
    p_cont_vs_ot_bonf         = pmin(p_cont_vs_ot_raw  * total_corr_strict, 1),
    p_const_vs_ot_bonf_within = pmin(p_const_vs_ot_raw * total_corr_within, 1),
    p_cont_vs_ot_bonf_within  = pmin(p_cont_vs_ot_raw  * total_corr_within, 1),
    pair_winner_const_vs_ot = w_const_vs_ot,
    pair_winner_cont_vs_ot  = w_cont_vs_ot,
    sig_const_vs_ot = !is.na(p_const_vs_ot_bonf) & p_const_vs_ot_bonf < SIG_LEVEL,
    sig_cont_vs_ot  = !is.na(p_cont_vs_ot_bonf)  & p_cont_vs_ot_bonf  < SIG_LEVEL,
    sig_const_vs_ot_within = !is.na(p_const_vs_ot_bonf_within) &
      p_const_vs_ot_bonf_within < SIG_LEVEL,
    sig_cont_vs_ot_within  = !is.na(p_cont_vs_ot_bonf_within) &
      p_cont_vs_ot_bonf_within  < SIG_LEVEL,
    OW_short = .detector_short(Observed_Winner),
    # Number of significant wins/losses among the Observed_Winner's DEFINED
    # pairs only (1 for a TA detector, 2 for Outbreak Threshold, 0 for
    # Alarm Threshold: DETECTOR_PAIRS deliberately does not test Alarm
    # Threshold head-to-head against anything, per the two comparisons
    # instructed, so an Alarm-Threshold-led region always falls through to
    # "contested" below rather than being silently mis-scored).
    OW_n_defined_pairs = dplyr::case_when(
      OW_short %in% c("ConstantTA", "ContinuousTA") ~ 1L,
      OW_short == "OutbreakThr"                     ~ 2L,
      OW_short == "AlarmThr"                        ~ 0L,
      TRUE                                          ~ NA_integer_
    ),
    OW_n_sig_wins = dplyr::case_when(
      OW_short == "ConstantTA"   ~ as.integer(sig_const_vs_ot & pair_winner_const_vs_ot == "ConstantTA"),
      OW_short == "ContinuousTA" ~ as.integer(sig_cont_vs_ot  & pair_winner_cont_vs_ot  == "ContinuousTA"),
      OW_short == "OutbreakThr"  ~ as.integer(sig_const_vs_ot & pair_winner_const_vs_ot == "OutbreakThr") +
                                   as.integer(sig_cont_vs_ot  & pair_winner_cont_vs_ot  == "OutbreakThr"),
      OW_short == "AlarmThr"     ~ 0L,
      TRUE ~ NA_integer_
    ),
    OW_n_sig_losses = dplyr::case_when(
      OW_short == "ConstantTA"   ~ as.integer(sig_const_vs_ot & pair_winner_const_vs_ot == "OutbreakThr"),
      OW_short == "ContinuousTA" ~ as.integer(sig_cont_vs_ot  & pair_winner_cont_vs_ot  == "OutbreakThr"),
      OW_short == "OutbreakThr"  ~ as.integer(sig_const_vs_ot & pair_winner_const_vs_ot == "ConstantTA") +
                                   as.integer(sig_cont_vs_ot  & pair_winner_cont_vs_ot  == "ContinuousTA"),
      OW_short == "AlarmThr"     ~ 0L,
      TRUE ~ NA_integer_
    ),
    consensus_R1_pass = !is.na(Observed_Winner) & nzchar(Observed_Winner),
    consensus_tier = dplyr::case_when(
      !consensus_R1_pass ~ "contested",
      OW_n_sig_losses > 0L ~ "contested",
      OW_n_sig_wins == OW_n_defined_pairs & OW_n_sig_wins > 0L ~ "strong",
      OW_n_sig_wins > 0L & OW_n_sig_wins < OW_n_defined_pairs  ~ "partial",
      OW_n_sig_wins == 0L & OW_n_sig_losses == 0L &
        !is.na(Dominance_Probability) & Dominance_Probability >= 0.50 ~ "lead_only",
      TRUE ~ "contested"
    ),
    # Explicit compatibility field retained for the downstream table builder.
    # In the six-detector implementation a consensus pass means the observed
    # winner has at least one statistically supported primary head-to-head win
    # and no significant loss (the strong or partial tiers).
    consensus_pass = consensus_tier %in% c("strong", "partial"),
    # Dominance-significance p-value used only for display/export. For a TA
    # winner this is its Bonferroni-adjusted comparison with Outbreak Threshold.
    # For an Outbreak-Threshold winner, strong consensus uses the weakest-link
    # (larger) p-value across its two required wins; partial consensus uses the
    # significant winning comparison. Alarm/WHO winners have no prespecified
    # primary head-to-head test, so this remains NA rather than inventing one.
    p_consensus = dplyr::case_when(
      OW_short == "ConstantTA" &
        pair_winner_const_vs_ot == "ConstantTA" ~ p_const_vs_ot_bonf,
      OW_short == "ContinuousTA" &
        pair_winner_cont_vs_ot == "ContinuousTA" ~ p_cont_vs_ot_bonf,
      OW_short == "OutbreakThr" & OW_n_sig_wins >= 2L ~
        pmax(p_const_vs_ot_bonf, p_cont_vs_ot_bonf),
      OW_short == "OutbreakThr" & OW_n_sig_wins == 1L &
        pair_winner_const_vs_ot == "OutbreakThr" ~ p_const_vs_ot_bonf,
      OW_short == "OutbreakThr" & OW_n_sig_wins == 1L &
        pair_winner_cont_vs_ot == "OutbreakThr" ~ p_cont_vs_ot_bonf,
      TRUE ~ NA_real_
    ),
    sig_star_dominance = dplyr::case_when(
      is.na(p_consensus)  ~ "na",
      p_consensus < 0.001 ~ "***",
      p_consensus < 0.01  ~ "**",
      p_consensus < 0.05  ~ "*",
      TRUE                ~ "ns"
    ),
    Leader_Label_short = dplyr::case_when(
      Observed_Winner == "Constant Transmission Acceleration"   ~ "Constant TA",
      Observed_Winner == "Continuous Transmission Acceleration" ~ "Continuous TA",
      Observed_Winner == "Outbreak Threshold"                   ~ "Outbreak Threshold",
      Observed_Winner == "Alarm Threshold"                      ~ "Alarm Threshold",
      Observed_Winner == "WHO 75th Percentile Threshold"          ~ "WHO 75th",
      Observed_Winner == "WHO 90th Percentile Threshold"          ~ "WHO 90th",
      TRUE                                                       ~ NA_character_
    ),
    consensus_winner = dplyr::case_when(
      consensus_tier == "strong"    ~ Observed_Winner,
      consensus_tier == "partial"   ~ Observed_Winner,
      consensus_tier == "lead_only" ~ "No consensus",
      TRUE                          ~ "No consensus"
    )
  ) |>
  dplyr::select(-OW_short)

n_strong    <- sum(winner_rows$consensus_tier == "strong",    na.rm = TRUE)
n_partial   <- sum(winner_rows$consensus_tier == "partial",   na.rm = TRUE)
n_lead_only <- sum(winner_rows$consensus_tier == "lead_only", na.rm = TRUE)
n_contested <- sum(winner_rows$consensus_tier == "contested", na.rm = TRUE)

n_strong_within <- sum(
  with(winner_rows, {
    OW_within_wins <- dplyr::case_when(
      Observed_Winner == "Constant Transmission Acceleration"   ~
        as.integer(sig_const_vs_ot_within & pair_winner_const_vs_ot == "ConstantTA"),
      Observed_Winner == "Continuous Transmission Acceleration" ~
        as.integer(sig_cont_vs_ot_within  & pair_winner_cont_vs_ot  == "ContinuousTA"),
      Observed_Winner == "Outbreak Threshold" ~
        as.integer(sig_const_vs_ot_within & pair_winner_const_vs_ot == "OutbreakThr") +
        as.integer(sig_cont_vs_ot_within  & pair_winner_cont_vs_ot  == "OutbreakThr"),
      TRUE ~ NA_integer_
    )
    OW_within_wins == OW_n_defined_pairs & OW_within_wins > 0L
  }),
  na.rm = TRUE
)

cat(sprintf(
  "\nConsensus tiers: strong=%d, partial=%d, lead_only=%d, contested=%d (of %d).\n",
  n_strong, n_partial, n_lead_only, n_contested, n_regions))
cat(sprintf(
  "Within-region-only Bonferroni sensitivity: %d regions meet 'strong' (vs %d under strict).\n",
  n_strong_within, n_strong))

# ------------------------------------------------------------------------------
# 21. PER-REGION PER-METRIC SIGNIFICANCE (four ordered pairs)
# ------------------------------------------------------------------------------
# REVISION: trimmed from three unordered comparisons (six ordered pairs) to
# the two instructed (Constant TA vs Outbreak Threshold, Continuous TA vs
# Outbreak Threshold; four ordered pairs). Each ordered pair X vs Y is its
# own one-sided test of H1: X > Y on the explicit signed difference X - Y.
# The two ordered directions of an unordered comparison give DIFFERENT
# p-values that sum to ~1.
#
# Test method: paired bootstrap on year-cluster replicates (rigorous).
# Bonferroni: factor 5 within (REGION x ordered_pair). Each ordered pair is
# its own family of 5 metrics; the four ordered pairs are four separate
# families.
# ------------------------------------------------------------------------------

# REVISION: trimmed from three unordered pairs (six ordered) to the two
# instructed comparisons (four ordered), both against Outbreak Threshold.
# Constant TA vs Continuous TA is no longer computed. See Item 8 of the
# Delphi review.
DETECTOR_PAIRS_FOR_MAP <- list(
  list(A = "Constant Transmission Acceleration",
       B = "Outbreak Threshold",
       short = "Constant TA vs OT",
       code  = "const_vs_ot"),
  list(A = "Continuous Transmission Acceleration",
       B = "Outbreak Threshold",
       short = "Continuous TA vs OT",
       code  = "cont_vs_ot")
)

ORDERED_PAIRS <- list(
  list(label = "Constant TA vs OT",
       X_full = "Constant Transmission Acceleration",
       Y_full = "Outbreak Threshold",
       X_short = "Constant TA", Y_short = "ET",
       code = "const_vs_ot", sym_code = "ot_vs_const"),
  list(label = "Continuous TA vs OT",
       X_full = "Continuous Transmission Acceleration",
       Y_full = "Outbreak Threshold",
       X_short = "Continuous TA", Y_short = "ET",
       code = "cont_vs_ot", sym_code = "ot_vs_cont"),
  list(label = "Outbreak Threshold vs Constant TA",
       X_full = "Outbreak Threshold",
       Y_full = "Constant Transmission Acceleration",
       X_short = "ET", Y_short = "Constant TA",
       code = "ot_vs_const", sym_code = "const_vs_ot"),
  list(label = "Outbreak Threshold vs Continuous TA",
       X_full = "Outbreak Threshold",
       Y_full = "Continuous Transmission Acceleration",
       X_short = "ET", Y_short = "Continuous TA",
       code = "ot_vs_cont", sym_code = "cont_vs_ot")
)

CI_LO_HI <- list(
  TAM            = c("TAM_lo",            "TAM_hi"),
  N_True_Alarms  = c("N_True_Alarms_lo",  "N_True_Alarms_hi"),
  Sensitivity    = c("Sens_lo",           "Sens_hi"),
  Mean_Lead_Time = c("MLT_lo",            "MLT_hi"),
  WP             = c("WP_lo",             "WP_hi")
)

# Build long-format replicates frame (with canonical region names) for the
# paired-bootstrap test.
all_reps_canonical <- dplyr::bind_rows(
  lapply(names(boot_results), function(reg) {
    reps <- boot_results[[reg]]$replicates
    reps$REGION <- canonical_region(reg)
    reps
  })
)

compute_level1_paired_bootstrap_ordered <- function(reps_df, regions,
                                                    ordered_pairs, metrics) {
  out <- list()
  for (reg in regions) {
    sub <- reps_df |> dplyr::filter(REGION == reg)
    if (nrow(sub) == 0L) next
    for (op in ordered_pairs) {
      x_df <- sub |> dplyr::filter(Method == op$X_full) |> dplyr::arrange(.rep)
      y_df <- sub |> dplyr::filter(Method == op$Y_full) |> dplyr::arrange(.rep)
      if (nrow(x_df) == 0L || nrow(y_df) == 0L) next
      n_common <- min(nrow(x_df), nrow(y_df))
      x_df <- x_df[seq_len(n_common), , drop = FALSE]
      y_df <- y_df[seq_len(n_common), , drop = FALSE]
      for (m in metrics) {
        diffs <- x_df[[m]] - y_df[[m]]
        diffs <- diffs[is.finite(diffs)]
        if (length(diffs) < 50L) {
          out[[length(out) + 1L]] <- data.frame(
            REGION = reg, Ordered_Pair = op$label, Pair_Code = op$code,
            Sym_Pair_Code = op$sym_code,
            Detector_X = op$X_full, Detector_Y = op$Y_full,
            Metric = m, n_used = length(diffs),
            prob_X_better = NA_real_, median_diff_XY = NA_real_,
            p_one_sided = NA_real_, p_two_sided = NA_real_,
            method = "bootstrap_paired", stringsAsFactors = FALSE)
          next
        }
        prX     <- mean(diffs > 0)
        med     <- stats::median(diffs)
        p_floor <- 1 / length(diffs)
        p_one   <- max(1 - prX, p_floor)
        p_two   <- max(2 * min(prX, 1 - prX), p_floor)
        out[[length(out) + 1L]] <- data.frame(
          REGION = reg, Ordered_Pair = op$label, Pair_Code = op$code,
          Sym_Pair_Code = op$sym_code,
          Detector_X = op$X_full, Detector_Y = op$Y_full,
          Metric = m, n_used = length(diffs),
          prob_X_better = prX, median_diff_XY = med,
          p_one_sided = p_one, p_two_sided = p_two,
          method = "bootstrap_paired", stringsAsFactors = FALSE)
      }
    }
  }
  dplyr::bind_rows(out)
}

regions_in_order <- unique(winner_rows$REGION)

cat("\nPer-region per-metric inference: paired bootstrap, 6 ordered pairs.\n")
level1_results_ordered <- compute_level1_paired_bootstrap_ordered(
  all_reps_canonical, regions_in_order, ORDERED_PAIRS, HH_METRICS_REGIONAL
)
level1_method <- "bootstrap_paired"

# Bonferroni adjustment within (REGION x Pair_Code).
if (nrow(level1_results_ordered) > 0L) {
  level1_results_ordered <- level1_results_ordered |>
    dplyr::group_by(REGION, Pair_Code) |>
    dplyr::mutate(
      p_bonf = pmin(p_one_sided * BONF_LEVEL1_FACTOR, 1),
      sig_star = dplyr::case_when(
        is.na(p_bonf)  ~ "ns",
        p_bonf < 0.001 ~ "***",
        p_bonf < 0.01  ~ "**",
        p_bonf < 0.05  ~ "*",
        TRUE           ~ "ns"
      )
    ) |>
    dplyr::ungroup()
}

cat("=== Per-region per-metric significance (6 ordered pairs x 5 metrics x ",
    n_regions, " regions = ", nrow(level1_results_ordered), " tests) ===\n",
    sep = "")
cat("Method:                     ", level1_method, "\n", sep = "")
cat("Test:                       one-sided H1: X > Y per ordered pair\n")
cat("Bonferroni factor:          ", BONF_LEVEL1_FACTOR,
    " (within REGION x ordered pair)\n", sep = "")
cat("Significant at p<0.05:      ",
    sum(level1_results_ordered$sig_star %in% c("*", "**", "***"), na.rm = TRUE),
    " of ", sum(!is.na(level1_results_ordered$p_bonf)), "\n", sep = "")

# ------------------------------------------------------------------------------
# 22. CROSS-REGION PER-METRIC WILCOXON (formatted from in-memory results)
# ------------------------------------------------------------------------------
level2_results <- all_wilcoxon_results |>
  dplyr::transmute(
    Metric, Comparison,
    Detector_A, Detector_B,
    n_pairs     = N_Regions_Paired,
    median_diff = Median_Diff_AminusB,
    p_value     = p_value
  ) |>
  dplyr::group_by(Metric) |>
  dplyr::mutate(
    p_bonf = pmin(p_value * BONF_LEVEL2_FACTOR, 1),
    sig_star = dplyr::case_when(
      is.na(p_bonf)  ~ "ns",
      p_bonf < 0.001 ~ "***",
      p_bonf < 0.01  ~ "**",
      p_bonf < 0.05  ~ "*",
      TRUE           ~ "ns"
    )
  ) |>
  dplyr::ungroup()

# Per-metric significance for the consensus winner (weakest-link rule among
# the two ordered pairs the winner participates in as X).
build_per_metric_sig_for_region <- function(reg, consensus_winner) {
  sub <- level1_results_ordered |> dplyr::filter(REGION == reg)
  if (nrow(sub) == 0L) {
    return(data.frame(REGION = reg, Metric = HH_METRICS_REGIONAL,
                      sig_star = "na", p_bonf = NA_real_,
                      reference_pair = NA_character_,
                      stringsAsFactors = FALSE))
  }
  rows <- list()
  for (m in HH_METRICS_REGIONAL) {
    sub_m <- sub |> dplyr::filter(Metric == m)
    if (nrow(sub_m) == 0L) {
      rows[[length(rows) + 1L]] <- data.frame(
        REGION = reg, Metric = m, sig_star = "na",
        p_bonf = NA_real_, reference_pair = NA_character_,
        stringsAsFactors = FALSE)
      next
    }
    if (consensus_winner == "Constant Transmission Acceleration") {
      ref_codes <- c("const_vs_ot")
    } else if (consensus_winner == "Continuous Transmission Acceleration") {
      ref_codes <- c("cont_vs_ot")
    } else if (consensus_winner == "Outbreak Threshold") {
      ref_codes <- c("ot_vs_const", "ot_vs_cont")
    } else {
      ref_codes <- unique(sub_m$Pair_Code)
    }
    sub_m_ref <- sub_m |> dplyr::filter(Pair_Code %in% ref_codes)
    if (nrow(sub_m_ref) == 0L) {
      rows[[length(rows) + 1L]] <- data.frame(
        REGION = reg, Metric = m, sig_star = "na",
        p_bonf = NA_real_, reference_pair = NA_character_,
        stringsAsFactors = FALSE)
      next
    }
    if (consensus_winner %in% c("Constant Transmission Acceleration",
                                "Continuous Transmission Acceleration",
                                "Outbreak Threshold")) {
      idx_pick <- which.max(sub_m_ref$p_bonf)[1]
    } else {
      idx_pick <- which.min(sub_m_ref$p_bonf)[1]
    }
    pick <- sub_m_ref[idx_pick, , drop = FALSE]
    rows[[length(rows) + 1L]] <- data.frame(
      REGION = reg, Metric = m,
      sig_star = pick$sig_star, p_bonf = pick$p_bonf,
      reference_pair = pick$Ordered_Pair, stringsAsFactors = FALSE)
  }
  dplyr::bind_rows(rows)
}

per_region_per_metric_sig <- dplyr::bind_rows(lapply(seq_len(nrow(winner_rows)),
  function(i) {
    build_per_metric_sig_for_region(
      winner_rows$REGION[i], winner_rows$consensus_winner[i]
    )
  }))

sig_wide <- per_region_per_metric_sig |>
  dplyr::select(REGION, Metric, sig_star) |>
  tidyr::pivot_wider(names_from = Metric, values_from = sig_star,
                     names_prefix = "sig_") |>
  as.data.frame()

winner_rows <- winner_rows |> dplyr::left_join(sig_wide, by = "REGION")

# Fail before map/table construction if a regional consensus column was lost
# during an upstream rewrite. These are the exact fields consumed by the Figure
# B5 map, display table, dominance panel, CSVs and legend builder.
.required_winner_cols <- c(
  "REGION", "Observed_Winner", "Dominance_Probability",
  "consensus_R1_pass", "consensus_pass", "consensus_tier",
  "consensus_winner", "p_consensus", "sig_star_dominance",
  "Leader_Label_short", "TAM", "N_True_Alarms", "Sensitivity",
  "Mean_Lead_Time", "WP",
  paste0("sig_", HH_METRICS_REGIONAL)
)
.missing_winner_cols <- setdiff(.required_winner_cols, names(winner_rows))
if (length(.missing_winner_cols) > 0L) {
  stop(
    "Regional Figure B5 schema error: winner_rows is missing required column(s): ",
    paste(.missing_winner_cols, collapse = ", "),
    call. = FALSE
  )
}

# ------------------------------------------------------------------------------
# 23. PHILIPPINE GEOMETRY
# ------------------------------------------------------------------------------
# Convert a serialized/local vector object to sf without relying on implicit
# method dispatch. `terra::saveRDS()` can serialize SpatVector objects as
# PackedSpatVector/Packed objects; older terra versions return that packed class
# from readRDS(), and sf::st_as_sf() has no method for it. Unpack first, then
# convert the resulting SpatVector explicitly.
as_sf_geometry_safe <- function(obj, source_label = "geometry object") {
  if (inherits(obj, "sf")) return(obj)

  if (inherits(obj, "PackedSpatVector") || inherits(obj, "Packed")) {
    if (!requireNamespace("terra", quietly = TRUE)) {
      stop("The bundled map geometry is a packed terra object, but package 'terra' is not installed. ",
           "Run source('R/install_dependencies.R') and re-run the pipeline.", call. = FALSE)
    }

    # terra serializes SpatVector objects as PackedSpatVector/Packed when saved
    # with saveRDS(). sf::st_as_sf() has no method for that packed wrapper, so
    # unpack it explicitly before conversion.
    unpacked <- tryCatch(terra::unwrap(obj), error = function(e) e)
    if (inherits(unpacked, "error")) {
      # Compatibility fallback for terra builds where vect() can reconstruct
      # the packed representation.
      unpacked <- tryCatch(terra::vect(obj), error = function(e) e)
    }
    if (inherits(unpacked, "error") || is.null(unpacked)) {
      msg <- if (inherits(unpacked, "error")) conditionMessage(unpacked) else "NULL result"
      stop("Failed to unpack ", source_label, " (class: ",
           paste(class(obj), collapse = "/"), "): ", msg, call. = FALSE)
    }
    obj <- unpacked
  }

  if (inherits(obj, "SpatVector")) {
    converted <- tryCatch(sf::st_as_sf(obj), error = function(e) e)
    if (inherits(converted, "error")) {
      stop("Failed to convert unpacked SpatVector to sf for ", source_label,
           ": ", conditionMessage(converted), call. = FALSE)
    }
    return(converted)
  }

  converted <- tryCatch(sf::st_as_sf(obj), error = function(e) e)
  if (inherits(converted, "error")) {
    stop("Unsupported map geometry class for ", source_label, ": ",
         paste(class(obj), collapse = "/"), ". ",
         conditionMessage(converted), call. = FALSE)
  }
  converted
}

get_ph_regions <- function() {
  g <- NULL

  # 1. Bundled RDS or local vector boundary (preferred: fully offline).
  rds_candidates <- list.files(DIR_GEOMETRY, pattern = "\\.rds$",
                               recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  if (length(rds_candidates) > 0L) {
    message("[map] Geometry source: bundled RDS -> ", basename(rds_candidates[1]))
    obj <- readRDS(rds_candidates[1])
    g <- as_sf_geometry_safe(obj, paste0("bundled RDS ", basename(rds_candidates[1])))
  } else if (!is.null(PH_SHAPEFILE) && file.exists(PH_SHAPEFILE)) {
    message("[map] Geometry source: local file -> ", basename(PH_SHAPEFILE))
    g <- sf::st_read(PH_SHAPEFILE, quiet = TRUE)

  # 2. GADM via geodata. Cached under data/geometry/ rather than tempdir(), so
  #    a second run does not re-download and the archive becomes reproducible.
  } else if (requireNamespace("geodata", quietly = TRUE)) {
    message("[map] Geometry source: GADM via geodata (cached in data/geometry/)")
    for (lvl in c(1L, 2L)) {
      gg <- try(geodata::gadm(country = "PHL", level = lvl,
                              path = DIR_GEOMETRY), silent = TRUE)
      if (!inherits(gg, "try-error") && !is.null(gg)) {
        g <- as_sf_geometry_safe(gg, paste0("GADM level ", lvl))
        break
      }
    }

  # 3. Natural Earth states as a last resort.
  } else if (requireNamespace("rnaturalearth", quietly = TRUE)) {
    message("[map] Geometry source: Natural Earth states")
    gg <- try(rnaturalearth::ne_states(country = "Philippines",
                                       returnclass = "sf"), silent = TRUE)
    if (!inherits(gg, "try-error") && !is.null(gg)) g <- gg
  }

  if (is.null(g)) {
    stop("No Philippines geometry available.\n",
         "  Either place a .shp/.gpkg/.geojson in data/geometry/,\n",
         "  or install 'geodata' (needs internet on first run).",
         call. = FALSE)
  }

  # --- CRS validation ---------------------------------------------------------
  # Sources vary: GADM ships EPSG:4326, some national shapefiles ship PRS92 or
  # a local grid, and a few carry no CRS at all. Normalising to WGS84 here means
  # the later projection to EPSG:3123 is well defined, so the scale bar measures
  # true ground distance rather than degrees.
  if (is.na(sf::st_crs(g))) {
    warning("[map] Geometry has no CRS; assuming EPSG:4326 (WGS84).",
            call. = FALSE)
    sf::st_crs(g) <- CRS_WGS84
  }
  g <- sf::st_transform(g, CRS_WGS84)

  # Repair invalid rings (self-intersections in coastline data are common and
  # otherwise cause st_union/summarise to fail).
  if (!all(sf::st_is_valid(g))) {
    message("[map] Repairing invalid geometries with st_make_valid().")
    g <- sf::st_make_valid(g)
  }
  g
}
ph_geom <- get_ph_regions()

name_candidates <- c("REGION", "REGION_NAM", "Region", "name", "NAME_1",
                     "name_1", "NAME_2", "name_2", "REGNAME",
                     "ADM1_EN", "ADM1_PCODE", "ADM1_NAME",
                     "region", "Region_1")
nm_col <- intersect(name_candidates, names(ph_geom))[1]
if (is.na(nm_col))
  stop("Cannot identify region/province column in geometry.", call. = FALSE)
ph_geom$REGION_RAW <- ph_geom[[nm_col]]
ph_geom$REGION     <- canonical_region(ph_geom$REGION_RAW)
if (length(unique(ph_geom$REGION)) > 17 || anyDuplicated(ph_geom$REGION) > 0) {
  ph_geom <- ph_geom |>
    dplyr::group_by(REGION) |>
    dplyr::summarise(.groups = "drop")
}

# Join audit
audit_df <- dplyr::full_join(
  ph_geom |> sf::st_drop_geometry() |> dplyr::distinct(REGION) |>
    dplyr::mutate(in_geometry = TRUE),
  winner_rows |> dplyr::distinct(REGION) |> dplyr::mutate(in_csv = TRUE),
  by = "REGION"
) |>
  dplyr::mutate(
    in_geometry = !is.na(in_geometry),
    in_csv      = !is.na(in_csv),
    status = dplyr::case_when(
      in_geometry &  in_csv  ~ "matched",
      in_geometry & !in_csv  ~ "geometry only",
      !in_geometry &  in_csv ~ "csv only",
      TRUE                   ~ "neither"
    )
  ) |>
  dplyr::arrange(status, REGION)

map_df <- ph_geom |>
  dplyr::filter(REGION %in% CANONICAL_17) |>
  dplyr::left_join(winner_rows, by = "REGION") |>
  dplyr::mutate(
    # The map encodes the observed composite winner among ALL SIX requested
    # detectors. Pairwise significance remains available in the exported
    # inference tables, but it does not suppress Alarm/WHO winners on the map.
    Detector_Label = dplyr::case_when(
      Observed_Winner == "Constant Transmission Acceleration"   ~ "Constant TA",
      Observed_Winner == "Continuous Transmission Acceleration" ~ "Continuous TA",
      Observed_Winner == "Outbreak Threshold"                   ~ "Outbreak Threshold",
      Observed_Winner == "Alarm Threshold"                      ~ "Alarm Threshold",
      Observed_Winner == "WHO 75th Percentile Threshold"        ~ "WHO 75th",
      Observed_Winner == "WHO 90th Percentile Threshold"        ~ "WHO 90th",
      TRUE                                                       ~ NA_character_
    ),
    Confidence_Tier = dplyr::case_when(
      Dominance_Probability >= CONFIDENCE_HI       ~ "Solid (Pr >= 0.90)",
      Dominance_Probability >= DOMINANCE_THRESHOLD ~ "Mid (Pr >= 0.75)",
      TRUE                                          ~ "Light (Pr < 0.75)"
    ),
    Tier_Fill_Group = Detector_Label,
    # Dashed region boundary flags a weaker (<0.75) bootstrap dominance
    # probability; fill still identifies the winning detector.
    Tier_Linetype = ifelse(Dominance_Probability >= DOMINANCE_THRESHOLD,
                           "solid", "dashed"),
    Tier_Linewidth = ifelse(Dominance_Probability >= DOMINANCE_THRESHOLD,
                            0.30, 0.55)
  )
map_df$Detector_Label <- factor(
  map_df$Detector_Label,
  levels = c("Constant TA", "Continuous TA", "Outbreak Threshold",
             "Alarm Threshold", "WHO 75th", "WHO 90th")
)
map_df$Confidence_Tier <- factor(
  map_df$Confidence_Tier,
  levels = c("Solid (Pr >= 0.90)", "Mid (Pr >= 0.75)", "Light (Pr < 0.75)")
)
map_df$Tier_Fill_Group <- factor(
  map_df$Tier_Fill_Group,
  levels = c("Constant TA", "Continuous TA", "Outbreak Threshold",
             "Alarm Threshold", "WHO 75th", "WHO 90th")
)
map_regions <- unique(as.character(map_df$REGION))
if (!setequal(map_regions, CANONICAL_17) || length(map_regions) != 17L) {
  stop("Regional map geometry/join must contain exactly all 17 canonical regions. Missing: ",
       paste(setdiff(CANONICAL_17, map_regions), collapse = ", "), call. = FALSE)
}
if (any(!is.finite(map_df$Dominance_Probability))) {
  bad <- unique(as.character(map_df$REGION[!is.finite(map_df$Dominance_Probability)]))
  stop("Dominance probability is missing for mapped region(s): ",
       paste(bad, collapse = ", "), call. = FALSE)
}

# ------------------------------------------------------------------------------
# 24. FIGURE 5 PALETTES
# ------------------------------------------------------------------------------
# Okabe-Ito colourblind-safe palette. The previous scheme paired #2CA02C
# (green) with #D62728 (red), which deuteranopic and protanopic readers cannot
# separate; these three hues remain distinct under all common forms of colour
# vision deficiency and stay separable in greyscale.
# Four target detectors, all Okabe-Ito so the set stays colourblind-safe.
detector_palette <- c(
  "Constant TA" = "#0072B2",
  "Continuous TA" = "#E69F00",
  "Outbreak Threshold" = "#009E73",
  "Alarm Threshold" = "#CC79A7",
  "WHO 75th" = "#56B4E9",
  "WHO 90th" = "#D55E00",
  "No consensus" = "#9A9A9A"
)
alpha_palette <- c(
  "Solid (Pr >= 0.90)" = 1.00,
  "Mid (Pr >= 0.75)"   = 0.82,
  "Light (Pr < 0.75)"  = 0.60
)
# Red used for the Continuous TA circles in Figure 5C only (see section 28).
PANEL_C_CONTINUOUS_TA_FILL <- "#D55E00"   # Okabe-Ito vermillion

tier_fill_palette <- c(
  "Constant TA"        = unname(detector_palette["Constant TA"]),
  "Continuous TA"      = unname(detector_palette["Continuous TA"]),
  "Outbreak Threshold" = unname(detector_palette["Outbreak Threshold"]),
  "Alarm Threshold"    = unname(detector_palette["Alarm Threshold"]),
  "WHO 75th"           = unname(detector_palette["WHO 75th"]),
  "WHO 90th"           = unname(detector_palette["WHO 90th"])
)

# ------------------------------------------------------------------------------
# 25. FIGURE 5 PANEL A - CHOROPLETH DETECTOR MAP
# ------------------------------------------------------------------------------
# Geometry is projected to EPSG:3123 (PRS92 / Philippines Zone III) BEFORE
# centroids are computed. Two reasons this matters:
#   * st_point_on_surface() on unprojected lon/lat degrees places labels
#     slightly off for elongated islands; on a projected grid it is correct.
#   * a scale bar over geographic degrees is meaningless, because a degree of
#     longitude shrinks with latitude. On a projected CRS the bar represents a
#     true, constant ground distance.
map_df_proj <- sf::st_transform(map_df, CRS_PH_PROJECTED)

suppressWarnings({
  centroids_sf <- sf::st_point_on_surface(map_df_proj)

  # Plain region names only -- no detector line, no boxes. Dominant detector is
  # already encoded by the polygon fill and its legend, so the second label line
  # was redundant; removing it lets 17 single-line names sit legibly on the map.
  centroids_sf$label <- centroids_sf$REGION
})

# RF-product-specific map: environmental RF qualification of the intended
# disease anchor, not a repetition of the six-detector case-based dominance map.
if (isTRUE(RF_GATE$active) && nrow(RF_TO_ANCHOR_SUMMARY) > 0L) {
  rf_map_df <- map_df_proj |>
    dplyr::select(REGION, geometry) |>
    dplyr::left_join(RF_TO_ANCHOR_SUMMARY, by = "REGION")
  if (nrow(rf_map_df) != 17L || anyDuplicated(rf_map_df$REGION)) {
    stop("RF-to-anchor map must contain exactly 17 unique regions.", call. = FALSE)
  }
  rf_map_title <- paste0(
    unique(RF_TO_ANCHOR_SUMMARY$RF_Representation)[1], " preceding ",
    unique(RF_TO_ANCHOR_SUMMARY$Disease_Anchor)[1])
  p_rf_anchor_map <- ggplot2::ggplot(rf_map_df) +
    ggplot2::geom_sf(ggplot2::aes(fill = RF_Supported_PPV), colour = "grey30", linewidth = 0.30) +
    ggrepel::geom_text_repel(
      data = sf::st_point_on_surface(rf_map_df),
      ggplot2::aes(geometry = geometry, label = REGION), stat = "sf_coordinates",
      size = pub_text_size(PUB_ANNOT), family = PUB_FAMILY, fontface = "bold",
      bg.colour = "white", bg.r = 0.12, max.overlaps = Inf,
      min.segment.length = 0, segment.color = "grey55", seed = GLOBAL_SEED,
      show.legend = FALSE) +
    ggplot2::scale_fill_viridis_c(
      option = "C", limits = c(0, 1), na.value = "grey92",
      name = "RF-supported\ntrue-alarm PPV",
      labels = scales::label_percent(accuracy = 1)) +
    ggplot2::labs(title = rf_map_title, x = "Longitude (°E)", y = "Latitude (°N)") +
    theme_pub_map() +
    ggplot2::theme(
      legend.position = "bottom", legend.direction = "horizontal",
      legend.box = "horizontal", legend.justification = "center",
      legend.key.width = grid::unit(40, "pt"),
      plot.margin = ggplot2::margin(10, 12, 18, 12))
  p_rf_anchor_map <- pub_map_frame(p_rf_anchor_map, crs = CRS_PH_PROJECTED)
  save_plot_file(
    paste0("FigureB7_RF_to_", if (RF_GATE$product %in% c(1L,3L)) "ConstantTA" else "OutbreakThreshold", "_Map"),
    p_rf_anchor_map, width = NC_W_DOUBLE, height = 7.2)
}

p_map <- ggplot2::ggplot() +
  ggplot2::geom_sf(
    data = map_df_proj,
    # alpha is not mapped here: fills are uniform, without
    # bootstrap-confidence shading, so there is nothing extra to explain in
    # the legend.
    ggplot2::aes(fill = Tier_Fill_Group),
    colour = "grey25",
    linetype = map_df_proj$Tier_Linetype,
    linewidth = pmax(map_df_proj$Tier_Linewidth, 0.30)  # keep boundaries above
                                                        # the press hairline floor
  ) +
  # Centroid markers, drawn from the sf geometry so they stay registered with
  # the polygons under any coord_sf transformation.
  ggplot2::geom_sf(
    data = centroids_sf, colour = "grey15", size = 0.6,
    show.legend = FALSE
  ) +
  # geom_TEXT_repel, not geom_label_repel: plain names drawn directly on the map
  # with no rectangle, fill or callout background. A soft white halo (bg.colour)
  # preserves contrast over dark choropleth fills without drawing a box.
  ggrepel::geom_text_repel(
    data = centroids_sf,
    ggplot2::aes(geometry = geometry, label = label),
    stat = "sf_coordinates",
    size = pub_text_size(PUB_ANNOT),   # points -> mm, matching the shared scale
    family = PUB_FAMILY,
    fontface = "bold", colour = "grey10",
    bg.colour = "white", bg.r = 0.14,
    box.padding = grid::unit(0.32, "lines"),
    point.padding = grid::unit(0.14, "lines"),
    force = 5, force_pull = 0.7, max.overlaps = Inf,
    min.segment.length = 0, segment.color = "grey45",
    segment.size = 0.28, seed = GLOBAL_SEED,
    show.legend = FALSE,
    # Confine repelled region names to the map's own bounding box so they
    # cannot be pushed off the canvas and cut.
    xlim = sf::st_bbox(map_df_proj)[c("xmin", "xmax")],
    ylim = sf::st_bbox(map_df_proj)[c("ymin", "ymax")]
  ) +
  # ONE legend only: the detector that dominates each region. Confidence is
  # carried by boundary style rather than a second fill/alpha legend, leaving a
  # single unambiguous detector key.
  ggplot2::scale_fill_manual(
    name   = "Detectors",
    values = tier_fill_palette,
    # Always display the complete six-detector key in canonical order. This
    # keeps the map legend stable across reruns even when a detector does not
    # happen to dominate any region in a particular analysis.
    limits = c("Constant TA", "Continuous TA", "Outbreak Threshold",
               "Alarm Threshold", "WHO 75th", "WHO 90th"),
    breaks = c("Constant TA", "Continuous TA", "Outbreak Threshold",
               "Alarm Threshold", "WHO 75th", "WHO 90th"),
    drop   = FALSE, na.translate = FALSE
  ) +
  ggplot2::guides(
    fill = ggplot2::guide_legend(order = 1, nrow = 2, byrow = TRUE,
                                 direction = "horizontal",
                                 override.aes = list(alpha = 0.9,
                                                     linewidth = 0.3))
  ) +
  ggplot2::labs(title = "a    Regional detector map",
               x = "Longitude (\u00b0E)", y = "Latitude (\u00b0N)") +
  theme_pub_map() +
  ggplot2::theme(
    # geom_sf()/coord_sf() default to literal "x" and "y" axis titles unless
    # overridden above; axis.title is re-enabled here (matching the bold,
    # PUB_AXIS_TIT house style used for the other Figure 5 panel) in case
    # theme_pub_map() blanks axis titles for other map panels.
    axis.title.x    = ggplot2::element_text(size = PUB_AXIS_TIT, face = "bold",
                                             colour = "grey25",
                                             margin = ggplot2::margin(t = 6)),
    axis.title.y    = ggplot2::element_text(size = PUB_AXIS_TIT, face = "bold",
                                             colour = "grey25",
                                             margin = ggplot2::margin(r = 6)),
    legend.position = "bottom",
    legend.box      = "horizontal",
    legend.box.just = "center",
    # EVEN SPACING between "Detectors", the first key and the second key.
    # The title-to-first-key gap is set by legend.title's right margin; the
    # key-to-key gap by legend.key.spacing.x plus the label's right margin.
    # They are matched at 22 pt so the three elements sit evenly apart.
    legend.title         = ggplot2::element_text(
                             size = PUB_LEG_TIT, face = "bold",
                             colour = "black", vjust = 0.5,
                             margin = ggplot2::margin(r = 7)),
    legend.text          = ggplot2::element_text(
                             size = PUB_LEG_TXT, colour = "black",
                             margin = ggplot2::margin(l = 2, r = 7)),
    legend.key.spacing.x = grid::unit(5, "pt"),
    legend.spacing.x     = grid::unit(5, "pt"),
    legend.spacing.y     = grid::unit(2, "pt"),
    legend.margin        = ggplot2::margin(3, 6, 3, 6),
    plot.margin          = ggplot2::margin(8, 12, 16, 12)
  )

# Cartographic furniture requested by the reviewer: scale bar (bottom-left),
# north arrow (top-right) and an explicit projected CRS. Added last so the
# annotations sit above the polygon layers.
p_map <- pub_map_frame(
  p_map,
  crs          = CRS_PH_PROJECTED,
  graticule    = TRUE,
  scalebar_loc = "bl",
  arrow_loc    = "tr"
)

# ------------------------------------------------------------------------------
# 26. PER-REGION METRIC TABLE DATA
# ------------------------------------------------------------------------------
# The dengue reference once built a graphical metric-table panel here. The
# revised leptospirosis Figure B5 contains only the requested detector map and
# six-detector dominance-probability panel, so constructing that unused ggplot
# adds failure surface without producing a deliverable. The complete per-region
# metrics, consensus fields and per-metric significance results are still
# exported below as CSV files from `winner_rows`, `level1_results_ordered` and
# `panel_d_data`.

# ------------------------------------------------------------------------------
# 27. FIGURE 5 PANEL C - PER-DETECTOR DOT PLOT
# ------------------------------------------------------------------------------
# All six target detectors are pivoted from TARGET_P_COLS so the
# dominance-probability figure stays synchronized with TARGET_DETECTORS.
# Maps a full detector name (as stored in Observed_Winner) to its short label.
WINNER_TO_SHORT <- c(
  "Constant Transmission Acceleration" = "Constant TA",
  "Continuous Transmission Acceleration" = "Continuous TA",
  "Outbreak Threshold" = "Outbreak Threshold",
  "Alarm Threshold" = "Alarm Threshold",
  "WHO 75th Percentile Threshold" = "WHO 75th",
  "WHO 90th Percentile Threshold" = "WHO 90th"
)

DETECTOR_SHORT_LABEL <- c(
  P_ConstantTA = "Constant TA",
  P_ContinuousTA = "Continuous TA",
  P_OutbreakThreshold = "Outbreak Threshold",
  P_AlarmThreshold = "Alarm Threshold",
  P_WHO75 = "WHO 75th",
  P_WHO90 = "WHO 90th"
)

panel_c_long <- dominance_df_canonical |>
  dplyr::select(REGION, dplyr::all_of(unname(TARGET_P_COLS))) |>
  tidyr::pivot_longer(cols = dplyr::all_of(unname(TARGET_P_COLS)),
                      names_to = "Detector_Code", values_to = "Dominance_P") |>
  dplyr::mutate(
    Detector_Label = unname(DETECTOR_SHORT_LABEL[Detector_Code]),
    Detector_Label = factor(Detector_Label,
                            levels = unname(DETECTOR_SHORT_LABEL[
                              unname(TARGET_P_COLS)]))
  ) |>
  dplyr::left_join(
    winner_rows |>
      dplyr::select(REGION, consensus_tier, Observed_Winner, Leader_Label_short),
    by = "REGION"
  ) |>
  dplyr::mutate(
    # A dot is the leader when its detector is the observed winner among
    # the six requested target detectors.
    is_leader_dot = !is.na(Observed_Winner) &
      Detector_Label == unname(WINNER_TO_SHORT[Observed_Winner]),
    border_colour = dplyr::case_when(
      is_leader_dot & consensus_tier == "strong"    ~ "#1F1F1F",
      is_leader_dot & consensus_tier == "partial"   ~ "#1F1F1F",
      is_leader_dot & consensus_tier == "lead_only" ~ "grey45",
      TRUE                                          ~ "grey55"
    ),
    border_stroke = dplyr::case_when(
      is_leader_dot & consensus_tier == "strong"    ~ 0.95,
      is_leader_dot & consensus_tier == "partial"   ~ 0.65,
      is_leader_dot & consensus_tier == "lead_only" ~ 0.55,
      TRUE                                          ~ 0.30
    )
  ) |>
  dplyr::arrange(REGION, Detector_Label)

pj_c <- ggplot2::position_jitter(width = 0.18, height = 0, seed = 12345)
chance_p <- 1 / length(TARGET_DETECTORS)

# Significance brackets for the dominance-probability panel: Constant TA versus
# each comparator, paired by region and tested with the same paired
# Wilcoxon signed-rank used everywhere else. ONLY SIGNIFICANT comparisons get a
# bracket and a p-value; non-significant pairs are simply not drawn. Every
# comparison remains in the exported Wilcoxon CSV.
# Full pairwise test table FIRST: every comparison, significant or not, with
# V, p, CI and a Drawn_On_Figure flag. This is the file that supports the
# significance brackets on the figure -- exported below so every p-value shown
# on the panel can be traced to a row.
dom_tests <- dominance_pairwise_tests(
  d         = as.data.frame(panel_c_long),
  unit_col  = "REGION",
  level_col = "Detector_Label",
  value_col = "Dominance_P",
  focal     = "Constant TA")

if (!is.null(dom_tests)) {
  utils::write.csv(dom_tests, file.path(OUTPUT_DIR, "FigureB5_DominanceProbability_Wilcoxon_Results.csv"),
                   row.names = FALSE)
  cat("  Saved: ", file.path(OUTPUT_DIR, "FigureB5_DominanceProbability_Wilcoxon_Results.csv"), "\n", sep = "")
}

# Geometry is built FROM that same table, so figure and CSV cannot disagree.
dom_brackets <- build_sig_brackets(
  d         = as.data.frame(panel_c_long),
  unit_col  = "REGION",
  level_col = "Detector_Label",
  value_col = "Dominance_P",
  focal     = "Constant TA",
  tests     = dom_tests)
dom_vrange <- diff(range(panel_c_long$Dominance_P, na.rm = TRUE))
if (!is.finite(dom_vrange) || dom_vrange <= 0) dom_vrange <- 1
# build_sig_brackets() returns NULL when nothing is significant; nrow(NULL) is
# NULL, not 0, so the count is normalised here rather than inline.
n_dom_brackets <- if (is.null(dom_brackets)) 0L else nrow(dom_brackets)

p_dot <- ggplot2::ggplot(panel_c_long,
                         ggplot2::aes(x = Detector_Label, y = Dominance_P)) +
  # Decisive / chance reference lines, their labels and the shaded decisive
  # band removed. DOMINANCE_THRESHOLD still drives the tier
  # classification and the "Dominance (score >= x)" column; it is simply no
  # longer drawn on this panel.
  ggplot2::geom_point(
    ggplot2::aes(fill = Detector_Label,
                 stroke = border_stroke,
                 colour = border_colour),
    shape = 21, position = pj_c,
    size = 3.0, alpha = 0.92, na.rm = TRUE, show.legend = FALSE
  ) +
  ggplot2::stat_summary(
    fun = stats::median, geom = "crossbar",
    width = 0.55, linewidth = 0.45, colour = "grey15",
    fatten = 0, fill = NA, na.rm = TRUE
  ) +
  # Dot-plot-only override: Continuous TA circles are drawn red, as requested.
  # PANEL_C_CONTINUOUS_TA_FILL is Okabe-Ito vermillion -- a true red that stays
  # separable from the blue and green under deuteranopia and protanopia, so the
  # panel keeps the colourblind-safe property of the rest of the figure.
  #
  # NOTE: this makes Continuous TA red in panel b (the dot plot) but orange in
  # panel a of the same figure. To use red everywhere instead, set
  # detector_palette["Continuous TA"] <- PANEL_C_CONTINUOUS_TA_FILL where the
  # palette is defined (section 24) and delete this override.
  # All six target detectors, built from detector_palette so the scale cannot
  # fall out of step with what is plotted, then the panel-specific red override
  # for Continuous TA is applied on top.
  ggplot2::scale_fill_manual(
    values = {
      .v <- detector_palette[unname(DETECTOR_SHORT_LABEL[unname(TARGET_P_COLS)])]
      .v["Continuous TA"] <- PANEL_C_CONTINUOUS_TA_FILL
      .v
    }
  ) +
  ggplot2::scale_colour_identity() +
  ggplot2::scale_y_continuous(
    limits = c(0, 1.05 + 0.16 * n_dom_brackets),
    breaks = c(0, 0.25, 0.50, 0.75, 1.00),
    labels = c("0", "0.25", "0.50", "0.75", "1.0"),
    expand = ggplot2::expansion(mult = c(0, 0))
  ) +
  ggplot2::scale_x_discrete(
    # All six detectors get an explicit two-line label.
    labels = c("Constant TA" = "Constant\nTA",
               "Continuous TA" = "Continuous\nTA",
               "Outbreak Threshold" = "Outbreak\nThreshold",
               "Alarm Threshold" = "Alarm\nThreshold",
               "WHO 75th" = "WHO 75th",
               "WHO 90th" = "WHO 90th"),
    expand = ggplot2::expansion(add = c(0.6, 0.6))
  ) +
  ggplot2::labs(
    title = "Per-detector regional dominance probabilities",
    x = NULL, y = "Dominance probability"
  ) +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::theme_minimal(base_size = PUB_BASE, base_family = PUB_FAMILY) +
  theme_bold_axes() +
  ggplot2::theme(
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor   = ggplot2::element_blank(),
        axis.text.x  = ggplot2::element_text(size = PUB_AXIS_TXT, face = "bold",
                                         lineheight = 0.85,
                                         hjust = 0.5, vjust = 1,
                                         margin = ggplot2::margin(t = 4)),
    axis.text.y  = ggplot2::element_text(size = PUB_AXIS_TXT),
    axis.title.y = ggplot2::element_text(size = PUB_AXIS_TIT, face = "bold",
                                         colour = "grey25",
                                         margin = ggplot2::margin(r = 6)),
    plot.title   = ggplot2::element_text(face = "bold", size = 13,
                                         margin = ggplot2::margin(b = 8)),
    plot.margin  = ggplot2::margin(20, 30, 14, 18)
  )

# ------------------------------------------------------------------------------
# 28. PER-REGION PER-METRIC HEATMAP DATA (saved as CSV; not a figure panel)
# ------------------------------------------------------------------------------
# REVISION: trimmed from six subgrids (three unordered pairs, both
# directions) to the four that exist now that Constant TA vs Continuous TA
# is no longer computed. See Item 8 of the Delphi review.
PANEL_D_SUBGRIDS <- list(
  list(label = "Constant TA vs Outbreak Threshold",
       pair_code = "const_vs_ot"),
  list(label = "Continuous TA vs Outbreak Threshold",
       pair_code = "cont_vs_ot"),
  list(label = "Outbreak vs Constant TA",
       pair_code = "ot_vs_const"),
  list(label = "Outbreak vs Continuous TA",
       pair_code = "ot_vs_cont")
)

build_panel_d_data <- function() {
  rows <- list()
  for (sg in PANEL_D_SUBGRIDS) {
    sub <- level1_results_ordered |>
      dplyr::filter(Pair_Code == sg$pair_code,
                    REGION %in% CANONICAL_17)
    if (nrow(sub) == 0L) next
    sub <- sub |>
      dplyr::mutate(
        direction = dplyr::case_when(
          is.na(p_bonf)                      ~ "na",
          sig_star %in% c("*", "**", "***")  ~ "win",
          TRUE                               ~ "ns"
        ),
        fill_key = dplyr::case_when(
          direction == "na"  ~ "na",
          direction == "ns"  ~ "ns",
          direction == "win" ~ paste0("win_", sig_star),
          TRUE               ~ "na"
        ),
        sig_label = ifelse(direction == "win", sig_star, ""),
        Subgrid_Label = sg$label,
        Subgrid_Pair_Code = sg$pair_code
      ) |>
      dplyr::select(REGION, Metric,
                    sig_star, p_bonf, p_one_sided, median_diff_XY,
                    direction, fill_key, sig_label,
                    Subgrid_Label, Subgrid_Pair_Code)
    rows[[length(rows) + 1L]] <- sub
  }
  if (length(rows) == 0L) return(data.frame())
  dplyr::bind_rows(rows)
}
panel_d_data <- build_panel_d_data()

# ------------------------------------------------------------------------------
# 29. COMPOSE AND SAVE FIGURE 5
# ------------------------------------------------------------------------------
# Panel B (the per-region metrics table) has been removed. Figure 5 is now two
# panels: the map (a) above the dominance-probability plot (b). Former panel c
# is therefore relabelled b throughout -- titles, filenames and captions.
fig5 <- (p_map / p_dot) +
  patchwork::plot_layout(heights = c(1.35, 1.00)) &
  # Padding between the map and dot panels so their titles, legends and axis
  # labels cannot collide in the composite.
  ggplot2::theme(plot.margin = ggplot2::margin(10, 12, 8, 8))

fig5 <- fig5 +
  patchwork::plot_annotation(
    title = "Bootstrap-supported regional dominance of six requested detectors",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = PUB_TITLE + 1,
                                         family = base_family_global,
                                         margin = ggplot2::margin(b = 6)),
      plot.tag = ggplot2::element_text(size = PUB_TAG, face = "bold",
                                       family = base_family_global,
                                       hjust = 0, vjust = 1),
      # Same left edge as the title: tags anchor to the plot, and
      # plot.title.position = "plot" moves the title off the panel edge to match.
      plot.tag.position     = "topleft",
      plot.title.position   = "plot"
    )
  )

# Figure 5 (map + dominance probability), authored at final print size.
FIG5_W <- NC_W_DOUBLE
FIG5_H <- 8.40
# Attach the significance brackets (none are drawn if nothing is significant).
p_dot <- add_sig_brackets(p_dot, dom_brackets, dom_vrange)

save_plot_file("FigureB5_Combined_RegionalDominance", fig5, FIG5_W, FIG5_H)

# ------------------------------------------------------------------------------
# 29b. SEPARATE PANEL EXPORTS (Figures 5A, 5B, 5C)
# ------------------------------------------------------------------------------
# Each panel is also exported on its own at full double-column width. Because a
# standalone panel has the whole width to itself rather than a share of the
# composite, titles are restated and the panels are given generous margins so
# no label, legend or annotation can collide.
panel_standalone_theme <- ggplot2::theme(
  plot.margin   = ggplot2::margin(10, 12, 8, 10),
  plot.title    = ggplot2::element_text(size = PUB_TITLE, face = "bold",
                                        family = base_family_global,
                                        margin = ggplot2::margin(b = 5)),
  legend.margin = ggplot2::margin(3, 3, 3, 3)
)

# Standalone titles carry no figure number, so each panel reads independently
# when reproduced on its own.
#
# NOTE ON WORDING: the requested examples were "National / Provincial /
# Regional Detector". Those do not describe these panels -- Stage 4 is entirely
# a regional analysis (17 Philippine administrative regions), and each
# panel differs by *what is shown* (map / dominance probability), not by
# administrative tier. There is no national or provincial stratification
# anywhere in this stage.
# Nature Communications style: lowercase letter, no trailing period.
# The per-region metrics table is not shown as a panel in Figure 5.
# The dominance-probability plot is panel b.
fig5a <- p_map + panel_standalone_theme +
  ggplot2::labs(title = "a    Regional detector map")
fig5b <- p_dot + panel_standalone_theme +
  ggplot2::labs(title = "b    Per-detector regional dominance probabilities")

# Heights are per-panel: the map is near-square, the table is tall and
# row-driven, the dot plot is short and wide.
save_plot_file("FigureB5_panel_a_DetectorMap",          fig5a, NC_W_DOUBLE, 5.60)
save_plot_file("FigureB5_panel_b_DominanceProbability", fig5b, NC_W_DOUBLE, 4.40)

# ------------------------------------------------------------------------------
# 30. WRITE FIGURE 5 CSV TABLES
# ------------------------------------------------------------------------------
readr::write_csv(winner_rows,
                 file.path(OUTPUT_DIR, "StageB_Detector_Map_RegionTable.csv"))
readr::write_csv(level1_results_ordered,
                 file.path(OUTPUT_DIR, "StageB_Detector_Map_PerMetricSignificance.csv"))
readr::write_csv(level2_results,
                 file.path(OUTPUT_DIR, "StageB_Detector_Map_CrossRegion_Wilcoxon.csv"))
readr::write_csv(
  panel_c_long |>
    dplyr::transmute(REGION, Detector = as.character(Detector_Label),
                     Dominance_P, is_leader_dot, consensus_tier,
                     Observed_Winner),
  file.path(OUTPUT_DIR, "StageB_Detector_Map_PanelC_DotPlot.csv")
)
readr::write_csv(panel_d_data,
                 file.path(OUTPUT_DIR, "StageB_Detector_Map_PanelD_HeatmapData.csv"))
readr::write_csv(audit_df,
                 file.path(OUTPUT_DIR, "StageB_Detector_Map_JoinAudit.csv"))

# ==============================================================================
# 31. LEGENDS
# ==============================================================================
# Two text files written to OUTPUT_DIR:
#   StageB_FigureB4_Legend.txt   (one block per metric + method summary)
#   StageB_FigureB5_Legend.txt    (single block for the detector map figure)
# ------------------------------------------------------------------------------

build_figure4_metric_legend <- function(metric_id) {
  metric_label <- hh_metric_full_name[[metric_id]]
  res <- all_wilcoxon_results |>
    dplyr::filter(Metric == metric_id) |>
    dplyr::arrange(Comparison)
  n_pairs <- max(res$N_Regions_Paired, na.rm = TRUE)
  sig_pairs <- res |>
    dplyr::filter(!is.na(Significant_005), Significant_005) |>
    dplyr::pull(Comparison) |> as.character()
  outcome <- if (length(sig_pairs) == 0L) {
    paste0("Neither of the two detector-pair comparisons (both against ",
           "Outbreak Threshold) reached the per-pairwise significance ",
           "threshold (p < ", sprintf("%.2f", HH_ALPHA),
           ") on the available evidence base of ", n_pairs,
           " paired regional observations.")
  } else if (length(sig_pairs) == 2L) {
    "Both detector-pair comparisons reached the per-pairwise significance threshold."
  } else {
    paste0(length(sig_pairs), " of 2 detector-pair comparisons reached the ",
           "per-pairwise significance threshold: ",
           paste(sig_pairs, collapse = "; "), ".")
  }
  paste0(
    "Figure 4 (", metric_label, ") | Regional generalisability of ",
    metric_label, " across the ", length(region_ppv_order),
    " evaluable Philippine regions, 2018-2025 (excluding 2020 and 2021; n = ", length(EVALUABLE_YEARS), " evaluable years per region; ",
    "regional coverage floor = ", MIN_OBSERVED_CASE_WEEKS_PER_YEAR,
    " observed case weeks/year and >= ", MIN_EVALUABLE_YEARS_PER_REGION,
    " evaluable years). This figure is one of the five per-metric panels ",
    "in the Figure 4 set (TAM, Number of True Alarms, Sensitivity, Mean Lead ",
    "Time, Warning Persistence). The early-warning timeliness metrics (Mean ",
    "Lead Time, Warning Persistence) use the same-denominator zero-coerced ",
    "scheme: years where the metric is not computable due to no qualifying ",
    "triggers contribute zero rather than being excluded. Sensitivity is A1-",
    "restricted (proportion of evaluable seasons with at least one True Alarm ",
    "in the Actionable Window).\n",
    "(a) Regional dominance matrix (method-centric layout). Rows = 6 ",
    "outbreak-detection methods, ordered top-to-bottom by sweep count ",
    "(descending; ties broken by descending mean score). Columns = ",
    length(region_ppv_order), " evaluable regions, ordered left-to-right by ",
    "region mean dominance score across the 6 methods (descending). Cell ",
    "shading encodes the normalized dominance score on this metric, computed ",
    "by min-max normalization within each region across detectors: 1.00 = ",
    "within-region leader on ", metric_label, "; 0.00 = within-region laggard. ",
    "Single-color blue ramp from white (0) through mid-blue (0.50) to dark ",
    "navy (1.00); a horizontal blue-intensity legend strip below the matrix ",
    "maps shade to score with reference labels at 0.00 (weak), 0.50 ",
    "(moderate), ", sprintf("%.2f", DOMINANCE_THRESHOLD),
    " (dominant), and 1.00 (sweep). The right-most annotation reports the ",
    "per-METHOD count of regions swept by that method on this metric (regions ",
    "where the method achieves dominance score >= ",
    sprintf("%.2f", DOMINANCE_THRESHOLD), "), out of ",
    length(region_ppv_order), " total. ",
    "(b) Per-detector regional values on ", metric_label,
    ". One dot per region per detector, jittered horizontally; mean across ",
    "regions shown as a horizontal crossbar; vertical lines indicate year-",
    "cluster bootstrap 95% CIs (B = ", BOOT_N_CI, " replicates) where ",
    "available. The pale-blue band marks the operationally favourable zone ",
    "where defined. ",
    "(c) Wilcoxon detector-paired significance bars. Two pairwise detector ",
    "comparisons, both against Outbreak Threshold: Constant TA vs Outbreak ",
    "Threshold; Continuous TA vs Outbreak Threshold. Wilcoxon signed-rank, ",
    "paired by region (n = ", n_pairs,
    "), two-sided, evaluated PER PAIRWISE COMPARISON at alpha = ",
    sprintf("%.2f", HH_ALPHA),
    " (no across-pair multiplicity correction). Bars show |median ",
    "difference|; navy bars (with adjacent '*') indicate p < ",
    sprintf("%.2f", HH_ALPHA), "; pale grey otherwise. ",
    "Outcome on this metric: ", outcome, "\n"
  )
}

build_figure4_summary_legend <- function() {
  # REVISION: trimmed from three pairs (a full round robin among Constant
  # TA, Continuous TA and Outbreak Threshold) to the two comparisons
  # instructed, both against Outbreak Threshold. Alarm Threshold is not
  # part of either pair.
  TARGETS_LOCAL <- c("Constant Transmission Acceleration",
                     "Continuous Transmission Acceleration",
                     "Outbreak Threshold")
  PAIRS_LOCAL <- list(
    list(a = "Constant Transmission Acceleration",
         b = "Outbreak Threshold",
         label = "Constant TA vs Outbreak Threshold"),
    list(a = "Continuous Transmission Acceleration",
         b = "Outbreak Threshold",
         label = "Continuous TA vs Outbreak Threshold")
  )
  pivot_score_lg <- regional_dominance_long |>
    dplyr::filter(as.character(Method) %in% TARGETS_LOCAL,
                  as.character(Metric) %in% HH_METRICS_REGIONAL) |>
    dplyr::select(REGION, Metric, Method, Dominance_Score) |>
    dplyr::mutate(Method = as.character(Method)) |>
    tidyr::pivot_wider(names_from = Method, values_from = Dominance_Score)
  sig_pairs_summary <- character(0)
  n_pairs_typical <- NA_integer_
  for (pr in PAIRS_LOCAL) {
    a_sc <- pivot_score_lg[[pr$a]]; b_sc <- pivot_score_lg[[pr$b]]
    keep <- !is.na(a_sc) & !is.na(b_sc)
    a_sc <- a_sc[keep]; b_sc <- b_sc[keep]
    if (length(a_sc) >= 2L && any(a_sc != b_sc)) {
      pv <- suppressWarnings(stats::wilcox.test(a_sc, b_sc, paired = TRUE,
                                                exact = FALSE))$p.value
      if (!is.na(pv) && pv < 0.05) {
        sig_pairs_summary <- c(sig_pairs_summary, pr$label)
      }
    }
    n_pairs_typical <- length(a_sc)
  }
  outcome <- if (length(sig_pairs_summary) == 0L) {
    paste0("Neither of the two detector-pair comparisons reached p < 0.05 on ",
           "the available evidence base of approximately ",
           n_pairs_typical, " region\u00D7metric paired observations.")
  } else if (length(sig_pairs_summary) == 2L) {
    "Both detector-pair comparisons reached p < 0.05."
  } else {
    paste0(length(sig_pairs_summary), " of 2 detector-pair comparisons ",
           "reached p < 0.05: ",
           paste(sig_pairs_summary, collapse = "; "), ".")
  }

  paste0(
    "Figure 4 (Method-Across-Metrics Summary) | Cross-metric summary of ",
    "detector performance across the 17 Philippine regions and the five ",
    "reported operational metrics (TAM, Number of True Alarms, Sensitivity, ",
    "Mean Lead Time, Warning Persistence), 2018-2025 (excluding 2020 and 2021; n = ", length(EVALUABLE_YEARS),
    " evaluable years per region; regional coverage floor = ",
    MIN_OBSERVED_CASE_WEEKS_PER_YEAR, " observed case weeks/year and >= ",
    MIN_EVALUABLE_YEARS_PER_REGION, " evaluable years). This composite ",
    "figure complements the five per-metric multipanels by aggregating ",
    "detector performance across all five metrics simultaneously, ",
    "answering: which detector achieves the best operational performance ",
    "across the broadest combination of metrics and regions?\n",
    "(a) Regional dominance matrix on five operational metrics. Rows = 6 ",
    "outbreak-detection methods, ordered top-to-bottom by mean sweep count ",
    "across the 5 metrics (descending; best-overall method at top). Columns ",
    "= the five reported metrics, ordered left-to-right as TAM, Number of ",
    "True Alarms, Sensitivity, Mean Lead Time, Warning Persistence. Cell ",
    "shading encodes the per-(method, metric) sweep count: the number of ",
    "regions (out of ", MAX_REGIONS_POSSIBLE,
    " evaluable) where that method achieves within-region dominance score ",
    ">= ", sprintf("%.2f", DOMINANCE_THRESHOLD),
    " on that metric. Single-color blue ramp from white (0 regions) through ",
    "mid-blue to dark navy (all ", MAX_REGIONS_POSSIBLE, " regions). ",
    "Numeric annotation inside each cell shows the raw sweep count (auto-",
    "contrast: dark text on light fills, white text on dark fills). Column ",
    "headers (TAM, Number of True Alarms, Sensitivity, Mean Lead Time, ",
    "Warning Persistence; and the rightmost Mean Sweep / z-score header) ",
    "are rendered in-cell as text inside neutral-grey header tiles ",
    "immediately above the top method row. The right-most data column ",
    "reports per-method aggregate values: the mean sweep count across the ",
    "5 metrics, followed by a parenthetical z-score computed across the 4 ",
    "methods. Classification thresholds: z >= +1.0 = 'Sweeper'; z <= -1.0 = ",
    "'Non-sweeper'; otherwise 'Average'. ",
    "(b) Per-method dot plot. One dot per metric per method (5 dots per row, ",
    "jittered horizontally for visibility) showing the raw sweep count on ",
    "that metric; method-mean across the 5 metrics shown as a horizontal ",
    "crossbar. ",
    "(c) Wilcoxon detector-paired significance bars (DUAL-STATISTIC). Two ",
    "pairwise detector contrasts, both against Outbreak Threshold: Constant ",
    "TA vs Outbreak Threshold; Continuous TA vs Outbreak Threshold. The bar ",
    "length per contrast is |median difference in regions-swept count across ",
    "the 5 metrics| (units: regions; computed from the 5 paired metric-level ",
    "counts per detector). The p-value per contrast is computed from a ",
    "separate properly-powered test: paired Wilcoxon signed-rank on within-",
    "region Dominance_Score, paired by (REGION x METRIC) cells, n_pairs ~ ",
    n_pairs_typical, " region\u00D7metric paired observations. Bars are ",
    "navy (with adjacent '*') if p < 0.05; pale grey otherwise. The dual ",
    "specification answers two distinct questions in one panel: is the ",
    "detector difference statistically real (n_pairs ~ ", n_pairs_typical,
    " powered test), and how big is the practical difference in operationally ",
    "interpretable units (n=5 sweep-count median across the 5 metrics)? ",
    "Outcome: ", outcome, "\n"
  )
}

# Combined Figure 4 legend (one block per metric + Method Summary).
fig4_legend_lines <- c(
  paste0("================================================================"),
  paste0("FIGURE 4  |  LEGENDS  (per-metric and method summary)"),
  paste0("================================================================")
)
for (m in HH_METRICS_REGIONAL) {
  fig4_legend_lines <- c(
    fig4_legend_lines, "",
    paste0("--- FigureB4_", hh_metric_slug[[m]], " ---"),
    build_figure4_metric_legend(m)
  )
}
fig4_legend_lines <- c(
  fig4_legend_lines, "",
  "--- FigureB4_Method_Summary ---",
  build_figure4_summary_legend()
)
fig4_legend_text <- paste(fig4_legend_lines, collapse = "\n")
writeLines(fig4_legend_text,
           file.path(OUTPUT_DIR, "StageB_FigureB4_Legend.txt"))
cat("\nSaved: ", file.path(OUTPUT_DIR, "StageB_FigureB4_Legend.txt"),
    "\n", sep = "")

# Figure 5 legend
build_figure5_legend <- function() {
  paste0(
    "Figure 5 | Regional detector dominance for leptospirosis. ",
    "Panel a maps the observed composite winner among the six requested ",
    "detectors: Constant TA, Continuous TA, Outbreak Threshold, Alarm ",
    "Threshold, WHO 75th and WHO 90th. Polygon fill identifies the winner; ",
    "a dashed boundary denotes dominance probability below ",
    sprintf("%.2f", DOMINANCE_THRESHOLD),
    ". Only inclusion-qualifying regions are mapped. Plain region names are ",
    "shown on the map. Panel b shows the bootstrap dominance probability for ",
    "all six detectors by region, with the chance reference equal to 1/6. ",
    "Dominance probabilities are based on ", N_BOOTS,
    " year-cluster bootstrap replicates per region. The per-region metrics, ",
    "consensus fields, pairwise significance results, geometry join audit, ",
    "and six-detector probability table are exported as CSV files in the same ",
    "output directory."
  )
}

fig5_legend_text <- build_figure5_legend()
writeLines(fig5_legend_text,
           file.path(OUTPUT_DIR, "StageB_FigureB5_Legend.txt"))
cat("Saved: ", file.path(OUTPUT_DIR, "StageB_FigureB5_Legend.txt"),
    "\n", sep = "")

# Echo legends to console.
cat("\n", fig4_legend_text, "\n", sep = "")
cat("\n", fig5_legend_text, "\n", sep = "")

# ==============================================================================
# 32. PARAMETER REPORT AND END-OF-RUN BANNER
# ==============================================================================
cat("\n============================================================\n")
cat("STAGE B - REGIONAL LEPTOSPIROSIS ANALYSIS  ", SCRIPT_VERSION, "\n", sep = "")
cat("============================================================\n")
cat("Anchor framework parameters (two-anchor; T = A1 union A2):\n")
cat("  A1 lead window         = [", A1_LEAD_MIN, ", ", A1_LEAD_MAX,
    "] weeks\n", sep = "")
cat("  A2 burden fraction     = ", A2_BURDEN_FRAC, "\n", sep = "")
cat("  True-Alarm rule        = t in A1 OR t in A2\n")
cat("\nDetector parameters:\n")
cat("  Continuous TA threshold = ", ETA_ON_CLASSIC, "\n", sep = "")
cat("  Constant TA   eta_ON    = ", ETA_ON,         "\n", sep = "")
cat("  Constant TA   eta_OFF   = ", ETA_OFF,        "\n", sep = "")
cat("  Constant TA STA/LTA/GUARD = ", STA_WIN_CONSTANT, " / ", LTA_WIN_CONSTANT,
    " / ", GUARD, "\n", sep = "")
cat("  Continuous TA STA/LTA     = ", STA_WIN_CONTINUOUS, " / ", LTA_WIN_CONTINUOUS,
    "\n", sep = "")
cat("\nRegional inclusion:\n")
cat("  MIN_OBSERVED_CASE_WEEKS_PER_YEAR        = ", MIN_OBSERVED_CASE_WEEKS_PER_YEAR,
    "\n", sep = "")
cat("  MIN_EVALUABLE_YEARS_PER_REGION = ", MIN_EVALUABLE_YEARS_PER_REGION,
    "\n", sep = "")
cat("  Excluded years                 = ",
    paste(EXCLUDED_YEARS, collapse = ", "), "\n", sep = "")
cat("\nBootstrap configuration:\n")
cat("  B = ", BOOT_N_CI, " year-cluster replicates per region\n", sep = "")
cat("  Wilcoxon detector-paired alpha (per pairwise comparison) = ",
    sprintf("%.2f", HH_ALPHA), "\n", sep = "")
cat("  Per-region per-metric Bonferroni factor                  = ",
    BONF_LEVEL1_FACTOR, " (within REGION x ordered pair)\n", sep = "")
cat("  Cross-region per-metric Bonferroni factor                = ",
    BONF_LEVEL2_FACTOR, " (within metric)\n", sep = "")

cat("\n=== Saved files (", OUTPUT_DIR, ") ===\n", sep = "")
cat("Figure 4 (per-metric, 5 figures + 1 method summary):\n")
cat("  FigureB4_TAM.{pdf,png}\n")
cat("  FigureB4_N_True_Alarms.{pdf,png}\n")
cat("  FigureB4_Sensitivity.{pdf,png}\n")
cat("  FigureB4_Mean_Lead_Time.{pdf,png}\n")
cat("  FigureB4_Warning_Persistence.{pdf,png}\n")
cat("  FigureB4_Method_Summary.{pdf,png}\n")
cat("Figure 5 (composite detector map):\n")
cat("  FigureB5_Combined_RegionalDominance.{pdf,png}\n")
cat("  FigureB5_panel_a_DetectorMap.{pdf,png}\n")
cat("  FigureB5_panel_b_DominanceProbability.{pdf,png}\n")
cat("Tables (Figure 4 supporting):\n")
cat("  StageB_Regional_Framework_Metrics.csv\n")
cat("  StageB_Regional_Framework_Metrics_with_CIs.csv\n")
cat("  StageB_Regional_8Metric_Summary.csv\n")
cat("  StageB_Regional_8Metric_Summary_with_CIs.csv\n")
cat("  StageB_Regional_Dominance_Matrix.csv\n")
cat("  StageB_Regional_Dominance_Probabilities.csv\n")
cat("  StageB_Regional_Wilcoxon_PerMetric.csv\n")
cat("  StageB_Regional_Wilcoxon_ConstantTA_vs_Comparators.csv\n")
cat("  StageB_Method_Summary_Long.csv\n")
cat("  StageB_Method_Summary_Aggregate.csv\n")
cat("  StageB_Regional_Bootstrap_Replicates.csv\n")
cat("Tables (Figure 5 supporting):\n")
cat("  StageB_Detector_Map_RegionTable.csv\n")
cat("  StageB_Detector_Map_PerMetricSignificance.csv\n")
cat("  StageB_Detector_Map_CrossRegion_Wilcoxon.csv\n")
cat("  StageB_Detector_Map_PanelC_DotPlot.csv\n")
cat("  StageB_Detector_Map_PanelD_HeatmapData.csv\n")
cat("  StageB_Detector_Map_JoinAudit.csv\n")
cat("Legends:\n")
cat("  StageB_FigureB4_Legend.txt\n")
cat("  StageB_FigureB5_Legend.txt\n")

cat("\nDone.\n")

# ==============================================================================
# END OF SCRIPT - STAGE 4 REGIONAL ANALYSIS
# ==============================================================================
