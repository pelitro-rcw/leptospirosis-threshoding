# =============================================================================
# STAGE 2 - NCR LEPTOSPIROSIS OUTBREAK-THRESHOLD DRIFT QC
# =============================================================================
# Reproduces the dengue threshold-drift QC logic for NCR leptospirosis using
# LC_DOH from the Regional Data sheet. Because leptospirosis has no observed
# rows in 2020-2021, a literal "pandemic included vs excluded" comparison is
# impossible without inventing data. This QC therefore compares:
#   A. calendar-window baseline: observed donors available in y-5..y-1; and
#   B. gap-aware baseline: the last up to five OBSERVED prior leptospirosis
#      years, requiring at least 3 donor years.
# Both use week-specific Alarm = mean + 1 SD and Outbreak = mean + 2 SD.
# =============================================================================

REQUIRED_PACKAGES <- c("readxl","dplyr","tidyr","ggplot2","ISOweek","patchwork","scales")
.bootstrap <- function() {
  root <- Sys.getenv("TA_PROJECT_ROOT", unset="")
  if (!nzchar(root) || !dir.exists(root)) {
    d <- normalizePath(getwd(), winslash="/")
    for (i in seq_len(6)) {
      if (dir.exists(file.path(d,"R")) && dir.exists(file.path(d,"scripts"))) {root <- d; break}
      p <- dirname(d); if (identical(p,d)) break; d <- p
    }
  }
  if (!nzchar(root)) stop("Set the working directory to the project root before sourcing.", call.=FALSE)
  source(file.path(root,"R","00_config.R"))
}
.bootstrap()
require_packages(REQUIRED_PACKAGES, "Stage 2 NCR threshold-drift QC")
invisible(lapply(REQUIRED_PACKAGES, function(p) suppressPackageStartupMessages(library(p, character.only=TRUE))))
source(file.path(DIR_R,"01_publication_theme.R"), local=TRUE)
set.seed(GLOBAL_SEED)

OUT_DIR <- file.path(DIR_OUTPUT,"Stage2_outbreak_threshold_drift_NCR_Lepto")
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

d <- readxl::read_excel(DATA_FILE, sheet=SHEET_REGIONAL)
need <- c("REGION","YR","WN","LC_DOH")
miss <- setdiff(need,names(d)); if(length(miss)) stop("Missing columns: ",paste(miss,collapse=", "))
d <- d |> dplyr::transmute(
  REGION=as.character(REGION), YR=as.integer(YR), WN=as.integer(WN),
  Cases=suppressWarnings(as.numeric(LC_DOH))) |>
  dplyr::filter(REGION=="NCR", !is.na(YR), !is.na(WN), WN>=1, WN<=53,
                YR>=2018, YR<=2025) |>
  dplyr::mutate(Date=ISOweek::ISOweek2date(sprintf("%d-W%02d-1",YR,WN))) |>
  dplyr::arrange(Date)
if(!nrow(d) || all(is.na(d$Cases))) stop("No NCR LC_DOH observations available.")
observed_years <- sort(unique(d$YR[is.finite(d$Cases)]))

make_map <- function(mode=c("calendar","gap_aware")) {
  mode <- match.arg(mode); out <- list()
  for(y in observed_years) {
    prior <- observed_years[observed_years < y]
    donors <- if(mode=="calendar") intersect(seq.int(y-5L,y-1L),prior) else tail(prior,5L)
    if(length(donors)<3L) donors <- integer(0)
    out[[as.character(y)]] <- donors
  }
  out
}

calc_thresholds <- function(map, specification) {
  rows <- list()
  for(nm in names(map)) {
    y <- as.integer(nm); donors <- map[[nm]]
    weeks <- sort(unique(d$WN[d$YR==y]))
    for(w in weeks) {
      vals <- d$Cases[d$YR %in% donors & d$WN==w]; vals <- vals[is.finite(vals)]
      n <- length(vals)
      mu <- if(n) mean(vals) else NA_real_
      sdv <- if(n>=2L) stats::sd(vals) else if(n==1L) 0 else NA_real_
      rows[[length(rows)+1L]] <- data.frame(
        Specification=specification,YR=y,WN=w,N_Donors=n,
        Donor_Years=if(length(donors)) paste(donors,collapse=",") else "",
        Mean=mu,SD=sdv,
        Alarm_Threshold=if(is.finite(mu)&&is.finite(sdv)) mu+sdv else NA_real_,
        Outbreak_Threshold=if(is.finite(mu)&&is.finite(sdv)) mu+2*sdv else NA_real_)
    }
  }
  dplyr::bind_rows(rows)
}

th_cal <- calc_thresholds(make_map("calendar"),"Calendar y-5..y-1 observed donors")
th_gap <- calc_thresholds(make_map("gap_aware"),"Last up to 5 observed prior years")
utils::write.csv(th_cal,file.path(OUT_DIR,"FigureA1_thresholds_calendar_window.csv"),row.names=FALSE)
utils::write.csv(th_gap,file.path(OUT_DIR,"FigureA1_thresholds_gap_aware.csv"),row.names=FALSE)

cmp <- th_cal |> dplyr::select(YR,WN,Outbreak_Calendar=Outbreak_Threshold,Alarm_Calendar=Alarm_Threshold) |>
  dplyr::full_join(th_gap |> dplyr::select(YR,WN,Outbreak_GapAware=Outbreak_Threshold,Alarm_GapAware=Alarm_Threshold), by=c("YR","WN")) |>
  dplyr::mutate(Outbreak_Drift=Outbreak_GapAware-Outbreak_Calendar,
                Alarm_Drift=Alarm_GapAware-Alarm_Calendar)
utils::write.csv(cmp,file.path(OUT_DIR,"FigureA1_threshold_drift_comparison.csv"),row.names=FALSE)

# Figure display is intentionally restricted to 2023-2025. Earlier observed
# years remain in `d` because they are legitimate donor years for the rolling
# threshold calculation, but they are not displayed in Figure A1.
FIGURE_YEARS <- 2023:2025

plot_dat <- d |>
  dplyr::filter(YR %in% FIGURE_YEARS, is.finite(Cases)) |>
  dplyr::left_join(
    th_cal |> dplyr::select(
      YR, WN,
      Outbreak_Calendar = Outbreak_Threshold,
      Alarm_Calendar = Alarm_Threshold
    ),
    by = c("YR", "WN")
  ) |>
  dplyr::arrange(Date)

if (!nrow(plot_dat)) {
  stop("FigureA1: no finite NCR leptospirosis observations in 2023-2025.", call. = FALSE)
}
missing_display_years <- setdiff(FIGURE_YEARS, sort(unique(plot_dat$YR)))
if (length(missing_display_years)) {
  stop("FigureA1: requested display year(s) have no NCR leptospirosis data: ",
       paste(missing_display_years, collapse = ", "), call. = FALSE)
}

# Single-panel NCR threshold QC restricted to 2023-2025. The plotted
# thresholds use the same rolling calendar-window baseline specification used
# by the primary NCR analysis. The gap-aware alternative remains available in
# the supporting drift tables but is intentionally not drawn here.
fig <- ggplot2::ggplot(plot_dat, ggplot2::aes(x = Date)) +
  ggplot2::geom_line(
    ggplot2::aes(y = Cases, colour = "Leptospirosis Cases",
                 linetype = "Leptospirosis Cases"),
    linewidth = 0.55, na.rm = TRUE) +
  ggplot2::geom_line(
    ggplot2::aes(y = Outbreak_Calendar,
                 colour = "Outbreak Threshold (mean + 2SD)",
                 linetype = "Outbreak Threshold (mean + 2SD)"),
    linewidth = 0.75, na.rm = TRUE) +
  ggplot2::geom_line(
    ggplot2::aes(y = Alarm_Calendar,
                 colour = "Alarm Threshold (mean + SD)",
                 linetype = "Alarm Threshold (mean + SD)"),
    linewidth = 0.75, na.rm = TRUE) +
  ggplot2::geom_vline(
    xintercept = as.numeric(as.Date(c("2024-01-01", "2025-01-01"))),
    linewidth = 0.35, linetype = "dotted", colour = "grey55",
    show.legend = FALSE) +
  ggplot2::scale_colour_manual(
    name = NULL,
    values = c(
      "Leptospirosis Cases" = "black",
      "Outbreak Threshold (mean + 2SD)" = "#D62728",
      "Alarm Threshold (mean + SD)" = "#377EB8"
    ),
    breaks = c(
      "Leptospirosis Cases",
      "Outbreak Threshold (mean + 2SD)",
      "Alarm Threshold (mean + SD)"
    ),
    guide = ggplot2::guide_legend(nrow = 2, byrow = TRUE, direction = "horizontal")) +
  ggplot2::scale_linetype_manual(
    name = NULL,
    values = c(
      "Leptospirosis Cases" = "solid",
      "Outbreak Threshold (mean + 2SD)" = "dashed",
      "Alarm Threshold (mean + SD)" = "dotdash"
    ),
    breaks = c(
      "Leptospirosis Cases",
      "Outbreak Threshold (mean + 2SD)",
      "Alarm Threshold (mean + SD)"
    ),
    guide = ggplot2::guide_legend(nrow = 2, byrow = TRUE, direction = "horizontal")) +
  ggplot2::scale_x_date(
    limits = c(as.Date("2023-01-01"), as.Date("2025-12-31")),
    date_breaks = "3 months", date_labels = "%b\n%Y",
    expand = ggplot2::expansion(mult = c(0.005, 0.005))) +
  ggplot2::scale_y_continuous(
    breaks = scales::pretty_breaks(n = 5),
    expand = ggplot2::expansion(mult = c(0, 0.06))) +
  ggplot2::labs(
    title = "NCR outbreak-threshold drift, 2023-2025",
    x = NULL, y = "Weekly leptospirosis cases") +
  theme_pub() +
  ggplot2::theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.justification = "center",
    legend.box.just = "center",
    legend.text = ggplot2::element_text(size = 6.5, lineheight = 0.92,
                                        margin = ggplot2::margin(l = 2, r = 6)),
    legend.key.width = grid::unit(11, "pt"),
    legend.key.height = grid::unit(8.5, "pt"),
    legend.spacing.x = grid::unit(5, "pt"),
    legend.spacing.y = grid::unit(2, "pt"),
    legend.margin = ggplot2::margin(3, 6, 3, 6),
    panel.grid.minor = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(t = 10, r = 16, b = 16, l = 16),
    axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 8)))

save_pub("FigureA1_NCR_OutbreakThreshold_Drift", fig,
         NC_W_DOUBLE, min(5.6, NC_H_SUPP), OUT_DIR,
         max_height = NC_H_SUPP)

# Figure-specific data export: only the years shown in the single-panel figure.
utils::write.csv(
  plot_dat |> dplyr::select(YR, WN, Date, Cases,
                             Outbreak_Calendar, Alarm_Calendar),
  file.path(OUT_DIR, "FigureA1_NCR_OutbreakThreshold_Drift_2023_2025.csv"),
  row.names = FALSE)

writeLines(c(
  "NCR leptospirosis outbreak-threshold drift QC.",
  "The source workbook contains no leptospirosis observations for 2020-2021, so a literal pandemic-included baseline cannot be estimated without imputation.",
  "The single displayed figure is restricted to 2023-2025; earlier observed years remain eligible as donor years for threshold estimation.",
  "The calendar specification uses observed donors inside y-5..y-1. The gap-aware specification uses the last up to five observed prior years. Both require at least three donors.",
  "Figure A1 displays exactly three legend series: Leptospirosis Cases, Outbreak Threshold (mean + 2SD), and Alarm Threshold (mean + SD). No other series are drawn in the figure.",
  "The figure uses a data-driven weekly leptospirosis y-axis with automatic headroom so observed peaks and threshold curves are not cropped.",
  "Alarm threshold = week-specific donor mean + 1 SD; outbreak threshold = donor mean + 2 SD. The supporting drift tables retain the gap-aware comparison for QC."
),file.path(OUT_DIR,"Stage2_Figure_Legend.txt"))
cat("Stage 2 complete: ",OUT_DIR,"\n",sep="")
