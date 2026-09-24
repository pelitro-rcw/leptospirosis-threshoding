
# Preflight parser, input-schema, and bundled-geometry audit.
cat("\n=== PRE-FLIGHT SELF-CHECK ===\n")
root <- Sys.getenv("TA_PROJECT_ROOT", unset="")
if(!nzchar(root)) root <- normalizePath(getwd(), winslash="/")
source(file.path(root,"R","00_config.R"))

# RF operator regression test. Current operational RF products are strictly
# upper-tail activations (>=).
.stopif_rf_operator <- function() {
  if (normalize_rf_operator(c(RF_OPERATOR_PRODUCT_TEST = ">=")) != ">=")
    stop("RF operator normalization self-check failed for >=.", call. = FALSE)
  x <- c(1, 2, 3, NA_real_)
  if (!identical(apply_rf_gate(x, 2, ">="), c(FALSE, TRUE, TRUE, FALSE)))
    stop("RF gate self-check failed for >=.", call. = FALSE)
  cls <- classify_rf_activation(c(1, 2, 3, 4), point = 3, lower = 2, upper = 4,
                                operator = ">=")
  exp <- c("No RF activation", "Precautionary RF trigger",
           "Formal RF trigger", "Strong RF trigger")
  if (!identical(cls, exp))
    stop("RF upper-tail activation classification self-check failed.", call. = FALSE)
}
.stopif_rf_operator()
rm(.stopif_rf_operator)

PREFLIGHT_REQUIRED_PACKAGES <- c(
  "readxl", "dplyr", "tidyr", "purrr", "tibble", "stringr", "readr",
  "zoo", "ISOweek", "scales", "rlang", "pROC", "MASS", "boot",
  "ggplot2", "ggrepel", "patchwork", "cowplot", "viridisLite",
  "sf", "terra", "ggspatial", "prettymapr", "here", "surveillance"
)
require_packages(PREFLIGHT_REQUIRED_PACKAGES, "full pipeline pre-flight self-check")

r_files <- c(list.files(file.path(root,"R"),pattern="\\.R$",full.names=TRUE),
             list.files(file.path(root,"scripts"),pattern="\\.R$",full.names=TRUE),
             file.path(root,"run_all.R"))
# Native-pipe reproducibility guard. The project requires R >= 4.1 and must
# not depend on an attached magrittr/dplyr pipe operator. This specifically
# prevents the historical missing legacy-pipe-operator preflight failure.
legacy_pipe_token <- paste0("%", ">", "%")
legacy_pipe_files <- r_files[vapply(r_files, function(f) {
  any(grepl(legacy_pipe_token, readLines(f, warn = FALSE), fixed = TRUE))
}, logical(1))]
if (length(legacy_pipe_files)) {
  stop("Legacy magrittr pipe detected in project R file(s):\n  ",
       paste(substring(legacy_pipe_files, nchar(root) + 2L), collapse = "\n  "),
       "\nUse the native |> pipe so execution does not depend on package attachment order.",
       call. = FALSE)
}

# Schema-safe row-binding guard. The historical crash came from
# do.call(rbind, ...) on heterogeneous detector/region branch outputs.
# Executable project code must use dplyr::bind_rows() instead.
base_rbind_token <- paste0("do.call(", "rbind")
base_rbind_files <- r_files[vapply(r_files, function(f) {
  lines <- trimws(readLines(f, warn = FALSE))
  code_lines <- lines[!startsWith(lines, "#")]
  any(grepl(base_rbind_token, code_lines, fixed = TRUE))
}, logical(1))]
if (length(base_rbind_files)) {
  stop("Schema-fragile base row-binding call detected in project R file(s):\n  ",
       paste(substring(base_rbind_files, nchar(root) + 2L), collapse = "\n  "),
       "\nUse dplyr::bind_rows() for heterogeneous analytical branches.",
       call. = FALSE)
}

# RF-figure regression contract. Products 1/2 must use an mm environmental
# axis and Products 3/4 a rainfall STA/LTA R(t) axis. RF-product figures are
# anchor-specific: P1/P3 retain Constant TA as the disease anchor; P2/P4 retain
# the Outbreak Threshold. The no-RF NCR run must preserve its earlier TA R(t)
# display structure.
stagea_src <- paste(readLines(file.path(root, "scripts", "StageA_NCR_Lepto_Analysis.R"),
                              warn = FALSE), collapse = "\n")
stageb_src <- paste(readLines(file.path(root, "scripts", "StageB_Regional_Lepto_Analysis.R"),
                              warn = FALSE), collapse = "\n")
rf_plot_tokens <- c(
  "Antecedent rainfall, t-4 to t-1 (mm)",
  "Antecedent rainfall STA/LTA R(t)",
  "Leptospirosis cases" = "#A6CEE3",
  "RF = ",
  "RF threshold lower/upper bounds",
  "RF R(t) lower/upper bounds",
  "make_detection_panel_lepto_only",
  "make_detection_panel_rf",
  "Constant TA trigger",
  "Outbreak threshold (mean + 2SD)",
  'legend.box = "horizontal"'
)
missing_rf_plot_tokens <- rf_plot_tokens[!vapply(rf_plot_tokens, grepl, logical(1),
                                                  x = stagea_src, fixed = TRUE)]
if (length(missing_rf_plot_tokens)) {
  stop("NCR RF/disease-only figure regression check failed. Missing plotting contract token(s): ",
       paste(missing_rf_plot_tokens, collapse = "; "), call. = FALSE)
}
# Regional bootstrap empty-trigger contract. RF gates may legitimately remove
# every trigger for a region/replicate; the bootstrap must still emit one row
# per detector so Method is always present downstream.
bootstrap_schema_tokens <- c(
  "trig_full_override = NULL",
  "Regional bootstrap zero-trigger regression: PASS",
  ".bootstrap_required_cols",
  "expected_boot_rows <- BOOT_N_CI * length(surge_defs)",
  "missing_boot_cols <- setdiff",
  "nrow(.zero_boot) != length(surge_defs)"
)
missing_bootstrap_schema <- bootstrap_schema_tokens[!vapply(
  bootstrap_schema_tokens, grepl, logical(1), x = stageb_src, fixed = TRUE)]
if (length(missing_bootstrap_schema)) {
  stop("Regional RF bootstrap schema regression check failed. Missing token(s): ",
       paste(missing_bootstrap_schema, collapse = "; "), call. = FALSE)
}

if (!grepl("FigureB6A_RF_Overlay_", stageb_src, fixed = TRUE) ||
    !grepl("Antecedent rainfall, t-4 to t-1 (mm)", stageb_src, fixed = TRUE) ||
    !grepl("Antecedent rainfall STA/LTA R(t)", stageb_src, fixed = TRUE) ||
    !grepl("Constant TA trigger", stageb_src, fixed = TRUE) ||
    !grepl('legend.box = "horizontal"', stageb_src, fixed = TRUE)) {
  stop("Regional RF figure regression check failed: anchor-specific mm/R(t) overlay or horizontal legend code is incomplete.",
       call. = FALSE)
}
stagec_src <- paste(readLines(file.path(root, "scripts", "StageC_RF_Threshold_Derivation.R"),
                              warn = FALSE), collapse = "\n")
config_src <- paste(readLines(file.path(root, "R", "00_config.R"),
                              warn = FALSE), collapse = "\n")
rf_bound_tokens <- c(
  "RF_Threshold_Operational_Bounds.csv",
  "RF_Operational_Threshold_Candidate_Grid.csv",
  "RF_Threshold_Bound_Sensitivity.csv",
  "RF_Threshold_HeadToHead_Overall.csv",
  "RF_Threshold_HeadToHead_ExactPairedTests.csv",
  "RF_Threshold_HeadToHead_LORO.csv",
  "RF_Threshold_HeadToHead_Decision.csv",
  "FigureC_NCR_RF_HeadToHead_ConstantTA",
  "FigureC_Regional_RF_HeadToHead_ConstantTA",
  "FigureC_NCR_RF_HeadToHead_OutbreakThreshold",
  "FigureC_Regional_RF_HeadToHead_OutbreakThreshold",
  "pROC::roc.test",
  "stagec_h2h_exact_signflip",
  "Threshold_Lower",
  "Threshold_Upper",
  "Precautionary_Activation_Bound",
  "Strong_Activation_Bound",
  "RF_THRESHOLD_LOWER_PRODUCT1_MM",
  "RF_THRESHOLD_UPPER_PRODUCT4_RT",
  "RF_Threshold_Derivation_Summary_NCR.csv",
  "RF_Threshold_LOYO_Stability_NCR.csv",
  "Higher rainfall activates (>=)",
  "False_Alarm_Rate",
  "RF_WARNING_HORIZON_WEEKS",
  "label_future_anchor_window",
  "RF_Operational_Threshold_Candidate_Grid.csv",
  "window is re-selected inside every LOYO/LORO stability fold",
  "preceding disease triggers"
)
missing_bound_tokens <- rf_bound_tokens[!vapply(
  rf_bound_tokens, function(tok) grepl(tok, stagec_src, fixed = TRUE) ||
    grepl(tok, config_src, fixed = TRUE), logical(1))]
if (length(missing_bound_tokens)) {
  stop("RF lower/point/upper threshold regression check failed. Missing token(s): ",
       paste(missing_bound_tokens, collapse = "; "), call. = FALSE)
}
# Stage C publication figures must keep x-axis labels horizontal; long product
# labels are wrapped with explicit line breaks rather than rotated. Figure
# titles/labels must explicitly state that rainfall precedes/precedes its anchor.
stagec_xaxis_lines <- grep("axis.text.x", strsplit(stagec_src, "\n", fixed = TRUE)[[1]],
                           value = TRUE, fixed = TRUE)
if (!length(stagec_xaxis_lines) || any(!grepl("angle = 0", stagec_xaxis_lines, fixed = TRUE))) {
  stop("Stage C figure layout regression: all x-axis labels must be horizontal (angle = 0).",
       call. = FALSE)
}
if (!grepl("preceding", stagec_src, fixed = TRUE) ||
    !grepl("precedes", stagec_src, fixed = TRUE)) {
  stop("Stage C figure wording regression: figures/methods must explicitly state preceding/precedes.",
       call. = FALSE)
}

# Stage C scope-separated layout regression. NCR and Regional threshold-product
# figures each contain two full-width panels in a 1 x 2 vertical layout; each
# scope-specific head-to-head contains three full-width panels in a 1 x 3
# vertical layout. Mixed-scope Stage C figure exports are not allowed.
stagec_layout_tokens <- c(
  'save_pub("FigureC_NCR_RF_Threshold_Products"',
  'save_pub("FigureC_Regional_RF_Threshold_Products"',
  'list(fig_products_ncr$mm, fig_products_ncr$rt), ncol = 1, nrow = 2',
  'list(fig_products_reg$mm, fig_products_reg$rt), ncol = 1, nrow = 2',
  'patchwork::wrap_plots(list(p_a, p_b, p_c), ncol = 1, nrow = 3',
  'save_pub("FigureC_NCR_RF_STA_LTA_Window_Selection"',
  'save_pub("FigureC_Regional_RF_STA_LTA_Window_Selection"',
  'save_pub("FigureC_NCR_RF_Operational_Threshold_Bounds"',
  'save_pub("FigureC_Regional_RF_Operational_Threshold_Bounds"',
  'save_pub("FigureC_NCR_RF_HeadToHead_ConstantTA"',
  'save_pub("FigureC_Regional_RF_HeadToHead_ConstantTA"',
  'save_pub("FigureC_NCR_RF_HeadToHead_OutbreakThreshold"',
  'save_pub("FigureC_Regional_RF_HeadToHead_OutbreakThreshold"'
)
missing_stagec_layout <- stagec_layout_tokens[!vapply(
  stagec_layout_tokens, grepl, logical(1), x = stagec_src, fixed = TRUE)]
if (length(missing_stagec_layout)) {
  stop("Stage C scope-separated panel-layout regression is incomplete. Missing token(s): ",
       paste(missing_stagec_layout, collapse = "; "), call. = FALSE)
}
for (old_name in c(
  'save_pub("FigureC_RF_Threshold_Products"',
  'save_pub("FigureC_RF_STA_LTA_Window_Selection"',
  'save_pub("FigureC_RF_Operational_Threshold_Bounds"',
  'save_pub("FigureC_RF_HeadToHead_ConstantTA"',
  'save_pub("FigureC_RF_HeadToHead_OutbreakThreshold"')) {
  if (grepl(old_name, stagec_src, fixed = TRUE)) {
    stop("Stage C mixed-scope figure regression: legacy combined export remains: ",
         old_name, call. = FALSE)
  }
}

if (!grepl("classify_rf_activation", config_src, fixed = TRUE) ||
    !grepl("Primary_Metrics_Threshold = \"Point estimate\"", stagea_src, fixed = TRUE) ||
    !grepl("Primary_Metrics_Threshold=\"Point estimate\"", stageb_src, fixed = TRUE)) {
  stop("RF operational-bound contract is incomplete: point estimate must remain primary and upper-tail activation classification must be available.",
       call. = FALSE)
}

if (!grepl("NCR_WEEKLY_CASE_AXIS_MAX <- Y_MAX * 1.10", stagea_src, fixed = TRUE) ||
    !grepl("case_upper <- NCR_WEEKLY_CASE_AXIS_MAX", stagea_src, fixed = TRUE) ||
    !grepl("ylim = c(0, NCR_WEEKLY_CASE_AXIS_MAX)", stagea_src, fixed = TRUE)) {
  stop("NCR y-axis regression: RF and non-RF weekly leptospirosis panels must share the same case-axis ceiling.",
       call. = FALSE)
}

# Legend integrity regression checks. Multi-year NCR figures must use one
# explicitly extracted shared legend rather than data-dependent patchwork guide
# collection; trigger categories must be trained independently of event presence.
legend_contract_tokens <- c(
  "compose_detection_grid_with_legend",
  "cowplot::get_legend",
  "legend_shape_df",
  'legend.box        = "horizontal"',
  "limits = shape_breaks"
)
missing_legend_contract <- legend_contract_tokens[!vapply(
  legend_contract_tokens, grepl, logical(1), x = stagea_src, fixed = TRUE)]
if (length(missing_legend_contract)) {
  stop("NCR shared-legend regression check failed. Missing token(s): ",
       paste(missing_legend_contract, collapse = "; "), call. = FALSE)
}
if (grepl('patchwork::plot_layout(guides = "collect")', stagea_src, fixed = TRUE)) {
  stop("NCR legend regression check failed: data-dependent patchwork guide collection was reintroduced.",
       call. = FALSE)
}
regional_legend_tokens <- c(
  'limits = c("Constant TA", "Continuous TA", "Outbreak Threshold"',
  'breaks = c("Constant TA", "Continuous TA", "Outbreak Threshold"',
  "nrow = 2, byrow = TRUE"
)
missing_regional_legend <- regional_legend_tokens[!vapply(
  regional_legend_tokens, grepl, logical(1), x = stageb_src, fixed = TRUE)]
if (length(missing_regional_legend)) {
  stop("Regional legend regression check failed. Missing stable six-detector legend token(s): ",
       paste(missing_regional_legend, collapse = "; "), call. = FALSE)
}

if (!grepl(".data[[r_col]] * SCALE", stagea_src, fixed = TRUE) ||
    !grepl("Activation threshold", stagea_src, fixed = TRUE)) {
  stop("Disease-only NCR plotting regression: the prior TA R(t)/activation-threshold structure was not restored.",
       call. = FALSE)
}
if (!grepl("plot.subtitle = ggplot2::element_blank()", stagea_src, fixed = TRUE) ||
    !grepl("plot.subtitle = ggplot2::element_blank()", stageb_src, fixed = TRUE)) {
  stop("RF-only export sanitation regression: RF panels must suppress subtitles/captions.",
       call. = FALSE)
}

# Legend-layout regression contract. All publication legends must default to
# bottom/horizontal placement with compact spacing; long legends must wrap
# rather than run beyond the 180-mm canvas.
theme_src <- paste(readLines(file.path(root, "R", "01_publication_theme.R"),
                             warn = FALSE), collapse = "\n")
drift_src <- paste(readLines(file.path(root, "scripts", "Stage2_outbreak_threshold_drift_NCR_Lepto.R"),
                             warn = FALSE), collapse = "\n")
legend_theme_tokens <- c(
  'legend.position      = "bottom"',
  'legend.direction     = "horizontal"',
  'legend.box           = "horizontal"',
  'legend.spacing.x     = grid::unit(5, "pt")'
)
missing_legend_tokens <- legend_theme_tokens[!vapply(
  legend_theme_tokens, grepl, logical(1), x = theme_src, fixed = TRUE)]
if (length(missing_legend_tokens)) {
  stop("Shared publication legend regression: missing compact horizontal legend token(s): ",
       paste(missing_legend_tokens, collapse = "; "), call. = FALSE)
}
if (!grepl('guide_legend(nrow = 2, byrow = TRUE, direction = "horizontal")',
           drift_src, fixed = TRUE)) {
  stop("NCR drift legend regression: the three-entry legend must wrap into horizontal rows.",
       call. = FALSE)
}
if (!grepl('fill = ggplot2::guide_legend(order = 1, nrow = 2',
           stageb_src, fixed = TRUE)) {
  stop("Regional map legend regression: detector keys must wrap to two horizontal rows.",
       call. = FALSE)
}
if (grepl('grid::unit(26, "pt")', stagea_src, fixed = TRUE)) {
  stop("NCR legend regression: oversized 26-pt legend spacing was reintroduced.",
       call. = FALSE)
}

parse_rows <- lapply(r_files,function(f){
  ok<-TRUE; msg<-"OK"
  tryCatch(parse(file=f, encoding="UTF-8"),error=function(e){ok<<-FALSE;msg<<-conditionMessage(e)})
  data.frame(File=substring(f,nchar(root)+2L),Parse_OK=ok,Message=msg)
})
parse_report <- dplyr::bind_rows(parse_rows)
if(any(!parse_report$Parse_OK)) {print(parse_report[!parse_report$Parse_OK,],row.names=FALSE);stop("R parser self-check failed.")}

sheets <- readxl::excel_sheets(DATA_FILE)
if(!SHEET_REGIONAL %in% sheets) stop("Regional Data sheet not found.")
head0 <- readxl::read_excel(DATA_FILE,sheet=SHEET_REGIONAL,n_max=1)
need <- c("REGION","YR","WN","LC_DOH","RF_HDX")
miss <- setdiff(need,names(head0)); if(length(miss))stop("Regional Data missing: ",paste(miss,collapse=", "))

# The supplied workbook must expose all 17 canonical regions for the regional map.
canon17 <- c("BARMM","CAR","MIMAROPA","NCR","REGION I","REGION II","REGION III",
             "REGION IV-A","REGION V","REGION VI","REGION VII","REGION VIII",
             "REGION IX","REGION X","REGION XI","REGION XII","REGION XIII")
reg_check <- readxl::read_excel(DATA_FILE, sheet=SHEET_REGIONAL,
                                col_types="text")
obs_regions <- sort(unique(as.character(reg_check$REGION[!is.na(reg_check$REGION)])))
if (!setequal(obs_regions, canon17) || length(obs_regions) != 17L) {
  stop("Regional Data must contain exactly the 17 canonical regions. Missing: ",
       paste(setdiff(canon17, obs_regions), collapse=", "), call.=FALSE)
}
# Coverage-based regional inclusion must retain every canonical region in each
# evaluable year; low incidence alone is not a reason for exclusion.
reg_cov <- reg_check |>
  dplyr::transmute(REGION=as.character(REGION),
                   YR=suppressWarnings(as.integer(YR)),
                   Cases=suppressWarnings(as.numeric(LC_DOH)),
                   RF=suppressWarnings(as.numeric(RF_HDX))) |>
  dplyr::filter(REGION %in% canon17, YR %in% c(2018L,2019L,2022L,2023L,2024L,2025L)) |>
  dplyr::group_by(REGION,YR) |>
  dplyr::summarise(n_case=sum(is.finite(Cases)), n_rf=sum(is.finite(RF)), .groups="drop")
expected_pairs <- expand.grid(REGION=canon17,
                              YR=c(2018L,2019L,2022L,2023L,2024L,2025L),
                              stringsAsFactors=FALSE)
reg_cov_full <- dplyr::left_join(expected_pairs,reg_cov,by=c("REGION","YR"))
if (any(is.na(reg_cov_full$n_case)) || any(reg_cov_full$n_case < 40L) ||
    any(is.na(reg_cov_full$n_rf)) || any(reg_cov_full$n_rf < 40L)) {
  bad <- reg_cov_full |> dplyr::filter(is.na(n_case) | n_case < 40L | is.na(n_rf) | n_rf < 40L)
  print(as.data.frame(bad),row.names=FALSE)
  stop("Coverage preflight failed: every canonical region/evaluable year must have at least 40 observed case and rainfall weeks.",call.=FALSE)
}
geom <- list.files(DIR_GEOMETRY,pattern="\\.(rds|shp|gpkg|geojson)$",recursive=TRUE,full.names=TRUE,ignore.case=TRUE)
if(!length(geom)) stop("No local Philippine geometry found; bundled geometry is required for an offline reproducible run.", call.=FALSE)

# Dry-run the exact bundled geometry conversion before any expensive analysis.
# This specifically guards against the PackedSpatVector -> st_as_sf failure.
geom_test <- geom[1]
ext <- tolower(tools::file_ext(geom_test))
if (ext == "rds") {
  obj <- readRDS(geom_test)
  if (inherits(obj, "PackedSpatVector") || inherits(obj, "Packed")) {
    obj <- tryCatch(terra::unwrap(obj), error=function(e) e)
    if (inherits(obj, "error")) {
      stop("Bundled geometry could not be unpacked with terra::unwrap(): ",
           conditionMessage(obj), call.=FALSE)
    }
  }
  if (inherits(obj, "SpatVector")) {
    sf_obj <- tryCatch(sf::st_as_sf(obj), error=function(e) e)
  } else if (inherits(obj, "sf")) {
    sf_obj <- obj
  } else {
    sf_obj <- tryCatch(sf::st_as_sf(obj), error=function(e) e)
  }
  if (inherits(sf_obj, "error")) {
    stop("Bundled geometry could not be converted to sf: ",
         conditionMessage(sf_obj), call.=FALSE)
  }
  if (!inherits(sf_obj, "sf") || nrow(sf_obj) < 1L) {
    stop("Bundled geometry conversion produced an empty/non-sf object.", call.=FALSE)
  }
} else {
  sf_obj <- tryCatch(sf::st_read(geom_test, quiet=TRUE), error=function(e) e)
  if (inherits(sf_obj, "error") || !inherits(sf_obj, "sf") || nrow(sf_obj) < 1L) {
    stop("Local geometry could not be read as sf.", call.=FALSE)
  }
}

# Stage H formal head-to-head contract: anchor-focused exact paired inference,
# dominance probability, threshold stability and data-driven preference.
stageh_src <- paste(readLines(file.path(root, "scripts", "StageH_RF_mm_vs_Rt_HeadToHead.R"),
                              warn = FALSE), collapse = "\n")
stageh_tokens <- c(
  "exact_signflip_p",
  "RF_HeadToHead_NCR_ExactPairedTests.csv",
  "RF_HeadToHead_Regional_ClusterPairedTests.csv",
  "RF_HeadToHead_DominanceProbability.csv",
  "RF_HeadToHead_ThresholdStability.csv",
  "RF_HeadToHead_BoundSensitivity.csv",
  "RF_HeadToHead_OverallDecision.csv",
  "No clear empirical superiority"
)
missing_stageh_tokens <- stageh_tokens[!vapply(stageh_tokens, grepl, logical(1),
                                               x = stageh_src, fixed = TRUE)]
if (length(missing_stageh_tokens)) {
  stop("Stage H formal RF head-to-head regression check failed. Missing token(s): ",
       paste(missing_stageh_tokens, collapse = "; "), call. = FALSE)
}

dir.create(file.path(DIR_OUTPUT,"_self_check"),recursive=TRUE,showWarnings=FALSE)
utils::write.csv(parse_report,file.path(DIR_OUTPUT,"_self_check","R_Parse_Report.csv"),row.names=FALSE)
utils::write.csv(data.frame(
  Check=c("Input workbook","Regional Data schema","17-region source coverage","17-region weekly coverage","Local geometry","Geometry conversion"),
  Status=c("OK","OK","OK","OK","OK","OK")),
  file.path(DIR_OUTPUT,"_self_check","Preflight_Report.csv"),row.names=FALSE)
cat("Parsed ",length(r_files)," R files successfully. Required schema and 17-region coverage present. Geometry conversion: OK.\n",sep="")

# RF truth-definition regression: Stage C must use the established A1/A2
# reference and downstream RF gates must evaluate the antecedent t-4..t-1
# feature directly at the disease-trigger week.
stagec_src <- paste(readLines(file.path(root, "scripts", "StageC_RF_Threshold_Derivation.R"),
                              warn = FALSE), collapse = "\n")
if (!grepl("Primary true anchor = first detector-active trigger week in A1 OR A2 per region-year", stagec_src, fixed = TRUE) ||
    !grepl("select_stagec_derivation_inventory", stagec_src, fixed = TRUE) ||
    !grepl("Primary true anchor: first trigger in A1 or A2", stagec_src, fixed = TRUE) ||
    !grepl("False-alarm comparator: trigger outside A1 and A2", stagec_src, fixed = TRUE) ||
    !grepl("Later true trigger: audit only", stagec_src, fixed = TRUE) ||
    !grepl("compute_stagec_a1", stagec_src, fixed = TRUE) ||
    !grepl("compute_stagec_a2", stagec_src, fixed = TRUE) ||
    !grepl("False_Alarm_Weeks_per_Evaluable_Unit", stagec_src, fixed = TRUE)) {
  stop("Stage C RF threshold regression: A1/A2 true/false alarm derivation is incomplete.", call. = FALSE)
}
# P1/P3 must share Constant-TA truth and P2/P4 must share Outbreak-Threshold
# truth. A previous implementation filtered anchor_truth_lookup by Product_ID 3/4,
# which yielded zero rows and caused P3/P4 empirical selection to fail.
rf_anchor_map_tokens <- c(
  "rf_anchor_truth_product_id <- function",
  "if (pid %in% c(1L, 3L)) 1L else 2L",
  "truth_pid <- rf_anchor_truth_product_id(product_id)",
  "Product_ID == truth_pid"
)
missing_rf_anchor_map <- rf_anchor_map_tokens[!vapply(
  rf_anchor_map_tokens, grepl, logical(1), x = stagec_src, fixed = TRUE)]
if (length(missing_rf_anchor_map)) {
  stop("Stage C RF anchor-product mapping regression: P1/P3 and P2/P4 anchor truth mapping is incomplete. Missing token(s): ",
       paste(missing_rf_anchor_map, collapse = "; "), call. = FALSE)
}
if (!grepl("make_direct_q <- function", stagea_src, fixed = TRUE) ||
    !grepl("make_direct_q <- function", stageb_src, fixed = TRUE)) {
  stop("RF gate regression: NCR/Regional reruns must evaluate the antecedent t-4..t-1 feature directly at the disease trigger.", call. = FALSE)
}

# -----------------------------------------------------------------------------
# Independent per-region replication + RF-to-anchor regression contracts
# -----------------------------------------------------------------------------
stagei_path <- file.path(root, "scripts", "StageI_PerRegion_Independent_Analysis.R")
if (!file.exists(stagei_path)) {
  stop("Independent per-region Stage I script is missing.", call. = FALSE)
}
stagei_src <- paste(readLines(stagei_path, warn = FALSE), collapse = "\n")
runall_src <- paste(readLines(file.path(root, "run_all.R"), warn = FALSE), collapse = "\n")
per_region_tokens <- c(
  "TA_ANALYSIS_REGION",
  "TA_ANALYSIS_MODE",
  "IndependentRegion",
  "StageI_PerRegion_Independent",
  "RF_Threshold_Derivation_Summary_PerRegion.csv",
  "17 regions x 5 analyses"
)
missing_per_region <- per_region_tokens[!vapply(
  per_region_tokens,
  function(tok) grepl(tok, stagea_src, fixed = TRUE) ||
    grepl(tok, stagec_src, fixed = TRUE) ||
    grepl(tok, stagei_src, fixed = TRUE) ||
    grepl(tok, config_src, fixed = TRUE),
  logical(1))]
if (length(missing_per_region)) {
  stop("Independent per-region replication regression check failed. Missing token(s): ",
       paste(missing_per_region, collapse = "; "), call. = FALSE)
}
if (!grepl('run_stage("scripts/StageI_PerRegion_Independent_Analysis.R")',
           runall_src, fixed = TRUE)) {
  stop("run_all.R does not execute the independent per-region Stage I branch.", call. = FALSE)
}
if (!grepl('REGION == ANALYSIS_REGION', stagea_src, fixed = TRUE)) {
  stop("Stage A is not parameterized to filter the requested independent region.", call. = FALSE)
}

rf_anchor_tokens <- c(
  "RF_to_Anchor_Trigger_Detail.csv",
  "RF_to_Anchor_Regional_Summary.csv",
  "FigureB7_RF_to_",
  "RF_Supported_PPV",
  "N_RF_Supported_True_Alarms",
  "N_RF_Supported_False_Alarms",
  '"Constant TA" else "Outbreak Threshold"'
)
missing_anchor <- rf_anchor_tokens[!vapply(rf_anchor_tokens, grepl, logical(1),
                                           x = stageb_src, fixed = TRUE)]
if (length(missing_anchor)) {
  stop("Regional RF-to-anchor regression check failed. Missing structural token(s): ",
       paste(missing_anchor, collapse = "; "), call. = FALSE)
}
# Do not regression-test human-facing legend prose: line wrapping is allowed.
# Test the actual RF-to-anchor data contract and product mapping instead.
if (!grepl("RF_GATE$product %in% c(1L, 3L)", stageb_src, fixed = TRUE) ||
    !grepl("ungated_surge_sta_lta_vaezi", stageb_src, fixed = TRUE) ||
    !grepl("ungated_surge_mean_2sd", stageb_src, fixed = TRUE) ||
    !grepl("RF_Supported & IsTrue_A1A2", stageb_src, fixed = TRUE) ||
    !grepl("RF_Supported & !IsTrue_A1A2", stageb_src, fixed = TRUE)) {
  stop("Regional RF-to-anchor regression check failed: P1/P3 -> Constant TA, P2/P4 -> Outbreak Threshold, and A1/A2 true/false RF-support fields must all be present.",
       call. = FALSE)
}
