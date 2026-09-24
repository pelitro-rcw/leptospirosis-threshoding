# =============================================================================
# STAGE I - INDEPENDENT NCR-STYLE ANALYSIS FOR EACH OF THE 17 REGIONS
# =============================================================================
# This branch is separate from the pooled Stage B Regional analysis.
# Each region receives DiseaseOnly plus RF Products 1-4 using its own Stage C
# region-specific RF threshold triplet.
# =============================================================================

.run_stage_i_per_region <- function() {
  root <- Sys.getenv("TA_PROJECT_ROOT", unset = normalizePath(getwd(), winslash = "/"))
  source(file.path(root, "R", "00_config.R"), local = TRUE)

  regions <- CANONICAL_17
  if (length(regions) != 17L) stop("Stage I requires the canonical 17-region list.", call. = FALSE)
  threshold_file <- file.path(DIR_OUTPUT, "StageC_RF_Threshold_Derivation",
                              "RF_Threshold_Derivation_Summary_PerRegion.csv")
  if (!file.exists(threshold_file)) {
    stop("Stage I requires independent per-region RF thresholds from Stage C: ",
         threshold_file, call. = FALSE)
  }
  thresholds <- utils::read.csv(threshold_file, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("REGION", "Product_ID", "Adopted_Threshold") %in% names(thresholds)) ||
      nrow(thresholds) != 68L || length(unique(thresholds$REGION)) != 17L ||
      any(table(thresholds$REGION) != 4L)) {
    stop("Stage I threshold handoff must contain exactly 4 products x 17 regions.", call. = FALSE)
  }

  env_names <- c("TA_ANALYSIS_REGION", "TA_ANALYSIS_MODE",
                 "TA_RF_GATE_PRODUCT", "TA_OUTPUT_SUFFIX")
  old_env <- Sys.getenv(env_names, unset = NA_character_)
  names(old_env) <- env_names
  on.exit({
    for (nm in env_names) {
      val <- old_env[[nm]]
      if (is.na(val)) {
        Sys.unsetenv(nm)
      } else {
        args <- setNames(list(val), nm)
        do.call(Sys.setenv, args)
      }
    }
  }, add = TRUE)

  manifest <- list()
  for (reg in regions) {
    safe_reg <- gsub("[^A-Za-z0-9]+", "_", reg)
    safe_reg <- gsub("^_+|_+$", "", safe_reg)
    for (product in 0:4) {
      Sys.setenv(TA_ANALYSIS_REGION = reg,
                 TA_ANALYSIS_MODE = "IndependentRegion",
                 TA_RF_GATE_PRODUCT = as.character(product),
                 TA_OUTPUT_SUFFIX = "")
      label <- if (product == 0L) "DiseaseOnly" else paste0("RF_Product", product)
      cat("\n[Stage I] ", reg, " | ", label, "\n", sep = "")
      tryCatch(
        sys.source(file.path(DIR_SCRIPTS, "StageA_NCR_Lepto_Analysis.R"),
                   envir = new.env(parent = globalenv())),
        error = function(e) {
          stop("Independent per-region analysis failed for ", reg, " / ", label,
               ": ", conditionMessage(e), call. = FALSE)
        }
      )
      out_dir <- file.path(DIR_OUTPUT, "StageI_PerRegion_Independent", safe_reg, label)
      if (!dir.exists(out_dir)) {
        stop("Stage I expected output directory was not created: ", out_dir, call. = FALSE)
      }
      .pair_exists <- function(stem) {
        all(file.exists(file.path(out_dir, paste0(stem, c(".pdf", ".png")))))
      }
      .required_files <- c(
        file.path(out_dir, "FigureA2_method_summary_primary.csv"),
        file.path(out_dir, "evaluation_framework", "Table2_Per_Year_Detail.csv")
      )
      .required_pairs <- c("FigureA2_panel_a_DominanceMatrix", "FigureA2_panel_b_Scatter")
      if (product == 0L || product %in% c(1L, 3L)) {
        .required_pairs <- c(.required_pairs, "FigureA3_Comparison_2019_2024_Constant_TA")
      } else {
        .required_pairs <- c(.required_pairs, "FigureA3_Comparison_2019_2024_Outbreak_Threshold")
      }
      if (product > 0L) {
        .required_files <- c(
          .required_files,
          file.path(out_dir, "RF_Gate_Metadata.csv"),
          file.path(out_dir, paste0("RF_Operational_Weekly_", safe_reg, ".csv"))
        )
      }
      .missing_files <- .required_files[!file.exists(.required_files)]
      .missing_pairs <- .required_pairs[!vapply(.required_pairs, .pair_exists, logical(1))]
      if (length(.missing_files) || length(.missing_pairs)) {
        stop("Stage I output contract failed for ", reg, " / ", label,
             ". Missing file(s): ", paste(basename(.missing_files), collapse = ", "),
             if (length(.missing_pairs)) paste0("; missing PDF/PNG pair(s): ", paste(.missing_pairs, collapse = ", ")) else "",
             call. = FALSE)
      }
      manifest[[length(manifest) + 1L]] <- data.frame(
        REGION = reg, Product = product, Analysis = label,
        Output_Directory = out_dir, Completed = TRUE,
        stringsAsFactors = FALSE)
    }
  }
  manifest_df <- dplyr::bind_rows(manifest)
  if (nrow(manifest_df) != 85L || any(!manifest_df$Completed)) {
    stop("Stage I must complete 17 regions x 5 analyses = 85 independent runs.", call. = FALSE)
  }
  manifest_dir <- file.path(DIR_OUTPUT, "StageI_PerRegion_Independent")
  readr::write_csv(manifest_df,
                   file.path(manifest_dir, "PerRegion_Independent_Run_Manifest.csv"), na = "")
  cat("\nStage I complete: 17 regions x 5 analyses (85 independent runs).\n")
  invisible(manifest_df)
}

.run_stage_i_per_region()
