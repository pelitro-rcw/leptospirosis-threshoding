# =============================================================================
# RUN ALL - LEPTOSPIROSIS ANALYTICAL PIPELINE
# =============================================================================
root <- normalizePath(getwd(), winslash="/")
if(!dir.exists(file.path(root,"R")) || !dir.exists(file.path(root,"scripts")))
  stop("Run source('run_all.R') from the project root.", call.=FALSE)

# Bootstrap syntax gate. This runs before any stage is sourced, so a syntax
# error inside 00_self_check.R (or any later script) cannot prevent the parser
# audit from running. Only base R is required here.
.parse_targets <- unique(c(
  file.path(root, "run_all.R"),
  list.files(file.path(root, "R"), pattern="\\.R$", full.names=TRUE, recursive=TRUE),
  list.files(file.path(root, "scripts"), pattern="\\.R$", full.names=TRUE, recursive=TRUE)
))
.parse_failures <- lapply(.parse_targets, function(f) {
  tryCatch({
    parse(file=f, encoding="UTF-8")
    NULL
  }, error=function(e) {
    list(File=substring(f, nchar(root) + 2L), Error=conditionMessage(e))
  })
})
.parse_failures <- Filter(Negate(is.null), .parse_failures)
if (length(.parse_failures)) {
  .parse_report <- data.frame(
    File=vapply(.parse_failures, function(x) x$File, character(1)),
    Error=vapply(.parse_failures, function(x) x$Error, character(1)),
    stringsAsFactors=FALSE
  )
  print(.parse_report, row.names=FALSE)
  stop("Bootstrap R syntax validation failed. No analytical stage was executed.",
       call.=FALSE)
}
rm(.parse_targets, .parse_failures)
cat("Bootstrap syntax validation: all project R files parsed successfully.\n")
Sys.setenv(TA_PROJECT_ROOT=root, TA_RF_GATE_PRODUCT="0", TA_OUTPUT_SUFFIX="", TA_ANALYSIS_REGION="NCR", TA_ANALYSIS_MODE="NCR")

# Start from a clean output tree by default so stale files from an interrupted
# earlier run cannot satisfy the final validator. Set TA_CLEAN_OUTPUTS=0 only
# when intentionally preserving prior outputs.
clean_outputs <- tolower(trimws(Sys.getenv("TA_CLEAN_OUTPUTS", unset="1"))) %in%
  c("1", "true", "yes", "y")
out_dir <- file.path(root, "outputs")
if (clean_outputs && dir.exists(out_dir)) {
  resolved_root <- normalizePath(root, winslash="/", mustWork=TRUE)
  resolved_out <- normalizePath(out_dir, winslash="/", mustWork=TRUE)
  if (!startsWith(resolved_out, paste0(resolved_root, "/"))) {
    stop("Refusing to clean outputs outside the project root.", call.=FALSE)
  }
  unlink(out_dir, recursive=TRUE, force=TRUE)
}
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

log_dir <- file.path(out_dir, "_run_logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
.warning_rows <- list()
.error_rows <- list()
.flush_run_logs <- function() {
  if (length(.warning_rows)) {
    w <- data.frame(
      Stage = vapply(.warning_rows, function(x) x$Stage, character(1)),
      Warning = vapply(.warning_rows, function(x) x$Warning, character(1)),
      stringsAsFactors = FALSE
    )
  } else {
    w <- data.frame(Stage=character(), Warning=character(), stringsAsFactors=FALSE)
  }
  if (length(.error_rows)) {
    e <- data.frame(
      Stage = vapply(.error_rows, function(x) x$Stage, character(1)),
      Error = vapply(.error_rows, function(x) x$Error, character(1)),
      stringsAsFactors = FALSE
    )
  } else {
    e <- data.frame(Stage=character(), Error=character(), stringsAsFactors=FALSE)
  }
  utils::write.csv(w, file.path(log_dir, "Warnings.csv"), row.names=FALSE)
  utils::write.csv(e, file.path(log_dir, "Errors.csv"), row.names=FALSE)
}

run_stage <- function(path) {
  cat("\n>>> ", path, "\n", sep="")
  stage_warnings <- character()
  tryCatch(
    withCallingHandlers(
      sys.source(file.path(root, path), envir=new.env(parent=globalenv())),
      warning=function(w) {
        stage_warnings <<- c(stage_warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error=function(e) {
      if (length(stage_warnings)) {
        for (msg in unique(stage_warnings)) {
          .warning_rows[[length(.warning_rows)+1L]] <<- list(Stage=path, Warning=msg)
        }
      }
      .error_rows[[length(.error_rows)+1L]] <<- list(
        Stage=path, Error=conditionMessage(e))
      .flush_run_logs()
      stop("Pipeline failed in ", path, ": ", conditionMessage(e), call.=FALSE)
    }
  )
  if (length(stage_warnings)) {
    for (msg in unique(stage_warnings)) {
      .warning_rows[[length(.warning_rows)+1L]] <<- list(Stage=path, Warning=msg)
    }
    cat("[warnings] ", length(unique(stage_warnings)),
        " unique warning(s) recorded in outputs/_run_logs/Warnings.csv\n", sep="")
  }
  .flush_run_logs()
  invisible(TRUE)
}
run_stage("scripts/00_self_check.R")
run_stage("scripts/Stage1_eta_threshold_derivation_NCR.R")
run_stage("scripts/Stage2_outbreak_threshold_drift_NCR_Lepto.R")
run_stage("scripts/StageA_NCR_Lepto_Analysis.R")
run_stage("scripts/StageB_Regional_Lepto_Analysis.R")
run_stage("scripts/StageC_RF_Threshold_Derivation.R")
run_stage("scripts/StageI_PerRegion_Independent_Analysis.R")
run_stage("scripts/StageD_RF_Product1_NCR_Regional.R")
run_stage("scripts/StageE_RF_Product2_NCR_Regional.R")
run_stage("scripts/StageF_RF_Product3_NCR_Regional.R")
run_stage("scripts/StageG_RF_Product4_NCR_Regional.R")
Sys.setenv(TA_RF_GATE_PRODUCT="0", TA_OUTPUT_SUFFIX="", TA_ANALYSIS_REGION="NCR", TA_ANALYSIS_MODE="NCR")
run_stage("scripts/StageH_RF_mm_vs_Rt_HeadToHead.R")
run_stage("scripts/99_validate_outputs.R")
 .flush_run_logs()
cat("\nAll scheduled stages completed and required outputs validated. See outputs/.\n")
