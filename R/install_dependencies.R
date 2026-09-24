# =============================================================================
# DEPENDENCY INSTALLER  --  run once on a fresh R installation
# -----------------------------------------------------------------------------
#   source("R/install_dependencies.R")
#
# Deliberately separated from the analysis scripts. The originals called
# install.packages() inline at run time, which silently mutated the user's
# library, stalled on offline machines and made non-interactive execution
# unreliable. Installation is now an explicit, auditable one-off step.
# =============================================================================

REPO <- getOption("repos")
if (is.null(REPO[["CRAN"]]) || REPO[["CRAN"]] == "@CRAN@") {
  REPO <- c(CRAN = "https://cloud.r-project.org")
}

# --- Required ----------------------------------------------------------------
pkgs_required <- c(
  # data handling
  "readxl", "dplyr", "tidyr", "purrr", "tibble", "stringr", "readr",
  "zoo", "ISOweek", "scales", "rlang",
  # statistics
  "pROC", "MASS", "boot", "surveillance",
  # graphics
  # ggrepel >= 0.9.0 required: Stage 4 map labels use bg.colour/bg.r halos
  "ggplot2", "ggrepel", "patchwork", "cowplot", "viridisLite", "grid",
  # spatial (Stage B map). terra is required because the bundled GADM RDS may
  # deserialize as PackedSpatVector and must be unpacked before sf conversion.
  "sf", "terra", "ggspatial", "prettymapr",
  # project utilities
  "here"
)

# --- Optional ----------------------------------------------------------------
# Geometry sources for the Stage 4 map. Only needed if you have not placed a
# local boundary file in data/geometry/.
pkgs_optional <- c("geodata", "rnaturalearth", "rnaturalearthdata")

install_if_missing <- function(pkgs, label) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) == 0) {
    message("[deps] All ", label, " packages already present.")
    return(invisible(character(0)))
  }
  message("[deps] Installing ", label, ": ", paste(missing, collapse = ", "))
  utils::install.packages(missing, repos = REPO, dependencies = TRUE)
  still <- missing[!vapply(missing, requireNamespace, logical(1), quietly = TRUE)]
  if (length(still) > 0) {
    warning("[deps] Failed to install: ", paste(still, collapse = ", "),
            call. = FALSE)
  }
  invisible(still)
}

message("R version: ", R.version.string)
if (getRversion() < "4.1.0") {
  stop("R >= 4.1 is required (the scripts use the native |> pipe).", call. = FALSE)
}

failed_required <- install_if_missing(pkgs_required, "required")
if (length(failed_required) > 0L) {
  stop("Required package installation failed: ",
       paste(failed_required, collapse = ", "),
       ". Resolve these packages before running run_all.R.", call. = FALSE)
}
install_if_missing(pkgs_optional, "optional (map geometry)")

# --- System library note -----------------------------------------------------
# 'sf' requires GDAL, GEOS and PROJ. If its installation fails:
#   Ubuntu/Debian : sudo apt install libgdal-dev libgeos-dev libproj-dev libudunits2-dev
#   macOS (brew)  : brew install gdal geos proj udunits
#   Windows       : binary CRAN packages bundle these; no action needed.
if ("sf" %in% failed_required) {
  message("\n[deps] 'sf' failed to build. Install the GDAL/GEOS/PROJ system ",
          "libraries listed at the foot of this script, then re-run.")
}

message("\n[deps] Done. Next: source(\"run_all.R\")")
