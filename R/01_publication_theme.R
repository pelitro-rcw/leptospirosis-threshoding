# =============================================================================
# SHARED PUBLICATION THEME  --  Nature Communications
# -----------------------------------------------------------------------------
# Single source of truth for typography, sizing, palettes, save helpers and
# map cartography across Stages 1-5. Sourced by every stage script.
#
# DESIGN PRINCIPLE -- "author at final print size"
# ---------------------------------------------------------------------------
# The original scripts exported figures at 11-20 inches wide. Nature
# Communications reproduces figures at a maximum of 180 mm (7.09 in), so those
# files were reduced by 0.35x-0.64x during production, shrinking every font by
# the same factor (a nominal 9 pt axis label printed at ~3-5 pt). That
# down-scaling -- not the nominal font settings -- is why the figures read as
# illegible.
#
# The fix is therefore twofold: author each figure at its true final width AND
# set type in real points. Nominal base_size drops (9 -> 8) while EFFECTIVE
# printed size rises substantially, because no reduction is applied.
#
# All sizes below are final printed points. Nature Communications requires a
# 5 pt minimum; nothing here falls below 6 pt.
# =============================================================================


# -----------------------------------------------------------------------------
# 1. CANVAS DIMENSIONS (Nature Communications)
# -----------------------------------------------------------------------------
MM_PER_IN <- 25.4

NC_W_SINGLE <- 88  / MM_PER_IN   # 3.46 in - single column
NC_W_MEDIUM <- 120 / MM_PER_IN   # 4.72 in - 1.5 column
NC_W_DOUBLE <- 180 / MM_PER_IN   # 7.09 in - double column (maximum width)
NC_H_MAX    <- 170 / MM_PER_IN   # 6.69 in - main-text figure height cap
# Supplementary / multi-page figures may run longer than a single text column.
# The 13-panel Figure 3 grid and the 8-metric head-to-head figures cannot be
# rendered legibly inside 170 mm: at that height each sub-panel gets ~30 mm,
# which is what forced the crowding. They are authored taller and are captioned
# as full-page/supplementary figures.
NC_H_SUPP   <- 240 / MM_PER_IN   # 9.45 in - dense multipanel / supplementary

#' Clamp a figure to the journal's printable area.
#' Guarantees no figure is ever silently down-scaled by production.
#' @param max_height Height cap. Defaults to the main-text limit; pass
#'   NC_H_SUPP for dense multipanel or supplementary figures.
nc_fit <- function(width, height, max_height = NC_H_MAX) {
  if (width > NC_W_DOUBLE) {
    height <- height * (NC_W_DOUBLE / width)
    width  <- NC_W_DOUBLE
  }
  if (height > max_height) {
    width  <- width * (max_height / height)
    height <- max_height
  }
  list(width = width, height = height)
}


# -----------------------------------------------------------------------------
# 2. TYPOGRAPHY
# -----------------------------------------------------------------------------
# TYPE_SCALE is the single tuning knob for the whole manuscript. If a reviewer
# asks for larger type again, raise this one value and re-run run_all.R; every
# figure in every stage responds consistently.
TYPE_SCALE <- getOption("ta.type_scale", 1.0)

PUB_BASE     <- 8.0 * TYPE_SCALE   # base_size for every theme
PUB_AXIS_TIT <- 8.0 * TYPE_SCALE   # axis titles
PUB_AXIS_TXT <- 7.0 * TYPE_SCALE   # axis tick labels
PUB_LEG_TIT  <- 8.0 * TYPE_SCALE   # legend titles
PUB_LEG_TXT  <- 7.0 * TYPE_SCALE   # legend text
PUB_STRIP    <- 8.0 * TYPE_SCALE   # facet strip labels
PUB_TITLE    <- 9.0 * TYPE_SCALE   # plot titles
PUB_SUBTITLE <- 7.5 * TYPE_SCALE   # plot subtitles
PUB_ANNOT    <- 6.5 * TYPE_SCALE   # in-panel annotations
PUB_TAG      <- 9.0 * TYPE_SCALE   # panel labels (a, b, c ...)
PUB_CAPTION  <- 6.5 * TYPE_SCALE   # footnotes / captions

# geom_text/geom_label take size in mm, not points. Always convert.
PT_TO_MM <- 1 / .pt
pub_text_size <- function(pt = PUB_ANNOT) pt * PT_TO_MM

# Font family: prefer Helvetica/Arial (NC house style), fall back gracefully.
.pub_font <- function() {
  fonts <- names(grDevices::pdfFonts())
  if ("Helvetica" %in% fonts) return("Helvetica")
  if ("Arial"     %in% fonts) return("Arial")
  "sans"
}
PUB_FAMILY <- .pub_font()


# -----------------------------------------------------------------------------
# 3. LINE AND POINT GEOMETRY
# -----------------------------------------------------------------------------
# Journal production cannot reliably reproduce strokes below ~0.25 pt. The
# originals used 0.25-0.35 linewidth on a canvas later reduced by up to 0.35x,
# putting effective strokes near 0.09 pt (hairlines that drop out in print).
PUB_LW_THIN   <- 0.30   # gridlines, panel borders
PUB_LW_BASE   <- 0.50   # standard data lines
PUB_LW_EMPH   <- 0.80   # emphasised / threshold lines
PUB_PT_SMALL  <- 0.90   # dense scatter
PUB_PT_BASE   <- 1.60   # standard points
PUB_PT_EMPH   <- 2.40   # highlighted points
PUB_STROKE    <- 0.40   # point outline


# -----------------------------------------------------------------------------
# 4. COLOURBLIND-SAFE PALETTES
# -----------------------------------------------------------------------------
# Okabe-Ito: the standard 8-colour qualitative palette, distinguishable under
# deuteranopia, protanopia and tritanopia, and separable in greyscale.
OKABE_ITO <- c(
  black      = "#000000",
  orange     = "#E69F00",
  sky_blue   = "#56B4E9",
  green      = "#009E73",
  yellow     = "#F0E442",
  blue       = "#0072B2",
  vermillion = "#D55E00",
  purple     = "#CC79A7"
)

# Detector palette - fixed assignment so colours mean the same thing in every
# figure and every stage. Chosen for maximum separation of the three detectors.
PAL_DETECTOR <- c(
  "Constant Transmission Acceleration"   = "#0072B2",  # blue
  "Continuous Transmission Acceleration" = "#E69F00",  # orange
  "Outbreak Threshold"                   = "#009E73",  # green
  "No consensus"                         = "#9A9A9A",  # neutral grey
  "Contested"                            = "#9A9A9A"
)
# Short aliases used by some panels
PAL_DETECTOR_SHORT <- c(
  "Constant TA"        = "#0072B2",
  "Continuous TA"      = "#E69F00",
  "Outbreak Threshold" = "#009E73",
  "Contested"          = "#9A9A9A",
  "No consensus"       = "#9A9A9A"
)

# Sequential ramp for continuous fills (viridis is colourblind-safe and
# perceptually uniform; used where a magnitude, not a category, is shown).
PAL_SEQ <- function(n) viridisLite::viridis(n, option = "D", begin = 0.08, end = 0.95)

# Diverging ramp for signed quantities, safe for red-green deficiency.
PAL_DIV <- c(low = "#B2182B", mid = "#F7F7F7", high = "#2166AC")

# Neutral greys used consistently for structure
GREY_TEXT   <- "grey20"
GREY_AXIS   <- "grey35"
GREY_GRID   <- "grey88"
GREY_BORDER <- "grey55"


# -----------------------------------------------------------------------------
# 5. CORE THEME
# -----------------------------------------------------------------------------
#' Unified publication theme.
#'
#' @param base_size Base font size in points at final print size.
#' @param grid      "both", "x", "y" or "none".
#' @param border    Draw a panel border.
theme_pub <- function(base_size = PUB_BASE,
                      base_family = PUB_FAMILY,
                      grid = "both",
                      border = TRUE) {

  s <- base_size / PUB_BASE   # proportional scaling of the type hierarchy

  th <- ggplot2::theme_minimal(base_size = base_size, base_family = base_family) %+replace%
    ggplot2::theme(
      # --- backgrounds -------------------------------------------------------
      # White (not #FAFAFA): grey panel fills lose contrast in CMYK print.
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background  = ggplot2::element_rect(fill = "white", colour = NA),

      # --- grid --------------------------------------------------------------
      # NO GRIDLINES ANYWHERE. Set in the base theme rather than patched per
      # figure, so a plot cannot reintroduce them by forgetting a modifier.
      panel.grid       = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),

      # --- axes --------------------------------------------------------------
      # Axis furniture is bold BLACK in the base theme for the same reason.
      axis.title   = ggplot2::element_text(size = PUB_AXIS_TIT * s,
                                           face = "bold", colour = "black"),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 4)),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 4),
                                           angle = 90),
      axis.text    = ggplot2::element_text(size = PUB_AXIS_TXT * s,
                                           face = "bold", colour = "black"),
      axis.text.x  = ggplot2::element_text(size = PUB_AXIS_TXT * s,
                                           face = "bold", colour = "black"),
      axis.text.y  = ggplot2::element_text(size = PUB_AXIS_TXT * s,
                                           face = "bold", colour = "black"),
      axis.text.x.top   = ggplot2::element_text(size = PUB_AXIS_TXT * s,
                                                face = "bold", colour = "black"),
      axis.text.y.right = ggplot2::element_text(size = PUB_AXIS_TXT * s,
                                                face = "bold", colour = "black"),
      axis.line    = ggplot2::element_line(colour = "black", linewidth = 0.6),
      axis.ticks   = ggplot2::element_line(colour = "black", linewidth = 0.6),
      axis.ticks.length = grid::unit(2.2, "pt"),

      # --- legend ------------------------------------------------------------
      # Publication-safe default: all legends sit beneath the plotting panel,
      # read left-to-right, and use compact keys/spacing. Figure-specific
      # guide_legend(nrow=...) calls may wrap long keys into two or three
      # horizontal rows without changing this common alignment.
      legend.position      = "bottom",
      legend.direction     = "horizontal",
      legend.box           = "horizontal",
      legend.justification = "center",
      legend.box.just      = "center",
      legend.background    = ggplot2::element_blank(),
      legend.key           = ggplot2::element_blank(),
      legend.title         = ggplot2::element_text(size = PUB_LEG_TIT * s,
                                                   face = "bold",
                                                   colour = GREY_TEXT,
                                                   margin = ggplot2::margin(r = 5)),
      legend.text          = ggplot2::element_text(size = PUB_LEG_TXT * s,
                                                   colour = GREY_TEXT,
                                                   lineheight = 0.92,
                                                   margin = ggplot2::margin(l = 2, r = 6)),
      legend.key.size      = grid::unit(8.5, "pt"),
      legend.key.width     = grid::unit(11, "pt"),
      legend.key.height    = grid::unit(8.5, "pt"),
      legend.spacing.x     = grid::unit(5, "pt"),
      legend.spacing.y     = grid::unit(2, "pt"),
      legend.margin        = ggplot2::margin(3, 6, 3, 6),
      legend.box.spacing   = grid::unit(4, "pt"),
      legend.text.align    = 0,

      # --- facets ------------------------------------------------------------
      strip.background = ggplot2::element_rect(fill = "grey94", colour = NA),
      strip.text       = ggplot2::element_text(size = PUB_STRIP * s,
                                               face = "bold", colour = GREY_TEXT,
                                               margin = ggplot2::margin(3, 3, 3, 3)),

      # --- titles ------------------------------------------------------------
      plot.title    = ggplot2::element_text(size = PUB_TITLE * s, face = "bold",
                                            hjust = 0, colour = GREY_TEXT,
                                            margin = ggplot2::margin(b = 4)),
      plot.subtitle = ggplot2::element_text(size = PUB_SUBTITLE * s, hjust = 0,
                                            colour = GREY_AXIS,
                                            margin = ggplot2::margin(b = 5)),
      plot.caption  = ggplot2::element_text(size = PUB_CAPTION * s, hjust = 0,
                                            colour = GREY_AXIS,
                                            margin = ggplot2::margin(t = 5)),
      plot.tag      = ggplot2::element_text(size = PUB_TAG * s, face = "bold",
                                            colour = "black",
                                            hjust = 0, vjust = 1,
                                            margin = ggplot2::margin(b = 2)),
      # ALIGNMENT. plot.tag.position = c(0, 1) anchors the letter to the outer
      # PLOT corner, but plot.title.position defaults to "panel", which indents
      # the title past the y-axis labels and axis title. The letter therefore
      # sat well to the left of, and above, the title it belongs to.
      # Setting both to the "plot" reference puts them on the same left edge.
      plot.tag.position   = "topleft",
      plot.title.position = "plot",
      plot.caption.position = "plot",

      # ANTI-CROPPING. Several panels use coord_cartesian(clip = "off") so that
      # arrows, trigger labels and panel letters can sit outside the panel. With
      # a tight plot margin those render clipped at the canvas edge. These
      # values give every side room for an out-of-panel annotation; the top is
      # largest because that is where the arrow rows and tags live.
      plot.margin = ggplot2::margin(t = 14, r = 16, b = 16, l = 12)
    )

  # Panel border suppressed: with bold black axis lines drawn, a border would
  # double them up. `border` is kept for call-site compatibility.
  th <- th %+replace% ggplot2::theme(panel.border = ggplot2::element_blank())

  # `grid` is retained as an argument for call-site compatibility but is now a
  # no-op: gridlines are removed unconditionally above.
  th
}

#' Theme for map panels: no axes, no grid, room for cartographic furniture.
theme_pub_map <- function(base_size = PUB_BASE, base_family = PUB_FAMILY) {
  theme_pub(base_size, base_family, grid = "none", border = FALSE) %+replace%
    ggplot2::theme(
      axis.text   = ggplot2::element_blank(),
      axis.title  = ggplot2::element_blank(),
      axis.ticks  = ggplot2::element_blank(),
      panel.grid  = ggplot2::element_blank(),
      # Map panel background is plain white per revision (was "#F4F8FB").
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.margin = ggplot2::margin(4, 4, 4, 4)
    )
}

#' Theme for table-like grob panels drawn with geom_tile/geom_text.
theme_pub_table <- function(base_size = PUB_BASE, base_family = PUB_FAMILY) {
  theme_pub(base_size, base_family, grid = "none", border = FALSE) %+replace%
    ggplot2::theme(
      axis.text  = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      axis.line  = ggplot2::element_blank(),
      # Table/heatmap panels carry in-cell text and often rotated row and column
      # labels that sit outside the plotting rectangle. Generous margins stop
      # those being cut at the canvas edge.
      plot.margin = ggplot2::margin(t = 14, r = 18, b = 14, l = 16)
    )
}

# Applied globally so any plot built without an explicit theme still complies.
ggplot2::theme_set(theme_pub())


# -----------------------------------------------------------------------------
# 6. PANEL TAGGING
# -----------------------------------------------------------------------------
# Nature Communications uses lower-case bold panel letters. Centralised here so
# tag style is identical across stages.
pub_tag_layout <- function(tag_levels = "a") {
  patchwork::plot_annotation(
    tag_levels = tag_levels,
    theme = ggplot2::theme(
      plot.tag = ggplot2::element_text(size = PUB_TAG, face = "bold",
                                       family = PUB_FAMILY, colour = "black",
                                       hjust = 0, vjust = 1,
                                       margin = ggplot2::margin(b = 2)),
      plot.tag.position     = "topleft",
      plot.title.position   = "plot",
      plot.caption.position = "plot"
    )
  )
}


# -----------------------------------------------------------------------------
# 7. DEVICE AND SAVE HELPERS
# -----------------------------------------------------------------------------
#' Cairo PDF where available (embeds fonts, honours transparency), else base pdf.
safe_pdf_device <- function() {
  ok <- tryCatch(isTRUE(capabilities("cairo")), error = function(e) FALSE)
  if (ok) grDevices::cairo_pdf else grDevices::pdf
}

#' Save one figure as both vector PDF (submission) and 600 dpi PNG (preview).
#'
#' Dimensions are clamped to the journal's printable area, so a figure can
#' never be reduced -- and its type never shrunk -- during production.
#'
#' @param stem  Filename without extension.
#' @param plot  ggplot/patchwork object.
#' @param width,height Inches, at final print size.
#' @param dir   Output directory.
#' @param dpi   Raster resolution (>= 300 required; 600 used throughout).
save_pub <- function(stem, plot, width, height, dir, dpi = 600,
                     max_height = NC_H_MAX) {
  fit <- nc_fit(width, height, max_height = max_height)
  if (abs(fit$width - width) > 1e-6 || abs(fit$height - height) > 1e-6) {
    message(sprintf(
      "  [size] '%s' clamped %.2fx%.2f -> %.2fx%.2f in to fit %.0f x %.0f mm",
      stem, width, height, fit$width, fit$height,
      NC_W_DOUBLE * MM_PER_IN, NC_H_MAX * MM_PER_IN))
  }
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)

  ggplot2::ggsave(file.path(dir, paste0(stem, ".pdf")), plot,
                  device = safe_pdf_device(),
                  width = fit$width, height = fit$height,
                  units = "in", bg = "white", limitsize = FALSE)

  ggplot2::ggsave(file.path(dir, paste0(stem, ".png")), plot,
                  device = "png", dpi = dpi,
                  width = fit$width, height = fit$height,
                  units = "in", bg = "white", limitsize = FALSE)

  message("  [saved] ", stem, ".pdf / .png  (",
          sprintf("%.2f x %.2f in @ %d dpi", fit$width, fit$height, dpi), ")")
  invisible(file.path(dir, paste0(stem, c(".pdf", ".png"))))
}


# -----------------------------------------------------------------------------
# 8. MAP CARTOGRAPHY HELPERS
# -----------------------------------------------------------------------------
# Addresses the reviewer's request for standard cartographic elements. Sizes
# are tuned for a double-column map at final print size.

#' Scale bar sized and styled for print.
pub_scalebar <- function(location = "bl",
                         width_hint = 0.22,
                         text_cex = PUB_ANNOT / 10) {
  ggspatial::annotation_scale(
    location   = location,
    width_hint = width_hint,
    height     = grid::unit(0.14, "cm"),
    bar_cols   = c("grey20", "white"),
    line_width = 0.5,
    line_col   = "grey20",
    text_col   = GREY_TEXT,
    text_cex   = text_cex,
    text_family = PUB_FAMILY,
    pad_x      = grid::unit(0.30, "cm"),
    pad_y      = grid::unit(0.30, "cm")
  )
}

#' North arrow sized and styled for print.
pub_north_arrow <- function(location = "tr",
                            which_north = "true") {
  ggspatial::annotation_north_arrow(
    location    = location,
    which_north = which_north,
    height      = grid::unit(0.80, "cm"),
    width       = grid::unit(0.80, "cm"),
    pad_x       = grid::unit(0.30, "cm"),
    pad_y       = grid::unit(0.30, "cm"),
    style = ggspatial::north_arrow_fancy_orienteering(
      line_width = 0.7,
      line_col   = "grey20",
      fill       = c("white", "grey20"),
      text_col   = GREY_TEXT,
      text_family = PUB_FAMILY,
      text_size  = PUB_ANNOT
    )
  )
}

# Projected CRS for the Philippines.
# EPSG:3123 (PRS92 / Philippines Zone III) is a transverse-Mercator system
# covering the main archipelago. Using a projected CRS -- rather than plotting
# unprojected lon/lat degrees -- means the scale bar represents a true, constant
# ground distance and areas are not latitude-distorted.
CRS_PH_PROJECTED <- 3123
CRS_WGS84        <- 4326

#' Attach CRS handling, graticule styling and cartographic furniture to a map.
#'
#' @param p        A ggplot with geom_sf layers.
#' @param crs      Target CRS (EPSG code) for coord_sf.
#' @param graticule Draw a light lon/lat graticule.
pub_map_frame <- function(p,
                          crs = CRS_PH_PROJECTED,
                          graticule = TRUE,
                          scalebar_loc = "bl",
                          arrow_loc = "tr") {
  p <- p +
    pub_scalebar(location = scalebar_loc) +
    pub_north_arrow(location = arrow_loc) +
    ggplot2::coord_sf(crs = sf::st_crs(crs), expand = TRUE, clip = "off")

  # Graticule suppressed: the revision removes gridlines from every plot, and a
  # lon/lat graticule is a gridline. `graticule` is kept as an argument for
  # call-site compatibility.
  p <- p + ggplot2::theme(panel.grid = ggplot2::element_blank())
  p
}


# -----------------------------------------------------------------------------
# 9. SCALE SHORTCUTS
# -----------------------------------------------------------------------------
scale_fill_detector <- function(name = "Dominant detector", ...) {
  ggplot2::scale_fill_manual(name = name, values = PAL_DETECTOR,
                             na.value = "grey85", drop = TRUE, ...)
}
scale_colour_detector <- function(name = "Detector", ...) {
  ggplot2::scale_colour_manual(name = name, values = PAL_DETECTOR,
                               na.value = "grey85", drop = TRUE, ...)
}
scale_color_detector <- scale_colour_detector

message("[theme] Publication theme loaded  |  base = ", PUB_BASE,
        " pt  |  TYPE_SCALE = ", TYPE_SCALE,
        "  |  family = ", PUB_FAMILY)


# -----------------------------------------------------------------------------
# 10. BOLD / GRIDLESS VARIANT
# -----------------------------------------------------------------------------
# Requested house style for the main manuscript figures: no gridlines, and all
# axis furniture -- titles, tick labels and the axis lines themselves -- in bold
# black rather than grey.
#
# Applied as a MODIFIER on top of an existing theme so a plot keeps its own
# legend placement, facet strips and margins; only the axis styling and grid are
# overridden. Use as:  p + theme_bold_axes()
theme_bold_axes <- function(base_size = PUB_BASE, axis_line_width = 0.6) {
  ggplot2::theme(
    # --- gridlines removed entirely -----------------------------------------
    panel.grid       = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),

    # --- axis furniture: bold black -----------------------------------------
    axis.title   = ggplot2::element_text(size = PUB_AXIS_TIT, face = "bold",
                                         colour = "black"),
    axis.title.x = ggplot2::element_text(size = PUB_AXIS_TIT, face = "bold",
                                         colour = "black",
                                         margin = ggplot2::margin(t = 5)),
    axis.title.y = ggplot2::element_text(size = PUB_AXIS_TIT, face = "bold",
                                         colour = "black", angle = 90,
                                         margin = ggplot2::margin(r = 5)),
    axis.text    = ggplot2::element_text(size = PUB_AXIS_TXT, face = "bold",
                                         colour = "black"),
    axis.text.x  = ggplot2::element_text(size = PUB_AXIS_TXT, face = "bold",
                                         colour = "black"),
    axis.text.y  = ggplot2::element_text(size = PUB_AXIS_TXT, face = "bold",
                                         colour = "black"),
    # Secondary axes must be named explicitly. ggplot2's add_theme() merges
    # element properties -- NULL fields of a later element inherit from the
    # earlier one -- so a downstream theme(axis.text.y = element_text(size = 6.8))
    # keeps face = "bold". But that inheritance only happens for elements this
    # theme actually SETS. axis.text.x.top and axis.text.y.right are used by the
    # dominance matrix and the secondary R(t) axis; without these two lines they
    # would silently stay non-bold.
    axis.text.x.top   = ggplot2::element_text(size = PUB_AXIS_TXT, face = "bold",
                                              colour = "black"),
    axis.text.y.right = ggplot2::element_text(size = PUB_AXIS_TXT, face = "bold",
                                              colour = "black"),
    axis.title.x.top   = ggplot2::element_text(size = PUB_AXIS_TIT, face = "bold",
                                               colour = "black"),
    axis.title.y.right = ggplot2::element_text(size = PUB_AXIS_TIT, face = "bold",
                                               colour = "black", angle = 90),
    axis.line    = ggplot2::element_line(colour = "black",
                                         linewidth = axis_line_width),
    axis.ticks   = ggplot2::element_line(colour = "black",
                                         linewidth = axis_line_width),

    # A panel border would duplicate the now-visible axis lines.
    panel.border     = ggplot2::element_blank(),
    panel.background = ggplot2::element_rect(fill = "white", colour = NA)
  )
}

#' Same, plus removal of the legend box/frame (Figure 6 request).
theme_bold_axes_nolegendbox <- function(...) {
  theme_bold_axes(...) +
    ggplot2::theme(
      legend.background = ggplot2::element_blank(),
      legend.box.background = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank()
    )
}


# -----------------------------------------------------------------------------
# 11. SIGNIFICANCE BRACKETS FOR CATEGORY-vs-CATEGORY PANELS
# -----------------------------------------------------------------------------
# Builds the bracket geometry for "focal vs each comparator" comparisons on a
# discrete-x panel, using a paired Wilcoxon signed-rank test paired by unit.
#
# ONLY SIGNIFICANT COMPARISONS ARE RETURNED. A non-significant pair gets no
# bracket and no p-value: the connector is what asserts a comparison, so drawing
# it for a null result invites over-reading. Brackets stack by position among
# those actually drawn, so an omitted pair leaves no vertical gap.
#
# @param d      long data: one row per unit per level, with `unit_col`,
#               `level_col` (a factor giving x order) and `value_col`.
# @param focal  the level compared against every other.
# @return data.frame(x, xend, y, label) -- possibly zero rows.
#' Conventional significance stars for a p-value.
#' *** p < 0.001, ** p < 0.01, * p < 0.05, "" otherwise.
sig_stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  ""
}

#' p-value with its stars, e.g. "p = 0.032 *" or "p < 0.001 ***".
p_label_starred <- function(p) {
  if (is.na(p)) return("")
  base <- if (p < 0.001) "p < 0.001" else paste0("p = ", sprintf("%.3f", p))
  paste0(base, " ", sig_stars(p))
}

#' Paired Wilcoxon of `focal` against every other level, paired by unit.
#'
#' Returns EVERY comparison, significant or not, with its full statistics. This
#' is the table that must be exported alongside any figure drawing significance
#' brackets, so the p-values on the panel can be traced to a file.
dominance_pairwise_tests <- function(d, unit_col, level_col, value_col, focal,
                                     alpha = 0.05, boot_B = 200L) {
  lv <- levels(d[[level_col]])
  if (is.null(lv)) lv <- sort(unique(as.character(d[[level_col]])))
  if (!(focal %in% lv)) return(NULL)

  fa <- d[as.character(d[[level_col]]) == focal, ]
  rows <- list()
  for (cm in setdiff(lv, focal)) {
    fb <- d[as.character(d[[level_col]]) == cm, ]
    common <- intersect(fa[[unit_col]], fb[[unit_col]])
    if (length(common) == 0L) next
    a <- fa[[value_col]][match(common, fa[[unit_col]])]
    b <- fb[[value_col]][match(common, fb[[unit_col]])]
    w <- hh_paired_wilcoxon(a, b, boot_B = boot_B, alpha = alpha)
    rows[[length(rows) + 1L]] <- cbind(
      data.frame(Comparison   = paste0(focal, " vs ", cm),
                 Detector_A   = focal,
                 Detector_B   = cm,
                 Pairing_Unit = unit_col,
                 Metric       = value_col,
                 Mean_A       = mean(a, na.rm = TRUE),
                 Mean_B       = mean(b, na.rm = TRUE),
                 stringsAsFactors = FALSE), w)
  }
  if (length(rows) == 0L) return(NULL)
  out <- dplyr::bind_rows(rows)
  out$Drawn_On_Figure <- !is.na(out$p_value) & out$Significant_005
  out
}

#' Bracket geometry for the significant comparisons only.
#'
#' @param tests Optional pre-computed dominance_pairwise_tests() table; supplied
#'   so the figure and the exported CSV are guaranteed to come from one
#'   computation rather than two independent ones that could drift apart.
build_sig_brackets <- function(d, unit_col, level_col, value_col, focal,
                               alpha = 0.05, step_frac = 0.13,
                               base_frac = 0.10, tests = NULL) {
  lv <- levels(d[[level_col]])
  if (is.null(lv)) lv <- sort(unique(as.character(d[[level_col]])))
  if (!(focal %in% lv)) return(NULL)

  if (is.null(tests)) {
    tests <- dominance_pairwise_tests(d, unit_col, level_col, value_col,
                                      focal, alpha = alpha)
  }
  if (is.null(tests) || nrow(tests) == 0L) return(NULL)

  vals <- d[[value_col]]
  vmax <- suppressWarnings(max(vals, na.rm = TRUE))
  vmin <- suppressWarnings(min(vals, na.rm = TRUE))
  if (!is.finite(vmax) || !is.finite(vmin)) return(NULL)
  vrange <- vmax - vmin
  if (vrange <= 0) vrange <- max(abs(vmax), 1)

  sig <- tests[isTRUE_vec(tests$Drawn_On_Figure), , drop = FALSE]
  if (nrow(sig) == 0L) return(NULL)

  out <- list()
  for (k in seq_len(nrow(sig))) {
    pv  <- sig$p_value[k]
    lab <- p_label_starred(pv)
    out[[k]] <- data.frame(
      x     = which(lv == focal),
      xend  = which(lv == sig$Detector_B[k]),
      y     = vmax + base_frac * vrange + (k - 1L) * step_frac * vrange,
      label = lab, stringsAsFactors = FALSE)
  }
  dplyr::bind_rows(out)
}

#' TRUE-only logical helper (NA-safe), used to select drawn comparisons.
isTRUE_vec <- function(x) !is.na(x) & x

#' Add bracket layers produced by build_sig_brackets() to a plot.
add_sig_brackets <- function(p, brack_df, vrange, family = PUB_FAMILY) {
  if (is.null(brack_df) || nrow(brack_df) == 0L) return(p)
  tick <- 0.030 * vrange
  p +
    ggplot2::geom_segment(data = brack_df,
      ggplot2::aes(x = x, xend = xend, y = y, yend = y),
      colour = "black", linewidth = 0.45, inherit.aes = FALSE) +
    ggplot2::geom_segment(data = brack_df,
      ggplot2::aes(x = x, xend = x, y = y, yend = y - tick),
      colour = "black", linewidth = 0.45, inherit.aes = FALSE) +
    ggplot2::geom_segment(data = brack_df,
      ggplot2::aes(x = xend, xend = xend, y = y, yend = y - tick),
      colour = "black", linewidth = 0.45, inherit.aes = FALSE) +
    ggplot2::geom_text(data = brack_df,
      ggplot2::aes(x = (x + xend) / 2, y = y, label = label),
      vjust = -0.45, colour = "black", fontface = "bold",
      size = pub_text_size(PUB_ANNOT - 0.5), family = family,
      inherit.aes = FALSE)
}
