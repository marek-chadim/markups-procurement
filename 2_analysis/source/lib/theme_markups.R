## theme_markups.R
## Shared ggplot2 theme and palette helper for markups-procurement figures.
## Source this file at the top of any R script that produces figures to get
## a consistent look across the paper's R/ggplot2 figures.
##
## Reference: Healy (2026), Data Visualization: A Practical Introduction,
## Princeton UP, 2nd ed., Chapter 8 (Refine your Plots).
##
## Provides:
##   - theme_set() with theme_minimal(base_size = 11) + light refinements
##   - Paul Tol's bright palette as a named vector `tol_bright`
##   - scale_color_markups() / scale_fill_markups() convenience scales
##   - ggsave_markups() with opinionated defaults (width, height, dpi, bg)

suppressPackageStartupMessages({
  library(ggplot2)
})

## ----- Paul Tol bright palette (colorblind-safe, up to 7 categories) -----
## Source: https://personal.sron.nl/~pault/  (colorblind-safe collection)
tol_bright <- c(
  blue   = "#4477AA",
  red    = "#EE6677",
  green  = "#228833",
  yellow = "#CCBB44",
  cyan   = "#66CCEE",
  purple = "#AA3377",
  grey   = "#BBBBBB"
)

## Named shortcuts matching the colors already used in lalonde_catt.R and
## lalonde_overlap.R (preserves visual continuity with existing figures).
markups_blue <- tol_bright[["blue"]]     # "#4477AA"
markups_pink <- tol_bright[["purple"]]   # "#AA3377"

## ----- Global theme: apply once per session ------------------------------
## Consistent base size, clean grid, bold facet strips, legend at bottom.
theme_set(
  theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor  = element_blank(),
      plot.title        = element_text(face = "bold"),
      plot.subtitle     = element_text(color = "grey40"),
      plot.caption      = element_text(color = "grey50", size = rel(0.8),
                                       hjust = 0),
      strip.background  = element_rect(fill = "grey92", color = NA),
      strip.text        = element_text(face = "bold"),
      legend.position   = "bottom"
    )
)

## ----- Discrete color / fill scales using Paul Tol bright ----------------
scale_color_markups <- function(...) {
  ggplot2::scale_color_manual(values = unname(tol_bright), ...)
}
scale_fill_markups <- function(...) {
  ggplot2::scale_fill_manual(values = unname(tol_bright), ...)
}

## ----- ggsave defaults (single-column journal figure, 300 DPI) -----------
MARKUPS_FIG_WIDTH  <- 6.5   # inches  (single-column journal width)
MARKUPS_FIG_HEIGHT <- 4.0   # inches  (~golden ratio)
MARKUPS_FIG_DPI    <- 300

ggsave_markups <- function(filename, plot = last_plot(),
                           width  = MARKUPS_FIG_WIDTH,
                           height = MARKUPS_FIG_HEIGHT,
                           dpi    = MARKUPS_FIG_DPI,
                           bg     = "white",
                           ...) {
  ggplot2::ggsave(
    filename = filename,
    plot     = plot,
    width    = width,
    height   = height,
    dpi      = dpi,
    bg       = bg,
    ...
  )
}
