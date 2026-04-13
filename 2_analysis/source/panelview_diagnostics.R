#===============================================================================
# panelview_diagnostics.R
#
# Panel data diagnostic visualizations for the Czech construction procurement
# panel using panelView (Mou, Liu & Xu), the companion visualization package
# to `fect` (Liu, Wang & Xu 2022).
#
# Produces four diagnostic figures that motivate the choice of panel data
# estimators (FEct, IFEct, MC, TROP) in the robustness stack:
#
#   1. Treatment status — staggered procurement entry pattern
#   2. Missing data    — unbalanced panel structure (~63% density)
#   3. Outcome dynamics — log markup evolution by treatment status
#   4. Cohort view     — log markup by entry-year cohort
#
# These visualizations surface identification-relevant patterns (selection into
# treatment cohorts, attrition, outcome trends) that should be described before
# applying weighted / factor-model estimators.
#
# Usage:
#   Rscript panelview_diagnostics.R
#
# Inputs:
#   ../input/data.dta                 (Czech firm-year panel)
#   ../output/data/paper_markups.dta  (wide-format markups)
#
# Outputs:
#   ../output/figures/panelview_treatment.pdf
#   ../output/figures/panelview_missing.pdf
#   ../output/figures/panelview_outcome.pdf
#   ../output/figures/panelview_cohort.pdf
#===============================================================================

suppressPackageStartupMessages({
  library(panelView)
  library(haven)
  library(ggplot2)
  library(dplyr)
})

## Source the shared theme helper (sets theme_minimal(base_size=11) as default
## and provides ggsave_markups with dpi=300 + bg=white). Note: panelView() uses
## theme_bw() internally, so the helper's theme_set only affects non-panelView
## ggsaves; the main benefit here is standardized ggsave parameters.
source("theme_markups.R")

input_dir  <- file.path("..", "input")
output_dir <- file.path("..", "output")
fig_dir    <- file.path(output_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(42)

# ------------------------------------------------------------------------------
# Data loading (mirrors fect_estimation.R and trop_estimation.R)
# ------------------------------------------------------------------------------
cat("[panelView] Loading data...\n")
df <- read_dta(file.path(input_dir, "data.dta"))
mk <- read_dta(file.path(output_dir, "data", "paper_markups.dta"))

if (inherits(mk$year, c("Date", "POSIXct", "POSIXt"))) {
  mk$year <- as.integer(format(mk$year, "%Y"))
} else {
  mk$year <- as.integer(mk$year)
}
df$year <- as.integer(df$year)

panel <- merge(
  df[, c("id", "year", "pp_dummy", "k", "cogs")],
  mk[, c("id", "year", "markup_A")],
  by = c("id", "year")
)
panel <- panel[!is.na(panel$markup_A) & panel$markup_A > 0, ]
panel$log_mu <- log(panel$markup_A)
panel$id     <- as.integer(panel$id)
panel$year   <- as.integer(panel$year)
panel$D      <- as.integer(panel$pp_dummy)

cat(sprintf("[panelView] Panel: %d obs, %d firms, %d years\n",
            nrow(panel),
            length(unique(panel$id)),
            length(unique(panel$year))))

# ------------------------------------------------------------------------------
# Sample for the firm-level figures (plotting 1298 firms is unreadable)
# ------------------------------------------------------------------------------
# Stratified: 50 ever-treated + 50 never-treated
ever_treated  <- unique(panel$id[panel$D == 1])
never_treated <- setdiff(unique(panel$id), ever_treated)
sample_ids <- c(
  sample(ever_treated,  min(50, length(ever_treated))),
  sample(never_treated, min(50, length(never_treated)))
)
panel_sample <- panel[panel$id %in% sample_ids, ]
cat(sprintf("[panelView] Visualization sample: %d firms (%d treated + %d control)\n",
            length(sample_ids),
            sum(unique(panel_sample$id) %in% ever_treated),
            sum(unique(panel_sample$id) %in% never_treated)))

# ------------------------------------------------------------------------------
# Helper: save any returned plot object to PDF
# panelView returns ggplot objects (it's built on ggplot2), but to be safe
# we handle both direct-return and side-effect-print styles.
# ------------------------------------------------------------------------------
save_plot <- function(p, path, width, height) {
  result <- tryCatch({
    if (inherits(p, "ggplot")) {
      ggsave_markups(path, p, width = width, height = height)
    } else {
      pdf(path, width = width, height = height, bg = "white")
      print(p)
      dev.off()
    }
    TRUE
  }, error = function(e) {
    cat("  [save error]", conditionMessage(e), "\n")
    FALSE
  })
  return(result)
}

# ------------------------------------------------------------------------------
# Figure 1: Treatment status — staggered adoption pattern
# ------------------------------------------------------------------------------
cat("[panelView] Figure 1/4: treatment status (by.timing)...\n")
p1 <- tryCatch(
  panelview(
    data    = panel_sample,
    Y       = "log_mu",
    D       = "D",
    index   = c("id", "year"),
    type    = "treat",
    by.timing = TRUE,
    main    = "Treatment Status: Procurement Entry Pattern",
    xlab    = "Year",
    ylab    = "Firms (sorted by treatment timing)",
    theme.bw = TRUE,
    legend.labs = c("Not in procurement", "In procurement")
  ),
  error = function(e) { cat("  [panelView error]", conditionMessage(e), "\n"); NULL }
)
if (!is.null(p1)) save_plot(p1, file.path(fig_dir, "panelview_treatment.pdf"), 9, 6)

# ------------------------------------------------------------------------------
# Figure 2: Missing data pattern — unbalanced panel visualization
# ------------------------------------------------------------------------------
cat("[panelView] Figure 2/4: missing data pattern...\n")
p2 <- tryCatch(
  panelview(
    data    = panel_sample,
    Y       = "log_mu",
    D       = "D",
    index   = c("id", "year"),
    type    = "missing",
    main    = "Missing Data Pattern (~63% density)",
    xlab    = "Year",
    ylab    = "Firms",
    theme.bw = TRUE,
    by.timing = TRUE
  ),
  error = function(e) { cat("  [panelView error]", conditionMessage(e), "\n"); NULL }
)
if (!is.null(p2)) save_plot(p2, file.path(fig_dir, "panelview_missing.pdf"), 9, 6)

# ------------------------------------------------------------------------------
# Figure 3: Outcome dynamics — log markup by treatment status
# ------------------------------------------------------------------------------
cat("[panelView] Figure 3/4: outcome dynamics (by.group)...\n")
p3 <- tryCatch(
  panelview(
    data    = panel_sample,
    Y       = "log_mu",
    D       = "D",
    index   = c("id", "year"),
    type    = "outcome",
    by.group = TRUE,
    main    = "Log Markup Dynamics by Treatment Status",
    xlab    = "Year",
    ylab    = "log(markup)",
    theme.bw = TRUE
  ),
  error = function(e) { cat("  [panelView error]", conditionMessage(e), "\n"); NULL }
)
if (!is.null(p3)) save_plot(p3, file.path(fig_dir, "panelview_outcome.pdf"), 9, 6)

# ------------------------------------------------------------------------------
# Figure 4: Cohort-stratified view — define absorbing treatment as first-entry
# panelView's by.cohort requires staggered adoption; Czech procurement has
# reversals, so we construct D_absorbing = 1 from the first treatment year
# onward. This shows what the cohort structure would look like if we treated
# procurement as absorbing — useful for comparing with BJS / SDID approaches
# that also treat entry as absorbing.
# ------------------------------------------------------------------------------
cat("[panelView] Figure 4/4: cohort view (absorbing-treatment variant)...\n")
panel_absorbing <- panel %>%
  group_by(id) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(
    first_treat_year = suppressWarnings(min(year[D == 1])),
    D_absorbing = as.integer(year >= first_treat_year & is.finite(first_treat_year))
  ) %>%
  ungroup() %>%
  as.data.frame()

p4 <- tryCatch(
  panelview(
    data    = panel_absorbing,
    Y       = "log_mu",
    D       = "D_absorbing",
    index   = c("id", "year"),
    type    = "outcome",
    by.cohort = TRUE,
    main    = "Log Markup by Entry Cohort (Absorbing Treatment)",
    xlab    = "Year",
    ylab    = "log(markup)",
    theme.bw = TRUE
  ),
  error = function(e) { cat("  [panelView error]", conditionMessage(e), "\n"); NULL }
)
if (!is.null(p4)) save_plot(p4, file.path(fig_dir, "panelview_cohort.pdf"), 10, 6)

# ------------------------------------------------------------------------------
# Report
# ------------------------------------------------------------------------------
cat("\n", strrep("=", 60), "\n", sep = "")
cat("panelView diagnostics complete\n")
cat(strrep("=", 60), "\n", sep = "")
cat("Output directory:", fig_dir, "\n")
for (f in c("panelview_treatment.pdf", "panelview_missing.pdf",
            "panelview_outcome.pdf", "panelview_cohort.pdf")) {
  full_path <- file.path(fig_dir, f)
  if (file.exists(full_path)) {
    size_kb <- round(file.info(full_path)$size / 1024, 1)
    cat(sprintf("  %s  (%.1f KB)\n", f, size_kb))
  } else {
    cat(sprintf("  %s  MISSING\n", f))
  }
}
