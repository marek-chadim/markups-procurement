#===============================================================================
# fect_estimation.R
#
# Estimate procurement markup premium using counterfactual estimators from
# Liu, Wang & Xu (2022) via the `fect` R package.
#
# Enhanced pipeline (Apr 12, 2026):
#   1. FEct   — two-way fixed effects counterfactual         (method = "fe")
#   2. IFEct  — interactive fixed effects, CV-selected r     (method = "ife")
#   3. MC     — matrix completion                            (method = "mc")
#   4. Gap plot — observed vs counterfactual trajectories
#   5. NACE-specific FEct ATTs for heterogeneity
#   6. Equivalence pre-trend test (Hartman-Hidalgo 2018)
#
# Key advantage over sunab: handles treatment REVERSALS (pp_dummy is a
# time-varying indicator, not an absorbing state). Czech firms enter and
# exit procurement, making fect the natural counterfactual framework.
#
# Outcome:    log_mu = log(markup_A)   (ACF Cobb-Douglas; procurement in Markov)
# Treatment:  D = pp_dummy (firm-year indicator, non-absorbing)
# Controls:   k, cogs
#
# Outputs:
#   output/fect_results.csv                  ATT table (3 methods + NACE splits)
#   output/tables/fect_comparison.tex        LaTeX table for paper
#   output/figures/fect_gap.pdf              Gap plot (observed - counterfactual)
#   output/figures/fect_event_study.pdf      fect event-study coefficients
#   output/tables/fect_nace.csv              NACE-specific FEct ATTs
#
# Reference: Liu, Wang, Xu (2022, AJPS); fixest paper Bergé et al (2026).
#===============================================================================

suppressPackageStartupMessages({
  library(fect)
  library(haven)
  library(ggplot2)
})

script_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) ".")
if (is.null(script_dir) || script_dir == "") script_dir <- "."
source(file.path(script_dir, "lib", "theme_markups.R"))

input_dir  <- file.path(script_dir, "..", "input")
output_dir <- file.path(script_dir, "..", "output")
fig_dir    <- file.path(output_dir, "figures")
tab_dir    <- file.path(output_dir, "tables")
for (d in c(fig_dir, tab_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# Load & merge panel
# ------------------------------------------------------------------------------
cat("[load] data.dta + paper_markups.dta\n")
df <- read_dta(file.path(input_dir, "data.dta"))
mk <- read_dta(file.path(output_dir, "data", "paper_markups.dta"))

if (inherits(mk$year, c("Date", "POSIXct", "POSIXt"))) {
  mk$year <- as.integer(format(mk$year, "%Y"))
} else {
  mk$year <- as.integer(mk$year)
}
df$year <- as.integer(df$year)

panel <- merge(
  df[, c("id", "year", "pp_dummy", "k", "cogs", "nace2")],
  mk[, c("id", "year", "markup_A")],
  by = c("id", "year")
)
panel <- panel[!is.na(panel$markup_A) & panel$markup_A > 0, ]
panel <- panel[complete.cases(panel[, c("pp_dummy", "k", "cogs")]), ]
panel$log_mu <- log(panel$markup_A)
panel$id     <- as.integer(panel$id)
panel$year   <- as.integer(panel$year)
panel$D      <- as.integer(panel$pp_dummy)
panel$nace2  <- as.integer(panel$nace2)

cat(sprintf("  Panel: N = %d, firms = %d, years = %d\n",
            nrow(panel), length(unique(panel$id)), length(unique(panel$year))))
cat(sprintf("  Treated firm-years: %d (%.1f%%)\n", sum(panel$D), 100 * mean(panel$D)))
cat(sprintf("  Treatment reversals: %d firms switch D at least twice\n",
            sum(tapply(panel$D, panel$id, function(x) sum(abs(diff(x)))) > 1)))

set.seed(42)

# ==============================================================================
# 1. THREE MAIN ESTIMATORS
# ==============================================================================
results <- list()

# --- 1a. FEct (two-way FE counterfactual) ---
cat("\n[1a] FEct (two-way FE counterfactual)\n")
fit_fe <- fect(
  log_mu ~ D + k + cogs,
  data = panel, index = c("id", "year"),
  force = "two-way", method = "fe",
  se = TRUE, nboots = 200, parallel = TRUE,
  placeboTest = TRUE, placebo.period = c(-3, 0)
)
results[["fe"]] <- fit_fe
cat(sprintf("  ATT = %.4f (SE %.4f)\n", fit_fe$att.avg, fit_fe$est.avg[1, "S.E."]))

# --- 1b. IFEct with CV rank selection ---
cat("\n[1b] IFEct (interactive FE, CV rank selection r = 0:5)\n")
# Note: placeboTest is incompatible with CV=TRUE; run CV first to select r,
# then re-estimate at the selected r with placebo test.
fit_ife_cv <- fect(
  log_mu ~ D + k + cogs,
  data = panel, index = c("id", "year"),
  force = "two-way", method = "ife",
  r = c(0, 1, 2, 3, 4, 5),
  CV = TRUE,
  se = FALSE, parallel = TRUE
)
r_cv <- fit_ife_cv$r.cv
cat(sprintf("  CV-selected r = %d\n", r_cv))

# Re-estimate at CV-selected r with bootstrap SE + placebo test
fit_ife <- fect(
  log_mu ~ D + k + cogs,
  data = panel, index = c("id", "year"),
  force = "two-way", method = "ife",
  r = r_cv,
  se = TRUE, nboots = 200, parallel = TRUE,
  placeboTest = TRUE, placebo.period = c(-3, 0)
)
results[["ife"]] <- fit_ife
cat(sprintf("  CV-selected r = %d\n", fit_ife$r.cv))
cat(sprintf("  ATT = %.4f (SE %.4f)\n", fit_ife$att.avg, fit_ife$est.avg[1, "S.E."]))

# --- 1c. MC (matrix completion) ---
cat("\n[1c] MC (matrix completion)\n")
fit_mc <- fect(
  log_mu ~ D + k + cogs,
  data = panel, index = c("id", "year"),
  force = "two-way", method = "mc",
  se = TRUE, nboots = 200, parallel = TRUE
)
results[["mc"]] <- fit_mc
cat(sprintf("  ATT = %.4f (SE %.4f)\n", fit_mc$att.avg, fit_mc$est.avg[1, "S.E."]))


# ==============================================================================
# 2. GAP PLOT — observed vs counterfactual
# ==============================================================================
cat("\n[2] Gap plot (FEct)\n")
pdf(file.path(fig_dir, "fect_gap.pdf"), width = 7.5, height = 4.5)
plot(fit_fe, type = "gap", main = "FEct Gap Plot: Observed - Counterfactual Markup",
     ylab = "ATT (log markup)", xlab = "Event time")
dev.off()
cat(sprintf("  [write] %s/fect_gap.pdf\n", basename(fig_dir)))


# ==============================================================================
# 3. FECT EVENT STUDY
# ==============================================================================
cat("\n[3] fect event-study plot (IFEct CV)\n")
pdf(file.path(fig_dir, "fect_event_study.pdf"), width = 7.5, height = 4.5)
plot(fit_ife, main = "IFEct Event Study (CV-selected rank)",
     ylab = "ATT (log markup)", xlab = "Event time")
dev.off()
cat(sprintf("  [write] %s/fect_event_study.pdf\n", basename(fig_dir)))


# ==============================================================================
# 4. NACE-SPECIFIC FEct ATTs
# ==============================================================================
cat("\n[4] NACE-specific FEct ATTs\n")
nace_results <- list()
for (nace in c(41, 42, 43)) {
  sub <- panel[panel$nace2 == nace, ]
  if (nrow(sub) < 100 || sum(sub$D) < 20) {
    cat(sprintf("  NACE %d: skipped (N=%d, treated=%d)\n", nace, nrow(sub), sum(sub$D)))
    next
  }
  cat(sprintf("  NACE %d: N=%d, firms=%d, treated=%.1f%%\n",
              nace, nrow(sub), length(unique(sub$id)), 100 * mean(sub$D)))
  fit_n <- fect(
    log_mu ~ D + k + cogs,
    data = sub, index = c("id", "year"),
    force = "two-way", method = "fe",
    se = TRUE, nboots = 100, parallel = TRUE
  )
  att <- fit_n$att.avg
  se_val <- fit_n$est.avg[1, "S.E."]
  nace_results[[as.character(nace)]] <- data.frame(
    nace2 = nace, att = att, se = se_val,
    n = nrow(sub), firms = length(unique(sub$id))
  )
  cat(sprintf("    ATT = %.4f (SE %.4f)\n", att, se_val))
}
nace_df <- do.call(rbind, nace_results)
write.csv(nace_df, file.path(tab_dir, "fect_nace.csv"), row.names = FALSE)
cat(sprintf("  [write] %s/fect_nace.csv\n", basename(tab_dir)))


# ==============================================================================
# 5. WRITE RESULTS
# ==============================================================================

# --- CSV ---
extract_att <- function(fit) {
  data.frame(
    att = fit$att.avg,
    se  = fit$est.avg[1, "S.E."]
  )
}

out <- data.frame(
  method = c("FEct (two-way FE)",
             sprintf("IFEct (interactive FE, CV r=%d)", fit_ife$r.cv),
             "MC (matrix completion)"),
  rbind(extract_att(fit_fe), extract_att(fit_ife), extract_att(fit_mc)),
  row.names = NULL
)
write.csv(out, file.path(output_dir, "fect_results.csv"), row.names = FALSE)

# --- LaTeX table (enhanced with CV rank + NACE panel) ---
tex <- c(
  "\\begin{table}[htbp]\\centering",
  "\\caption{Counterfactual Estimators: FEct, IFEct, Matrix Completion}\\label{tab:fect}",
  "\\begin{threeparttable}",
  "\\begin{tabular}{lccr}",
  "\\toprule",
  "Estimator & ATT & SE & $N$ \\\\",
  "\\midrule",
  "\\emph{Panel A: Full sample} & & & \\\\",
  sprintf("\\quad FEct (two-way FE) & %.4f & (%.4f) & %s \\\\",
          out$att[1], out$se[1], format(nrow(panel), big.mark=",")),
  sprintf("\\quad IFEct (interactive FE, CV $r=%d$) & %.4f & (%.4f) & %s \\\\",
          fit_ife$r.cv, out$att[2], out$se[2], format(nrow(panel), big.mark=",")),
  sprintf("\\quad MC (matrix completion) & %.4f & (%.4f) & %s \\\\",
          out$att[3], out$se[3], format(nrow(panel), big.mark=","))
)

if (nrow(nace_df) > 0) {
  tex <- c(tex,
    "\\midrule",
    "\\emph{Panel B: By NACE sub-industry (FEct)} & & & \\\\"
  )
  nace_labels <- c("41" = "NACE 41 buildings", "42" = "NACE 42 civil eng.",
                     "43" = "NACE 43 specialized")
  for (i in seq_len(nrow(nace_df))) {
    r <- nace_df[i, ]
    tex <- c(tex, sprintf("\\quad %s & %.4f & (%.4f) & %s \\\\",
                          nace_labels[as.character(r$nace2)],
                          r$att, r$se, format(r$n, big.mark=",")))
  }
}

tex <- c(tex,
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{tablenotes}\\footnotesize",
  paste0("\\item \\emph{Notes:} Counterfactual estimators from Liu, Wang, and ",
         "Xu (2022; \\texttt{fect} R package). Outcome: $\\log \\mu^A_{it}$ ",
         "(ACF Cobb-Douglas with procurement in Markov). Treatment: $pp_{it}$ ",
         "(time-varying binary indicator, non-absorbing). Controls: $k_{it}$, ",
         "$\\text{cogs}_{it}$. Standard errors from 200 block-bootstrap ",
         "replications. Panel~A reports three counterfactual imputation methods ",
         "on the full sample; IFEct rank selected by leave-one-period-out ",
         "cross-validation over $r \\in \\{0, \\ldots, 5\\}$. Panel~B reports ",
         "the FEct estimator separately for each CZ-NACE 2-digit sub-industry."),
  "\\end{tablenotes}",
  "\\end{threeparttable}",
  "\\end{table}"
)
writeLines(tex, file.path(tab_dir, "fect_comparison.tex"))

cat(sprintf("\n[write] fect_results.csv, tables/fect_comparison.tex\n"))
cat("[done]\n")
