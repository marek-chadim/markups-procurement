#===============================================================================
# fect_estimation.R
#
# Estimate procurement markup premium using counterfactual estimators from
# Liu, Wang & Xu (2022) via the `fect` R package.
#
# Specifications:
#   1. FEct  — two-way fixed effects counterfactual  (method = "fe")
#   2. IFEct — interactive fixed effects, r = 2       (method = "ife")
#   3. MC    — matrix completion                      (method = "mc")
#
# Outcome:    log_mu = log(markup_A)   (ACF Cobb-Douglas; procurement in Markov)
# Treatment:  D = pp_dummy (firm-year indicator)
# Controls:   k, cogs
#
# Inputs:
#   ../input/data.dta                 (Czech panel; firm-year)
#   ../output/data/paper_markups.dta  (wide-format markups; markup_A)
#
# Outputs:
#   ../output/fect_results.csv
#   ../output/tables/fect_comparison.tex
#===============================================================================

suppressPackageStartupMessages({
  library(fect)
  library(haven)
  library(dplyr)
})

input_dir  <- file.path("..", "input")
output_dir <- file.path("..", "output")

# ------------------------------------------------------------------------------
# Load & merge panel
# ------------------------------------------------------------------------------
df <- read_dta(file.path(input_dir, "data.dta"))
mk <- read_dta(file.path(output_dir, "data", "paper_markups.dta"))

# year is numeric in both; coerce defensively (handles Stata datetime if present)
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
panel <- panel[complete.cases(panel[, c("pp_dummy", "k", "cogs")]), ]
panel$log_mu <- log(panel$markup_A)
panel$id     <- as.integer(panel$id)
panel$year   <- as.integer(panel$year)
panel$D      <- as.integer(panel$pp_dummy)

cat("Panel: N =", nrow(panel),
    ", firms =", length(unique(panel$id)),
    ", years =", length(unique(panel$year)), "\n")
cat("Treated firm-years:", sum(panel$D),
    " (", round(100 * mean(panel$D), 1), "%)\n", sep = "")

# ------------------------------------------------------------------------------
# Run three specifications
# ------------------------------------------------------------------------------
set.seed(42)
results <- list()
method_labels <- c(
  fe  = "FEct (two-way FE)",
  ife = "IFEct (interactive FE, r=2)",
  mc  = "MC (matrix completion)"
)

for (method in c("fe", "ife", "mc")) {
  r_val <- if (method == "ife") 2 else 0
  cat("\n=== Running ", method, " ===\n", sep = "")

  # MC uses cross-validation to pick lambda by default, which is incompatible
  # with placeboTest=TRUE. Run MC with CV but without placebo; FE/IFE get placebo.
  do_placebo <- (method != "mc")

  fit <- fect(
    log_mu ~ D + k + cogs,
    data           = panel,
    index          = c("id", "year"),
    force          = "two-way",
    method         = method,
    r              = r_val,
    se             = TRUE,
    nboots         = 100,
    parallel       = TRUE,
    placeboTest    = do_placebo,
    placebo.period = if (do_placebo) c(-2, 0) else NULL
  )

  att_avg <- as.numeric(fit$att.avg)
  est_avg <- fit$est.avg
  if (is.matrix(est_avg) || is.data.frame(est_avg)) {
    se_val <- as.numeric(est_avg[1, "S.E."])
  } else {
    se_val <- as.numeric(est_avg["S.E."])
  }

  results[[method]] <- list(
    att_avg = att_avg,
    se      = se_val,
    method  = method
  )
  cat("  ATT avg: ", round(att_avg, 4),
      "  SE: ", round(se_val, 4), "\n", sep = "")
}

# ------------------------------------------------------------------------------
# Write CSV
# ------------------------------------------------------------------------------
out <- data.frame(
  method = unname(method_labels[c("fe", "ife", "mc")]),
  att    = sapply(results[c("fe", "ife", "mc")], function(x) x$att_avg),
  se     = sapply(results[c("fe", "ife", "mc")], function(x) x$se),
  row.names = NULL,
  stringsAsFactors = FALSE
)
write.csv(out, file.path(output_dir, "fect_results.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# Write LaTeX table
# ------------------------------------------------------------------------------
tex <- c(
  "\\begin{table}[htbp]\\centering",
  "\\caption{Counterfactual Estimators: FEct, IFEct, Matrix Completion}\\label{tab:fect}",
  "\\begin{threeparttable}",
  "\\begin{tabular}{lcc}",
  "\\toprule",
  "Estimator & ATT & SE \\\\",
  "\\midrule",
  sprintf("FEct (two-way FE) & %.4f & (%.4f) \\\\", out$att[1], out$se[1]),
  sprintf("IFEct (interactive FE, $r=2$) & %.4f & (%.4f) \\\\", out$att[2], out$se[2]),
  sprintf("MC (matrix completion) & %.4f & (%.4f) \\\\", out$att[3], out$se[3]),
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{tablenotes}\\footnotesize",
  "\\item \\textit{Notes:} Counterfactual estimators from Liu, Wang, and Xu (2022; \\texttt{fect} R package). Outcome: $\\log \\mu^A_{it}$ (ACF Cobb-Douglas with procurement in Markov). Treatment: $pp_{it}$. Controls: $k_{it}$, $\\text{cogs}_{it}$. Standard errors from 100 bootstrap replications. Two-way fixed effects (firm + year).",
  "\\end{tablenotes}",
  "\\end{threeparttable}",
  "\\end{table}"
)
writeLines(tex, file.path(output_dir, "tables", "fect_comparison.tex"))

cat("\nSaved: fect_results.csv, tables/fect_comparison.tex\n")
