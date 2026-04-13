## lalonde_sens.R
## Ports lalonde7_sens.R from Imbens & Xu (2025, JEP 39(4)) to the Czech
## construction panel. Runs Cinelli & Hazlett (2020) sensemakr sensitivity
## analysis on the OLS regression adjustment with three benchmark covariates:
## cogs (variable-input intensity, strongest linear predictor of markup),
## empl_mid (firm size), and mktshare. Produces robustness values and
## benchmark bounds at kd = 1, 2, 3.

suppressPackageStartupMessages({
  library(sensemakr)
})

source("lalonde_functions_helpers.R")

input_dir  <- file.path("..", "input")
output_dir <- file.path("..", "output")

dir.create(file.path(output_dir, "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "tables"),  showWarnings = FALSE, recursive = TRUE)

## ---- Fit linear regression adjustment ----
panel <- load_lalonde_panel(input_dir, output_dir)
cat(sprintf("Panel: N=%d, treated=%d\n", nrow(panel), sum(panel$D)))

fml <- as.formula(paste("log_mu ~ D +", paste(LALONDE_COVARIATES, collapse = " + ")))
mod <- lm(fml, data = panel)
coef_row <- summary(mod)$coefficients["D", ]
cat("OLS regression adjustment:\n")
print(round(coef_row, 5))

## ---- sensemakr with three benchmarks ----
benchmarks <- c("cogs", "empl_mid", "mktshare")
kd_values  <- 1:3

sens_list <- lapply(benchmarks, function(bm) {
  sensemakr(
    model                = mod,
    treatment            = "D",
    benchmark_covariates = bm,
    kd                   = kd_values,
    sensitivity.of       = "t-value"
  )
})
names(sens_list) <- benchmarks

## ---- Diagnostic contour plot using cogs benchmark ----
pdf(file.path(output_dir, "figures", "lalonde_sensemakr_cogs.pdf"),
    width = 7, height = 5.5)
plot(sens_list[["cogs"]])
dev.off()
cat("Saved: figures/lalonde_sensemakr_cogs.pdf\n")

## ---- Extract robustness values and bounds into a table ----
## sensemakr stores results directly on the object: $sensitivity_stats is a
## one-row data frame with estimate/se/t_statistic/r2yd.x/rv_q/rv_qa/f2yd.x;
## $bounds is an ovb_bounds data frame keyed by kd (one row per kd).
sens_stats <- sens_list[["cogs"]]$sensitivity_stats
rv_q  <- sens_stats$rv_q    # confounder strength to drive ATT to zero
rv_qa <- sens_stats$rv_qa   # confounder strength to cross significance

cat(sprintf("Point estimate (reg): %.5f\n", coef_row["Estimate"]))
cat(sprintf("t-value:              %.3f\n", coef_row["t value"]))
cat(sprintf("Robustness value (zero):         %.4f\n", rv_q))
cat(sprintf("Robustness value (significance): %.4f\n", rv_qa))

adjusted_row <- function(sens_obj, kd) {
  bnd <- sens_obj$bounds
  row <- bnd[kd, ]
  c(est = as.numeric(row$adjusted_estimate),
    tv  = as.numeric(row$adjusted_t),
    lb  = as.numeric(row$adjusted_lower_CI),
    ub  = as.numeric(row$adjusted_upper_CI))
}

bound_table <- list()
for (bm in benchmarks) {
  for (kd in kd_values) {
    r <- adjusted_row(sens_list[[bm]], kd)
    bound_table[[length(bound_table) + 1]] <- data.frame(
      benchmark = bm, kd = kd,
      adj_estimate = r["est"], adj_t = r["tv"],
      adj_lower = r["lb"], adj_upper = r["ub"]
    )
  }
}
bound_df <- do.call(rbind, bound_table)
rownames(bound_df) <- NULL
write.csv(bound_df, file.path(output_dir, "lalonde_sensemakr_bounds.csv"),
          row.names = FALSE)
cat("Saved: lalonde_sensemakr_bounds.csv\n")
print(bound_df)

## ---- LaTeX table ----
bm_label <- c(
  cogs     = "Variable input (cogs)",
  empl_mid = "Employment",
  mktshare = "Market share"
)

tex <- c(
  "\\begin{table}[htbp]\\centering",
  "\\caption{Cinelli--Hazlett (2020) Sensitivity Analysis via \\texttt{sensemakr}}\\label{tab:lalonde_sens}",
  "\\begin{threeparttable}",
  "\\begin{tabular}{lccc}",
  "\\toprule",
  "& Point est.\\ & $t$-stat & 95\\% CI \\\\",
  "\\midrule",
  sprintf("Unadjusted OLS (reg) & %.4f & %.2f & [%.4f, %.4f] \\\\",
          coef_row["Estimate"], coef_row["t value"],
          coef_row["Estimate"] - 1.96 * coef_row["Std. Error"],
          coef_row["Estimate"] + 1.96 * coef_row["Std. Error"]),
  "\\midrule",
  "\\multicolumn{4}{l}{\\textit{Robustness values}} \\\\",
  sprintf("~~RV to drive ATT to zero & \\multicolumn{3}{c}{%.4f} \\\\", rv_q),
  sprintf("~~RV to erase significance ($q = 0.05$) & \\multicolumn{3}{c}{%.4f} \\\\", rv_qa),
  "\\midrule",
  "\\multicolumn{4}{l}{\\textit{Bounds under confounder $k_D \\times$ stronger than benchmark}} \\\\"
)
for (bm in benchmarks) {
  tex <- c(tex, sprintf("~~\\textit{%s} & & & \\\\", bm_label[[bm]]))
  for (kd in kd_values) {
    row <- bound_df[bound_df$benchmark == bm & bound_df$kd == kd, ]
    tex <- c(tex, sprintf(
      "~~~~$k_D = %d$ & %.4f & %.2f & [%.4f, %.4f] \\\\",
      kd, row$adj_estimate, row$adj_t, row$adj_lower, row$adj_upper))
  }
}
tex <- c(tex,
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{tablenotes}\\footnotesize",
  "\\item \\textit{Notes:} Cinelli and Hazlett (2020) sensitivity analysis via the \\texttt{sensemakr} R package, applied to the linear regression adjustment $\\log \\mu^A_{it} = \\beta pp_{it} + X_{it}'\\gamma + \\varepsilon_{it}$ on the full Czech panel. The robustness value (RV) is the partial $R^2$ an unobserved confounder would need to share with both the treatment and the outcome to drive the estimated ATT to zero (first row) or to erase statistical significance at the 5\\% level (second row). Benchmark bounds show what the adjusted point estimate, $t$-statistic, and 95\\% confidence interval would become if an unobserved confounder were $k_D$ times as strong as the named benchmark covariate in explaining variation in the treatment. Following Imbens and Xu (2025, \\textit{JEP} 39(4)).",
  "\\end{tablenotes}",
  "\\end{threeparttable}",
  "\\end{table}"
)
writeLines(tex, file.path(output_dir, "tables", "lalonde_sensemakr.tex"))
cat("Saved: tables/lalonde_sensemakr.tex\n")
