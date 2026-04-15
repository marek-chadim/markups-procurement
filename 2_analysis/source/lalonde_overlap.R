## lalonde_overlap.R
## Ports lalonde2_trim.R + lalonde3_overlap.R from Imbens & Xu (2025, JEP 39(4))
## to the Czech construction panel. Runs the 11-method CIA batch on both the
## full sample and the Crump-Hotz-Imbens-Mitnik (2009) symmetric [alpha, 1-alpha]
## trimmed sample, writes a paired Full/Trimmed LaTeX table, and exports a
## propensity-score density diagnostic.

suppressPackageStartupMessages({
  library(ggplot2)
})

source("lib/theme_markups.R")
source("lib/lalonde_functions_helpers.R")

input_dir    <- file.path("..", "input")
output_dir   <- file.path("..", "output")
lalonde_code <- "/Users/marek/Desktop/io/project/thesis/observables/code"

source(file.path(lalonde_code, "functions_est.R"))
diff_est <- diff

dir.create(file.path(output_dir, "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "tables"),  showWarnings = FALSE, recursive = TRUE)

## ---- Load data & compute propensity ----
panel   <- load_lalonde_panel(input_dir, output_dir)
ps      <- compute_propensity(panel)
trimmed <- trim_propensity(panel, ps, alpha = 0.1)

cat(sprintf("Full:    N=%d, treated=%d, control=%d\n",
            nrow(panel),   sum(panel$D),   sum(1 - panel$D)))
cat(sprintf("Trimmed: N=%d, treated=%d, control=%d (alpha=0.1)\n",
            nrow(trimmed), sum(trimmed$D), sum(1 - trimmed$D)))

## ---- Overlap diagnostic figures ----
plot_df <- data.frame(
  ps    = ps,
  treat = factor(panel$D, levels = 0:1, labels = c("Control", "Treated"))
)
p <- ggplot(plot_df, aes(x = ps, fill = treat)) +
  geom_density(alpha = 0.45, color = NA) +
  geom_vline(xintercept = c(0.1, 0.9), linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c("Control" = markups_blue, "Treated" = markups_pink)) +
  labs(x = "GRF propensity score (probability_forest, num.trees = 4000)",
       y = "Density", fill = NULL,
       title = "Propensity score overlap, Czech construction panel",
       subtitle = sprintf("Full N=%d; dashed lines mark alpha=0.1 symmetric trim",
                          nrow(panel)))
ggsave_markups(file.path(output_dir, "figures", "lalonde_ps_overlap.pdf"), p)
cat("Saved: figures/lalonde_ps_overlap.pdf\n")

## Imbens-Wooldridge (2007) Figure 3-6 style: 2x2 grid of histograms,
## treated vs control in columns, full vs trimmed in rows. Separate histograms
## (not overlapping densities) make the mass-at-zero overlap problem visible
## in the style of the NBER lecture notes.
ps_trim <- attr(trimmed, "ps_kept")
hist_df <- rbind(
  data.frame(
    ps     = ps,
    treat  = ifelse(panel$D == 1, "Treated", "Control"),
    sample = "Full sample",
    stringsAsFactors = FALSE
  ),
  data.frame(
    ps     = ps_trim,
    treat  = ifelse(trimmed$D == 1, "Treated", "Control"),
    sample = "Trimmed ($[0.1, 0.9]$)",
    stringsAsFactors = FALSE
  )
)
hist_df$treat  <- factor(hist_df$treat,  levels = c("Control", "Treated"))
hist_df$sample <- factor(hist_df$sample, levels = c("Full sample", "Trimmed ($[0.1, 0.9]$)"))

p_hist <- ggplot(hist_df, aes(x = ps)) +
  geom_histogram(binwidth = 0.025, boundary = 0, fill = markups_blue, color = "white") +
  geom_vline(xintercept = c(0.1, 0.9), linetype = "dashed", color = "grey40") +
  facet_grid(sample ~ treat, scales = "free_y") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(x = "GRF propensity score", y = "Count",
       title = "Propensity-score histograms before and after trim",
       subtitle = "Imbens-Wooldridge (2007) Figure 3-6 style")
ggsave_markups(file.path(output_dir, "figures", "lalonde_ps_histograms.pdf"),
               p_hist, width = 8, height = 6)
cat("Saved: figures/lalonde_ps_histograms.pdf\n")

## ---- Run 11-method batch on both samples ----
methods <- c("diff", "reg", "om.reg", "om.grf", "matching", "psm",
             "ipw", "cbps", "ebal", "dml", "aipw_grf")
covariates <- LALONDE_COVARIATES

run_one <- function(name, fn, data, needs_covar = TRUE) {
  cat(sprintf("  %-10s ", name))
  out <- tryCatch({
    r <- if (needs_covar) fn(data, "log_mu", "D", covariates)
         else             fn(data, "log_mu", "D")
    as.numeric(r)[1:4]
  }, error = function(e) {
    cat("FAILED: ", conditionMessage(e), "\n")
    return(rep(NA_real_, 4))
  })
  if (!any(is.na(out))) cat(sprintf("ATT=%.4f SE=%.4f\n", out[1], out[2]))
  out
}

run_batch <- function(data, label) {
  cat(sprintf("\n[%s] Running %d methods on N=%d\n", label, length(methods), nrow(data)))
  res <- as.data.frame(matrix(NA_real_, length(methods), 4))
  rownames(res) <- methods
  colnames(res) <- c("Estimate", "SE", "CI_lower", "CI_upper")
  set.seed(1234)
  res["diff", ]     <- run_one("diff",     diff_est, data, needs_covar = FALSE)
  res["reg", ]      <- run_one("reg",      reg,      data)
  res["om.reg", ]   <- run_one("om.reg",   om.reg,   data)
  res["om.grf", ]   <- run_one("om.grf",   om.grf,   data)
  res["matching", ] <- run_one("matching", matching, data)
  res["psm", ]      <- run_one("psm",      psm,      data)
  res["ipw", ]      <- run_one("ipw",      ipw,      data)
  res["cbps", ]     <- run_one("cbps",     cbps,     data)
  res["ebal", ]     <- run_one("ebal",     function(d,Y,Tr,X) quiet(ebal(d,Y,Tr,X)), data)
  res["dml", ]      <- run_one("dml",      dml,      data)
  res["aipw_grf", ] <- run_one("aipw_grf", aipw,     data)
  res
}

full_res    <- run_batch(panel,   "Full")
trimmed_res <- run_batch(trimmed, "Trimmed")

## ---- Combine and write CSV ----
combined <- cbind(
  full_est    = full_res$Estimate,
  full_se     = full_res$SE,
  trim_est    = trimmed_res$Estimate,
  trim_se     = trimmed_res$SE,
  delta       = trimmed_res$Estimate - full_res$Estimate
)
rownames(combined) <- methods
write.csv(combined, file.path(output_dir, "lalonde_overlap_results.csv"))
cat("\nSaved: lalonde_overlap_results.csv\n")

## ---- Paired LaTeX table ----
label_map <- c(
  diff      = "Difference in Means",
  reg       = "Regression Adjustment",
  om.reg    = "Outcome Model (OLS)",
  om.grf    = "Outcome Model (GRF)",
  matching  = "Nearest-Neighbor Matching",
  psm       = "Propensity Score Matching",
  ipw       = "IPW",
  cbps      = "CBPS",
  ebal      = "Entropy Balancing",
  dml       = "Double/Debiased ML",
  aipw_grf  = "AIPW-GRF"
)

tex <- c(
  "\\begin{table}[htbp]\\centering",
  "\\caption{Overlap Trimming: Full vs. Crump et al.\\ (2009) Trimmed Sample}\\label{tab:lalonde_overlap}",
  "\\begin{threeparttable}",
  "\\begin{tabular}{lcccc}",
  "\\toprule",
  "& \\multicolumn{2}{c}{Full} & \\multicolumn{2}{c}{Trimmed} \\\\",
  "\\cmidrule(lr){2-3} \\cmidrule(lr){4-5}",
  "Estimator & ATT & SE & ATT & SE \\\\",
  "\\midrule"
)
for (m in methods) {
  lab <- label_map[[m]]
  fe  <- full_res[m, "Estimate"];    fs <- full_res[m, "SE"]
  te  <- trimmed_res[m, "Estimate"]; ts <- trimmed_res[m, "SE"]
  tex <- c(tex, sprintf(
    "%s & %.4f & (%.4f) & %.4f & (%.4f) \\\\",
    lab, fe, fs, te, ts))
}
tex <- c(tex,
  "\\midrule",
  sprintf("$N$ & \\multicolumn{2}{c}{%d} & \\multicolumn{2}{c}{%d} \\\\",
          nrow(panel), nrow(trimmed)),
  sprintf("Treated & \\multicolumn{2}{c}{%d} & \\multicolumn{2}{c}{%d} \\\\",
          sum(panel$D), sum(trimmed$D)),
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{tablenotes}\\footnotesize",
  "\\item \\textit{Notes:} Eleven selection-on-observables estimators from Imbens and Xu (2025, \\textit{JEP} 39(4)) applied to the full Czech construction panel and to a symmetric propensity-trimmed subsample. Trimming follows Crump, Hotz, Imbens, and Mitnik (2009): keep observations with $\\hat{\\pi}(X_{it}) \\in [0.1, 0.9]$, where $\\hat{\\pi}$ is a generalized random forest propensity score (4{,}000 trees, seed 1234). Outcome: $\\log \\mu^A_{it}$. Treatment: $pp_{it}$. Covariates as in Table~\\ref{tab:lalonde}.",
  "\\end{tablenotes}",
  "\\end{threeparttable}",
  "\\end{table}"
)
writeLines(tex, file.path(output_dir, "tables", "lalonde_overlap.tex"))
cat("Saved: tables/lalonde_overlap.tex\n")
