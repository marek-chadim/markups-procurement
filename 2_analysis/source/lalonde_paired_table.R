## lalonde_paired_table.R
## Consolidates the main-outcome (lalonde_results.csv) and pre-treatment
## placebo (lalonde_placebo_results.csv) eleven-method batches into a single
## paired LaTeX table following Imbens and Wooldridge (2007, NBER Summer
## Institute Lecture 1, Table 5): side-by-side placebo vs main outcome. The
## placebo column should hug zero; the main column is the procurement premium.
## The visual contrast between the two columns is the strongest empirical
## evidence that the CIA paradigm recovers the right answer on this panel.
##
## Outputs:
##   output/tables/lalonde_paired.tex

input_dir  <- file.path("..", "input")
output_dir <- file.path("..", "output")

main    <- read.csv(file.path(output_dir, "lalonde_results.csv"),
                    row.names = 1, stringsAsFactors = FALSE)
placebo <- read.csv(file.path(output_dir, "lalonde_placebo_results.csv"),
                    row.names = 1, stringsAsFactors = FALSE)

methods <- c("diff", "reg", "om.reg", "om.grf", "matching", "psm",
             "ipw", "cbps", "ebal", "dml", "aipw_grf")
stopifnot(all(methods %in% rownames(main)))
stopifnot(all(methods %in% rownames(placebo)))

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
  "\\caption{Paired Placebo and Main-Outcome Estimates, Imbens--Wooldridge (2007) Style}\\label{tab:lalonde_paired}",
  "\\begin{threeparttable}",
  "\\begin{tabular}{lcc@{\\hspace{2em}}cc}",
  "\\toprule",
  "& \\multicolumn{2}{c}{Placebo: pre-treatment $\\log \\mu$} & \\multicolumn{2}{c}{Main: $\\log \\mu^A_{it}$} \\\\",
  "\\cmidrule(lr){2-3} \\cmidrule(lr){4-5}",
  "Estimator & ATT & SE & ATT & SE \\\\",
  "\\midrule"
)
for (m in methods) {
  lab <- label_map[[m]]
  pe  <- placebo[m, "Estimate"]; ps <- placebo[m, "SE"]
  me  <- main[m,    "Estimate"]; ms <- main[m,    "SE"]
  tex <- c(tex, sprintf(
    "%s & %.4f & (%.4f) & %.4f & (%.4f) \\\\",
    lab, pe, ps, me, ms))
}
tex <- c(tex,
  "\\midrule",
  "$N$ & \\multicolumn{2}{c}{1{,}025 firms} & \\multicolumn{2}{c}{7{,}666 firm-years} \\\\",
  "Treated & \\multicolumn{2}{c}{483} & \\multicolumn{2}{c}{3{,}443} \\\\",
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{tablenotes}\\footnotesize",
  "\\item \\textit{Notes:} Eleven selection-on-observables estimators from Imbens and Xu (2025, \\textit{JEP} 39(4)) run twice. Placebo column: each firm's earliest panel observation with $pp_{it} = 0$ is treated as a cross-sectional observation, the log markup at that year is the outcome, and an indicator for ``ever later treated'' is the placebo treatment. Under the conditional independence assumption, the placebo ATT should be zero because future procurement cannot causally raise present markups. Main column: full pooled firm-year cross-section with the contemporaneous procurement indicator. The pairing follows Imbens and Wooldridge (2007, NBER Summer Institute Lecture Notes 1, Table 5), which uses 1975 earnings as the placebo outcome for the LaLonde 1978 earnings main outcome.",
  "\\end{tablenotes}",
  "\\end{threeparttable}",
  "\\end{table}"
)
writeLines(tex, file.path(output_dir, "tables", "lalonde_paired.tex"))
cat("Saved: tables/lalonde_paired.tex\n")
