## lalonde_balance.R
## Writes a normalized-difference balance table comparing the full panel to
## the Crump et al. alpha=0.1 trimmed subsample, following Imbens and
## Wooldridge (2007, NBER Summer Institute Lecture 1, Tables 1 and 4).
##
## Imbens and Wooldridge argue (pp. 33-34) that the normalized difference
##     (mean_T - mean_C) / sqrt((var_T + var_C) / 2)
## is the correct overlap diagnostic, *not* the t-statistic, because the
## t-statistic mechanically grows with sample size. Their rule of thumb is
## that |diff/sd| > 0.25 indicates a substantial imbalance.
##
## Outputs:
##   output/lalonde_balance.csv
##   output/tables/lalonde_balance.tex

source("lalonde_functions_helpers.R")

input_dir  <- file.path("..", "input")
output_dir <- file.path("..", "output")

dir.create(file.path(output_dir, "tables"), showWarnings = FALSE, recursive = TRUE)

## ---- Load panel and compute balance for full vs trimmed ----
panel   <- load_lalonde_panel(input_dir, output_dir)
ps      <- compute_propensity(panel)
trimmed <- trim_propensity(panel, ps, alpha = 0.1)

bal_full <- normalized_balance(panel)
bal_trim <- normalized_balance(trimmed)

combined <- data.frame(
  covariate    = bal_full$covariate,
  mean_T_full  = bal_full$mean_T,
  mean_C_full  = bal_full$mean_C,
  diff_full    = bal_full$diff_over_sd,
  mean_T_trim  = bal_trim$mean_T,
  mean_C_trim  = bal_trim$mean_C,
  diff_trim    = bal_trim$diff_over_sd
)

cat("Imbens-Wooldridge (2007) normalized balance:\n")
print(combined, digits = 3)

write.csv(combined, file.path(output_dir, "lalonde_balance.csv"),
          row.names = FALSE)
cat("Saved: lalonde_balance.csv\n")

## ---- LaTeX formatting ----
cov_label <- c(
  k        = "Capital $k_{it}$",
  cogs     = "Variable input $\\text{cogs}_{it}$",
  empl_mid = "Employment",
  foreign  = "Foreign ownership",
  mktshare = "Market share",
  nace41   = "NACE 41 (buildings)",
  nace42   = "NACE 42 (civil eng.)"
)

## Highlight cells with |diff/sd| > 0.25 using \mathbf{}
fmt_diff <- function(x) {
  if (is.na(x)) return("---")
  if (abs(x) > 0.25) sprintf("$\\mathbf{%+.3f}$", x)
  else               sprintf("$%+.3f$", x)
}

tex <- c(
  "\\begin{table}[htbp]\\centering",
  "\\caption{Normalized Covariate Balance: Full vs.\\ Trimmed Sample}\\label{tab:lalonde_balance}",
  "\\begin{threeparttable}",
  "\\begin{tabular}{lccc@{\\hspace{1em}}ccc}",
  "\\toprule",
  "& \\multicolumn{3}{c}{Full ($N = 7{,}666$)} & \\multicolumn{3}{c}{Trimmed ($N = 7{,}376$)} \\\\",
  "\\cmidrule(lr){2-4} \\cmidrule(lr){5-7}",
  "Covariate & Treated & Control & diff$/$sd & Treated & Control & diff$/$sd \\\\",
  "\\midrule"
)
for (i in seq_len(nrow(combined))) {
  cv  <- combined$covariate[i]
  lab <- cov_label[[cv]]
  tex <- c(tex, sprintf(
    "%s & %.3f & %.3f & %s & %.3f & %.3f & %s \\\\",
    lab,
    combined$mean_T_full[i], combined$mean_C_full[i], fmt_diff(combined$diff_full[i]),
    combined$mean_T_trim[i], combined$mean_C_trim[i], fmt_diff(combined$diff_trim[i])))
}
tex <- c(tex,
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{tablenotes}\\footnotesize",
  "\\item \\textit{Notes:} Normalized mean differences follow Imbens and Wooldridge (2007, NBER Summer Institute, Lecture 1): $\\text{diff}/\\text{sd} = (\\bar{X}_T - \\bar{X}_C) / \\sqrt{(s^2_T + s^2_C)/2}$. Values in bold exceed the 0.25 threshold that Imbens and Wooldridge (p.\\ 34) describe as substantial imbalance. Trimming uses the symmetric Crump, Hotz, Imbens, and Mitnik (2009) rule $\\hat{\\pi}(X_{it}) \\in [0.1, 0.9]$. The $t$-statistic for each difference is not reported because, as Imbens and Wooldridge note, it confounds covariate imbalance with sample size: a larger $t$-statistic for the difference between covariate means does not indicate a harder identification problem, whereas a larger normalized difference does.",
  "\\end{tablenotes}",
  "\\end{threeparttable}",
  "\\end{table}"
)
writeLines(tex, file.path(output_dir, "tables", "lalonde_balance.tex"))
cat("Saved: tables/lalonde_balance.tex\n")
