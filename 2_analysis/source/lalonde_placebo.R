## lalonde_placebo.R
## Structural adaptation of the Imbens & Xu (2025, JEP 39(4)) re75 placebo to
## Czech construction panel. For each firm, take the log markup at the first
## untreated firm-year as a sham outcome and define the placebo treatment as
## "was ever treated at any later year in the sample". Under the conditional
## independence assumption, an unconfounded estimator should return ATT ~ 0
## because pre-treatment markups are not caused by future procurement. Any
## systematic positive premium on the placebo flags selection-on-unobservables
## or a violation of the parallel-trends assumption implicit in the CIA.

suppressPackageStartupMessages({})

source("lalonde_functions_helpers.R")

input_dir    <- file.path("..", "input")
output_dir   <- file.path("..", "output")
lalonde_code <- "/Users/marek/Desktop/io/project/thesis/observables/code"

source(file.path(lalonde_code, "functions_est.R"))
diff_est <- diff

dir.create(file.path(output_dir, "tables"), showWarnings = FALSE, recursive = TRUE)

## ---- Build pre-treatment cross-section ----
panel <- load_lalonde_panel(input_dir, output_dir)
cat(sprintf("Panel: N=%d firm-years, firms=%d\n",
            nrow(panel), length(unique(panel$id))))

panel <- panel[order(panel$id, panel$year), ]

## For each firm: earliest observation where D = 0 (not yet treated).
## Firms whose very first observation has D = 1 are dropped because their
## pre-treatment state is unobserved in the panel.
first_untreated <- do.call(rbind, lapply(split(panel, panel$id), function(g) {
  untreated <- g[g$D == 0, ]
  if (nrow(untreated) == 0) return(NULL)
  first_row   <- untreated[which.min(untreated$year), ]
  ever_later  <- any(g$year > first_row$year & g$D == 1)
  first_row$placebo_D      <- as.integer(ever_later)
  first_row$first_year     <- first_row$year
  first_row
}))
cross <- as.data.frame(first_untreated)
cat(sprintf("Pre-treatment cross-section: N=%d firms, eventual winners=%d, never winners=%d\n",
            nrow(cross), sum(cross$placebo_D), sum(1 - cross$placebo_D)))
cat("Distribution of first untreated year:\n")
print(table(cross$first_year))

## Overwrite D with placebo_D for the batch estimator
cross$log_mu_placebo <- cross$log_mu
cross$D              <- cross$placebo_D

## ---- Run 11-method batch on the placebo cross-section ----
methods <- c("diff", "reg", "om.reg", "om.grf", "matching", "psm",
             "ipw", "cbps", "ebal", "dml", "aipw_grf")
covariates <- LALONDE_COVARIATES

run_one <- function(name, fn, data, needs_covar = TRUE) {
  cat(sprintf("  %-10s ", name))
  out <- tryCatch({
    r <- if (needs_covar) fn(data, "log_mu_placebo", "D", covariates)
         else             fn(data, "log_mu_placebo", "D")
    as.numeric(r)[1:4]
  }, error = function(e) {
    cat("FAILED: ", conditionMessage(e), "\n")
    return(rep(NA_real_, 4))
  })
  if (!any(is.na(out))) cat(sprintf("ATT=%.4f SE=%.4f\n", out[1], out[2]))
  out
}

res <- as.data.frame(matrix(NA_real_, length(methods), 4))
rownames(res) <- methods
colnames(res) <- c("Estimate", "SE", "CI_lower", "CI_upper")

cat("\n[Placebo] Running 11 methods on N =", nrow(cross), "\n")
set.seed(1234)
res["diff", ]     <- run_one("diff",     diff_est, cross, needs_covar = FALSE)
res["reg", ]      <- run_one("reg",      reg,      cross)
res["om.reg", ]   <- run_one("om.reg",   om.reg,   cross)
res["om.grf", ]   <- run_one("om.grf",   om.grf,   cross)
res["matching", ] <- run_one("matching", matching, cross)
res["psm", ]      <- run_one("psm",      psm,      cross)
res["ipw", ]      <- run_one("ipw",      ipw,      cross)
res["cbps", ]     <- run_one("cbps",     cbps,     cross)
res["ebal", ]     <- run_one("ebal",     function(d,Y,Tr,X) quiet(ebal(d,Y,Tr,X)), cross)
res["dml", ]      <- run_one("dml",      dml,      cross)
res["aipw_grf", ] <- run_one("aipw_grf", aipw,     cross)

cat("\n--- Placebo Results ---\n")
print(round(res, 4))

write.csv(res, file.path(output_dir, "lalonde_placebo_results.csv"),
          row.names = TRUE)
cat("\nSaved: lalonde_placebo_results.csv\n")

## ---- LaTeX table ----
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
  "\\caption{Placebo: Pre-Treatment Markup vs.\\ Future Procurement}\\label{tab:lalonde_placebo}",
  "\\begin{threeparttable}",
  "\\begin{tabular}{lccc}",
  "\\toprule",
  "Estimator & Placebo ATT & SE & 95\\% CI \\\\",
  "\\midrule"
)
for (m in methods) {
  row <- res[m, ]
  lab <- label_map[[m]]
  if (any(is.na(row))) {
    tex <- c(tex, sprintf("%s & --- & --- & --- \\\\", lab))
  } else {
    tex <- c(tex, sprintf("%s & %.4f & (%.4f) & [%.4f, %.4f] \\\\",
                          lab, row[[1]], row[[2]], row[[3]], row[[4]]))
  }
}
tex <- c(tex,
  "\\midrule",
  sprintf("$N$ firms & \\multicolumn{3}{c}{%d} \\\\", nrow(cross)),
  sprintf("Eventual winners & \\multicolumn{3}{c}{%d} \\\\", sum(cross$D)),
  sprintf("Never winners & \\multicolumn{3}{c}{%d} \\\\", sum(1 - cross$D)),
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{tablenotes}\\footnotesize",
  "\\item \\textit{Notes:} Placebo test paralleling the Imbens and Xu (2025, \\textit{JEP} 39(4)) \\texttt{re75} placebo. Each firm contributes a single cross-sectional observation at its first year in the panel with $pp_{it} = 0$ (``pre-treatment'' state). The outcome is the log markup at that firm-year; the placebo treatment is an indicator for whether the firm was ever treated ($pp_{it} = 1$) in any later year in the panel. Covariates are those in Table~\\ref{tab:lalonde}, taken at the pre-treatment year. Under the conditional independence assumption, placebo ATTs should lie near zero: a positive estimate flags selection-on-unobservables because future procurement cannot causally raise current markups.",
  "\\end{tablenotes}",
  "\\end{threeparttable}",
  "\\end{table}"
)
writeLines(tex, file.path(output_dir, "tables", "lalonde_placebo.tex"))
cat("Saved: tables/lalonde_placebo.tex\n")
