## lalonde_estimation.R
## Estimates procurement markup premium using the 11 selection-on-observables
## estimators from Imbens & Xu (2025, JEP 39(4)) LaLonde replication package.
## Outcome: log(markup_A). Treatment: pp_dummy. Pooled cross-section.

suppressPackageStartupMessages({
  library(haven)
  library(dplyr)
})

input_dir    <- file.path("..", "input")
output_dir   <- file.path("..", "output")
lalonde_code <- "/Users/marek/Desktop/io/project/thesis/observables/code"

## Source lalonde estimator functions
source(file.path(lalonde_code, "functions_est.R"))

## The sourced `diff` masks base::diff; give it a unique name so our wrapper
## calls it directly rather than dispatching to the base generic.
diff_est <- diff

## ----- Load & merge data -----
df <- read_dta(file.path(input_dir, "data.dta"))
mk <- read_dta(file.path(output_dir, "data", "paper_markups.dta"))

## year is already integer in both files (not Stata date)
df$year <- as.integer(df$year)
mk$year <- as.integer(mk$year)

keep_df <- c("id", "year", "pp_dummy", "k", "cogs", "empl_mid",
             "foreign", "mktshare", "nace2")
keep_mk <- c("id", "year", "markup_A")

panel <- merge(df[, keep_df], mk[, keep_mk],
               by = c("id", "year"))
panel <- panel[!is.na(panel$markup_A) & panel$markup_A > 0, ]
panel$log_mu <- log(panel$markup_A)
panel$D      <- as.integer(panel$pp_dummy)
panel$foreign <- as.integer(panel$foreign)

## Impute empl_mid median for missing
med_empl <- median(panel$empl_mid, na.rm = TRUE)
panel$empl_mid[is.na(panel$empl_mid)] <- med_empl

## NACE dummies (nace2_41 omitted as base to avoid collinearity in some
## estimators; include 41 and 42 per task spec, drop 43 as reference)
panel$nace41 <- as.integer(panel$nace2 == 41)
panel$nace42 <- as.integer(panel$nace2 == 42)
panel$nace43 <- as.integer(panel$nace2 == 43)

## Drop any residual NAs on required covariates
req <- c("log_mu", "D", "k", "cogs", "empl_mid", "foreign", "mktshare",
         "nace41", "nace42", "nace43")
panel <- panel[complete.cases(panel[, req]), ]
panel <- as.data.frame(panel)

cat(sprintf("Sample: N=%d, treated=%d, control=%d\n",
            nrow(panel), sum(panel$D), sum(1 - panel$D)))

covariates <- c("k", "cogs", "empl_mid", "foreign", "mktshare",
                "nace41", "nace42")

## ----- Run estimators one-by-one with error handling -----
methods <- c("diff", "reg", "om.reg", "om.grf", "matching", "psm",
             "ipw", "cbps", "ebal", "dml", "aipw_grf")

res <- as.data.frame(matrix(NA_real_, length(methods), 4))
rownames(res) <- methods
colnames(res) <- c("Estimate", "SE", "CI_lower", "CI_upper")

run_one <- function(name, fn, needs_covar = TRUE) {
  cat(sprintf("Running %-10s ... ", name))
  out <- tryCatch({
    r <- if (needs_covar) {
      fn(panel, "log_mu", "D", covariates)
    } else {
      fn(panel, "log_mu", "D")
    }
    as.numeric(r)[1:4]
  }, error = function(e) {
    cat("FAILED: ", conditionMessage(e), "\n")
    return(rep(NA_real_, 4))
  })
  if (!any(is.na(out))) {
    cat(sprintf("ATT=%.4f SE=%.4f\n", out[1], out[2]))
  }
  out
}

set.seed(1234)
res["diff", ]      <- run_one("diff",     diff_est, needs_covar = FALSE)
res["reg", ]       <- run_one("reg",      reg)
res["om.reg", ]    <- run_one("om.reg",   om.reg)
res["om.grf", ]    <- run_one("om.grf",   om.grf)
res["matching", ]  <- run_one("matching", matching)
res["psm", ]       <- run_one("psm",      psm)
res["ipw", ]       <- run_one("ipw",      ipw)
res["cbps", ]      <- run_one("cbps",     cbps)
res["ebal", ]      <- run_one("ebal",     function(d,Y,Tr,X) quiet(ebal(d,Y,Tr,X)))
res["dml", ]       <- run_one("dml",      dml)
res["aipw_grf", ]  <- run_one("aipw_grf", aipw)

cat("\n--- Results ---\n")
print(round(res, 4))

## ----- Write outputs -----
write.csv(res, file.path(output_dir, "lalonde_results.csv"), row.names = TRUE)

## Pretty method labels
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

dir.create(file.path(output_dir, "tables"), showWarnings = FALSE,
           recursive = TRUE)

tex <- c(
  "\\begin{table}[htbp]\\centering",
  "\\caption{Selection on Observables: Eleven Estimators}\\label{tab:lalonde}",
  "\\begin{threeparttable}",
  "\\begin{tabular}{lccc}",
  "\\toprule",
  "Estimator & ATT & SE & 95\\% CI \\\\",
  "\\midrule"
)
for (m in methods) {
  row <- res[m, ]
  lab <- label_map[[m]]
  if (any(is.na(row))) {
    tex <- c(tex, sprintf("%s & --- & --- & --- \\\\", lab))
  } else {
    tex <- c(tex, sprintf("%s & %.4f & (%.4f) & [%.3f, %.3f] \\\\",
                          lab, row[[1]], row[[2]], row[[3]], row[[4]]))
  }
}
tex <- c(tex,
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{tablenotes}\\footnotesize",
  "\\item \\textit{Notes:} Eleven selection-on-observables estimators from Imbens and Xu (2025, \\textit{JEP} 39(4)). Outcome: $\\log \\mu^A_{it}$. Treatment: $pp_{it}$. Covariates: $k_{it}$, cogs$_{it}$, employment (midpoint), foreign-ownership dummy, market share, NACE-2 dummies. Pooled cross-section.",
  "\\end{tablenotes}",
  "\\end{threeparttable}",
  "\\end{table}"
)
writeLines(tex, file.path(output_dir, "tables", "lalonde_comparison.tex"))

cat("\nSaved: lalonde_results.csv, tables/lalonde_comparison.tex\n")
