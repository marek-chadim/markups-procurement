################################################################################
# nnr_propensity.R — NNRPanel (Zeleneev & Zhang 2026) logit-IFE applied to
#                    the procurement-propensity model on the Czech balanced
#                    subpanel
#
# Applies the Zeleneev-Zhang (2026) two-step nonlinear panel estimator with
# interactive fixed effects to the selection equation
#
#   P(pp_dummy_it = 1 | cogs_it, IFE)
#     = Lambda(beta_c * cogs_it + alpha_i + gamma_t + lambda_i' f_t)
#
# Only one covariate (standardized log cogs) is included because NNRPanel's
# data-driven tuning-parameter selection and MLE step fail to converge
# cleanly when two collinear covariates (log k and log cogs) are included
# simultaneously on the narrow Czech subpanel.
#
# Panel: 58 switcher firms x 11 years (2010-2020). The balanced-subpanel
# restriction, combined with the conditional-logit-style dropping of firms
# whose pp_dummy is constant across the window, leaves a much smaller
# effective sample than the N=T=200 benchmark used in the NNRPanel README.
# Results should be interpreted as illustrative rather than load-bearing.
#
# Output: $output/tables/nnr_propensity.tex
#         $output/tables/nnr_propensity.csv
#
# Reference: Zeleneev and Zhang (2026) "Tractable Estimation of Nonlinear
# Panels with Interactive Fixed Effects" arXiv:2511.15427v2. Package:
# https://github.com/wszhang-econ/NNRPanel
################################################################################

suppressPackageStartupMessages({
  library(haven)
  library(NNRPanel)
  library(fixest)
})

cat("--- nnr_propensity.R ---\n")

args <- commandArgs(trailingOnly = TRUE)
data_path <- if (length(args) >= 1) args[1] else {
  "/Users/marek/Desktop/io/markups-procurement/2_analysis/output/stata/markups_panel.dta"
}
out_tables <- if (length(args) >= 2) args[2] else {
  "/Users/marek/Desktop/io/markups-procurement/2_analysis/output/stata/tables"
}
dir.create(out_tables, showWarnings = FALSE, recursive = TRUE)

# ---- Build balanced subpanel of switchers -------------------------------
df <- read_dta(data_path)
df <- df[!is.na(df$mu_A) & df$mu_A > 0 & !is.na(df$k) &
         !is.na(df$cogs) & !is.na(df$pp_dummy), ]

win_start <- 2010L
win_end   <- 2020L
df <- df[df$year >= win_start & df$year <= win_end, ]
n_years <- length(unique(df$year))
firm_nobs <- table(df$id)
bal_ids <- as.integer(names(firm_nobs[firm_nobs == n_years]))
df <- df[df$id %in% bal_ids, ]

firm_var <- tapply(df$pp_dummy, df$id, function(v) length(unique(v)) > 1)
switcher_ids <- as.integer(names(firm_var[firm_var]))
df <- df[df$id %in% switcher_ids, ]
df <- df[order(df$id, df$year), ]

N <- length(switcher_ids)
T <- n_years
cat(sprintf("  Balanced switcher subpanel: N = %d, T = %d, obs = %d (window %d-%d)\n",
            N, T, nrow(df), win_start, win_end))
stopifnot(nrow(df) == N * T)

df$id_nnr   <- as.integer(factor(df$id))
df$t_nnr    <- df$year - (win_start - 1L)
df$cogs_std <- as.numeric(scale(df$cogs))

# ---- NNRPanel: logit-IFE on pp_dummy ~ cogs_std + IFE ------------------
data_frame <- cbind(df$id_nnr, df$t_nnr, df$pp_dummy, df$cogs_std)
set.seed(42)
cat("\n  Running NNRPanel_estimate(func='logit', R_max=2, R_true=1)...\n")
est <- tryCatch(
  NNRPanel_estimate(data_frame, func = "logit",
                    R_max = 2, R_true = 1, delta = 0.05, iter_max = 500),
  error = function(e) { cat("  ERROR: ", conditionMessage(e), "\n"); NULL }
)

if (is.null(est)) {
  cat("  NNRPanel estimation failed; writing placeholder table.\n")
  est <- list(beta_nnr = NA_real_, beta_fe = NA_real_,
              beta_corr = NA_real_, beta_corr_sp = NA_real_,
              std_beta_corr = NA_real_, num_factor_est = NA_integer_)
}

cat("\n  --- NNRPanel logit-IFE results (std log cogs) ---\n")
cat(sprintf("  beta_nnr (step 1, nuclear-norm):       %.4f\n", as.numeric(est$beta_nnr)))
cat(sprintf("  beta_fe  (step 2, local MLE):          %.4f\n", as.numeric(est$beta_fe)))
cat(sprintf("  beta_corr (analytical bias-corrected): %.4f\n", as.numeric(est$beta_corr)))
cat(sprintf("  std_beta_corr:                          %.4f\n", as.numeric(est$std_beta_corr)))
cat(sprintf("  num_factor_est:                         %d\n", as.integer(est$num_factor_est)))

# ---- Additive-FE logit baseline ----------------------------------------
df$firm_f <- factor(df$id_nnr)
df$year_f <- factor(df$t_nnr)
fe_fit <- tryCatch(
  feglm(pp_dummy ~ cogs_std | firm_f + year_f,
        data = df, family = "logit"),
  error = function(e) { cat("  feglm error: ", conditionMessage(e), "\n"); NULL }
)
if (!is.null(fe_fit)) {
  fe_beta <- coef(fe_fit)["cogs_std"]
  fe_se   <- sqrt(diag(vcov(fe_fit)))["cogs_std"]
} else {
  fe_beta <- NA_real_; fe_se <- NA_real_
}
cat(sprintf("\n  Additive-FE logit: beta_cogs = %.4f (SE %.4f)\n", fe_beta, fe_se))

# ---- Save outputs ------------------------------------------------------
results <- data.frame(
  estimator     = c("Additive FE logit (feglm)",
                    "NNR step 1",
                    "NNR step 2 (local MLE)",
                    "Bias-corrected (analytical)"),
  beta          = c(as.numeric(fe_beta),
                    as.numeric(est$beta_nnr),
                    as.numeric(est$beta_fe),
                    as.numeric(est$beta_corr)),
  se            = c(as.numeric(fe_se),
                    NA_real_,
                    NA_real_,
                    as.numeric(est$std_beta_corr))
)
cat("\n  --- Summary ---\n")
for (i in seq_len(nrow(results))) {
  cat(sprintf("  %-32s beta = %s  se = %s\n",
              results$estimator[i],
              if (is.na(results$beta[i])) "NA" else sprintf("%.4f", results$beta[i]),
              if (is.na(results$se[i]))   "NA" else sprintf("%.4f", results$se[i])))
}

csv_path <- file.path(out_tables, "nnr_propensity.csv")
write.csv(results, csv_path, row.names = FALSE)
cat("  Saved: ", csv_path, "\n", sep = "")

tex_path <- file.path(out_tables, "nnr_propensity.tex")
sink(tex_path)
cat("\\begin{table}[htbp]\\centering\n")
cat("\\caption{NNRPanel Logit-IFE Robustness for the Procurement Propensity Model}\n")
cat("\\label{tab:nnr_propensity}\n")
cat("\\begin{threeparttable}\n")
cat("\\begin{tabular}{lcc}\n")
cat("\\toprule\n")
cat("Estimator & $\\hat\\beta_{\\text{cogs}}$ & SE \\\\\n")
cat("\\midrule\n")
for (i in seq_len(nrow(results))) {
  r <- results[i, ]
  beta_str <- if (is.na(r$beta)) "--" else sprintf("%.4f", r$beta)
  se_str   <- if (is.na(r$se))   "--" else sprintf("%.4f", r$se)
  cat(sprintf("%s & %s & %s \\\\\n", r$estimator, beta_str, se_str))
}
cat("\\midrule\n")
cat(sprintf("\\multicolumn{3}{l}{Balanced switcher subpanel: $N = %d$ firms, $T = %d$ years, window %d--%d} \\\\\n",
            N, T, win_start, win_end))
if (!is.na(est$num_factor_est)) {
  cat(sprintf("\\multicolumn{3}{l}{Data-driven number of factors: $\\hat R = %d$} \\\\\n",
              as.integer(est$num_factor_est)))
}
cat("\\bottomrule\n")
cat("\\end{tabular}\n")
cat("\\begin{tablenotes}\\footnotesize\n")
cat("\\item \\emph{Notes:} Dependent variable is the firm-year procurement-participation indicator \\texttt{pp\\_dummy}; the single covariate is standardized log cogs. ")
cat("The first row is a conditional logit with additive firm and year fixed effects via \\texttt{fixest::feglm}, which is the propensity-model baseline used in the LaLonde battery of Section~\\ref{sec:unconfound}. ")
cat("The next three rows are the two-step Zeleneev-Zhang \\cite{ZeleneevZhang2026} nonlinear IFE estimator via the \\texttt{NNRPanel} R package: ``NNR step 1'' is the first-step nuclear-norm regularized estimator, ``NNR step 2'' is the second-step local MLE initialized at the NNR estimate, and the bias-corrected row applies the analytical bias correction of Chen, Fern\\'{a}ndez-Val, and Weidner \\cite{ChenFernandezValWeidner2021}. ")
cat("The panel is restricted to firms present in every year of the window and with at least one procurement switch, leaving $N = 58$ switcher firms. ")
cat("Only a single covariate is included because the data-driven tuning-parameter selection of NNRPanel does not converge cleanly with two potentially collinear covariates on a panel of this size; the method's asymptotics require both $N$ and $T$ to be large in the sense of the README benchmark ($N = T = 200$), and the Czech balanced switcher subpanel is well below that regime. ")
cat("The result should be read as an illustrative application of the nonlinear IFE estimator to the selection equation rather than as a load-bearing robustness check.\n")
cat("\\end{tablenotes}\n")
cat("\\end{threeparttable}\n")
cat("\\end{table}\n")
sink()
cat("  Saved: ", tex_path, "\n", sep = "")

cat("\n--- nnr_propensity.R DONE ---\n")
