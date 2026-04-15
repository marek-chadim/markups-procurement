################################################################################
# awz_panel_ife.R — Armstrong, Weidner, Zeleneev (2025) robust interactive
#                   fixed-effects estimation of the 2012 reform DiD
#
# Replaces the prior staggered-DiD battery (Callaway-Sant'Anna, Sun-Abraham,
# Borusyak-Jaravel-Spiess, Goodman-Bacon) and the de Chaisemartin-D'Haultfoeuille
# cumulative-effects result with a single methodologically unified panel
# robustness check: the AWZ 2025 debiased estimator with bias-aware CI, valid
# uniformly across weak/strong/no-factor DGPs.
#
# Package: PanelIFE (chenweihsiang/PanelIFE on GitHub), implements
#   - ls_factor()            Bai (2009) / Moon-Weidner (2017) LS estimator
#   - honest_weak_factors()  Armstrong-Weidner-Zeleneev (2025) debiased + CI
#
# Panel restriction: PanelIFE requires a balanced panel (N x T matrices).
# The Czech construction panel is unbalanced (9,164 obs over 1,521 firms,
# 17 years). This script restricts to the 2010-2020 window, keeping firms
# present in all 11 years (approximately 108 firms, 1,188 observations).
# The window covers all three transparency reforms (2012, 2016, 2017).
#
# Pipeline:
#   1. Load markups_panel.dta via haven
#   2. Restrict to balanced subpanel (2010-2020)
#   3. Partial out additive firm FE + year FE + log k + log cogs
#      from log_markup and tp2012 = ever_pp x post2012
#   4. Reshape to N x T matrices
#   5. Run ls_factor() and honest_weak_factors() for R in {0, 1, 2, 3}
#   6. Write a summary .tex table and a CSV with the numbers
#
# All numerical results are also printed to stdout so prose claims can cite
# exact values without reading the figure.
#
# Input:  $data/markups_panel.dta
# Output: $output/tables/awz_panel_ife.tex
#         $output/tables/awz_panel_ife.csv
#         $output/data/awz_panel_ife.RData   (full R results for audit)
################################################################################

suppressPackageStartupMessages({
  library(haven)
  library(PanelIFE)
  library(Matrix)
})

cat("--- awz_panel_ife.R ---\n")

# ---- Paths --------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
data_path <- if (length(args) >= 1) args[1] else {
  "/Users/marek/Desktop/io/markups-procurement/2_analysis/output/stata/markups_panel.dta"
}
out_tables <- if (length(args) >= 2) args[2] else {
  "/Users/marek/Desktop/io/markups-procurement/2_analysis/output/stata/tables"
}
out_data <- if (length(args) >= 3) args[3] else {
  "/Users/marek/Desktop/io/markups-procurement/2_analysis/output/stata/data"
}

dir.create(out_tables, showWarnings = FALSE, recursive = TRUE)
dir.create(out_data,   showWarnings = FALSE, recursive = TRUE)

# ---- Load data ----------------------------------------------------------
cat("  Loading ", data_path, "\n", sep = "")
df <- read_dta(data_path)

# Keep variables we need
df <- df[, c("id", "year", "mu_A", "pp_dummy", "k", "cogs", "nace2")]
df <- df[!is.na(df$mu_A) & df$mu_A > 0, ]
df <- df[!is.na(df$k) & !is.na(df$cogs), ]
df$log_mu <- log(df$mu_A)

# Ever-PP indicator (firm-level)
ever_pp <- tapply(df$pp_dummy, df$id, max, na.rm = TRUE)
df$ever_pp <- ever_pp[as.character(df$id)]

# Post-2012 indicator
df$post2012 <- as.integer(df$year >= 2012)
df$tp2012   <- df$ever_pp * df$post2012

# ---- Restrict to balanced subpanel 2010-2020 -----------------------------
win_start <- 2010L
win_end   <- 2020L
df_win <- df[df$year >= win_start & df$year <= win_end, ]
n_years <- length(unique(df_win$year))
stopifnot(n_years == (win_end - win_start + 1))
cat("  Window: ", win_start, "-", win_end, " (", n_years, " years)\n", sep = "")

# Keep only firms present in all years
firm_nobs <- table(df_win$id)
balanced_ids <- as.integer(names(firm_nobs[firm_nobs == n_years]))
df_bal <- df_win[df_win$id %in% balanced_ids, ]
N_firms <- length(balanced_ids)
T_years <- n_years
cat("  Balanced firms: ", N_firms, "\n", sep = "")
cat("  Balanced obs:   ", nrow(df_bal), " (should equal ", N_firms * T_years, ")\n", sep = "")
stopifnot(nrow(df_bal) == N_firms * T_years)

# ---- Partial out additive FEs and controls from log_mu and tp2012 -------
# Fit two auxiliary regressions with firm FE + year FE + k + cogs and take
# residuals. The resulting Y_perp and X_perp are the inputs to PanelIFE
# (which then estimates the residual interactive factor structure).
df_bal$firm <- factor(df_bal$id)
df_bal$yr   <- factor(df_bal$year)

lm_y <- lm(log_mu ~ firm + yr + k + cogs, data = df_bal)
lm_x <- lm(tp2012 ~ firm + yr + k + cogs, data = df_bal)

df_bal$y_perp <- residuals(lm_y)
df_bal$x_perp <- residuals(lm_x)

# ---- Reshape to N x T matrices ------------------------------------------
df_bal <- df_bal[order(df_bal$id, df_bal$year), ]
Y <- matrix(df_bal$y_perp, nrow = N_firms, ncol = T_years, byrow = TRUE)
X_mat <- matrix(df_bal$x_perp, nrow = N_firms, ncol = T_years, byrow = TRUE)
X <- array(X_mat, dim = c(1, N_firms, T_years))   # K x N x T tensor with K=1
cat("  Y dim: ", dim(Y), "\n")
cat("  X dim: ", dim(X), "\n")

# ---- Run estimators for R in {0, 1, 2, 3} -------------------------------
R_values <- c(0, 1, 2, 3)
# Each row reports:
#   LS (Bai 2009) point + CI
#   AWZ bias-aware CI under R_w = R (worst case: all R factors could be weak)
#   AWZ no-adjustment CI under R_w = 0 (all factors treated as strong)
results <- data.frame(
  R             = integer(),
  ls_beta       = double(),
  ls_se         = double(),
  ls_lo         = double(),
  ls_hi         = double(),
  awz_beta      = double(),
  awz_se        = double(),
  awz_biasaware_lo = double(),
  awz_biasaware_hi = double(),
  awz_noadj_lo  = double(),
  awz_noadj_hi  = double()
)

set.seed(42)

for (R in R_values) {
  cat("\n  --- R = ", R, " ---\n", sep = "")

  # Least-squares estimator (Bai 2009 / Moon-Weidner)
  if (R == 0) {
    # R=0 means no factors: the LS estimator reduces to pooled OLS on partialed-out
    # variables, which equals the additive-FE TWFE coefficient by FWL theorem.
    lm_pooled <- lm(y_perp ~ x_perp - 1, data = df_bal)
    ls_beta <- coef(lm_pooled)[1]
    ls_se   <- summary(lm_pooled)$coefficients[1, "Std. Error"]
  } else {
    ls_res <- tryCatch(
      ls_factor(Y = Y, X = X, R = R, report = "silent",
                precision_beta = 1e-8, method = "m1",
                start = c(0), repMIN = 3, repMAX = 10, M1 = 2, M2 = 2),
      error = function(e) { cat("    ls_factor error: ", conditionMessage(e), "\n"); NULL }
    )
    if (is.null(ls_res)) {
      ls_beta <- NA_real_
      ls_se   <- NA_real_
    } else {
      ls_beta <- as.numeric(ls_res$beta)
      ls_se   <- sqrt(as.numeric(ls_res$Vbeta2))  # heteroskedastic SE
    }
  }
  ls_lo <- ls_beta - 1.96 * ls_se
  ls_hi <- ls_beta + 1.96 * ls_se
  cat("    LS   beta = ", sprintf("%.4f", ls_beta),
      " SE = ", sprintf("%.4f", ls_se),
      " CI = [", sprintf("%.4f", ls_lo), ", ", sprintf("%.4f", ls_hi), "]\n", sep = "")

  # Armstrong-Weidner-Zeleneev debiased estimator with bias-aware CI.
  # The return has $beta, $se, and $LB/$UB as data.frames with two rows:
  #   R_w = 0 : "CI without adjustment" (treats all R factors as strong)
  #   R_w > 0 : bias-aware CI allowing that many factors to be weak
  if (R == 0) {
    # R=0 means no factors: AWZ reduces to additive-FE TWFE by FWL.
    awz_beta <- ls_beta
    awz_se   <- ls_se
    awz_biasaware_lo <- ls_lo
    awz_biasaware_hi <- ls_hi
    awz_noadj_lo <- ls_lo
    awz_noadj_hi <- ls_hi
  } else {
    awz_res <- tryCatch(
      honest_weak_factors(Y = Y, X = X, R = R,
                          Gamma_LS = NULL, alpha = 0.05,
                          clustered_se = TRUE),
      error = function(e) { cat("    honest_weak_factors error: ", conditionMessage(e), "\n"); NULL }
    )
    if (is.null(awz_res)) {
      awz_beta <- NA_real_
      awz_se   <- NA_real_
      awz_biasaware_lo <- NA_real_
      awz_biasaware_hi <- NA_real_
      awz_noadj_lo <- NA_real_
      awz_noadj_hi <- NA_real_
    } else {
      awz_beta <- as.numeric(awz_res$beta)
      awz_se   <- as.numeric(awz_res$se)
      # LB / UB are data.frames with NumberOfWeakFactors and Est1 columns.
      # Row 1: R_w=0 (no adjustment); Row 2: R_w = max (bias-aware worst case).
      lb_df <- awz_res$LB
      ub_df <- awz_res$UB
      awz_noadj_lo <- as.numeric(lb_df$Est1[lb_df$NumberOfWeakFactors == 0])
      awz_noadj_hi <- as.numeric(ub_df$Est1[ub_df$NumberOfWeakFactors == 0])
      # Use the largest reported weak-factor count as the bias-aware worst case
      max_rw <- max(lb_df$NumberOfWeakFactors)
      awz_biasaware_lo <- as.numeric(lb_df$Est1[lb_df$NumberOfWeakFactors == max_rw])
      awz_biasaware_hi <- as.numeric(ub_df$Est1[ub_df$NumberOfWeakFactors == max_rw])
    }
  }
  cat("    AWZ  beta = ", sprintf("%.4f", awz_beta),
      " SE = ", sprintf("%.4f", awz_se), "\n", sep = "")
  cat("       noadj CI    = [", sprintf("%.4f", awz_noadj_lo), ", ",
      sprintf("%.4f", awz_noadj_hi), "]\n", sep = "")
  cat("       biasaware CI = [", sprintf("%.4f", awz_biasaware_lo), ", ",
      sprintf("%.4f", awz_biasaware_hi), "]\n", sep = "")

  results <- rbind(results, data.frame(
    R = R, ls_beta = ls_beta, ls_se = ls_se, ls_lo = ls_lo, ls_hi = ls_hi,
    awz_beta = awz_beta, awz_se = awz_se,
    awz_biasaware_lo = awz_biasaware_lo, awz_biasaware_hi = awz_biasaware_hi,
    awz_noadj_lo = awz_noadj_lo, awz_noadj_hi = awz_noadj_hi
  ))
}

cat("\n  --- Summary ---\n")
print(round(results, 4))

# ---- Write outputs ------------------------------------------------------
csv_path <- file.path(out_tables, "awz_panel_ife.csv")
write.csv(results, csv_path, row.names = FALSE)
cat("  Saved: ", csv_path, "\n", sep = "")

save(results, Y, X, file = file.path(out_data, "awz_panel_ife.RData"))
cat("  Saved: awz_panel_ife.RData\n")

# ---- Write LaTeX table --------------------------------------------------
tex_path <- file.path(out_tables, "awz_panel_ife.tex")
sink(tex_path)
cat("\\begin{table}[htbp]\\centering\n")
cat("\\caption{Interactive Fixed Effects Robustness for the 2012 Reform DiD}\n")
cat("\\label{tab:awz_panel_ife}\n")
cat("\\begin{threeparttable}\n")
cat("\\begin{tabular}{lccccc}\n")
cat("\\toprule\n")
cat(" & \\multicolumn{2}{c}{Least-squares IFE} & \\multicolumn{3}{c}{AWZ bias-aware} \\\\\n")
cat("\\cmidrule(lr){2-3}\\cmidrule(lr){4-6}\n")
cat("$R$ & $\\hat\\beta_{\\text{LS}}$ & 95\\% CI & $\\hat\\beta_{\\text{AWZ}}$ & CI ($R_w=0$) & CI ($R_w=R$) \\\\\n")
cat("\\midrule\n")
for (i in seq_len(nrow(results))) {
  r <- results[i, ]
  cat(sprintf("%d & %.4f & [%.4f, %.4f] & %.4f & [%.4f, %.4f] & [%.4f, %.4f] \\\\\n",
              r$R, r$ls_beta, r$ls_lo, r$ls_hi, r$awz_beta,
              r$awz_noadj_lo, r$awz_noadj_hi,
              r$awz_biasaware_lo, r$awz_biasaware_hi))
}
cat("\\midrule\n")
cat(sprintf("\\multicolumn{6}{l}{Balanced panel: $N = %d$ firms, $T = %d$ years, window %d--%d} \\\\\n",
            N_firms, T_years, win_start, win_end))
cat("\\bottomrule\n")
cat("\\end{tabular}\n")
cat("\\begin{tablenotes}\\footnotesize\n")
cat("\\item \\emph{Notes:} Dependent variable is log ACF markup; treatment is $\\text{ever\\_pp} \\times \\text{post2012}$. ")
cat("Additive firm and year fixed effects plus log capital and log cogs controls are partialed out of the outcome and treatment before applying PanelIFE, so $R$ counts \\emph{additional} interactive factors beyond the additive terms. ")
cat("Least-squares IFE: Bai \\cite{Bai2009} / Moon and Weidner \\cite{MoonWeidner2015} LS estimator via \\texttt{PanelIFE::ls\\_factor}; CI from heteroskedasticity-robust variance. ")
cat("AWZ bias-aware: Armstrong, Weidner, and Zeleneev \\cite{ArmstrongWeidnerZeleneev2024} debiased estimator via \\texttt{PanelIFE::honest\\_weak\\_factors}. ")
cat("Two AWZ CIs are reported: CI\\,($R_w=0$) imposes that all $R$ factors are strong (reduces to a Bai-style bias-aware CI); CI\\,($R_w=R$) is the \\emph{bias-aware} worst-case interval allowing that up to $R$ of the factors could be weak, which is uniformly valid regardless of factor strength. ")
cat("$R = 0$ reduces to the additive-FE TWFE estimator by Frisch-Waugh. ")
cat("Panel is restricted to firms present in every year of the window to satisfy PanelIFE's balanced-panel requirement.\n")
cat("\\end{tablenotes}\n")
cat("\\end{threeparttable}\n")
cat("\\end{table}\n")
sink()
cat("  Saved: ", tex_path, "\n", sep = "")

cat("\n--- awz_panel_ife.R DONE ---\n")
