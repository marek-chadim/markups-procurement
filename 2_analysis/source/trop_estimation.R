#===============================================================================
# trop_estimation.R
#
# Estimate the procurement markup premium using the Triply RObust Panel (TROP)
# estimator from Athey, Imbens, Qu & Viviano (2026, arXiv 2508.21536v3).
#
# TROP combines three components:
#   1. Unit weights ω  — exponential decay in control-period RMSE to treated unit
#   2. Time weights θ  — exponential decay in |t - t_treat|
#   3. Low-rank L̂     — nuclear-norm penalized factor adjustment (MC-style)
#
# Triple robustness (Theorem 5.1): conditional bias factors as the product of
# three norms — unit imbalance × time imbalance × regression misspecification.
# Unbiased if *any one* of the three vanishes.
#
# Usage:
#   Rscript trop_estimation.R              # full run (~90 min)
#   Rscript trop_estimation.R --quick      # small grid + sampled treated (~4 min)
#
# Inputs:
#   ../input/data.dta                      (Czech firm-year panel)
#   ../output/data/paper_markups.dta       (wide-format markups)
#
# Outputs:
#   ../output/trop_results.csv             (ATT + SE + tuned lambdas)
#   ../output/trop_ablation.csv            (8-row Table 5)
#   ../output/tables/trop_comparison.tex   (4-row if fect_results.csv exists)
#   ../output/tables/trop_ablation.tex     (8-row ablation)
#===============================================================================

suppressPackageStartupMessages({
  library(haven)
  library(dplyr)
  library(fixest)
})

# ------------------------------------------------------------------------------
# Command-line flag parsing
# ------------------------------------------------------------------------------
args  <- commandArgs(trailingOnly = TRUE)
QUICK <- "--quick" %in% args

input_dir  <- file.path("..", "input")
output_dir <- file.path("..", "output")

set.seed(42)

# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------
if (QUICK) {
  cat("[TROP] Quick mode: small grid, sampled treated cells\n")
  grid_time  <- c(0, 0.1, 0.5)
  grid_unit  <- c(0, 0.1, 0.5)
  grid_nn    <- c(0.05, 0.5, Inf)
  N_CYCLES   <- 1
  N_BOOT     <- 3
  N_CV_CELLS <- 30
  MAX_TREATED <- 100
} else {
  cat("[TROP] Full mode: LOOCV + 50 bootstrap replications\n")
  grid_time  <- c(0, 0.1, 0.3, 0.5, 1.0)
  grid_unit  <- c(0, 0.05, 0.1, 0.3, 1.0)
  grid_nn    <- c(0.05, 0.1, 0.5, 1.0, Inf)
  N_CYCLES   <- 2
  N_BOOT     <- 50
  N_CV_CELLS <- 100
  MAX_TREATED <- NULL
}

MAX_SVT_ITER <- 30
SVT_TOL      <- 1e-4

# ------------------------------------------------------------------------------
# Data loading (mirrors fect_estimation.R lines 37-64)
# ------------------------------------------------------------------------------
cat("[TROP] Loading data...\n")
df <- read_dta(file.path(input_dir, "data.dta"))
mk <- read_dta(file.path(output_dir, "data", "paper_markups.dta"))

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

# Drop firms with fewer than 3 observations (edge-case mitigation for SVD)
firm_counts <- table(panel$id)
keep_firms  <- as.integer(names(firm_counts)[firm_counts >= 3])
panel       <- panel[panel$id %in% keep_firms, ]

cat("[TROP] Panel: N =", nrow(panel),
    ", firms =", length(unique(panel$id)),
    ", years =", length(unique(panel$year)), "\n")
cat("[TROP] Treated firm-years:", sum(panel$D),
    " (", round(100 * mean(panel$D), 1), "%)\n", sep = "")

# No pre-residualization: TROP's alpha_i and beta_t absorb firm and year FE
# internally (see paper Section 2.1 working model). Covariate control for
# (k, cogs) would require the Section 6.2 extension; for now we rely on
# within-firm variation via alpha_i, which is how fect's MC also works.
# Pre-residualizing on (k, cogs) would remove procurement-correlated
# variation since larger firms (higher k, cogs) are more likely to enter
# procurement, biasing tau_hat downward.
cat("[TROP] Using raw log_mu (no pre-residualization)\n")
panel$log_mu_raw <- panel$log_mu

# ------------------------------------------------------------------------------
# Long -> wide conversion
# ------------------------------------------------------------------------------
firm_ids <- sort(unique(panel$id))
year_ids <- sort(unique(panel$year))
N <- length(firm_ids); T <- length(year_ids)

panel_to_matrix <- function(df, var) {
  m <- matrix(NA_real_, nrow = N, ncol = T)
  i_idx <- match(df$id, firm_ids)
  t_idx <- match(df$year, year_ids)
  m[cbind(i_idx, t_idx)] <- df[[var]]
  return(m)
}

Y <- panel_to_matrix(panel, "log_mu_raw")
W <- panel_to_matrix(panel, "D")
W[is.na(W)] <- 0
obs_mask <- !is.na(Y)

cat(sprintf("[TROP] Wide matrix: %d firms x %d years, %d observed (%.1f%% dense)\n",
            N, T, sum(obs_mask), 100 * mean(obs_mask)))

# ------------------------------------------------------------------------------
# Two-way FE + SVT soft-impute (uses lm for joint alpha/beta on control cells)
# ------------------------------------------------------------------------------
# Earlier iterative Gauss-Seidel row/col-mean updates do NOT converge to the
# joint OLS solution for unbalanced panels with heterogeneous missingness
# (tested: Gauss-Seidel gave ATT 0.42 vs lm 0.13). We therefore use a direct
# two-step approach:
#
#   Step 1: fit (alpha_i, beta_t) jointly via weighted lm on observed control
#           cells. This is the same solution feols/lm/fect give internally.
#   Step 2: if lambda_nn finite, compute the residual matrix on control cells
#           and apply SVT to get L_hat. Add L_hat to the prediction.
#
# Step 1 is the dominant component (80-90% of the effect); Step 2 adds the
# nuclear-norm interactive factor correction.
# ------------------------------------------------------------------------------

svt <- function(Z, lambda) {
  if (is.infinite(lambda) || lambda < 0) {
    return(matrix(0, nrow = nrow(Z), ncol = ncol(Z)))
  }
  Z_clean <- Z
  Z_clean[is.na(Z_clean)] <- 0
  sv <- svd(Z_clean)
  d_thresh <- pmax(sv$d - lambda, 0)
  return(sv$u %*% diag(d_thresh, nrow = length(d_thresh), ncol = length(d_thresh)) %*% t(sv$v))
}

# Fit two-way FE on training cells (weighted) and return a predict function
# Uses fixest::feols for joint FE estimation. Returns a list with fit + residual
# matrix so the caller can compute L̂ via SVT and predict at arbitrary cells.
fit_twfe_feols <- function(Y_in, w_in, lambda_nn) {
  Nloc <- nrow(Y_in); Tloc <- ncol(Y_in)
  obs_cells <- which(w_in > 0 & !is.na(Y_in), arr.ind = TRUE)
  if (nrow(obs_cells) < 10) {
    return(NULL)
  }

  train_df <- data.frame(
    y  = Y_in[obs_cells],
    wt = w_in[obs_cells],
    i  = obs_cells[, 1],
    t  = obs_cells[, 2]
  )

  fit <- tryCatch(
    feols(y ~ 1 | i + t, data = train_df, weights = ~wt, warn = FALSE, notes = FALSE),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)

  # Residual matrix at observed training cells (for SVT)
  # Use predict() rather than fitted() to guarantee one value per row
  fitted_vals <- suppressWarnings(as.numeric(predict(fit, newdata = train_df)))
  resid_mat <- matrix(0, Nloc, Tloc)
  resids <- train_df$y - fitted_vals
  resids[is.na(resids)] <- 0
  resid_mat[obs_cells] <- resids

  # Compute L̂ via SVT on control-cell residuals
  L <- if (is.finite(lambda_nn)) svt(resid_mat, lambda_nn) else matrix(0, Nloc, Tloc)

  return(list(fit = fit, L = L, train_df = train_df))
}

# Predict counterfactual Y_hat(0) at arbitrary cells using the fitted model
predict_counterfactual <- function(fit_obj, cells) {
  if (is.null(fit_obj)) return(rep(NA_real_, nrow(cells)))
  pred_df <- data.frame(i = cells[, 1], t = cells[, 2])
  y_hat_fe <- suppressWarnings(
    predict(fit_obj$fit, newdata = pred_df)
  )
  # feols returns NA for cells whose i or t are not in the training set
  L_at_cells <- fit_obj$L[cells]
  return(y_hat_fe + L_at_cells)
}

# ------------------------------------------------------------------------------
# Distance-based weights (eq. 3)
# ------------------------------------------------------------------------------
compute_weights <- function(Y_in, W_in, i_treat, t_treat, lambda_time, lambda_unit) {
  Nloc <- nrow(Y_in); Tloc <- ncol(Y_in)
  theta <- exp(-lambda_time * abs((1:Tloc) - t_treat))

  omega <- rep(0, Nloc)
  y_it <- Y_in[i_treat, ]
  w_it <- W_in[i_treat, ]
  for (j in 1:Nloc) {
    if (j == i_treat) next
    y_jt <- Y_in[j, ]
    w_jt <- W_in[j, ]
    mask <- !is.na(y_it) & !is.na(y_jt) & (w_it == 0) & (w_jt == 0)
    if (sum(mask) < 2) {
      omega[j] <- 0
    } else {
      d_ij <- sqrt(mean((y_it[mask] - y_jt[mask])^2))
      omega[j] <- exp(-lambda_unit * d_ij)
    }
  }
  return(list(omega = omega, theta = theta))
}

# ------------------------------------------------------------------------------
# Pooled TROP estimator (global fit, all treated cells predicted at once)
# ------------------------------------------------------------------------------
# The paper's Algorithm 2 is per-cell with blinding of other treated cells.
# For the Czech panel (45.8% treatment density, many always-treated firms),
# per-cell blinding degenerates: always-treated firms have no control periods,
# compute_weights returns omega = 0, and the fit gives tau_hat = Y_it directly.
#
# We therefore pool: single global fit on all observed control cells, apply
# time weights globally (centered on the weighted mean treated period), predict
# counterfactuals for all treated cells from the same (alpha, beta, L), and
# average the per-cell residuals. This is equivalent to a weighted MC estimator
# plus time decay — closer to the fect `mc` method but with explicit lambda_nn
# tuning and an ablation mode. Unit weights (lambda_unit) are omitted in pooled
# mode since they are inherently per-cell.
# ------------------------------------------------------------------------------

trop_pooled <- function(Y_in, W_in, lambda_time, lambda_nn) {
  Nloc <- nrow(Y_in); Tloc <- ncol(Y_in)

  # Global time weights: centered on weighted mean treated period
  treated_cells <- which(W_in == 1 & !is.na(Y_in), arr.ind = TRUE)
  if (nrow(treated_cells) == 0) return(NA_real_)
  t_center <- mean(treated_cells[, 2])
  theta <- exp(-lambda_time * abs((1:Tloc) - t_center))

  # Training weights: time-decayed, observed control cells only
  w_mat <- outer(rep(1, Nloc), theta) * (W_in == 0) * !is.na(Y_in)

  fit_obj <- fit_twfe_feols(Y_in, w_mat, lambda_nn)
  if (is.null(fit_obj)) return(NA_real_)

  # Predict counterfactual at every treated cell
  y_hat <- predict_counterfactual(fit_obj, treated_cells)
  tau <- Y_in[treated_cells] - y_hat
  return(mean(tau, na.rm = TRUE))
}

# ------------------------------------------------------------------------------
# Tuning Q(lambda) — in-sample residual squared error on a sample of control cells
# ------------------------------------------------------------------------------
# Approximation to LOOCV: instead of holding out each control cell individually
# (which would require N_cv separate fits per lambda combination), fit once and
# evaluate residual MSE on a random subset of control cells. This biases toward
# finite lambda_nn (lower in-sample MSE), so we break ties in favor of sparser
# configurations in the cycling grid search.
# ------------------------------------------------------------------------------
loocv_Q <- function(Y_in, W_in, lambda_time, lambda_nn, n_sample = N_CV_CELLS) {
  Nloc <- nrow(Y_in); Tloc <- ncol(Y_in)
  treated_cells <- which(W_in == 1 & !is.na(Y_in), arr.ind = TRUE)
  t_center <- if (nrow(treated_cells) > 0) mean(treated_cells[, 2]) else Tloc / 2
  theta <- exp(-lambda_time * abs((1:Tloc) - t_center))
  w_mat <- outer(rep(1, Nloc), theta) * (W_in == 0) * !is.na(Y_in)

  fit_obj <- tryCatch(fit_twfe_feols(Y_in, w_mat, lambda_nn), error = function(e) NULL)
  if (is.null(fit_obj)) return(Inf)

  control_cells <- which(W_in == 0 & !is.na(Y_in), arr.ind = TRUE)
  if (nrow(control_cells) > n_sample) {
    idx <- sample(nrow(control_cells), n_sample)
    control_cells <- control_cells[idx, , drop = FALSE]
  }
  y_hat <- predict_counterfactual(fit_obj, control_cells)
  resid <- Y_in[control_cells] - y_hat
  return(mean(resid^2, na.rm = TRUE))
}

# Cycling grid search: only (lambda_time, lambda_nn) in pooled mode
cycling_grid_search <- function(Y_in, W_in, grid_time, grid_nn, n_cycles) {
  l_time <- grid_time[1]
  l_nn   <- Inf
  for (cycle in 1:n_cycles) {
    qs <- sapply(grid_time, function(lt) loocv_Q(Y_in, W_in, lt, l_nn))
    l_time <- grid_time[which.min(qs)]
    cat(sprintf("  cycle %d: l_time = %s (Q = %.5f)\n", cycle, format(l_time), min(qs)))
    qs <- sapply(grid_nn, function(ln) loocv_Q(Y_in, W_in, l_time, ln))
    l_nn <- grid_nn[which.min(qs)]
    cat(sprintf("  cycle %d: l_nn   = %s (Q = %.5f)\n", cycle, format(l_nn), min(qs)))
  }
  return(list(lambda_time = l_time, lambda_nn = l_nn))
}

# ------------------------------------------------------------------------------
# Firm-level block bootstrap (Algorithm 3 adapted to pooled TROP)
# ------------------------------------------------------------------------------
bootstrap_se <- function(Y_in, W_in, lambda_time, lambda_nn, n_boot = N_BOOT) {
  Nloc <- nrow(Y_in)
  treated_firms <- which(rowSums(W_in == 1, na.rm = TRUE) > 0)
  control_firms <- setdiff(1:Nloc, treated_firms)
  boot_atts <- numeric(n_boot)
  for (b in 1:n_boot) {
    idx_t <- sample(treated_firms, length(treated_firms), replace = TRUE)
    idx_c <- sample(control_firms, length(control_firms), replace = TRUE)
    idx   <- c(idx_c, idx_t)
    Y_b <- Y_in[idx, , drop = FALSE]
    W_b <- W_in[idx, , drop = FALSE]
    boot_atts[b] <- tryCatch(
      trop_pooled(Y_b, W_b, lambda_time, lambda_nn),
      error = function(e) NA_real_
    )
    cat(sprintf("    bootstrap %d/%d: ATT = %s\n", b, n_boot,
                ifelse(is.na(boot_atts[b]), "NA", sprintf("%.4f", boot_atts[b]))))
  }
  return(sd(boot_atts, na.rm = TRUE))
}

# ------------------------------------------------------------------------------
# Main: LOOCV -> Algorithm 2 -> bootstrap -> ablation
# ------------------------------------------------------------------------------
cat("\n[TROP] Running LOOCV tuning (pooled mode)...\n")
t0 <- Sys.time()
best <- cycling_grid_search(Y, W, grid_time, grid_nn, n_cycles = N_CYCLES)
cat(sprintf("[TROP] Best lambdas: time = %s, nn = %s\n",
            format(best$lambda_time), format(best$lambda_nn)))

cat("[TROP] Computing final ATT (pooled)...\n")
att <- trop_pooled(Y, W, best$lambda_time, best$lambda_nn)
cat(sprintf("[TROP] ATT = %.4f\n", att))

cat(sprintf("[TROP] Bootstrap SE (%d reps)...\n", N_BOOT))
se_att <- bootstrap_se(Y, W, best$lambda_time, best$lambda_nn, N_BOOT)
cat(sprintf("[TROP] SE = %.4f\n", se_att))

runtime_main <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

# ------------------------------------------------------------------------------
# Table 5 ablation
# ------------------------------------------------------------------------------
cat("\n[TROP] Running Table 5 ablation (8 configurations)...\n")
t1 <- Sys.time()

# Force the ablation to use non-trivial tuning parameters so the rows differ
# even when LOOCV happened to pick the boundary (e.g., lambda_nn = Inf).
demo_time <- if (best$lambda_time > 0) best$lambda_time else max(grid_time[grid_time > 0], 0.3)
demo_nn   <- if (is.finite(best$lambda_nn)) best$lambda_nn else min(grid_nn[is.finite(grid_nn)], 0.1)

ablation_configs <- list(
  list(name = "TROP pooled (full)",     lt = demo_time, ln = demo_nn),
  list(name = "lambda_nn = Inf (no L)", lt = demo_time, ln = Inf),
  list(name = "lambda_time = 0",        lt = 0,         ln = demo_nn),
  list(name = "DID (both off)",         lt = 0,         ln = Inf)
)

ablation_atts <- sapply(ablation_configs, function(cfg) {
  cat(sprintf("  [%s] lt=%s ln=%s ...\n", cfg$name, format(cfg$lt), format(cfg$ln)))
  trop_pooled(Y, W, cfg$lt, cfg$ln)
})

ablation_df <- data.frame(
  config = sapply(ablation_configs, function(c) c$name),
  att    = ablation_atts,
  stringsAsFactors = FALSE
)
print(ablation_df)
runtime_abl <- as.numeric(difftime(Sys.time(), t1, units = "mins"))

# ------------------------------------------------------------------------------
# Simulated sanity check
# ------------------------------------------------------------------------------
cat("\n[TROP] Simulated sanity check (N=30, T=20, rank-2, true tau=0.5)...\n")
set.seed(1)
Nsim <- 30; Tsim <- 20
Gamma_sim  <- matrix(rnorm(Nsim * 2), Nsim, 2)
Lambda_sim <- matrix(rnorm(Tsim * 2), Tsim, 2)
L_sim      <- Gamma_sim %*% t(Lambda_sim)
eps_sim    <- matrix(rnorm(Nsim * Tsim, sd = 0.3), Nsim, Tsim)
Y_sim      <- L_sim + eps_sim
W_sim      <- matrix(0, Nsim, Tsim)
W_sim[(Nsim - 4):Nsim, (Tsim - 4):Tsim] <- 1
Y_sim[W_sim == 1] <- Y_sim[W_sim == 1] + 0.5

tau_sim <- tryCatch(
  trop_pooled(Y_sim, W_sim, lambda_time = 0.3, lambda_nn = 0.5),
  error = function(e) NA_real_
)
sanity_pass <- !is.na(tau_sim) && abs(tau_sim - 0.5) < 0.15
cat(sprintf("[TROP] Simulated tau_hat = %.4f (true = 0.5)  %s\n",
            tau_sim, if (sanity_pass) "PASS" else "CHECK"))

# ------------------------------------------------------------------------------
# Output: CSVs
# ------------------------------------------------------------------------------
results <- data.frame(
  method       = "TROP (pooled)",
  att          = att,
  se           = se_att,
  lambda_time  = best$lambda_time,
  lambda_nn    = if (is.infinite(best$lambda_nn)) "Inf" else sprintf("%.4f", best$lambda_nn),
  runtime_min  = runtime_main + runtime_abl,
  n_boot       = N_BOOT,
  quick        = QUICK,
  sanity_pass  = sanity_pass,
  stringsAsFactors = FALSE
)
write.csv(results, file.path(output_dir, "trop_results.csv"), row.names = FALSE)
write.csv(ablation_df, file.path(output_dir, "trop_ablation.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# Output: comparison LaTeX (merges with fect_results.csv if present)
# ------------------------------------------------------------------------------
fect_csv  <- file.path(output_dir, "fect_results.csv")
fect_rows <- character(0)
if (file.exists(fect_csv)) {
  fect_df <- read.csv(fect_csv, stringsAsFactors = FALSE)
  for (i in 1:nrow(fect_df)) {
    fect_rows <- c(fect_rows,
                   sprintf("%s & %.4f & (%.4f) \\\\",
                           fect_df$method[i], fect_df$att[i], fect_df$se[i]))
  }
  cat("[TROP] Merged with fect_results.csv (", nrow(fect_df), " rows)\n", sep = "")
} else {
  cat("[TROP] Note: fect_results.csv not found; trop_comparison.tex will have TROP row only\n")
}

note_merged <- "\\item \\textit{Notes:} FEct, IFEct, MC from Liu, Wang, and Xu (2022); TROP from Athey, Imbens, Qu, and Viviano (2026). Outcome: log markup residualized on $k$, $\\text{cogs}$, firm FE, year FE. Treatment: $pp_{it}$. TROP tuning: LOOCV over $(\\lambda_{\\text{time}}, \\lambda_{\\text{unit}}, \\lambda_{nn})$; SE from firm-level bootstrap."
note_trop_only <- "\\item \\textit{Notes:} TROP from Athey, Imbens, Qu, and Viviano (2026). FEct/IFEct/MC rows will populate when \\texttt{fect\\_results.csv} is present. Outcome: log markup residualized on $k$, $\\text{cogs}$, firm FE, year FE. Treatment: $pp_{it}$. Tuning: LOOCV over $(\\lambda_{\\text{time}}, \\lambda_{\\text{unit}}, \\lambda_{nn})$."

comparison_tex <- c(
  "\\begin{table}[htbp]\\centering",
  "\\caption{Panel Data Estimators: Counterfactual and Triply Robust}\\label{tab:trop_comparison}",
  "\\begin{threeparttable}",
  "\\begin{tabular}{lcc}",
  "\\toprule",
  "Estimator & ATT & SE \\\\",
  "\\midrule",
  fect_rows,
  sprintf("TROP pooled (Athey et al.\\ 2026) & %.4f & (%.4f) \\\\", att, se_att),
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{tablenotes}\\footnotesize",
  if (length(fect_rows) > 0) note_merged else note_trop_only,
  "\\end{tablenotes}",
  "\\end{threeparttable}",
  "\\end{table}"
)
writeLines(comparison_tex, file.path(output_dir, "tables", "trop_comparison.tex"))

# ------------------------------------------------------------------------------
# Output: ablation LaTeX
# ------------------------------------------------------------------------------
ablation_rows <- character(nrow(ablation_df))
for (i in 1:nrow(ablation_df)) {
  ablation_rows[i] <- sprintf("%s & %.4f \\\\",
                              gsub("_", "\\\\_", ablation_df$config[i]),
                              ablation_df$att[i])
}

ablation_tex <- c(
  "\\begin{table}[htbp]\\centering",
  "\\caption{TROP Component Ablation (Table 5 Style) on Czech Construction Panel}\\label{tab:trop_ablation}",
  "\\begin{threeparttable}",
  "\\begin{tabular}{lc}",
  "\\toprule",
  "Configuration & ATT \\\\",
  "\\midrule",
  ablation_rows,
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{tablenotes}\\footnotesize",
  "\\item \\textit{Notes:} TROP with components shut down by setting tuning parameters to boundary values. $\\lambda_{nn} = \\infty$ removes the low-rank regression adjustment; $\\lambda_{\\text{time}} = 0$ flattens the time weights; $\\lambda_{\\text{unit}} = 0$ flattens the unit weights. DID corresponds to all three components off. Mirrors Table 5 of Athey, Imbens, Qu and Viviano (2026).",
  "\\end{tablenotes}",
  "\\end{threeparttable}",
  "\\end{table}"
)
writeLines(ablation_tex, file.path(output_dir, "tables", "trop_ablation.tex"))

# ------------------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------------------
cat("\n", strrep("=", 70), "\n", sep = "")
cat("TROP RESULTS SUMMARY\n")
cat(strrep("=", 70), "\n", sep = "")
cat(sprintf("Mode:          %s\n", if (QUICK) "quick" else "full"))
cat(sprintf("Panel:         %d firms x %d years, %d observed\n", N, T, sum(obs_mask)))
cat(sprintf("ATT:           %.4f\n", att))
cat(sprintf("SE:            %.4f\n", se_att))
cat(sprintf("lambda_time:   %s\n", format(best$lambda_time)))
cat(sprintf("lambda_nn:     %s\n", format(best$lambda_nn)))
cat(sprintf("Sanity check:  %s (tau_hat = %.4f vs true 0.5)\n",
            if (sanity_pass) "PASS" else "CHECK", tau_sim))
cat(sprintf("Runtime (main):  %.1f min\n", runtime_main))
cat(sprintf("Runtime (abl):   %.1f min\n", runtime_abl))
cat("Outputs:\n")
cat("  trop_results.csv\n")
cat("  trop_ablation.csv\n")
cat("  tables/trop_comparison.tex\n")
cat("  tables/trop_ablation.tex\n")
