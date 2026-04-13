## panel_treatment_weidner.R
## Cross-implementation replication of the interactive-fixed-effects estimator
## from Moon & Weidner (2015, Econometrica 83(4)) and the weak-factor-robust
## estimator of Armstrong, Weidner & Zeleneev (2024, IFS Working Paper CWP28/24;
## JPE forthcoming). Both are provided by the PanelIFE R package of Hsiang,
## Weidner, and coauthors.
##
## The existing fect_estimation.R already produces an Interactive-FE result via
## the Liu, Wang, and Xu (2022) fect package on the full unbalanced panel. This
## script provides a parallel canonical implementation on the same balanced
## sub-panel used by the SDID and Bacon decomposition scripts, so the cross-
## implementation comparison is apples-to-apples.
##
## Outputs:
##   ../output/panel_weidner_results.csv
##   ../output/tables/panel_weidner.tex

suppressPackageStartupMessages({
  library(haven)
  library(dplyr)
  library(PanelIFE)
  library(fect)
})

input_dir  <- file.path("..", "input")
output_dir <- file.path("..", "output")
dir.create(file.path(output_dir, "tables"), showWarnings = FALSE, recursive = TRUE)

## ---- Load panel (matches fect_estimation.R / lalonde_functions_helpers.R) ----
df <- read_dta(file.path(input_dir, "data.dta"))
mk <- read_dta(file.path(output_dir, "data", "paper_markups.dta"))
df$year <- as.integer(df$year)
mk$year <- as.integer(mk$year)
panel <- merge(
  df[, c("id", "year", "pp_dummy", "k", "cogs")],
  mk[, c("id", "year", "markup_A")],
  by = c("id", "year")
)
panel <- panel[!is.na(panel$markup_A) & panel$markup_A > 0, ]
panel <- panel[complete.cases(panel[, c("pp_dummy", "k", "cogs")]), ]
panel$log_mu <- log(panel$markup_A)
panel$D      <- as.integer(panel$pp_dummy)
panel$id     <- as.integer(panel$id)
panel$year   <- as.integer(panel$year)

cat(sprintf("Full panel: N = %d firm-years, %d firms\n",
            nrow(panel), length(unique(panel$id))))

## ---- Build balanced sub-panel ----
## Keep firms observed in every year (same construction as the Bacon script).
n_years_per_firm <- panel %>% count(id)
max_years <- max(n_years_per_firm$n)
balanced_ids <- n_years_per_firm$id[n_years_per_firm$n == max_years]
panel_bal <- panel[panel$id %in% balanced_ids, ]
panel_bal <- panel_bal[order(panel_bal$id, panel_bal$year), ]

n_firms <- length(unique(panel_bal$id))
n_years <- length(unique(panel_bal$year))
cat(sprintf("Balanced sub-panel: N = %d firm-years, %d firms x %d years\n",
            nrow(panel_bal), n_firms, n_years))

## ---- Reshape to Y (N x T) and X (K x N x T) ----
firms <- sort(unique(panel_bal$id))
years <- sort(unique(panel_bal$year))

mat_Y <- matrix(NA_real_, nrow = n_firms, ncol = n_years,
                dimnames = list(as.character(firms), as.character(years)))
mat_D <- matrix(NA_real_, nrow = n_firms, ncol = n_years,
                dimnames = list(as.character(firms), as.character(years)))
mat_K <- matrix(NA_real_, nrow = n_firms, ncol = n_years)
mat_C <- matrix(NA_real_, nrow = n_firms, ncol = n_years)

for (i in seq_along(firms)) {
  sub_i <- panel_bal[panel_bal$id == firms[i], ]
  for (t in seq_along(years)) {
    row_it <- sub_i[sub_i$year == years[t], ]
    if (nrow(row_it) == 1) {
      mat_Y[i, t] <- row_it$log_mu
      mat_D[i, t] <- row_it$D
      mat_K[i, t] <- row_it$k
      mat_C[i, t] <- row_it$cogs
    }
  }
}

stopifnot(all(!is.na(mat_Y)))  # balanced panel should have no NAs

X_arr <- array(NA_real_, dim = c(3, n_firms, n_years))
X_arr[1, , ] <- mat_D  # treatment indicator
X_arr[2, , ] <- mat_K  # capital
X_arr[3, , ] <- mat_C  # variable input

cat(sprintf("Y: %d x %d; X: 3 x %d x %d\n",
            nrow(mat_Y), ncol(mat_Y), n_firms, n_years))

## ---- Run Moon-Weidner (2015) LS estimator ----
cat("\n=== Moon-Weidner (2015) LS estimator, R = 2 ===\n")
set.seed(42)
ls_out <- ls_factor(
  Y = mat_Y, X = X_arr, R = 2,
  method = "m1", repMIN = 10, repMAX = 30,
  start = rep(0, 3), report = "silent",
  M1 = 2, M2 = 2
)
cat(sprintf("  exitflag: %d (1 = success)\n", ls_out$exitflag))
cat(sprintf("  beta (D, k, cogs): %s\n",
            paste(sprintf("%.4f", ls_out$beta), collapse = ", ")))
## The three variance estimates correspond to different assumptions on the
## errors: Vbeta1 = no serial correlation, Vbeta2 = no cross-section correlation,
## Vbeta3 = general case (see Christensen Lecture 15, section 15.3).
## For a 3-regressor model each Vbeta is a 3x3 matrix; extract the SE on D.
ls_beta_D <- as.numeric(ls_out$beta[1])
ls_se1    <- sqrt(as.numeric(ls_out$Vbeta1[1, 1]))
ls_se2    <- sqrt(as.numeric(ls_out$Vbeta2[1, 1]))
ls_se3    <- sqrt(as.numeric(ls_out$Vbeta3[1, 1]))
## Bias corrections (the sqrt(T/N)*B and sqrt(N/T)*C components). Apply the
## general-case correction to the D coefficient.
ls_bcorr_D <- as.numeric(ls_out$bcorr3[1])
ls_beta_D_bc <- ls_beta_D - ls_bcorr_D
cat(sprintf("  beta_D = %.4f, bias correction = %+.4f, beta_D_bc = %.4f\n",
            ls_beta_D, ls_bcorr_D, ls_beta_D_bc))
cat(sprintf("  SE estimates: V1 = %.4f, V2 = %.4f, V3 = %.4f\n",
            ls_se1, ls_se2, ls_se3))

## ---- Run Armstrong-Weidner-Zeleneev (2024) honest weak-factors estimator ----
cat("\n=== Armstrong-Weidner-Zeleneev (2024) honest weak factors, R = 2 ===\n")
hwf_out <- tryCatch(
  honest_weak_factors(
    Y = mat_Y, X = X_arr, R = 2, alpha = 0.05, clustered_se = TRUE
  ),
  error = function(e) {
    cat("honest_weak_factors failed:", conditionMessage(e), "\n")
    NULL
  }
)

if (!is.null(hwf_out)) {
  hwf_beta <- as.numeric(hwf_out$beta[1])
  hwf_se   <- as.numeric(hwf_out$se[1])
  ## LB and UB are data.frames with columns NumberOfWeakFactors and Est1.
  ## Row 1 (0 weak factors) is the standard CI; subsequent rows widen the CI
  ## assuming progressively more factors are weak. Report the standard CI
  ## (row 1) and the worst-case CI (last row) separately.
  lb_vec <- hwf_out$LB$Est1
  ub_vec <- hwf_out$UB$Est1
  nwf_vec <- hwf_out$LB$NumberOfWeakFactors
  hwf_lb_std   <- lb_vec[1]
  hwf_ub_std   <- ub_vec[1]
  hwf_lb_worst <- lb_vec[length(lb_vec)]
  hwf_ub_worst <- ub_vec[length(ub_vec)]
  cat(sprintf("  beta_D = %.4f, SE = %.4f\n", hwf_beta, hwf_se))
  cat(sprintf("  Standard (no weak factors) 95%% CI: [%.4f, %.4f]\n",
              hwf_lb_std, hwf_ub_std))
  cat(sprintf("  Worst-case (%d weak factors) 95%% CI: [%.4f, %.4f]\n",
              nwf_vec[length(nwf_vec)], hwf_lb_worst, hwf_ub_worst))
  cat("  Full robust CI table:\n")
  full_tab <- data.frame(
    n_weak_factors = nwf_vec,
    lb             = lb_vec,
    ub             = ub_vec
  )
  print(round(full_tab, 4), row.names = FALSE)
  hwf_lb <- hwf_lb_worst
  hwf_ub <- hwf_ub_worst
} else {
  hwf_beta <- NA; hwf_se <- NA; hwf_lb <- NA; hwf_ub <- NA
  hwf_lb_std <- NA; hwf_ub_std <- NA
}

## ---- fect IFEct on the same balanced sub-panel for apples-to-apples check ----
cat("\n=== fect IFEct on the balanced sub-panel (same sample) ===\n")
fect_fit <- tryCatch(
  fect(
    log_mu ~ D + k + cogs,
    data     = panel_bal,
    index    = c("id", "year"),
    force    = "two-way",
    method   = "ife",
    r        = 2,
    se       = TRUE,
    nboots   = 100,
    parallel = TRUE
  ),
  error = function(e) {
    cat("fect failed:", conditionMessage(e), "\n")
    NULL
  }
)
if (!is.null(fect_fit)) {
  att_fect <- as.numeric(fect_fit$att.avg)
  est_tab <- fect_fit$est.avg
  se_fect <- if (is.matrix(est_tab) || is.data.frame(est_tab))
    as.numeric(est_tab[1, "S.E."]) else as.numeric(est_tab["S.E."])
  cat(sprintf("  fect IFEct (balanced sub-panel): ATT = %.4f, SE = %.4f\n",
              att_fect, se_fect))
} else {
  att_fect <- NA; se_fect <- NA
}

## ---- Consolidate results ----
res <- data.frame(
  estimator = c(
    "Moon-Weidner (2015) LS, uncorrected",
    "Moon-Weidner (2015) LS, bias-corrected (V3)",
    "Armstrong-Weidner-Zeleneev (2024) honest weak factors",
    "fect IFEct (same balanced sub-panel)"
  ),
  att = c(ls_beta_D, ls_beta_D_bc, hwf_beta, att_fect),
  se  = c(ls_se3,    ls_se3,       hwf_se,   se_fect),
  ci_lb = c(ls_beta_D - 1.96 * ls_se3,
            ls_beta_D_bc - 1.96 * ls_se3,
            hwf_lb,
            att_fect - 1.96 * se_fect),
  ci_ub = c(ls_beta_D + 1.96 * ls_se3,
            ls_beta_D_bc + 1.96 * ls_se3,
            hwf_ub,
            att_fect + 1.96 * se_fect),
  stringsAsFactors = FALSE
)
cat("\n--- Interactive-FE cross-implementation results ---\n")
res_print <- res
num_cols <- sapply(res_print, is.numeric)
res_print[num_cols] <- lapply(res_print[num_cols], round, 4)
print(res_print, row.names = FALSE)

write.csv(res, file.path(output_dir, "panel_weidner_results.csv"),
          row.names = FALSE)
cat("\nSaved: panel_weidner_results.csv\n")

## ---- LaTeX table ----
fmt_est <- function(e, s) {
  if (is.na(e)) "--- & ---"
  else sprintf("%.4f & (%.4f)", e, s)
}
fmt_ci <- function(lb, ub) {
  if (is.na(lb) || is.na(ub)) "---"
  else sprintf("[%.4f, %.4f]", lb, ub)
}

tex <- c(
  "\\begin{table}[htbp]\\centering",
  "\\caption{Interactive Fixed Effects: Canonical Cross-Implementation Check}\\label{tab:panel_weidner}",
  "\\begin{threeparttable}",
  "\\begin{tabular}{p{6.5cm}ccc}",
  "\\toprule",
  "Estimator & $\\hat{\\beta}_D$ & SE & 95\\% CI \\\\",
  "\\midrule",
  sprintf("Moon--Weidner (2015) LS, uncorrected & %s & %s \\\\",
          fmt_est(res$att[1], res$se[1]),
          fmt_ci(res$ci_lb[1], res$ci_ub[1])),
  sprintf("Moon--Weidner (2015) LS, bias-corrected & %s & %s \\\\",
          fmt_est(res$att[2], res$se[2]),
          fmt_ci(res$ci_lb[2], res$ci_ub[2])),
  sprintf("Armstrong--Weidner--Zeleneev (2024) weak-factor robust & %s & %s \\\\",
          fmt_est(res$att[3], res$se[3]),
          fmt_ci(res$ci_lb[3], res$ci_ub[3])),
  "\\midrule",
  sprintf("fect IFEct (same sub-panel, for comparison) & %s & %s \\\\",
          fmt_est(res$att[4], res$se[4]),
          fmt_ci(res$ci_lb[4], res$ci_ub[4])),
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{tablenotes}\\footnotesize\\raggedright",
  sprintf("\\item \\textit{Notes:} Cross-implementation replication of the interactive fixed-effects estimator on a balanced sub-panel of %d Czech construction firms observed in every panel year (%d firm-years). Model: $\\log \\mu^A_{it} = \\beta_D pp_{it} + \\beta_k k_{it} + \\beta_c \\text{cogs}_{it} + \\lambda'_i F_t + \\varepsilon_{it}$ with $R = 2$ factors. The first three rows use the \\texttt{PanelIFE} R package of Hsiang, Weidner, and coauthors (the canonical author-written implementation of Moon and Weidner 2015 and Armstrong, Weidner, and Zeleneev 2024). The Moon-Weidner bias correction applies the general-case ($B$, $C$) expansion from Bai (2009); the weak-factor-robust estimator gives bias-aware confidence intervals that are valid uniformly in factor strength. The fourth row re-runs the \\texttt{fect} IFEct estimator of Liu, Wang, and Xu (2022) on the same sub-panel for an apples-to-apples comparison.",
          n_firms, nrow(panel_bal)),
  "\\end{tablenotes}",
  "\\end{threeparttable}",
  "\\end{table}"
)
writeLines(tex, file.path(output_dir, "tables", "panel_weidner.tex"))
cat("Saved: tables/panel_weidner.tex\n")
