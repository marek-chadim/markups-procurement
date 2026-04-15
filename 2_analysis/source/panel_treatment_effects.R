## panel_treatment_effects.R
## Runs multiple panel treatment-effect estimators on the Czech construction
## panel under the absorbing-treatment reformulation. Complements the existing
## fect_estimation.R (counterfactual estimators) by adding the modern
## staggered-DiD estimators from Callaway & Sant'Anna (2021), Sun & Abraham
## (2021), Borusyak, Jaravel & Spiess (2024), plus a naive TWFE benchmark
## and an absorbing-treatment TWFE for reference.
##
## Absorbing treatment reformulation: a firm's first_treat year is the
## earliest year with pp_dummy = 1, and all subsequent firm-years are treated
## regardless of contemporaneous contract status. Never-treated firms have
## first_treat = 0 (the CS/BJS convention for "never treated" cohort).
##
## Outputs:
##   output/panel_treatment_effects_results.csv
##   output/tables/panel_treatment_effects.tex

suppressPackageStartupMessages({
  library(haven)
  library(dplyr)
  library(fixest)
  library(did)
  library(didimputation)
})

source("lib/lalonde_functions_helpers.R")

input_dir  <- file.path("..", "input")
output_dir <- file.path("..", "output")

dir.create(file.path(output_dir, "tables"), showWarnings = FALSE, recursive = TRUE)

## ---- Load & prepare panel ----
panel <- load_lalonde_panel(input_dir, output_dir)
panel <- panel[order(panel$id, panel$year), ]
panel$id_num <- as.integer(panel$id)

## Compute first treated year per firm (0 = never treated)
firm_first_treat <- panel %>%
  group_by(id_num) %>%
  summarize(first_treat = if (any(D == 1)) min(year[D == 1]) else 0L,
            .groups = "drop")
panel <- merge(panel, firm_first_treat, by = "id_num")

## Absorbing treatment indicator
panel$treat_absorb <- as.integer(panel$first_treat > 0 & panel$year >= panel$first_treat)

n_ever   <- sum(firm_first_treat$first_treat > 0)
n_never  <- sum(firm_first_treat$first_treat == 0)
cat(sprintf("Panel: N = %d firm-years, %d firms (ever-treated = %d, never-treated = %d)\n",
            nrow(panel), nrow(firm_first_treat), n_ever, n_never))
cat(sprintf("First-treat year distribution (top): "))
print(head(sort(table(firm_first_treat$first_treat[firm_first_treat$first_treat > 0]),
                decreasing = TRUE), 10))

## ---- (1) Naive TWFE with pp_dummy (non-absorbing) ----
cat("\n(1) Naive TWFE with pp_dummy (contemporaneous treatment):\n")
m_naive <- feols(log_mu ~ D | id_num + year, data = panel, cluster = "id_num")
naive_est <- as.numeric(coef(m_naive)["D"])
naive_se  <- as.numeric(se(m_naive)["D"])
cat(sprintf("  ATT = %.4f (SE %.4f)\n", naive_est, naive_se))

## ---- (2) Absorbing TWFE (ever-treated after first win) ----
cat("\n(2) Absorbing TWFE (ever-treated after first win):\n")
m_absorb <- feols(log_mu ~ treat_absorb | id_num + year,
                  data = panel, cluster = "id_num")
absorb_est <- as.numeric(coef(m_absorb)["treat_absorb"])
absorb_se  <- as.numeric(se(m_absorb)["treat_absorb"])
cat(sprintf("  ATT = %.4f (SE %.4f)\n", absorb_est, absorb_se))

## ---- (3) Callaway-Sant'Anna group-time ATT ----
cat("\n(3) Callaway-Sant'Anna (did::att_gt) with never-treated control:\n")
cs_out <- att_gt(
  yname              = "log_mu",
  gname              = "first_treat",
  tname              = "year",
  idname             = "id_num",
  xformla            = ~ 1,
  data               = panel,
  control_group      = "nevertreated",
  allow_unbalanced_panel = TRUE,
  panel              = TRUE,
  bstrap             = TRUE,
  cband              = FALSE,
  print_details      = FALSE
)
cs_simple <- aggte(cs_out, type = "simple", na.rm = TRUE)
cs_dyn    <- aggte(cs_out, type = "dynamic", na.rm = TRUE)
cs_est    <- as.numeric(cs_simple$overall.att)
cs_se     <- as.numeric(cs_simple$overall.se)
cat(sprintf("  Simple aggregate ATT = %.4f (SE %.4f)\n", cs_est, cs_se))

dyn_window <- data.frame(
  event_time = cs_dyn$egt,
  att        = cs_dyn$att.egt,
  se         = cs_dyn$se.egt
)
cat("  Dynamic event-time effects (subset):\n")
print(round(dyn_window[dyn_window$event_time >= -3 & dyn_window$event_time <= 5, ], 4),
      row.names = FALSE)

## ---- (4) Sun-Abraham via fixest::sunab() ----
cat("\n(4) Sun-Abraham interaction-weighted event study (fixest::sunab):\n")
## sunab() expects a cohort variable with a sentinel for never-treated.
## fixest recommends using a very large year (10000) for never-treated.
panel$sa_cohort <- ifelse(panel$first_treat > 0, panel$first_treat, 10000L)
m_sa <- feols(log_mu ~ sunab(sa_cohort, year) | id_num + year,
              data = panel, cluster = "id_num")

## summary(m, agg = "att") aggregates the full set of cohort x event-time
## interactions into a single post-treatment ATT row labeled "ATT".
sa_row <- summary(m_sa, agg = "att")$coeftable["ATT", ]
sa_est <- as.numeric(sa_row["Estimate"])
sa_se  <- as.numeric(sa_row["Std. Error"])
cat(sprintf("  Aggregate post-treatment ATT = %.4f (SE %.4f)\n", sa_est, sa_se))

## ---- (5) BJS imputation (didimputation::did_imputation) ----
cat("\n(5) BJS imputation estimator (didimputation):\n")
panel_bjs <- panel[, c("id_num", "year", "log_mu", "first_treat")]
bjs_out <- did_imputation(
  data    = panel_bjs,
  yname   = "log_mu",
  gname   = "first_treat",
  tname   = "year",
  idname  = "id_num",
  horizon = FALSE
)
bjs_est <- as.numeric(bjs_out$estimate[1])
bjs_se  <- as.numeric(bjs_out$std.error[1])
cat(sprintf("  Overall ATT = %.4f (SE %.4f)\n", bjs_est, bjs_se))

## ---- Consolidate results ----
res <- data.frame(
  estimator = c("Naive TWFE (pp_dummy)",
                "Absorbing TWFE (ever-treated)",
                "Callaway-Sant'Anna (2021)",
                "Sun-Abraham (2021)",
                "BJS imputation (2024)"),
  att = c(naive_est, absorb_est, cs_est, sa_est, bjs_est),
  se  = c(naive_se,  absorb_se,  cs_se,  sa_se,  bjs_se),
  stringsAsFactors = FALSE
)
cat("\n--- Panel treatment-effect estimates ---\n")
res_print <- res
res_print$att <- round(res_print$att, 4)
res_print$se  <- round(res_print$se,  4)
print(res_print, row.names = FALSE)

write.csv(res, file.path(output_dir, "panel_treatment_effects_results.csv"),
          row.names = FALSE)
cat("\nSaved: panel_treatment_effects_results.csv\n")

## ---- LaTeX table ----
tex <- c(
  "\\begin{table}[htbp]\\centering",
  "\\caption{Panel Treatment-Effect Estimators: Czech Construction Panel}\\label{tab:panel_treatment}",
  "\\begin{threeparttable}",
  "\\begin{tabular}{lcc}",
  "\\toprule",
  "Estimator & ATT & SE \\\\",
  "\\midrule",
  sprintf("Naive TWFE (pp$_{it}$) & %.4f & (%.4f) \\\\", naive_est, naive_se),
  sprintf("Absorbing TWFE (ever-treated) & %.4f & (%.4f) \\\\", absorb_est, absorb_se),
  sprintf("Callaway and Sant'Anna (2021) & %.4f & (%.4f) \\\\", cs_est, cs_se),
  sprintf("Sun and Abraham (2021) & %.4f & (%.4f) \\\\", sa_est, sa_se),
  sprintf("Borusyak, Jaravel, and Spiess (2024) & %.4f & (%.4f) \\\\", bjs_est, bjs_se),
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{tablenotes}\\footnotesize",
  sprintf("\\item \\textit{Notes:} Panel treatment-effect estimators applied to the full unbalanced Czech construction panel ($N = %d$ firm-years, %d firms, of which %d are ever treated and %d are never treated). Naive TWFE uses the contemporaneous year-specific procurement indicator $pp_{it}$. The other four estimators use the absorbing reformulation: a firm's cohort is the earliest year with $pp_{it} = 1$, and all subsequent firm-years are classified as treated. The Callaway-Sant'Anna overall ATT aggregates group-time ATTs with never-treated firms as the control group; the Sun-Abraham estimator uses \\texttt{fixest::sunab} with interaction-weighted event-time aggregation; the Borusyak-Jaravel-Spiess imputation estimator uses \\texttt{didimputation::did\\_imputation} with firm and year fixed effects. All standard errors are firm-clustered.",
          nrow(panel), nrow(firm_first_treat), n_ever, n_never),
  "\\end{tablenotes}",
  "\\end{threeparttable}",
  "\\end{table}"
)
writeLines(tex, file.path(output_dir, "tables", "panel_treatment_effects.tex"))
cat("Saved: tables/panel_treatment_effects.tex\n")
