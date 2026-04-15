#!/usr/bin/env Rscript
# kitagawa_iv_test.R — Kitagawa (2015) IV Validity Test
#
# Tests whether a BH-recentered reform-exposure instrument satisfies
# LATE IV validity using the Kitagawa (2015, Econometrica) specification test.
#
# Instrument: 2012 single-bid ban exposure
#   Z_it = pre_pp_intensity_i × post2012_t  (shift-share)
#   Recenter: z_tilde = Z_it - E[Z_it | X_i, counterfactual reform years]
#   Binarize: Z_binary = 1(z_tilde > 0)
#
# Test: variance-weighted KS statistic on the density nesting conditions
#   H0: P(B,1) - Q(B,1) >= 0 and Q(B,0) - P(B,0) >= 0 for all B
#   where P = Pr(Y,D | Z=1), Q = Pr(Y,D | Z=0)
#
# References:
#   Kitagawa (2015) "A Test for Instrument Validity" Econometrica 83(5)
#   Borusyak & Hull (2023) "Non-Random Exposure to Exogenous Shocks"
#   Borusyak & Hull (2026) "Optimal Formula Instruments"

library(haven)
library(dplyr)

# Paths
SCRIPT_DIR <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) ".")
INPUT_DIR  <- file.path(SCRIPT_DIR, "..", "input")
OUTPUT_DIR <- file.path(SCRIPT_DIR, "..", "output")

# Source Kitagawa functions
source(file.path(SCRIPT_DIR, "Kitagawa2015_functions.R"))

cat("=== Kitagawa (2015) IV Validity Test ===\n\n")

# ── Load data ──────────────────────────────────────────────────────────
df <- read_dta(file.path(INPUT_DIR, "data.dta"))
cat(sprintf("Loaded: %d obs, %d firms\n", nrow(df), length(unique(df$id))))

# Need: log markup, pp_dummy, year, id, nace2, single_bid_share
df <- df %>%
  filter(!is.na(go), !is.na(k), !is.na(cogs), !is.na(pp_dummy)) %>%
  mutate(
    # Log markup (ACF CD baseline — use expenditure share as proxy if markup_A not available)
    alpha = exp(cogs) / exp(go),
    markup_cs = 1 / alpha,  # cost-share markup (theta=1 assumption)
    lmu = log(markup_cs)
  )

cat(sprintf("After cleaning: %d obs\n", nrow(df)))

# ── Step 1: Construct shift-share instrument ───────────────────────────
# "Share": firm i's pre-reform procurement intensity
# "Shift": post-2012 indicator (single-bid ban)

reform_year <- 2012

# Pre-reform procurement intensity (firm-level average of pp_dummy before 2012)
pre_intensity <- df %>%
  filter(year < reform_year) %>%
  group_by(id) %>%
  summarise(
    pre_pp_intensity = mean(pp_dummy, na.rm = TRUE),
    pre_years = n(),
    .groups = "drop"
  ) %>%
  filter(pre_years >= 2)  # require >=2 pre-reform years for stable estimate

df <- df %>%
  left_join(pre_intensity, by = "id") %>%
  filter(!is.na(pre_pp_intensity))

cat(sprintf("Firms with >=2 pre-2012 obs: %d (%d obs)\n",
            length(unique(df$id)), nrow(df)))

# Shift-share instrument
df$post_reform <- as.numeric(df$year >= reform_year)
df$Z_raw <- df$pre_pp_intensity * df$post_reform

# ── Step 2: BH Recentering ─────────────────────────────────────────────
# Recenter by averaging Z_raw over counterfactual reform years
# Counterfactual: what if the reform hit in year t' instead of 2012?
# For each firm, compute E[Z_it | X_i, counterfactual reform year]

counterfactual_years <- c(2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016)
n_cf <- length(counterfactual_years)

cat(sprintf("Counterfactual reform years: %s\n",
            paste(counterfactual_years, collapse = ", ")))

# For each counterfactual reform year, compute Z^(cf) for each (id, year)
Z_cf_sum <- rep(0, nrow(df))
for (cf_year in counterfactual_years) {
  Z_cf_sum <- Z_cf_sum + df$pre_pp_intensity * as.numeric(df$year >= cf_year)
}
df$Z_expected <- Z_cf_sum / n_cf

# Recentered instrument: actual - expected
df$z_tilde <- df$Z_raw - df$Z_expected

cat(sprintf("z_tilde: mean=%.4f (should be ~0), sd=%.4f\n",
            mean(df$z_tilde), sd(df$z_tilde)))
cat(sprintf("  Positive: %d (%.1f%%), Negative: %d, Zero: %d\n",
            sum(df$z_tilde > 0), 100*mean(df$z_tilde > 0),
            sum(df$z_tilde < 0), sum(df$z_tilde == 0)))

# ── Step 3: Binarize ──────────────────────────────────────────────────
df$Z_binary <- as.numeric(df$z_tilde > 0)

cat(sprintf("\nBinarized instrument: Z=1: %d (%.1f%%), Z=0: %d\n",
            sum(df$Z_binary == 1), 100*mean(df$Z_binary),
            sum(df$Z_binary == 0)))

# First-stage check
fs <- lm(pp_dummy ~ Z_binary, data = df)
cat(sprintf("First stage: coef=%.4f (SE=%.4f), F=%.1f\n",
            coef(fs)["Z_binary"],
            summary(fs)$coefficients["Z_binary", "Std. Error"],
            summary(fs)$fstatistic[1]))

# ── Step 4: Run Kitagawa (2015) test ──────────────────────────────────
cat("\n--- Running Kitagawa (2015) IV validity test ---\n")

Y <- df$lmu
D <- as.integer(df$pp_dummy)
Z <- as.integer(df$Z_binary)

# Drop any NA
valid <- !is.na(Y) & !is.na(D) & !is.na(Z)
Y <- Y[valid]; D <- D[valid]; Z <- Z[valid]

cat(sprintf("Test sample: N=%d, D=1: %d, Z=1: %d\n",
            length(Y), sum(D), sum(Z)))
cat(sprintf("  Pr(D=1|Z=1)=%.3f, Pr(D=1|Z=0)=%.3f\n",
            mean(D[Z==1]), mean(D[Z==0])))

# Run the test with multiple xi values (trimming constants)
xis <- c(0.01, 0.05, 0.1)
set.seed(42)
result <- mergedZtest(Y = Y, D = D, Z = Z,
                      Z_order = c(0, 1),
                      xis = xis,
                      B = 1000,
                      alpha = c(0.10, 0.05, 0.01))

cat("\n=== KITAGAWA TEST RESULTS ===\n\n")
cat(sprintf("  %-8s  %12s  %12s  %12s  %12s\n",
            "xi", "T_N (stat)", "p-value", "cv_10%", "cv_5%"))
cat("  ", rep("-", 60), "\n", sep = "")
for (i in seq_along(xis)) {
  cat(sprintf("  %-8.2f  %12.4f  %12.4f  %12.4f  %12.4f\n",
              xis[i],
              result$teststat[i],
              result$pvals[i],
              result$cvs[1, i],
              result$cvs[2, i]))
}

cat("\n  Interpretation:\n")
cat("  H0: IV validity (exclusion + independence + monotonicity)\n")
if (all(result$pvals > 0.05)) {
  cat("  => CANNOT REJECT H0 at 5% for any xi.\n")
  cat("  => The recentered reform instrument passes the Kitagawa test.\n")
} else {
  cat("  => REJECT H0 at 5% for some xi.\n")
  cat("  => Evidence against IV validity for this instrument.\n")
}

# ── Step 5: Also test with 2016 reform ─────────────────────────────────
cat("\n--- Robustness: 2016 MEAT reform instrument ---\n")

reform_year_2016 <- 2016
cf_years_2016 <- c(2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020)
n_cf2 <- length(cf_years_2016)

df$post_2016 <- as.numeric(df$year >= reform_year_2016)
df$Z_raw_2016 <- df$pre_pp_intensity * df$post_2016

Z_cf_sum2 <- rep(0, nrow(df))
for (cf_year in cf_years_2016) {
  Z_cf_sum2 <- Z_cf_sum2 + df$pre_pp_intensity * as.numeric(df$year >= cf_year)
}
df$Z_expected_2016 <- Z_cf_sum2 / n_cf2
df$z_tilde_2016 <- df$Z_raw_2016 - df$Z_expected_2016
df$Z_binary_2016 <- as.numeric(df$z_tilde_2016 > 0)

Y2 <- df$lmu[valid]; D2 <- D; Z2 <- as.integer(df$Z_binary_2016[valid])

set.seed(42)
result2 <- mergedZtest(Y = Y2, D = D2, Z = Z2,
                       Z_order = c(0, 1),
                       xis = xis,
                       B = 1000,
                       alpha = c(0.10, 0.05, 0.01))

cat(sprintf("  %-8s  %12s  %12s\n", "xi", "T_N", "p-value"))
cat("  ", rep("-", 35), "\n", sep = "")
for (i in seq_along(xis)) {
  cat(sprintf("  %-8.2f  %12.4f  %12.4f\n",
              xis[i], result2$teststat[i], result2$pvals[i]))
}

# ── Save results ──────────────────────────────────────────────────────
results_df <- data.frame(
  reform = rep(c("2012 single-bid ban", "2016 MEAT criteria"), each = length(xis)),
  xi = rep(xis, 2),
  T_N = c(result$teststat, result2$teststat),
  p_value = c(result$pvals, result2$pvals),
  cv_10 = c(result$cvs[1,], result2$cvs[1,]),
  cv_05 = c(result$cvs[2,], result2$cvs[2,])
)

write.csv(results_df,
          file.path(OUTPUT_DIR, "data", "kitagawa_iv_test.csv"),
          row.names = FALSE)

# LaTeX table
tex_file <- file.path(OUTPUT_DIR, "tables", "kitagawa_iv_test.tex")
sink(tex_file)
cat("\\begin{table}[htbp]\\centering\n")
cat("\\caption{Kitagawa (2015) IV Validity Test}\n")
cat("\\label{tab:kitagawa}\n")
cat("\\begin{threeparttable}\n")
cat("\\begin{tabular}{llccc}\n")
cat("\\toprule\n")
cat("Reform instrument & $\\xi$ & $T_N$ & $p$-value & Reject 5\\%? \\\\\n")
cat("\\midrule\n")
for (i in 1:nrow(results_df)) {
  reject <- ifelse(results_df$p_value[i] < 0.05, "Yes", "No")
  cat(sprintf("%s & %.2f & %.3f & %.3f & %s \\\\\n",
              results_df$reform[i],
              results_df$xi[i],
              results_df$T_N[i],
              results_df$p_value[i],
              reject))
}
cat("\\bottomrule\n")
cat("\\end{tabular}\n")
cat("\\begin{tablenotes}\\footnotesize\n")
cat("\\item \\textit{Notes:} Kitagawa (2015) variance-weighted KS test of LATE IV validity. $H_0$: instrument exclusion, independence, and monotonicity jointly hold. Instrument: BH-recentered (Borusyak and Hull 2023) shift-share $\\tilde{z}_{it} = \\text{pre\\_pp}_i \\times \\text{post\\_reform}_t - E[\\cdot | X_i, \\text{counterfactual reform timing}]$, binarized at $\\tilde{z} > 0$. Counterfactual reform years: $\\pm 4$ years around actual reform. Bootstrap $p$-values (1,000 replications). $\\xi$: trimming constant for variance weighting (Kitagawa recommends 0.05--0.10).\n")
cat("\\end{tablenotes}\n")
cat("\\end{threeparttable}\n")
cat("\\end{table}\n")
sink()

cat(sprintf("\nSaved: %s\n", tex_file))
cat(sprintf("Saved: %s\n", file.path(OUTPUT_DIR, "data", "kitagawa_iv_test.csv")))
cat("\n=== Done ===\n")
