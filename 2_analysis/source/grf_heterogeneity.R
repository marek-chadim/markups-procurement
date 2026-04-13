## grf_heterogeneity.R — Full causal-forest heterogeneity pipeline.
##
## Adapts grf-labs/grf/experiments/ijmpr/analysis.R (Sverdrup, Petukhova, and
## Wager 2024, "Estimating Treatment Effect Heterogeneity in Psychiatry:
## A Review and Tutorial with Causal Forests") to the Czech procurement
## markup setting.
##
## Pipeline:
##   1. ATE via causal_forest AIPW (doubly robust)
##   2. CATE estimation (train/test split) + quartile heterogeneity check
##   3. TOC / AUTOC for heterogeneity evaluation
##   4. QINI curves via maq for policy targeting
##   5. Variable importance + covariate profiling of high/low CATE groups
##   6. Best linear projection (BLP) of CATEs onto NACE × reform
##   7. Risk vs CATE targeting comparison
##
## Outputs:
##   output/figures/grf_cate_quartiles.pdf
##   output/figures/grf_toc.pdf
##   output/figures/grf_qini.pdf
##   output/figures/grf_covariate_profiles.pdf
##   output/tables/grf_heterogeneity_summary.csv
##   output/tables/grf_blp.csv
##
## Reference: Sverdrup, Petukhova, and Wager (2024, IJMPR);
##            Wager and Athey (2018 JASA); Wager (2025 textbook Ch 4-5).

library(grf)
library(maq)
library(ggplot2)
library(haven)

# ---- Paths + style -----------------------------------------------------------
script_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) ".")
if (is.null(script_dir) || script_dir == "") script_dir <- "."
source(file.path(script_dir, "theme_markups.R"))

input_dir  <- file.path(script_dir, "..", "input")
output_dir <- file.path(script_dir, "..", "output")
fig_dir    <- file.path(output_dir, "figures")
tab_dir    <- file.path(output_dir, "tables")
dat_dir    <- file.path(output_dir, "data")
for (d in c(fig_dir, tab_dir, dat_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

set.seed(42)

# ---- Load data ---------------------------------------------------------------
cat("[load] paper_markups.dta + data_rebuilt.dta\n")
mk <- read_dta(file.path(dat_dir, "paper_markups.dta"))
rb <- read_dta(file.path(input_dir, "data_rebuilt.dta"))

# Merge on (id, year) — inner join on non-null markup_A
mk <- mk[!is.na(mk$markup_A), ]
df <- merge(mk[, c("id", "year", "nace2", "markup_A", "omega_A", "pp_dummy", "k", "cogs", "go")],
            rb[, c("id", "year", "mktshare", "foreign", "l", "w", "entry", "exit")],
            by = c("id", "year"))

df$log_mu <- log(df$markup_A)
cat(sprintf("  merged: %d firm-years, %d firms\n", nrow(df), length(unique(df$id))))

# ---- Pre-residualize by firm FE (panel DML approach) -------------------------
# The paper's 14% premium is identified from within-firm variation (firm + year FE).
# Running causal_forest on raw data would identify from cross-sectional variation,
# giving a much lower ATE (0.008). To make the forest comparable to the paper's
# headline, we pre-absorb firm fixed effects by within-firm demeaning of Y and W,
# following Chernozhukov-Newey-Singh (2022) panel DML with firm FE — the same
# strategy as dml_core.py's firm_demean().
#
# Year is kept as a covariate in X so the forest can estimate year-specific CATEs.
# NACE is NOT demeaned (it's time-invariant → collinear with firm FE after demeaning).

# ---- Absorb firm + year FE via explicit nuisance predictions -----------------
# The grf causal_forest accepts explicit Y.hat = E[Y|X] and W.hat = E[W|X].
# Internally, the forest residualizes: Y - Y.hat and W - W.hat, then estimates
# CATEs from the residuals. By passing Y.hat and W.hat from two-way FE
# regressions, we absorb firm + year effects while keeping the ORIGINAL binary
# W — so all grf functions (RATE, QINI, maq) that require binary treatment work.
#
# This is the grf-recommended panel approach; see Athey-Tibshirani-Wager (2019)
# Section 6.1 on "supplying external estimates".

cat("[residualize] Computing firm + year FE predictions for Y.hat and W.hat\n")
Y <- df$log_mu
W <- as.integer(df$pp_dummy)

# Two-way FE fitted values: Y.hat_it = firm_mean_i + year_mean_t - grand_mean
firm_mean_Y <- ave(Y, df$id, FUN = mean)
year_mean_Y <- ave(Y, df$year, FUN = mean)
grand_mean_Y <- mean(Y)
Y.hat.fe <- firm_mean_Y + year_mean_Y - grand_mean_Y

firm_mean_W <- ave(as.numeric(W), df$id, FUN = mean)
year_mean_W <- ave(as.numeric(W), df$year, FUN = mean)
grand_mean_W <- mean(W)
W.hat.fe <- firm_mean_W + year_mean_W - grand_mean_W
# Clip W.hat to avoid 0/1 boundaries for AIPW
W.hat.fe <- pmax(0.01, pmin(0.99, W.hat.fe))

cat(sprintf("  Y residual (Y - Y.hat) sd = %.4f (vs raw sd %.4f)\n", sd(Y - Y.hat.fe), sd(Y)))
cat(sprintf("  W residual (W - W.hat) sd = %.4f (vs raw sd %.4f)\n", sd(W - W.hat.fe), sd(W)))

# Covariates: exogenous controls + year. NOT k/cogs (mechanical link to markup_A).
X <- data.frame(
  mktshare  = df$mktshare,
  foreign   = as.numeric(df$foreign),
  year      = as.numeric(df$year)
)
# Replace NAs with column medians for stability
for (j in seq_len(ncol(X))) {
  na_mask <- is.na(X[, j])
  if (any(na_mask)) X[na_mask, j] <- median(X[!na_mask, j])
}


# ==========================================================================
# 1. ATE via causal_forest (doubly robust AIPW)
# ==========================================================================
cat("\n[1] ATE via causal_forest AIPW\n")
cf.full <- causal_forest(X, Y, W, Y.hat = Y.hat.fe, W.hat = W.hat.fe,
                          num.trees = 4000, seed = 42)
ate <- average_treatment_effect(cf.full, target.sample = "overlap")
cat(sprintf("  ATE = %.4f (SE %.4f)\n", ate[1], ate[2]))
cat(sprintf("  Compare: OLS firm+year FE premium ~ 0.132\n"))


# ==========================================================================
# 2. CATE estimation — train/test split + quartile heterogeneity
# ==========================================================================
cat("\n[2] CATE quartile heterogeneity check\n")
train <- sample(nrow(X), 0.6 * nrow(X))
test  <- setdiff(seq_len(nrow(X)), train)

cate.forest <- causal_forest(X[train, ], Y[train], W[train],
                              Y.hat = Y.hat.fe[train], W.hat = W.hat.fe[train],
                              num.trees = 4000, seed = 42)
tau.hat.test <- predict(cate.forest, X[test, ])$predictions

# Quartile groups
num.groups <- 4
quartile <- cut(tau.hat.test,
                quantile(tau.hat.test, seq(0, 1, by = 1 / num.groups)),
                labels = 1:num.groups, include.lowest = TRUE)
samples.by.quartile <- split(seq_along(quartile), quartile)

# Evaluation forest on test set — clip propensity scores to avoid 0/1 boundary.
# First estimate propensities from a regression forest, then clip and pass explicitly.
eval.forest <- causal_forest(X[test, ], Y[test], W[test],
                              Y.hat = Y.hat.fe[test], W.hat = W.hat.fe[test],
                              num.trees = 4000, seed = 42)

# ATEs by quartile
ate.by.quartile <- lapply(samples.by.quartile, function(s) {
  average_treatment_effect(eval.forest, subset = s, target.sample = "overlap")
})
df.plot.ate <- data.frame(
  estimate = sapply(ate.by.quartile, `[`, 1),
  std.err  = sapply(ate.by.quartile, `[`, 2),
  group    = 1:num.groups
)
cat("  Quartile ATEs:\n")
print(df.plot.ate)

# Plot
p_quartile <- ggplot(df.plot.ate, aes(x = group, y = estimate)) +
  geom_point(color = markups_blue, size = 3) +
  geom_errorbar(aes(ymin = estimate - 1.96 * std.err,
                     ymax = estimate + 1.96 * std.err),
                width = 0.2, color = markups_blue) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  labs(x = "Estimated CATE quartile", y = "Average treatment effect",
       title = "CATE Heterogeneity by Quartile (Causal Forest)") +
  theme_minimal(base_size = 11)
ggsave_markups(file.path(fig_dir, "grf_cate_quartiles.pdf"), plot = p_quartile, width = 6, height = 4)


# ==========================================================================
# 3. TOC / AUTOC for heterogeneity evaluation
# ==========================================================================
cat("\n[3] TOC/AUTOC heterogeneity evaluation\n")
rate.cate <- rank_average_treatment_effect(
  eval.forest, tau.hat.test,
  q = seq(0.05, 1, length.out = 100)
)
cat(sprintf("  AUTOC estimate: %.4f (SE %.4f)\n", rate.cate$estimate, rate.cate$std.err))
p_val_autoc <- 2 * pnorm(-abs(rate.cate$estimate / rate.cate$std.err))
cat(sprintf("  AUTOC p-value: %.4e\n", p_val_autoc))

# Plot TOC
pdf(file.path(fig_dir, "grf_toc.pdf"), width = 7, height = 4.5)
plot(rate.cate, main = "TOC: Targeting Operating Characteristic (Causal Forest)")
dev.off()
cat(sprintf("  [write] %s/grf_toc.pdf\n", basename(fig_dir)))


# ==========================================================================
# 4. QINI curves via maq
# ==========================================================================
cat("\n[4] QINI curves via maq\n")
cost <- 1
scores <- get_scores(eval.forest)
qini <- maq(tau.hat.test, cost, scores, R = 200)
qini.baseline <- maq(tau.hat.test, cost, scores, R = 200,
                      target.with.covariates = FALSE)

# Scale to real-world deployment: how many firms to "audit" (top N by predicted CATE)
max.firms <- length(test)

pdf(file.path(fig_dir, "grf_qini.pdf"), width = 7, height = 4.5)
plot(scale_maq(qini, max.firms),
     ylab = "Expected markup reduction (log points)",
     xlab = "Number of procurement firms targeted",
     main = "QINI Curve: Targeted vs. Uniform Procurement Monitoring")
plot(scale_maq(qini.baseline, max.firms), add = TRUE, ci.args = NULL)
legend("bottomright", c("CATE-targeted", "Uniform"), col = c("black", "black"),
       lty = c(1, 2), bty = "n")
dev.off()
cat(sprintf("  [write] %s/grf_qini.pdf\n", basename(fig_dir)))

# Report the gain at 25%, 50% targeting
for (frac in c(0.25, 0.50, 0.75)) {
  n_target <- round(max.firms * frac)
  gain <- average_gain(scale_maq(qini, max.firms), n_target)
  cat(sprintf("  QINI at %d%% targeted (%d firms): gain = %.4f\n",
              as.integer(frac * 100), n_target, gain[1]))
}


# ==========================================================================
# 5. Variable importance + covariate profiling
# ==========================================================================
cat("\n[5] Variable importance\n")
varimp <- variable_importance(cate.forest)
ranked <- order(varimp, decreasing = TRUE)
n_top <- min(ncol(X), 4)
top.varnames <- colnames(X)[ranked[1:n_top]]
cat("  Top 4 variables:", paste(top.varnames, collapse = ", "), "\n")
cat("  Importances:", paste(round(varimp[ranked[1:n_top]], 4), collapse = ", "), "\n")

# Covariate profiles of low/high CATE groups
low  <- samples.by.quartile[[1]]
high <- samples.by.quartile[[num.groups]]

df.lo <- data.frame(
  value = unlist(as.vector(X[test, top.varnames][low, ])),
  name  = rep(top.varnames, each = length(low)),
  group = "Low CATE (Q1)"
)
df.hi <- data.frame(
  value = unlist(as.vector(X[test, top.varnames][high, ])),
  name  = rep(top.varnames, each = length(high)),
  group = "High CATE (Q4)"
)
df.hist <- rbind(df.lo, df.hi)

p_profiles <- ggplot(df.hist, aes(x = value, fill = group)) +
  geom_histogram(alpha = 0.7, position = "identity", bins = 30) +
  scale_fill_manual(values = c(markups_blue, markups_pink)) +
  facet_wrap(~ name, scales = "free", ncol = 2) +
  labs(x = "Covariate value", fill = NULL,
       title = "Covariate Profiles: Low vs. High Predicted CATE") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")
ggsave_markups(file.path(fig_dir, "grf_covariate_profiles.pdf"), plot = p_profiles,
               width = 8, height = 6)


# ==========================================================================
# 6. Best linear projection (BLP) of CATEs onto NACE × reform
# ==========================================================================
cat("\n[6] Best Linear Projection (BLP)\n")
blp.vars <- data.frame(
  nace42     = as.numeric(df$nace2 == 42),
  nace43     = as.numeric(df$nace2 == 43),
  post2012   = as.numeric(df$year >= 2012),
  mktshare   = df$mktshare
)
# Replace NAs
for (j in seq_len(ncol(blp.vars))) {
  na_mask <- is.na(blp.vars[, j])
  if (any(na_mask)) blp.vars[na_mask, j] <- median(blp.vars[!na_mask, j])
}

blp <- best_linear_projection(cf.full, as.matrix(blp.vars))
cat("  BLP coefficients:\n")
print(blp)

# Save BLP to CSV
blp_df <- data.frame(
  variable  = rownames(blp),
  estimate  = blp[, "Estimate"],
  std.error = blp[, "Std. Error"],
  t.value   = blp[, "t value"],
  p.value   = blp[, "Pr(>|t|)"]
)
write.csv(blp_df, file.path(tab_dir, "grf_blp.csv"), row.names = FALSE)
cat(sprintf("  [write] %s/grf_blp.csv\n", basename(tab_dir)))


# ==========================================================================
# 7. Risk vs CATE targeting comparison
# ==========================================================================
cat("\n[7] Risk vs CATE targeting\n")

# "Risk" model: predict high markup from observables among non-procurement firms
train.ctrl <- train[W[train] == 0]
if (length(train.ctrl) > 50) {
  rf.risk <- regression_forest(X[train.ctrl, ], Y[train.ctrl] - Y.hat.fe[train.ctrl],
                                num.trees = 2000, seed = 42)
  risk.hat.test <- predict(rf.risk, X[test, ])$predictions

  rate.compare <- rank_average_treatment_effect(
    eval.forest,
    cbind(CATE = tau.hat.test, Risk = risk.hat.test)
  )
  cat("  AUTOC comparison (CATE vs Risk):\n")
  print(rate.compare)

  ci <- rate.compare$estimate + data.frame(
    lower = -1.96 * rate.compare$std.err,
    upper =  1.96 * rate.compare$std.err,
    row.names = rate.compare$target
  )
  cat("  95% CIs:\n")
  print(ci)
} else {
  cat("  Skipped: too few control-group training observations\n")
}


# ==========================================================================
# Summary table
# ==========================================================================
cat("\n[summary] Writing grf_heterogeneity_summary.csv\n")
summary_df <- data.frame(
  metric = c("ATE (full sample)", "ATE SE",
             "AUTOC", "AUTOC SE", "AUTOC p-value",
             "Q1 ATE", "Q2 ATE", "Q3 ATE", "Q4 ATE",
             paste0("Top variable ", seq_along(top.varnames))),
  value = c(ate[1], ate[2],
            rate.cate$estimate, rate.cate$std.err, p_val_autoc,
            df.plot.ate$estimate,
            top.varnames)
)
write.csv(summary_df, file.path(tab_dir, "grf_heterogeneity_summary.csv"), row.names = FALSE)
cat(sprintf("  [write] %s/grf_heterogeneity_summary.csv\n", basename(tab_dir)))

cat("\n[done]\n")
