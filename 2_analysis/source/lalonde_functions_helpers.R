## lalonde_functions_helpers.R
## Shared utilities for the Imbens & Xu (2025, JEP 39(4)) supplementary scripts:
## lalonde_overlap.R, lalonde_catt.R, lalonde_sens.R, lalonde_placebo.R.
##
## Provides: load_lalonde_panel(), compute_propensity(), trim_propensity(),
## and the canonical lalonde_covariates vector.
##
## Source via: source("lalonde_functions_helpers.R") from 2_analysis/source/.

suppressPackageStartupMessages({
  library(haven)
  library(dplyr)
  library(grf)
})

LALONDE_COVARIATES <- c("k", "cogs", "empl_mid", "foreign", "mktshare",
                        "nace41", "nace42")

## Load and clean the firm-year panel used by lalonde_estimation.R.
## Returns a data.frame with log_mu, D (treatment), the 7 covariates, plus
## id/year/nace2 for downstream joins. Row-count must match the baseline
## run (N=7666, treated=3443) for any Phase-3 script to be valid.
load_lalonde_panel <- function(input_dir = file.path("..", "input"),
                               output_dir = file.path("..", "output")) {
  df <- read_dta(file.path(input_dir, "data.dta"))
  mk <- read_dta(file.path(output_dir, "data", "paper_markups.dta"))

  df$year <- as.integer(df$year)
  mk$year <- as.integer(mk$year)

  keep_df <- c("id", "year", "pp_dummy", "k", "cogs", "empl_mid",
               "foreign", "mktshare", "nace2")
  keep_mk <- c("id", "year", "markup_A")

  panel <- merge(df[, keep_df], mk[, keep_mk], by = c("id", "year"))
  panel <- panel[!is.na(panel$markup_A) & panel$markup_A > 0, ]
  panel$log_mu  <- log(panel$markup_A)
  panel$D       <- as.integer(panel$pp_dummy)
  panel$foreign <- as.integer(panel$foreign)

  med_empl <- median(panel$empl_mid, na.rm = TRUE)
  panel$empl_mid[is.na(panel$empl_mid)] <- med_empl

  panel$nace41 <- as.integer(panel$nace2 == 41)
  panel$nace42 <- as.integer(panel$nace2 == 42)
  panel$nace43 <- as.integer(panel$nace2 == 43)

  req <- c("log_mu", "D", LALONDE_COVARIATES)
  panel <- panel[complete.cases(panel[, req]), ]
  as.data.frame(panel)
}

## GRF-based propensity score. Mirrors Imbens & Xu (2025) lalonde2_trim.R:
## probability_forest with 4000 trees and seed 1234, followed by floor at 1e-7
## to avoid log(0) in downstream IPW-style estimators.
compute_propensity <- function(data, covariates = LALONDE_COVARIATES,
                               num_trees = 4000, seed = 1234) {
  pf <- probability_forest(X = data[, covariates, drop = FALSE],
                           Y = as.factor(data$D),
                           seed = seed, num.trees = num_trees)
  ps <- pf$predictions[, 2]
  ps[abs(ps) <= 1e-7] <- 1e-7
  ps
}

## Crump, Hotz, Imbens, Mitnik (2009) symmetric trim ps in [alpha, 1-alpha].
## Default alpha=0.1 is the principled default when both tails are populous,
## which is the Czech case (45% treatment rate) unlike LDW (1% treatment).
## The Imbens & Xu (2025) package uses one-sided thresholds of 0.8-0.9 for
## LaLonde; symmetric is the right translation to a balanced panel.
trim_propensity <- function(data, ps, alpha = 0.1) {
  keep <- ps >= alpha & ps <= (1 - alpha)
  trimmed <- data[keep, , drop = FALSE]
  attr(trimmed, "n_full")    <- nrow(data)
  attr(trimmed, "n_trimmed") <- nrow(trimmed)
  attr(trimmed, "alpha")     <- alpha
  attr(trimmed, "ps_full")   <- ps
  attr(trimmed, "ps_kept")   <- ps[keep]
  trimmed
}

## Imbens & Wooldridge (2007, NBER Summer Institute Lecture 1, p.34):
## "A larger normalized difference ... unambiguously indicates a more severe
## overlap problem." They use this metric rather than the t-statistic because
## the t-statistic confounds imbalance with sample size. The rule of thumb
## from page 34 is that |diff/sd| > 0.25 is substantial.
##
## For each covariate, returns:
##   mean_T, mean_C, sd_T, sd_C,  diff_over_sd = (mean_T - mean_C) / sqrt((var_T + var_C)/2)
normalized_balance <- function(data, covariates = LALONDE_COVARIATES,
                               treatment = "D") {
  treat <- data[[treatment]]
  out <- lapply(covariates, function(v) {
    xt <- data[[v]][treat == 1]
    xc <- data[[v]][treat == 0]
    mt <- mean(xt, na.rm = TRUE)
    mc <- mean(xc, na.rm = TRUE)
    vt <- var(xt, na.rm = TRUE)
    vc <- var(xc, na.rm = TRUE)
    denom <- sqrt((vt + vc) / 2)
    data.frame(
      covariate   = v,
      mean_T      = mt,
      sd_T        = sqrt(vt),
      mean_C      = mc,
      sd_C        = sqrt(vc),
      diff_over_sd = if (is.finite(denom) && denom > 0)
                       (mt - mc) / denom else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}
