## sunab_event_study.R — Sun-Abraham (2021) heterogeneity-robust event study.
##
## Implements the Sun and Abraham (2021) estimator via fixest::sunab() as a
## robustness check for the TWFE event study in §5.3 of the paper. The SA
## estimator is robust to heterogeneous treatment effects across cohorts
## (firms that entered procurement in different years), which TWFE is not
## (Goodman-Bacon 2021).
##
## Output:
##   output/figures/sunab_event_study.pdf  — SA vs TWFE event-study plot
##   output/tables/sunab_summary.csv       — SA ATT + per-period coefficients
##
## Reference: Sun and Abraham (2021, JoE); Bergé, Butts, McDermott (2026, fixest).

library(fixest)
library(haven)
library(ggplot2)

script_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) ".")
if (is.null(script_dir) || script_dir == "") script_dir <- "."
source(file.path(script_dir, "theme_markups.R"))

input_dir  <- file.path(script_dir, "..", "input")
output_dir <- file.path(script_dir, "..", "output")
fig_dir    <- file.path(output_dir, "figures")
tab_dir    <- file.path(output_dir, "tables")
for (d in c(fig_dir, tab_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

set.seed(42)

# ---- Load and merge ----------------------------------------------------------
cat("[load] paper_markups.dta + data_rebuilt.dta\n")
mk <- read_dta(file.path(output_dir, "data", "paper_markups.dta"))
rb <- read_dta(file.path(input_dir, "data_rebuilt.dta"))
mk <- mk[!is.na(mk$markup_A), ]
df <- merge(mk[, c("id", "year", "markup_A", "pp_dummy", "nace2")],
            rb[, c("id", "year", "pp_entry_year")],
            by = c("id", "year"))
df$log_mu <- log(df$markup_A)
cat(sprintf("  N = %d, firms = %d\n", nrow(df), length(unique(df$id))))

# ---- Sun-Abraham event study -------------------------------------------------
cat("\n[sunab] Sun-Abraham (2021) heterogeneity-robust event study\n")
est_sa <- feols(log_mu ~ sunab(pp_entry_year, year) | id + year,
                data = df, vcov = ~id)

# Aggregate ATT
att_sa <- summary(est_sa, agg = "ATT")
cat(sprintf("  SA ATT = %.4f (SE %.4f, t = %.2f)\n",
            coef(att_sa), se(att_sa), tstat(att_sa)))

# ---- TWFE comparison on the same sample --------------------------------------
cat("\n[twfe] Standard TWFE on same sample\n")
# Create event-time variable
df$event_time <- df$year - df$pp_entry_year
# Standard TWFE with leads/lags (reference period: t = -1)
est_twfe <- feols(log_mu ~ i(event_time, ref = -1) | id + year,
                  data = df, vcov = ~id)
cat(sprintf("  TWFE sample N = %d (same as SA after NA/singleton drops)\n", est_twfe$nobs))

# ---- Event-study figure (SA vs TWFE side by side) ----------------------------
cat("\n[figure] SA vs TWFE event-study comparison\n")

# Extract SA coefficients by event time
sa_coefs <- summary(est_sa, agg = "cohort")
# The SA estimates are aggregated by relative period
sa_df <- data.frame(
  period = as.numeric(gsub(".*::", "", names(coef(est_sa)))),
  estimate = coef(est_sa),
  se = se(est_sa),
  method = "Sun-Abraham (2021)"
)
sa_df <- sa_df[!is.na(sa_df$period), ]

# Extract TWFE coefficients
twfe_names <- names(coef(est_twfe))
twfe_periods <- as.numeric(gsub("event_time::", "", twfe_names))
twfe_df <- data.frame(
  period = twfe_periods,
  estimate = coef(est_twfe),
  se = se(est_twfe),
  method = "TWFE"
)
twfe_df <- twfe_df[!is.na(twfe_df$period), ]

# Combine and trim to reasonable window
both <- rbind(sa_df, twfe_df)
both <- both[both$period >= -5 & both$period <= 8, ]
# Add reference period (t = -1, coefficient = 0)
ref <- data.frame(period = -1, estimate = 0, se = 0,
                   method = c("Sun-Abraham (2021)", "TWFE"))
both <- rbind(both, ref)

p <- ggplot(both, aes(x = period, y = estimate, color = method, shape = method)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
  geom_vline(xintercept = -0.5, linetype = "dashed", color = "grey70") +
  geom_point(position = position_dodge(width = 0.4), size = 2.5) +
  geom_errorbar(aes(ymin = estimate - 1.96 * se, ymax = estimate + 1.96 * se),
                position = position_dodge(width = 0.4), width = 0.3) +
  scale_color_manual(values = c("Sun-Abraham (2021)" = markups_blue,
                                 "TWFE" = markups_pink)) +
  labs(x = "Event time (years relative to first procurement)",
       y = "Coefficient on procurement entry",
       color = NULL, shape = NULL,
       title = "Event Study: Sun-Abraham (2021) vs. TWFE") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

ggsave_markups(file.path(fig_dir, "sunab_event_study.pdf"),
               plot = p, width = 7.5, height = 4.5)

# ---- Firm-specific trends robustness (fixest varying slopes) -----------------
cat("\n[8] Firm-specific linear trends (fixest varying slopes)\n")
est_plain <- fixest::feols(log_mu ~ pp_dummy | id + year, data = df, vcov = ~id)
est_trend <- fixest::feols(log_mu ~ pp_dummy | id[year] + year, data = df, vcov = ~id)
cat(sprintf("  Plain firm+year FE:      %.4f (SE %.4f)\n", coef(est_plain), se(est_plain)))
cat(sprintf("  Firm trends + year FE:   %.4f (SE %.4f)\n", coef(est_trend), se(est_trend)))

# ---- Summary table -----------------------------------------------------------
sa_att_row <- data.frame(
  metric = "SA ATT (aggregated)",
  estimate = coef(att_sa),
  se = se(att_sa),
  stringsAsFactors = FALSE
)

# Also report SA-specific per-period (post only)
sa_post <- sa_df[sa_df$period >= 0, ]
sa_post$metric <- paste0("SA t=", sa_post$period)
sa_rows <- sa_post[, c("metric", "estimate", "se")]

summary_df <- rbind(sa_att_row, sa_rows)
write.csv(summary_df, file.path(tab_dir, "sunab_summary.csv"), row.names = FALSE)
cat(sprintf("[write] %s/sunab_summary.csv\n", basename(tab_dir)))
cat(sprintf("[write] %s/sunab_event_study.pdf\n", basename(fig_dir)))

cat("\n[done]\n")
