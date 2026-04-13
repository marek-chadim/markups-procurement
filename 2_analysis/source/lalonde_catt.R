## lalonde_catt.R
## Ports lalonde5_catt.R from Imbens & Xu (2025, JEP 39(4)) to the Czech
## construction panel. Estimates firm-year level CATTs via grf::causal_forest
## and exposes heterogeneity along three dimensions: firm size (empl_mid),
## market share (mktshare), NACE-2 subsector (41/42/43). The Imbens-Xu
## package's native catt() function from functions_est.R is reused directly.

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
})

source("theme_markups.R")
source("lalonde_functions_helpers.R")

input_dir    <- file.path("..", "input")
output_dir   <- file.path("..", "output")
lalonde_code <- "/Users/marek/Desktop/io/project/thesis/observables/code"

source(file.path(lalonde_code, "functions_est.R"))

dir.create(file.path(output_dir, "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "tables"),  showWarnings = FALSE, recursive = TRUE)

## ---- Load panel and fit causal forest via Imbens-Xu's catt() ----
panel <- load_lalonde_panel(input_dir, output_dir)
cat(sprintf("Panel: N=%d, treated=%d\n", nrow(panel), sum(panel$D)))

set.seed(1234)
catt_out <- catt(panel, Y = "log_mu", treat = "D", covar = LALONDE_COVARIATES)

att_avg  <- as.numeric(catt_out$att["estimate"])
att_se   <- as.numeric(catt_out$att["std.err"])
catt_vec <- catt_out$catt
cat(sprintf("AIPW ATT = %.4f (SE = %.4f)\n", att_avg, att_se))
cat(sprintf("CATT on treated: N=%d, mean=%.4f, median=%.4f, sd=%.4f\n",
            length(catt_vec), mean(catt_vec), median(catt_vec), sd(catt_vec)))
cat(sprintf("CATT quantiles: q10=%.4f q25=%.4f q75=%.4f q90=%.4f\n",
            quantile(catt_vec, 0.10), quantile(catt_vec, 0.25),
            quantile(catt_vec, 0.75), quantile(catt_vec, 0.90)))

## Attach CATT to treated rows and write CSV of (firm, year, CATT, covariates)
treated <- panel[panel$D == 1, ]
treated$catt <- catt_vec
write.csv(
  treated[, c("id", "year", "nace2", "empl_mid", "mktshare", "log_mu", "catt")],
  file.path(output_dir, "lalonde_catt_results.csv"),
  row.names = FALSE
)
cat("Saved: lalonde_catt_results.csv\n")

## ---- Heterogeneity plot (4 panels) ----
plot_df <- treated
plot_df$nace2_lab <- factor(
  plot_df$nace2,
  levels = c(41, 42, 43),
  labels = c("41: Buildings", "42: Civil engineering", "43: Specialized")
)

## Top 5 and bottom 5 CATT firms — labeled on the scatter panels via ggrepel
## (Healy 2026, Ch 5: the subset-for-labels idiom).
outlier_order <- order(plot_df$catt)
labels_df <- plot_df[c(head(outlier_order, 5), tail(outlier_order, 5)), ]

p1 <- ggplot(plot_df, aes(x = catt)) +
  geom_density(fill = markups_blue, alpha = 0.5, color = NA) +
  geom_vline(xintercept = att_avg, linetype = "dashed", color = markups_pink) +
  labs(x = "CATT", y = "Density",
       title = "Distribution of conditional ATTs",
       subtitle = sprintf("Average (dashed) = %.3f", att_avg))

p2 <- ggplot(plot_df, aes(x = empl_mid, y = catt)) +
  geom_point(alpha = 0.15, color = markups_blue, size = 0.6) +
  geom_smooth(method = "loess", span = 0.8, color = markups_pink, se = TRUE) +
  geom_text_repel(
    data = labels_df,
    aes(label = id),
    size = 2.5,
    max.overlaps = Inf,
    min.segment.length = 0.1,
    color = "grey30",
    segment.color = "grey70"
  ) +
  scale_x_log10() +
  labs(x = "Employment (midpoint, log scale)", y = "CATT",
       title = "CATT by firm size")

p3 <- ggplot(plot_df, aes(x = mktshare, y = catt)) +
  geom_point(alpha = 0.15, color = markups_blue, size = 0.6) +
  geom_smooth(method = "loess", span = 0.8, color = markups_pink, se = TRUE) +
  geom_text_repel(
    data = labels_df,
    aes(label = id),
    size = 2.5,
    max.overlaps = Inf,
    min.segment.length = 0.1,
    color = "grey30",
    segment.color = "grey70"
  ) +
  labs(x = "Market share", y = "CATT",
       title = "CATT by market share")

p4 <- ggplot(plot_df, aes(x = nace2_lab, y = catt)) +
  geom_boxplot(fill = markups_blue, alpha = 0.5, outlier.size = 0.4) +
  geom_hline(yintercept = att_avg, linetype = "dashed", color = markups_pink) +
  labs(x = NULL, y = "CATT",
       title = "CATT by NACE-2 subsector") +
  theme(axis.text.x = element_text(angle = 15, hjust = 1))

combined <- ((p1 | p2) / (p3 | p4)) +
  plot_annotation(
    title      = "Heterogeneous Treatment Effects: Procurement Markup Premium",
    subtitle   = sprintf("Conditional ATT by firm size, market share, and NACE-2 subsector (N = %d treated firm-years)", nrow(plot_df)),
    caption    = "Notes: CATTs from grf::causal_forest (4,000 trees). Panels B and C label the top/bottom 5 firms by CATT. Dashed line = average ATT.",
    tag_levels = "A"
  )

ggsave_markups(
  file.path(output_dir, "figures", "lalonde_catt_heterogeneity.pdf"),
  combined,
  width  = 8,
  height = 6.5
)
cat("Saved: figures/lalonde_catt_heterogeneity.pdf\n")

## ---- Heterogeneity table: CATT means + F-test across NACE2 ----
means_by_nace <- aggregate(catt ~ nace2, data = plot_df, FUN = function(x)
  c(mean = mean(x), sd = sd(x), n = length(x)))
means_by_nace <- do.call(data.frame, means_by_nace)
names(means_by_nace) <- c("nace2", "mean", "sd", "n")
print(means_by_nace)

## F-test: is CATT mean different across NACE2?
fit_anova <- lm(catt ~ factor(nace2), data = plot_df)
aov_tab   <- anova(fit_anova)
f_stat    <- aov_tab$`F value`[1]
p_val     <- aov_tab$`Pr(>F)`[1]
cat(sprintf("F-test CATT mean across NACE2: F=%.3f p=%.4g\n", f_stat, p_val))

## LaTeX summary
tex <- c(
  "\\begin{table}[htbp]\\centering",
  "\\caption{Conditional ATT Heterogeneity (GRF Causal Forest)}\\label{tab:lalonde_catt}",
  "\\begin{threeparttable}",
  "\\begin{tabular}{lccc}",
  "\\toprule",
  "Stratum & Mean CATT & SD & $N$ \\\\",
  "\\midrule",
  sprintf("All treated firm-years & %.4f & %.4f & %d \\\\",
          mean(catt_vec), sd(catt_vec), length(catt_vec)),
  "\\midrule",
  sprintf("NACE 41 (buildings) & %.4f & %.4f & %d \\\\",
          means_by_nace$mean[means_by_nace$nace2 == 41],
          means_by_nace$sd[means_by_nace$nace2 == 41],
          means_by_nace$n[means_by_nace$nace2 == 41]),
  sprintf("NACE 42 (civil engineering) & %.4f & %.4f & %d \\\\",
          means_by_nace$mean[means_by_nace$nace2 == 42],
          means_by_nace$sd[means_by_nace$nace2 == 42],
          means_by_nace$n[means_by_nace$nace2 == 42]),
  sprintf("NACE 43 (specialized) & %.4f & %.4f & %d \\\\",
          means_by_nace$mean[means_by_nace$nace2 == 43],
          means_by_nace$sd[means_by_nace$nace2 == 43],
          means_by_nace$n[means_by_nace$nace2 == 43]),
  "\\midrule",
  sprintf("F-test across NACE2 & \\multicolumn{3}{c}{$F = %.2f$, $p = %.4g$} \\\\", f_stat, p_val),
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{tablenotes}\\footnotesize",
  "\\item \\textit{Notes:} Conditional average treatment effects on the treated estimated via generalized random forest causal forest (4{,}000 trees) following Imbens and Xu (2025, \\textit{JEP} 39(4)). Outcome: $\\log \\mu^A_{it}$. Treatment: $pp_{it}$. Covariates as in Table~\\ref{tab:lalonde}. Each treated firm-year receives a predicted CATT.",
  "\\end{tablenotes}",
  "\\end{threeparttable}",
  "\\end{table}"
)
writeLines(tex, file.path(output_dir, "tables", "lalonde_catt.tex"))
cat("Saved: tables/lalonde_catt.tex\n")
