# Markups and Public Procurement: Evidence from Czech Construction

**Author:** Marek Chadim (Pre-Doctoral Fellow, Tobin Center for Economic Policy, Yale)
**Predoc supervisors:** Seth Zimmerman and Barbara Biasi

Reproducible pipeline estimating firm-level markups in Czech construction and measuring the causal effect of public procurement on markups.

## Headline Result

Procurement firms charge **14% higher markups** than otherwise-similar private-sector firms. The premium declined from 30% in 2006 to 10% in 2021, with Act 55/2012 (single-bid ban) and Act 134/2016 (MEAT criteria) DiD estimates of −2.6 and −2.2 log points respectively.

The paper implements a **9-layer Andrews (2017, 2025) transparency framework** covering Oster sensitivity, ABGRS strong-exclusion residualization (baseline + Chamberlain-optimal instruments), Borusyak-Hull randomization inference, AGS numerical sensitivity, Andrews two-step identification-robust CIs, a 33-estimate specification curve, and a Double/Debiased Machine Learning layer (Chernozhukov et al. 2018) with DML-PLR, DML sensitivity, DML-CATE, and DML-IV diagnostics.

**Fiscal welfare** (§7): capping procurement contracts at engineer estimates would save ~10 bn CZK/year (3.6%); pricing at marginal cost implies ~70 bn CZK/year (12.3%). The top 10% of contracts account for 72% of the implied fiscal transfer, concentrated in NACE 42 civil engineering.

## Running the Pipeline

```bash
# Full end-to-end reproduction (raw data → paper PDF, ~11 min)
bash run_all.sh
```

Requires: Python 3.10+ (anaconda), R 4.2+ (with `fect`, `haven`, `dplyr`, `grf`, `CBPS`, `hbal`, `DoubleML`), Stata 17+, LaTeX (with `latexmk`).

The pipeline produces:
- `4_paper/output/markups_procurement.pdf` — 67-page working paper (JMP readiness ~92/100 per `4_paper/jmp_readiness_report.md`)
- `3_slides/output/slides.pdf` — beamer presentation

## Repository Structure

```
markups-procurement/
├── 0_raw/                       # Raw data (read-only, symlinked external)
│   ├── magnusweb/               # Czech firm financials (MagnusWeb)
│   ├── datlab/                  # Czech procurement register (Datlab)
│   ├── orbis/                   # Orbis all-industry panel (WRDS/BvD)
│   └── figures/                 # Static figures (external R-produced)
├── 1_data/                      # Data construction
│   └── source/
│       ├── rebuild_data.py      # MagnusWeb → data_rebuilt.dta (9,164 obs)
│       └── build_orbis_panel.py # Orbis → orbis_panel_construction.dta
├── 2_analysis/                  # Estimation + robustness (24 Python + 10 R + Stata replication branch)
│   └── source/
│       ├── acf_estimator.py     # ACF Cobb-Douglas + translog, analytical SEs
│       ├── paper_results.py     # Core estimation (6 specs)
│       ├── panel_treatment.py   # Firm FE + Oster bounds
│       ├── dlw_treatment_eval.py  # DLW entry/exit decomposition
│       ├── bmy_czech_analysis.py  # BMY sample sensitivity + percentiles
│       ├── dls_markup_comparison.py  # DLS 9-method comparison
│       ├── cwdl_robustness.py   # CWDL Markov variations
│       ├── raval_test.py        # Raval overidentification
│       ├── klms_analysis.py     # KLMS double market power (appendix)
│       ├── favoritism_decomposition.py  # Reform DiDs + Titl benchmark
│       ├── aggregate_markup_trends.py   # DLEU aggregation + OP decomposition
│       ├── premium_timeseries.py        # Annual premium with reform lines
│       ├── orbis_acf_estimation.py      # Cross-industry validation
│       ├── summary_stats.py     # 15-variable balance table
│       ├── specification_sensitivity_table.py  # Markov-spec sensitivity
│       ├── strong_exclusion_diagnostic.py   # ABGRS partial R² of instruments
│       ├── acf_strong_exclusion.py          # ABGRS residualized ACF
│       ├── borusyak_hull_randomization.py   # BH RI (1000 permutations)
│       ├── ags_twostep_identification.py    # AGS Λ + Andrews two-step CS
│       ├── specification_curve.py           # 33-estimate aggregation
│       ├── style_markups.py     # Python matplotlib style helper (Paul Tol + STIX)
│       │   # ─── Welfare (§7) ───
│       ├── fiscal_welfare_tenders.py     # Aggregate Kaldor-Hicks benchmarks + time-series figure
│       ├── contract_level_welfare.py     # Pass-through regression + heterogeneity + top-20 incidence
│       │   # ─── DML Suite (§6.8) ───
│       ├── dml_core.py          # Shared cross-fitting infrastructure (Chernozhukov et al 2018)
│       ├── dml_premium.py       # DML partially linear premium (Lasso/RF/GB)
│       ├── dml_sensitivity.py   # Chernozhukov-Cinelli sensitivity (RV = 0.54)
│       ├── dml_cate.py          # DML CATE by NACE and reform era
│       ├── dml_iv.py            # DML partially linear IV (shift-share instrument)
│       ├── dml_strong_exclusion.py  # Orthogonal moments for ABGRS strong exclusion
│       │   # ─── Panel / Counterfactual ───
│       ├── fect_estimation.R    # FEct / IFEct (CV r=1) / MC counterfactuals + NACE-specific ATTs
│       ├── sunab_event_study.R  # Sun-Abraham heterogeneity-robust event study (fixest)
│       ├── grf_heterogeneity.R  # Sverdrup-Petukhova-Wager IJMPR pipeline: AUTOC + QINI + BLP
│       │   # ─── Selection on Observables ───
│       ├── lalonde_estimation.R # 11 selection-on-observables estimators
│       ├── lalonde_catt.R       # GRF causal forest + CATT heterogeneity
│       ├── lalonde_overlap.R    # Propensity score overlap + 11-method CIA batch
│       │   # ─── Style helpers ───
│       ├── theme_markups.R      # R ggplot2 style helper (Paul Tol palette)
│       ├── stata/launcher.do    # Stata replication pipeline (end-to-end, 8 sec)
│       ├── stata/code/graph_markups.do  # Stata graph style helper (Paul Tol macros)
│       └── paper_tables.do      # Stata table formatting
├── 3_slides/                    # Beamer presentation
├── 4_paper/                     # Working paper source
│   ├── source/
│   │   ├── markups_procurement.tex   # 67 pages, ~60 bibitems
│   │   └── backups/             # Prior revision snapshots
│   ├── audit_reply.tex          # MAD audit response document
│   └── mad_audit_report.md      # Full audit findings
├── lib/shell/                   # Runner scripts (from Gentzkow template)
├── run_all.sh                   # Pipeline orchestrator
└── setup.sh                     # Initial configuration
```

## Data

- **MagnusWeb** (Czech firm financials, CZ-NACE F): 1,521 construction firms, 9,164 firm-year observations, 2005–2021
- **Datlab** (Czech procurement register): contract values, bidder identities, tender types
- **Orbis** (Bureau van Dijk, WRDS): 275K firms, 1.3M firm-year observations, all industries, 2006–2023

Raw data is NOT committed to this repository. To reproduce, populate `0_raw/` with:
- `0_raw/magnusweb/` — MagnusWeb export (CSVs) and deflator series
- `0_raw/datlab/` — Datlab tender-level CSV
- `0_raw/orbis/` — Orbis all-industry financials, NACE codes, and identifiers

## Methodology

### Production function estimation
**ACF (Ackerberg, Caves & Frazer 2015)** Cobb-Douglas, with lagged procurement indicator in the Markov productivity transition — the specification required for internal consistency per De Loecker (2013) and De Loecker-Syverson (2021, Handbook of IO). Translog functional form as robustness.

### Causal identification (five independent strategies)
1. **Panel firm FE**: exploits 474 switcher firms
2. **DLW entry/exit**: symmetric ±11.5%
3. **fect counterfactual** (Liu, Wang & Xu 2022): FEct 12.5%, IFEct 11.9%, MC 12.2%
4. **Lalonde 11 estimators** (Imbens-Xu 2024): median 15%, range [14%, 19%]
5. **Policy DiDs**: 2012 and 2016 reforms

### 9-Layer Robustness Stack
1. **Oster (2019)** bounds: δ* = −6.05
2. **ABGRS (2025, QJE)** strong exclusion: partial R² < 0.10 for all ACF instruments; residualized premium 0.132
3. **ABGRS optimal instruments** (Chamberlain 1987 sieve): premium 0.134, SE = 0.426 (2.4× tighter than baseline)
4. **Borusyak-Hull (2021, 2026)** randomization inference: p < 0.001
5. **AGS (2017, QJE)** numerical Λ: lagged COGS is load-bearing moment
6. **Andrews (2017, ECMA)** two-step CS: S-stat 2.6–5× wider than Wald
7. **Specification curve**: 33 estimates, median 0.138, range [0.003, 0.186]
8. **DML partially linear** (Chernozhukov et al. 2018 EJ): premium 0.132–0.135 across Lasso/RF/GB nuisance
9. **DML sensitivity** (Chernozhukov-Cinelli et al. 2019): robustness value RV = 0.54 (confounder would need >54% partial R² in BOTH equations to flip the sign)

### Complementary causal-forest heterogeneity (Appendix B.10–B.11)
- **grf::causal_forest** with firm+year FE via supplied Y.hat/W.hat: ATE = 0.134 (matches OLS 0.132), AUTOC = 0.008 (p = 0.035)
- **BLP**: NACE 42 = −0.061 (p < 10^{-14}), mktshare = +0.560 (p < 0.001)
- **QINI targeting curve** via `maq`: top 25% of firms → 60% of heterogeneity gain
- **Sun-Abraham (2021)** event study via `fixest::sunab()`: heterogeneity-robust ATT 0.074 (lower than TWFE due to early-cohort weighting)
- **fect counterfactual** (Liu-Wang-Xu 2022): CV-selected IFEct r=1, NACE 42 ATT ≈ 0 (null under counterfactual imputation)

## Figure Style (Apr 2026)

All paper figures use the colorblind-safe **Paul Tol bright** palette via three language-specific helpers, following Healy (2026, *Data Visualization: A Practical Introduction*, 2nd ed., Ch 8):

- **R / ggplot2:** `2_analysis/source/theme_markups.R` — `theme_set(theme_minimal(base_size=11))` + `tol_bright` named vector + `ggsave_markups()`
- **Python / matplotlib:** `2_analysis/source/style_markups.py` — `apply_markups_style()` + STIX Two Text serif + Computer Modern math (`mathtext.fontset='cm'`)
- **Stata:** `2_analysis/source/stata/code/graph_markups.do` — `${tol_*_hex}` globals + `${markups_gropts}` plot-region options

Palette: `#4477AA` blue, `#EE6677` red, `#228833` green, `#CCBB44` yellow, `#66CCEE` cyan, `#AA3377` purple, `#BBBBBB` grey.

All three helpers define the identical 7-color palette so R, Python, and Stata figures render with matching colors in the final paper. Source the appropriate helper at the top of any new figure-producing script for automatic cross-language consistency. Full study notes: `references/notes/healy_socviz_2026_study.md`.

## Key References

| Reference | Role in paper |
|-|-|
| Ackerberg, Caves & Frazer (2015, Econometrica) | Core production function identification |
| De Loecker & Warzynski (2012, AER) | Markup formula |
| De Loecker & Syverson (2021, Handbook of IO) | Markov-with-pp specification |
| Gandhi, Navarro & Rivers (2020, JPE) | Non-identification of gross-output PF |
| Ackerberg & De Loecker (2024) | Imperfect competition correction |
| De Ridder, Grassi & Morzenti (2026, Econometrica) | DGM cross-validation |
| Kroft, Luo, Mogstad & Setzler (2025, AER) | Double market power framework |
| Baránek & Titl (2024, JLE) | Favoritism channel benchmark |
| Benkard, Miller & Yurukoglu (2026) | Compositional decomposition |
| De Loecker, Eeckhout & Unger (2020, QJE) | Aggregate markup trends |
| Oster (2019, JBES) | Selection bounds |
| Andrews, Gentzkow & Shapiro (2017, QJE) | Λ sensitivity measure |
| Andrews (2017, Econometrica) | Two-step identification-robust CS |
| Andrews et al. (2025, QJE) | Strong-exclusion misspecification theory |
| Borusyak, Hull & Jaravel (2025, ECTJ) | Formula instruments review |
| Raval (2023, REStud) | Overidentification test |
| Chernozhukov et al. (2018, ECTJ) | DML foundation — partially linear model |
| Chernozhukov et al. (2019, NBER) | DML sensitivity analysis — partial R² bounds |
| Wager & Athey (2018, JASA) | Causal forest — heterogeneous treatment effects |
| Sun & Abraham (2021, JoE) | Heterogeneity-robust event study |
| Liu, Wang & Xu (2022, AJPS) | `fect` counterfactual estimators |
| Hendren & Sprung-Keyser (2020, QJE) | MVPF welfare framework |
| Finkelstein & Hendren (2020, JEP) | Welfare analysis meets causal inference |
| Bergé, Butts & McDermott (2026) | `fixest` fast fixed-effects framework |

## Acknowledgments

Repository structure adapted from [GentzkowLabTemplate](https://github.com/gentzkowlab/GentzkowLabTemplate) (© 2024 Matthew Gentzkow, MIT license). Build infrastructure based on the template's module/source/get_inputs/make.sh pattern.

MAD audit protocol adapted from [tjhavranek/research-audit-duel-protocol](https://github.com/tjhavranek/research-audit-duel-protocol).

## License

Code: MIT (inheriting from GentzkowLabTemplate). Working paper: all rights reserved by the author.
