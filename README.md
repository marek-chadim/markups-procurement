# Markups and Public Procurement: Evidence from Czech Construction

**Author:** Marek Chadim (Pre-Doctoral Fellow, Tobin Center for Economic Policy, Yale)
**Collaboration:** Chad Syverson (University of Chicago)

Reproducible pipeline estimating firm-level markups in Czech construction and measuring the causal effect of public procurement on markups.

## Headline Result

Procurement firms charge **14% higher markups** than otherwise-similar private-sector firms. The premium declined from 30% in 2006 to 10% in 2021, with Act 55/2012 (single-bid ban) and Act 134/2016 (MEAT criteria) DiD estimates of −2.6 and −2.2 log points respectively.

The paper implements a **6-layer Andrews (2017, 2025) transparency framework** covering Oster sensitivity, ABGRS strong-exclusion residualization, Borusyak-Hull randomization inference, AGS numerical sensitivity, Andrews two-step identification-robust CIs, and a 33-estimate specification curve.

## Running the Pipeline

```bash
# Full end-to-end reproduction (raw data → paper PDF, ~11 min)
bash run_all.sh
```

Requires: Python 3.10+ (anaconda), R 4.2+ (with `fect`, `haven`, `dplyr`, `grf`, `CBPS`, `hbal`, `DoubleML`), Stata 17+, LaTeX (with `latexmk`).

The pipeline produces:
- `4_paper/output/markups_procurement.pdf` — 44-page working paper
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
├── 2_analysis/                  # Estimation + robustness (16 Python + 2 R + 1 Stata)
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
│       ├── fect_estimation.R    # FEct / IFEct / MC counterfactuals
│       ├── lalonde_estimation.R # 11 selection-on-observables estimators
│       └── paper_tables.do      # Stata table formatting
├── 3_slides/                    # Beamer presentation
├── 4_paper/                     # Working paper source
│   ├── source/
│   │   ├── markups_procurement.tex   # 44 pages, ~50 bibitems
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

### 6-Layer Robustness Stack
1. **Oster (2019)** bounds: δ* = −6.05
2. **ABGRS (2025, QJE)** strong exclusion: Δ = 0.0002 log points
3. **Borusyak-Hull (2021, 2026)** randomization inference: p < 0.001
4. **AGS (2017, QJE)** numerical Λ: lagged COGS is load-bearing moment
5. **Andrews (2017, ECMA)** two-step CS: S-stat 2.6–5× wider than Wald
6. **Specification curve**: 33 estimates, median 0.138, range [0.003, 0.186]

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

## Acknowledgments

Repository structure adapted from [GentzkowLabTemplate](https://github.com/gentzkowlab/GentzkowLabTemplate) (© 2024 Matthew Gentzkow, MIT license). Build infrastructure based on the template's module/source/get_inputs/make.sh pattern.

MAD audit protocol adapted from [tjhavranek/research-audit-duel-protocol](https://github.com/tjhavranek/research-audit-duel-protocol).

## License

Code: MIT (inheriting from GentzkowLabTemplate). Working paper: all rights reserved by the author.
