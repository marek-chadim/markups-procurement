# Markups and Public Procurement: Evidence from Czech Construction

**Author:** Marek Chadim (Pre-Doctoral Fellow, Tobin Center for Economic Policy, Yale)

Reproducible pipeline estimating firm-level markups in Czech construction and measuring the causal effect of public procurement on markups.

## Headline Result

Procurement firms charge **14% higher markups** than otherwise-similar private-sector firms. The premium declined from 30% in 2006 to 10% in 2021, with Act 55/2012 (single-bid ban) and Act 134/2016 (MEAT criteria) DiD estimates of −2.6 and −2.2 log points respectively.

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
| De Loecker, Eeckhout & Unger (2020, QJE) | Aggregate markup trends |
| De Loecker & Syverson (2021, Handbook of IO) | Markov-with-pp specification |
| Ackerberg & De Loecker (2024) | Imperfect competition correction |
| De Ridder, Grassi & Morzenti (2026, Econometrica) | DGM cross-validation |
| Raval (2023, REStud) | Overidentification test |
| Kroft, Luo, Mogstad & Setzler (2025, AER) | Double market power framework |
| Benkard, Miller & Yurukoglu (2026) | Compositional decomposition |
| Hendren & Sprung-Keyser (2020, QJE) | MVPF welfare framework |
| Finkelstein & Hendren (2020, JEP) | Welfare analysis meets causal inference |
| Andrews, Gentzkow & Shapiro (2017, QJE) | Λ sensitivity measure |
| Andrews (2017, Econometrica) | Two-step identification-robust CS |
| Andrews et al. (2025, QJE) | Strong-exclusion misspecification theory |
| Baránek & Titl (2024, JLE) | Favoritism channel benchmark |

## Acknowledgments

Repository structure adapted from [GentzkowLabTemplate](https://github.com/gentzkowlab/GentzkowLabTemplate) (© 2024 Matthew Gentzkow, MIT license). Build infrastructure based on the template's module/source/get_inputs/make.sh pattern.

MAD audit protocol adapted from [tjhavranek/research-audit-duel-protocol](https://github.com/tjhavranek/research-audit-duel-protocol).

## License

Code: MIT (inheriting from GentzkowLabTemplate). Working paper: all rights reserved by the author.
