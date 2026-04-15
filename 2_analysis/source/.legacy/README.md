# .legacy/ — historical record of orphaned scripts

This directory holds R and Python scripts that were part of the
`markups-procurement/2_analysis/source/` pipeline at some earlier
point but are no longer called by `make.sh` or referenced by the
paper. They are kept on disk as historical record; none of them
are executed in the active pipeline.

## Moved on Apr 14 tidy pass

**Python scripts (10):**
- `klms_analysis.py` — KLMS double-market-power analysis. Superseded
  by Stata `table_klms_benchmark.do` which produces `klms_benchmark.tex`
  referenced in Appendix `app:klms`.
- `raval_test.py` — Raval overidentification test. The `app:raval`
  appendix subsection is disabled via `\iffalse` in the paper.
- `premium_timeseries.py` — disabled in make.sh; premium-by-year
  appendix subsection disabled in the paper.
- `specification_sensitivity_table.py` — superseded by Stata
  `table_spec_curve.do`.
- `acf_specification_tests.py` — superseded by Stata
  `table_acf_specification.do` (itself since moved to .legacy).
- `dls_table2_replication.py` — superseded by Stata
  `table_dls_comparison.do`.
- `acf_strong_exclusion.py` — merged into Stata
  `table_strong_exclusion.do`.
- `strong_exclusion_diagnostic.py` — merged into Stata
  `table_strong_exclusion.do`.
- `orbis_acf_estimation.py` — cross-industry Orbis validation
  disabled via `\iffalse` in the paper because Orbis data is not
  available.
- `check_python_counts.py` — one-off diagnostic script.

**R scripts (4):**
- `trop_estimation.R` — Athey-Imbens-Qu-Viviano (2026) triply
  robust panel (TROP) estimator. Disabled in make.sh; not cited
  in the paper.
- `panelview_diagnostics.R` — treatment-pattern visualization; never
  called by the active pipeline.
- `kitagawa_iv_test.R` — Kitagawa (2015) weak-IV test; disabled in
  make.sh.
- `Kitagawa2015_functions.R` — helper functions for the above.

## Why kept on disk

These files may be needed if a referee asks for a specific robustness
check that the paper previously reported, or if a future research
project builds on the same data. They are fully reversible: a single
`mv` brings any file back into the active `source/` directory.

See `/Users/marek/.claude/projects/-Users-marek-Desktop-io/memory/matej_mechanism_pivot.md`
for the full restructuring history.
