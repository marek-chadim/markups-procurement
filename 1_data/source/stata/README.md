# 1_data/source/stata — Stata data-build scripts

This directory holds Stata do-files that read raw sources from `0_raw/`
and write cleaned firm-year panels to `1_data/output/`. Together with the
Python builders in `1_data/source/`, these scripts constitute the
`1_data` module of the Gentzkow Lab template: they produce the inputs
that downstream `2_analysis` consumes.

## Scripts

- **`build_magnusweb_profit.do`** — imports three `0_raw/magnusweb/ratios*.csv`
  exports and writes a firm-year panel of Czech accounting profitability
  ratios (Contribution Margin III, wage share of value added, wage share
  of sales, labour productivity) to `$data/magnusweb_profit.dta`, winsorized
  at the 2nd and 98th percentiles. CM III is the primary accounting-profit
  outcome used in the Section 5.3 reform-mechanism DiD.

## Calling convention

Every script in this directory uses the `$data` global to resolve its
output directory. The caller (either `1_data/make.sh` for top-down
pipeline runs, or `2_analysis/source/stata/launcher.do` for the
Stata-only alternative entry point) sets `$data` before calling.

If the script is invoked directly without a wrapper, it defaults
`$data` to `../../output` relative to this directory, which resolves to
`1_data/output/`.

## Adding a new Stata data-build script

1. Add the `.do` file here with a header comment documenting input,
   output, and caller conventions.
2. Add a `run_stata <file>.do` line in `1_data/make.sh` inside the
   `(cd stata ; ...)` subshell block that sources `run_stata.sh`.
3. Update `2_analysis/source/stata/launcher.do` if the Stata-only
   pipeline needs to call it cross-module.
4. Update `../get_inputs.sh` in `2_analysis/` so downstream analysis
   can symlink the new output as input.

## Why Stata here at all

The Python builders in `1_data/source/` cover most data-cleaning tasks.
Stata is used specifically for the MagnusWeb ratios import because the
Stata `import delimited` with `case(lower)` + `destring` handles the
Czech CSV encoding and the extreme-value percentile winsorization
idiomatically, and because the `_pctile` command is the canonical
Stata percentile computation that matches the behavior downstream
analysis scripts in `2_analysis/source/stata/` rely on.
