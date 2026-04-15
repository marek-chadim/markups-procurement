*===============================================================================
* master_code.do — Run all tables and figures
*
* Analog of DGM master_code.do
*===============================================================================

dis _newline "--- master_code.do ---"

* ---- Core (existing) ----
do "$code/table_pf.do"
do "$code/table_premium.do"
do "$code/table_moments.do"
do "$code/table_correlations.do"
do "$code/table_dlw_eval.do"
do "$code/table_cwdl_acf.do"
do "$code/table_bh_ri.do"

* ---- Beer 2024 Table C2 (Python-Stata convergence) ----
cap noisily do "$code/table_beer_c2.do"

* ---- Phase 1.1 Extensions ----
cap noisily do "$code/table_strong_exclusion.do"
cap noisily do "$code/table_spec_curve.do"
cap noisily do "$code/table_ags_twostep.do"
cap noisily do "$code/table_adl_instruments.do"
* DEREGISTERED Apr 14 tidy (output raval_overid.tex only referenced in \iffalse):
*   table_raval_overid.do → moved to .legacy/

* ---- Phase 1.2 DML suite ----
cap noisily do "$code/table_dml_premium.do"
cap noisily do "$code/table_dml_iv.do"
cap noisily do "$code/table_dml_cate.do"
cap noisily do "$code/table_dml_sensitivity.do"

* ---- Phase 1.3 Welfare ----
cap noisily do "$code/table_welfare_aggregates.do"
cap noisily do "$code/table_welfare_contract.do"
cap noisily do "$code/table_welfare_incidence.do"

* ---- Phase 1.4 Descriptive + aggregates ----
cap noisily do "$code/table_summary_stats.do"
cap noisily do "$code/table_bmy_decomp.do"
cap noisily do "$code/table_favoritism.do"
cap noisily do "$code/table_klms_benchmark.do"
* DEREGISTERED Apr 14 tidy (aggregate_markup_trends.tex not in paper):
*   table_aggregate_trends.do → moved to .legacy/

* ---- Phase 1.5 Panel IFE (matrix completion + linear AWZ on baseline/reform + nonlinear NNRPanel) ----
* Apr 14 restructuring: replaced the old staggered-DiD battery (CS/SA/BJS/Goodman-Bacon/SDID)
* with PanelIFE (linear IFE, AWZ 2024) + NNRPanel (nonlinear IFE, Zeleneev-Zhang 2026).
* The orphaned Stata files table_panel_didid.do and the legacy R scripts
* panel_treatment.py, panel_treatment_effects.R, panel_treatment_weidner.R, sunab_event_study.R,
* kitagawa_iv_test.R, trop_estimation.R are kept on disk as historical record but not registered.
cap noisily do "$code/table_lalonde_cia.do"

* ---- Phase 2 Supplementary ----
cap noisily do "$code/table_misspec_diagnostics.do"
cap noisily do "$code/table_dls_comparison.do"
cap noisily do "$code/table_cwdl_robustness.do"
* DEREGISTERED Apr 14 tidy (acf_specification.tex not in paper):
*   table_acf_specification.do → moved to .legacy/

* ---- Matej mechanism pivot (Apr 14): joint markup + profitability DiD
* around exogenous Czech procurement reforms
cap noisily do "$code/table_reforms_mechanism.do"

* ---- Figures ----
cap noisily do "$code/figure_timeseries.do"
cap noisily do "$code/figure_binscatter.do"
cap noisily do "$code/figure_event_study_2012.do"
cap noisily do "$code/figure_ri_reform_2012.do"
cap noisily do "$code/figure_bmy_decomposition.do"
