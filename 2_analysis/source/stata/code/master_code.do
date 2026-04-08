*===============================================================================
* master_code.do — Run all tables and figures
*
* Analog of DGM master_code.do
*===============================================================================

dis _newline "--- master_code.do ---"

* Tables
do "$code/table_pf.do"
do "$code/table_premium.do"
do "$code/table_moments.do"
do "$code/table_correlations.do"

* Figures (require Stata graph capabilities)
cap do "$code/figure_timeseries.do"
cap do "$code/figure_binscatter.do"
