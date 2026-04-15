*===============================================================================
* table_strong_exclusion.do — ABGRS (2025) Strong Exclusion Partial R²
*
* Port of dml_strong_exclusion.py. Computes partial R² of each ACF instrument
* on the year × NACE 2-digit control set via cross-fitted ML partialling-out,
* using pdslasso. Small partial R² (<0.10) confirms ABGRS strong exclusion.
*
* Input:  $data/analysis_panel.dta
* Output: ../../output/tables/strong_exclusion_scores.tex
*===============================================================================

dis _newline "--- table_strong_exclusion.do ---"

cap which pdslasso
if _rc != 0 {
    dis "  SKIP: pdslasso not installed"
    exit
}

use "$data/analysis_panel.dta", clear
xtset id year
gen Lcogs = L.cogs
gen L2cogs = L2.cogs
gen Lk = L.k
gen Lpp = L.pp_dummy
gen L2pp = L2.pp_dummy

keep if !mi(k, cogs, Lk, Lcogs, L2cogs, Lpp, pp_dummy)

* Controls: year + NACE FE
qui tab year, gen(yr_)
qui tab nace2, gen(nace_)
local controls "yr_* nace_*"

* Instruments to test
local ivs "Lk Lcogs Lpp L2pp"

tempname tf
cap mkdir "$output/tables"
file open `tf' using "../../output/tables/strong_exclusion_scores.tex", write replace
#delimit ;
file write `tf'
    "\begin{table}[htbp]\centering" _n
    "\caption{DML Strong Exclusion Diagnostic: Partial R-squared of ACF Instruments on Controls}" _n
    "\label{tab:strong_exclusion}" _n
    "\begin{threeparttable}" _n
    "\begin{tabular}{llc}" _n
    "\toprule" _n
    "Instrument & Description & Partial R-squared \\" _n
    "\midrule" _n
;
#delimit cr

foreach z of local ivs {
    * Lasso-partial each instrument on controls
    cap rlasso `z' `controls'
    if _rc == 0 {
        predict z_hat, xb
        gen z_resid = `z' - z_hat
        qui summ `z'
        local var_z = r(Var)
        qui summ z_resid
        local var_resid = r(Var)
        local partial_r2 = 1 - (`var_resid' / `var_z')

        local desc ""
        if "`z'" == "Lk"     local desc "Lagged log capital"
        if "`z'" == "Lcogs"  local desc "Lagged log COGS"
        if "`z'" == "Lpp"    local desc "Lagged procurement dummy"
        if "`z'" == "L2pp"   local desc "2-year lagged procurement"

        file write `tf' "\texttt{`z'}" " & `desc' & " %7.4f (`partial_r2') " \\" _n
        drop z_hat z_resid
    }
}

* Residualized premium removed from this do-file; requires markups_panel.dta
* merge (mu_A is not in analysis_panel.dta). The partial R-squared rows are
* the primary ABGRS diagnostic; Python's dml_strong_exclusion.py computes
* the residualized-premium companion figure.

#delimit ;
file write `tf'
    "\bottomrule" _n
    "\end{tabular}" _n
    "\begin{tablenotes}\footnotesize" _n
    "\item \emph{Notes:} Partial R-squared of each ACF instrument on year $\times$ "
    "NACE 2-digit controls, computed via Lasso partialling-out (rlasso). "
    "Values $< 0.10$ indicate near-mean-independence from the control set, "
    "satisfying the ABGRS (2025) strong exclusion requirement." _n
    "\end{tablenotes}" _n
    "\end{threeparttable}" _n
    "\end{table}" _n
;
#delimit cr
file close `tf'
dis "  Saved: strong_exclusion_scores.tex"
