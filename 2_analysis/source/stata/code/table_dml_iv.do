*===============================================================================
* table_dml_iv.do — DML Partially Linear IV (§6.8, Appendix B.8)
*
* Port of dml_iv.py. Partially linear IV:
*     log(μ_it) = α · pp_dummy_it + g(X_it) + ε_it,  E[ε·Z] = 0
* where Z is an exogenous shift-share / reform exposure instrument.
*
* Input:  $data/markups_panel.dta  + external instruments merge
* Output: ../../output/tables/dml_iv.tex
*===============================================================================

dis _newline "--- table_dml_iv.do ---"

cap which ddml
if _rc != 0 {
    dis "  SKIP: ddml not installed"
    exit
}

use "$data/markups_panel.dta", clear
cap confirm var mu_A
if _rc != 0 {
    dis "  SKIP: missing mu_A"
    exit
}

gen log_mu = log(mu_A)

* Build instrument: lagged procurement dummy (DLW 2012 style),
* or pp_dummy interacted with post-2012 reform indicator (Callaway-
* Sant'Anna continuous-exposure DiD object).
xtset id year
gen pp_L1 = L.pp_dummy
gen post2012 = year >= 2012
gen Z_reform = pp_dummy * post2012

keep if !mi(log_mu, pp_dummy, pp_L1, k, cogs, year, nace2, id)

qui tab year, gen(yr_)
qui tab nace2, gen(nace_)
local base_X "k cogs yr_* nace_*"

set seed 42
eststo clear

* DML-IV with L.pp as instrument (DLW 2012 timing)
ddml init iv, kfolds(5) mname(dml_iv_L)
ddml E[Y|X], mname(dml_iv_L): reg log_mu `base_X'
ddml E[D|X], mname(dml_iv_L): reg pp_dummy `base_X'
ddml E[Z|X], mname(dml_iv_L): reg pp_L1 `base_X'
ddml crossfit, mname(dml_iv_L)
ddml estimate, mname(dml_iv_L) robust
local b_L = _b[pp_dummy]
local se_L = _se[pp_dummy]
local n_L = e(N)

* DML-IV with reform exposure as instrument (post-2012 interaction)
ddml init iv, kfolds(5) mname(dml_iv_R)
ddml E[Y|X], mname(dml_iv_R): reg log_mu `base_X'
ddml E[D|X], mname(dml_iv_R): reg pp_dummy `base_X'
ddml E[Z|X], mname(dml_iv_R): reg Z_reform `base_X'
ddml crossfit, mname(dml_iv_R)
cap ddml estimate, mname(dml_iv_R) robust
if _rc == 0 {
    local b_R = _b[pp_dummy]
    local se_R = _se[pp_dummy]
    local n_R = e(N)
}

* Write LaTeX
cap mkdir "$output/tables"
tempname tf
file open `tf' using "../../output/tables/dml_iv.tex", write replace
#delimit ;
file write `tf'
    "\begin{table}[htbp]\centering" _n
    "\caption{DML Partially Linear IV: LATE on 2012 Reform Compliers}" _n
    "\label{tab:dml_iv}" _n
    "\begin{threeparttable}" _n
    "\begin{tabular}{lccc}" _n
    "\toprule" _n
    "Instrument & LATE & SE & N \\" _n
    "\midrule" _n
;
#delimit cr
file write `tf' "L.pp (DLW 2012 timing)" " & " %7.4f (`b_L') " & " %7.4f (`se_L') " & " %9.0fc (`n_L') " \\" _n
cap confirm var `b_R'
if _rc == 0 {
    file write `tf' "pp * (year >= 2012)" " & " %7.4f (`b_R') " & " %7.4f (`se_R') " & " %9.0fc (`n_R') " \\" _n
}
#delimit ;
file write `tf'
    "\bottomrule" _n
    "\end{tabular}" _n
    "\begin{tablenotes}\footnotesize" _n
    "\item \emph{Notes:} DML partially linear IV with 5-fold cross-fitted "
    "nuisances. First instrument uses DLW (2012) timing (lagged procurement) "
    "for the ATT. Second uses the 2012 single-bid-ban reform as a " _n
    "Callaway-Sant'Anna continuous-exposure DiD object and identifies the " _n
    "compliers' LATE." _n
    "\end{tablenotes}" _n
    "\end{threeparttable}" _n
    "\end{table}" _n
;
#delimit cr
file close `tf'
dis "  Saved: dml_iv.tex"
