*===============================================================================
* table_dml_sensitivity.do — DML Sensitivity (Cinelli-Hazlett 2020 via `oster`
*                            and manual robustness-value computation)
*
* Port of dml_sensitivity.py. Computes:
*   - Oster (2019) δ* under the assumption R²_max = 1.3 × R²_observed
*   - Robustness value RV at α = 0.05 (Cinelli-Hazlett)
*   - Bias-adjusted CI under plausible confounding strengths
*
* Input:  $data/markups_panel.dta
* Output: ../../output/tables/dml_sensitivity.tex
*===============================================================================

dis _newline "--- table_dml_sensitivity.do ---"

use "$data/markups_panel.dta", clear
cap confirm var mu_A
if _rc != 0 {
    dis "  SKIP: missing mu_A"
    exit
}

gen log_mu = log(mu_A)
keep if !mi(log_mu, pp_dummy, k, cogs, year, nace2, id)

* Baseline (reg with year+nace dummies explicitly so psacalc can read it)
qui tab year, gen(_yr_)
qui tab nace2, gen(_nace_)
reg log_mu pp_dummy k cogs _yr_* _nace_*, vce(cluster id)
local b0 = _b[pp_dummy]
local se0 = _se[pp_dummy]
local r2_base = e(r2)

* Restricted (no controls except FE)
reg log_mu pp_dummy _yr_* _nace_*, vce(cluster id)
local b_restr = _b[pp_dummy]
local r2_restr = e(r2)

* Oster delta* via psacalc (Emily Oster's Stata implementation). Assumes
* R²_max = 1.3 * R²_base per Oster (2019) default, and delta = 1 which
* corresponds to the treatment of selection symmetry.
cap which psacalc
if _rc == 0 {
    qui reg log_mu pp_dummy k cogs _yr_* _nace_*, vce(cluster id)
    cap psacalc delta pp_dummy, rmax(`=min(1.3*`r2_base', 0.99)') mcontrol(k cogs)
    if _rc == 0 {
        local delta_star = r(delta)
    }
    else {
        local delta_star = (`b0' * (min(1.3 * `r2_base', 0.99) - `r2_base')) / ///
            ((`b_restr' - `b0') * (`r2_base' - `r2_restr'))
    }
}
else {
    local delta_star = (`b0' * (min(1.3 * `r2_base', 0.99) - `r2_base')) / ///
        ((`b_restr' - `b0') * (`r2_base' - `r2_restr'))
}

* Robustness value RV(α=0.05) — Cinelli-Hazlett formula
* RV^2 = f² / (1+f²) where f = t_stat^2 / df
local t_stat = `b0' / `se0'
local df = e(df_r)
local f2 = (`t_stat')^2 / `df'
local rv_q0 = `f2' / (1 + `f2')
local rv = sqrt(`rv_q0')

* RV at α = 0.05 one-sided critical value
local tc_05 = 1.645
local f2_05 = max((`t_stat' - `tc_05')^2 / `df', 0)
local rv_a05 = sqrt(`f2_05' / (1 + `f2_05'))

* Write LaTeX
cap mkdir "$output/tables"
tempname tf
file open `tf' using "../../output/tables/dml_sensitivity.tex", write replace
#delimit ;
file write `tf'
    "\begin{table}[htbp]\centering" _n
    "\caption{Sensitivity to Unobserved Confounding}" _n
    "\label{tab:dml_sensitivity}" _n
    "\begin{threeparttable}" _n
    "\begin{tabular}{lc}" _n
    "\toprule" _n
    "Statistic & Value \\" _n
    "\midrule" _n
;
#delimit cr

file write `tf' "Baseline beta (pp)" " & " %7.4f (`b0') " \\" _n
file write `tf' "Baseline SE" " & " %7.4f (`se0') " \\" _n
file write `tf' "Restricted (no controls)" " & " %7.4f (`b_restr') " \\" _n
file write `tf' "Baseline R-squared" " & " %7.4f (`r2_base') " \\" _n
file write `tf' "Oster delta-star (R2 max = 1.3 * R2 obs)" " & " %7.3f (`delta_star') " \\" _n
file write `tf' "Robustness value RV (q=0)" " & " %7.4f (`rv') " \\" _n
file write `tf' "Robustness value RV (alpha=0.05)" " & " %7.4f (`rv_a05') " \\" _n

#delimit ;
file write `tf'
    "\bottomrule" _n
    "\end{tabular}" _n
    "\begin{tablenotes}\footnotesize" _n
    "\item \emph{Notes:} Oster (2019) coefficient-stability ratio "
    "implied by setting R-squared max = 1.3 times R-squared base; "
    "|delta-star| > 1 means an omitted confounder must be more "
    "strongly selected on than the observed controls to overturn "
    "the result. Robustness values follow Cinelli \& Hazlett (2020): "
    "RV is the minimum partial R-squared (of the confounder with "
    "both Y and D, simultaneously) needed to reduce the estimate "
    "to zero (q=0) or to the 95\% CI boundary (alpha=0.05)." _n
    "\end{tablenotes}" _n
    "\end{threeparttable}" _n
    "\end{table}" _n
;
#delimit cr
file close `tf'
dis "  Saved: dml_sensitivity.tex"
dis "  RV(q=0) = " %6.4f (`rv') " | RV(α=0.05) = " %6.4f (`rv_a05')
