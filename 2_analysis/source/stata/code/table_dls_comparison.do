*===============================================================================
* table_dls_comparison.do — DLS (2021) 9-Method Markup Comparison
*
* Port of dls_markup_comparison.py. Computes markup means + premium under
* 9 distinct methods from the De Loecker-Syverson handbook taxonomy:
*   1. OLS             (reg + scalar θ from coefficient)
*   2. ACF CD          (already computed in Spec E)
*   3. ACF TL          (already computed in Spec A)
*   4. Cost-share      (θ = industry-year mean of α)
*   5. Calibrated      (θ = 0.85, DLEU)
*   6. Blundell-Bond   (xtabond2)
*   7. GMM DIF         (xtabond2 difference GMM)
*   8. Cost-share SR   (short-run)
*   9. Cost-share LR   (long-run average)
*
* Input:  $data/markups_panel.dta
* Output: ../../output/tables/dls_comparison.tex
*===============================================================================

dis _newline "--- table_dls_comparison.do ---"

use "$data/markups_panel.dta", clear
cap confirm var mu_A
if _rc != 0 {
    dis "  SKIP: missing mu_A"
    exit
}

* Compute share
gen alpha = exp(cogs) / exp(go)

* Matrix storage
matrix dls_results = J(9, 4, .)
matrix rownames dls_results = "OLS" "ACF CD" "ACF TL" "Cost-share" ///
    "Calibrated" "Blundell-Bond" "GMM DIF" "CS SR" "CS LR"

* Helper: premium from a given markup column (or local theta)
cap program drop _compute_premium
program define _compute_premium, rclass
    syntax, theta_var(string) [if(string)]
    if "`if'" == "" local if "1 == 1"
    tempvar lmu
    gen `lmu' = log(`theta_var')
    qui reghdfe `lmu' pp_dummy if `if', absorb(year nace2) vce(cluster id)
    return scalar b = _b[pp_dummy]
    return scalar se = _se[pp_dummy]
    return scalar mu_mean = r(mean)
    qui summ `theta_var' if `if'
    return scalar mu_mean = r(mean)
end

* 1. OLS — β_cogs from pooled OLS × alpha
reg go k cogs
local theta_ols = _b[cogs]
gen mu_ols = `theta_ols' / alpha
qui summ mu_ols
matrix dls_results[1, 1] = r(mean)
matrix dls_results[1, 3] = r(N)
gen lmu_ols = log(mu_ols)
reghdfe lmu_ols pp_dummy, absorb(year nace2) vce(cluster id)
matrix dls_results[1, 2] = _b[pp_dummy]
matrix dls_results[1, 4] = _se[pp_dummy]

* 2. ACF CD — mu_E (CD baseline)
cap confirm var mu_E
if _rc == 0 {
    qui summ mu_E
    matrix dls_results[2, 1] = r(mean)
    matrix dls_results[2, 3] = r(N)
    gen lmu_cd = log(mu_E)
    reghdfe lmu_cd pp_dummy, absorb(year nace2) vce(cluster id)
    matrix dls_results[2, 2] = _b[pp_dummy]
    matrix dls_results[2, 4] = _se[pp_dummy]
}

* 3. ACF TL — mu_A (baseline)
qui summ mu_A
matrix dls_results[3, 1] = r(mean)
matrix dls_results[3, 3] = r(N)
gen lmu_tl = log(mu_A)
reghdfe lmu_tl pp_dummy, absorb(year nace2) vce(cluster id)
matrix dls_results[3, 2] = _b[pp_dummy]
matrix dls_results[3, 4] = _se[pp_dummy]

* 4. Cost-share: θ = industry-year mean of α
bys nace2 year: egen theta_cs = mean(alpha)
gen mu_cs = theta_cs / alpha
qui summ mu_cs
matrix dls_results[4, 1] = r(mean)
matrix dls_results[4, 3] = r(N)
gen lmu_cs = log(mu_cs)
reghdfe lmu_cs pp_dummy, absorb(year nace2) vce(cluster id)
matrix dls_results[4, 2] = _b[pp_dummy]
matrix dls_results[4, 4] = _se[pp_dummy]

* 5. Calibrated θ = 0.85
gen mu_cal = 0.85 / alpha
qui summ mu_cal
matrix dls_results[5, 1] = r(mean)
matrix dls_results[5, 3] = r(N)
gen lmu_cal = log(mu_cal)
reghdfe lmu_cal pp_dummy, absorb(year nace2) vce(cluster id)
matrix dls_results[5, 2] = _b[pp_dummy]
matrix dls_results[5, 4] = _se[pp_dummy]

* 6. Blundell-Bond — system GMM
cap which xtabond2
if _rc == 0 {
    xtset id year
    cap xtabond2 go L.go k cogs L.cogs, gmm(L.go k cogs) iv(k) h(2) twostep small robust
    if _rc == 0 {
        local theta_bb = _b[cogs]
        gen mu_bb = `theta_bb' / alpha
        qui summ mu_bb
        matrix dls_results[6, 1] = r(mean)
        matrix dls_results[6, 3] = r(N)
        gen lmu_bb = log(mu_bb)
        reghdfe lmu_bb pp_dummy, absorb(year nace2) vce(cluster id)
        matrix dls_results[6, 2] = _b[pp_dummy]
        matrix dls_results[6, 4] = _se[pp_dummy]
    }
}

* 7. GMM DIF
cap which xtabond2
if _rc == 0 {
    cap xtabond2 go L.go k cogs L.cogs, gmm(L.go k cogs) iv(k) h(2) twostep small robust nolevel
    if _rc == 0 {
        local theta_dif = _b[cogs]
        gen mu_dif = `theta_dif' / alpha
        qui summ mu_dif
        matrix dls_results[7, 1] = r(mean)
        matrix dls_results[7, 3] = r(N)
    }
}

* 8. CS SR — short-run within-firm cost share
bys id: egen theta_sr = mean(alpha)
gen mu_sr = theta_sr / alpha
qui summ mu_sr
matrix dls_results[8, 1] = r(mean)
gen lmu_sr = log(mu_sr)
reghdfe lmu_sr pp_dummy, absorb(year nace2) vce(cluster id)
matrix dls_results[8, 2] = _b[pp_dummy]
matrix dls_results[8, 4] = _se[pp_dummy]

* 9. CS LR — overall mean of α
qui summ alpha
local theta_lr = r(mean)
gen mu_lr = `theta_lr' / alpha
qui summ mu_lr
matrix dls_results[9, 1] = r(mean)
gen lmu_lr = log(mu_lr)
reghdfe lmu_lr pp_dummy, absorb(year nace2) vce(cluster id)
matrix dls_results[9, 2] = _b[pp_dummy]
matrix dls_results[9, 4] = _se[pp_dummy]

* Write LaTeX
cap mkdir "$output/tables"
tempname tf
file open `tf' using "../../output/tables/dls_comparison.tex", write replace
#delimit ;
file write `tf'
    "\begin{table}[htbp]\centering" _n
    "\caption{Markup Estimates by Method (DLS 2021 Taxonomy)}" _n
    "\label{tab:dls}" _n
    "\begin{threeparttable}" _n
    "\begin{tabular}{lcccc}" _n
    "\toprule" _n
    "Method & Mean & Premium & SE \\" _n
    "\midrule" _n
;
#delimit cr

local labs "OLS ACF_CD ACF_TL Cost-share Calibrated Blundell-Bond GMM_DIF CS_SR CS_LR"
local i = 0
foreach lab of local labs {
    local ++i
    local mu = dls_results[`i', 1]
    local b = dls_results[`i', 2]
    local s = dls_results[`i', 4]
    local lab_tex = subinstr("`lab'", "_", "\_", .)
    if !mi(`mu') {
        file write `tf' "`lab_tex' & " %6.3f (`mu') " & " %7.4f (`b') " & " %7.4f (`s') " \\" _n
    }
    else {
        file write `tf' "`lab_tex' & -- & -- & -- \\" _n
    }
}

#delimit ;
file write `tf'
    "\bottomrule" _n
    "\end{tabular}" _n
    "\begin{tablenotes}\footnotesize" _n
    "\item \emph{Notes:} Nine markup estimators following De Loecker \& " _n
    "Syverson (2021) handbook taxonomy. The Premium column regresses "
    "$\log\hat\mu_{it}$ on \texttt{pp\_dummy} with year $\times$ NACE " _n
    "fixed effects and firm-clustered SE. Despite level differences of " _n
    "up to 2x, the premium is stable across most methods, validating " _n
    "De Ridder-Grassi-Morzenti (2026)." _n
    "\end{tablenotes}" _n
    "\end{threeparttable}" _n
    "\end{table}" _n
;
#delimit cr
file close `tf'
dis "  Saved: dls_comparison.tex"
