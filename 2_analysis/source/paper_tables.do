*===============================================================================
* paper_tables.do
* Generate LaTeX tables and cross-validate with Python
* Reads estimation output from paper_results.do
*===============================================================================

clear all
set more off

local scriptdir = "`c(pwd)'"
global outpath "`scriptdir'/../output"

cap log close
log using "$outpath/paper_tables.log", text replace

*===============================================================================
* 1. COEFFICIENTS TABLE
*===============================================================================

use "$outpath/paper_coefficients.dta", clear

dis _newline "============================================================"
dis "  TABLE 1: Production Function Estimates"
dis "============================================================"

dis _newline "  Spec A (Base: survival + pp in Markov):"
list nace2 b_k_A se_k_A b_cogs_A se_cogs_A N_obs, noobs

dis _newline "  Spec B (No survival):"
list nace2 b_k_B se_k_B b_cogs_B se_cogs_B N_obs, noobs

dis _newline "  Spec C (No pp in Markov):"
list nace2 b_k_C se_k_C b_cogs_C se_cogs_C N_obs, noobs

dis _newline "  Spec D (Plain ACF):"
list nace2 b_k_D se_k_D b_cogs_D se_cogs_D N_obs, noobs

dis _newline "  Spec E (Translog):"
list nace2 b_k_E b_cogs_E b_k2_E b_cogs2_E b_kcogs_E, noobs

dis _newline "  OLS:"
list nace2 b_k_OLS b_cogs_OLS, noobs

* RTS
foreach spec in A B C D OLS {
    gen rts_`spec' = b_k_`spec' + b_cogs_`spec'
}

dis _newline "  Returns to Scale:"
list nace2 rts_A rts_B rts_C rts_D rts_OLS, noobs

*===============================================================================
* 2. COMBINE MARKUPS AND COMPUTE PREMIUMS
*===============================================================================

use "$outpath/temp/paper_markups_41.dta", clear
foreach n in 42 43 {
    append using "$outpath/temp/paper_markups_`n'.dta"
}
save "$outpath/paper_markups.dta", replace

* Treatment definitions
sort id year
by id: gen pp_L1 = pp_dummy[_n-1]
by id: gen pp_L2 = pp_dummy[_n-2]
gen pp_ever_3y = (pp_dummy == 1 | pp_L1 == 1 | pp_L2 == 1)
replace pp_ever_3y = . if mi(pp_dummy)

dis _newline(2) "============================================================"
dis "  TABLE 2: PROCUREMENT PREMIUM"
dis "============================================================"

* Markup distributions
dis _newline "--- Markup Distributions ---"
foreach spec in A B C D E OLS {
    dis _newline "  Spec `spec':"
    tabstat markup_`spec', by(nace2) stat(mean sd p10 p50 p90 N) format(%9.3f)
}

* Log markups
foreach spec in A B C D E OLS {
    gen lmu_`spec' = ln(markup_`spec') if markup_`spec' > 0
}

dis _newline "--- Raw Premium (unconditional mean difference in log markups) ---"
dis "  {Spec}: pp_dummy premium (SE)  |  pp_ever_3y premium (SE)"
dis "  ---------------------------------------------------------------"

foreach spec in A B C D E OLS {
    * pp_dummy
    qui ttest lmu_`spec', by(pp_dummy)
    local raw = r(mu_2) - r(mu_1)
    local se  = r(se)

    * pp_ever_3y
    qui ttest lmu_`spec', by(pp_ever_3y)
    local raw3 = r(mu_2) - r(mu_1)
    local se3  = r(se)

    dis "  {`spec'}: " %7.4f `raw' " (" %5.4f `se' ")    |  " ///
        %7.4f `raw3' " (" %5.4f `se3' ")"
}

dis _newline "--- Regression Premium (+ k, cogs, year*nace2 FE, cluster(id)) ---"
dis "  {Spec}: pp_dummy coef (SE) [R2]  |  pp_ever_3y coef (SE)"
dis "  ---------------------------------------------------------------"

foreach spec in A B C D E OLS {
    * pp_dummy
    xi: qui reg lmu_`spec' pp_dummy k cogs i.year*i.nace2, cluster(id)
    local reg_pp = _b[pp_dummy]
    local se_pp  = _se[pp_dummy]
    local r2_pp  = e(r2)
    local N_pp   = e(N)

    * pp_ever_3y
    xi: qui reg lmu_`spec' pp_ever_3y k cogs i.year*i.nace2, cluster(id)
    local reg_3y = _b[pp_ever_3y]
    local se_3y  = _se[pp_ever_3y]

    dis "  {`spec'}: " %7.4f `reg_pp' " (" %5.4f `se_pp' ") [" %5.3f `r2_pp' "]" ///
        "  |  " %7.4f `reg_3y' " (" %5.4f `se_3y' ")  N=" `N_pp'
}

*===============================================================================
* 3. LaTeX TABLE 1: PF Estimates
*===============================================================================

use "$outpath/paper_coefficients.dta", clear

cap file close tab1
file open tab1 using "$outpath/table_pf_estimates.tex", write replace

#delimit ;
file write tab1
    "\begin{table}[htbp]\centering" _newline
    "\caption{Production Function Estimates by Industry}" _newline
    "\label{tab:pf_estimates}" _newline
    "\begin{tabular}{l*{3}{c}}" _newline
    "\toprule" _newline
    " & NACE 41 & NACE 42 & NACE 43 \\" _newline
    " & (Buildings) & (Civil Eng.) & (Specialized) \\" _newline
    "\midrule" _newline
    "\multicolumn{4}{l}{\textit{Panel A: Baseline (survival + pp in Markov)}} \\" _newline
;
#delimit cr

* Beta_k row
local bk41 = b_k_A[1]
local bk42 = b_k_A[2]
local bk43 = b_k_A[3]
local sk41 = se_k_A[1]
local sk42 = se_k_A[2]
local sk43 = se_k_A[3]
file write tab1 "$\hat{\beta}_k$"
file write tab1 " & " %6.4f (`bk41') " & " %6.4f (`bk42') " & " %6.4f (`bk43') " \\" _newline
file write tab1 " & (" %6.4f (`sk41') ") & (" %6.4f (`sk42') ") & (" %6.4f (`sk43') ") \\" _newline

* Beta_cogs row
local bc41 = b_cogs_A[1]
local bc42 = b_cogs_A[2]
local bc43 = b_cogs_A[3]
local sc41 = se_cogs_A[1]
local sc42 = se_cogs_A[2]
local sc43 = se_cogs_A[3]
file write tab1 "$\hat{\beta}_{\text{cogs}}$"
file write tab1 " & " %6.4f (`bc41') " & " %6.4f (`bc42') " & " %6.4f (`bc43') " \\" _newline
file write tab1 " & (" %6.4f (`sc41') ") & (" %6.4f (`sc42') ") & (" %6.4f (`sc43') ") \\" _newline

* RTS
file write tab1 "RTS"
forvalues i = 1/3 {
    local rts = b_k_A[`i'] + b_cogs_A[`i']
    file write tab1 " & " %5.3f (`rts')
}
file write tab1 " \\" _newline

* N
file write tab1 "$N$"
forvalues i = 1/3 {
    local nn = N_obs[`i']
    file write tab1 " & " %6.0fc (`nn')
}
file write tab1 " \\" _newline

* Panel B header
#delimit ;
file write tab1
    "\midrule" _newline
    "\multicolumn{4}{l}{\textit{Panel B: Robustness --- $\hat{\beta}_{\text{cogs}}$ (SE)}} \\" _newline
;
#delimit cr

* Robustness rows
foreach spec in B C D OLS {
    if "`spec'" == "B" local slabel "No survival"
    if "`spec'" == "C" local slabel "No pp in Markov"
    if "`spec'" == "D" local slabel "Plain ACF"
    if "`spec'" == "OLS" local slabel "OLS"

    file write tab1 "`slabel'"
    forvalues i = 1/3 {
        local b = b_cogs_`spec'[`i']
        file write tab1 " & " %6.4f (`b')
    }
    file write tab1 " \\" _newline

    if "`spec'" != "OLS" {
        file write tab1 " "
        forvalues i = 1/3 {
            local s = se_cogs_`spec'[`i']
            file write tab1 " & (" %6.4f (`s') ")"
        }
        file write tab1 " \\" _newline
    }
}

* Translog row
file write tab1 "Translog"
forvalues i = 1/3 {
    local b = b_cogs_E[`i']
    file write tab1 " & " %6.4f (`b')
}
file write tab1 " \\" _newline
file write tab1 " "
forvalues i = 1/3 {
    local s = se_cogs_E[`i']
    file write tab1 " & (" %6.4f (`s') ")"
}
file write tab1 " \\" _newline

#delimit ;
file write tab1
    "\bottomrule" _newline
    "\multicolumn{4}{p{10cm}}{\footnotesize " _newline
    "Notes: Two-step ACF estimator (Ackerberg et al.\ 2015). "
    "Baseline Markov: $\omega_t = \rho\omega_{t-1} + \gamma pp_{t-1} + \delta\hat{p}_{t-1} + \xi_t$. "
    "Instruments: $(1, k_t, \text{cogs}_{t-1})$. "
    "Analytical SEs (ACH 2012), clustered by firm.} \\" _newline
    "\end{tabular}" _newline
    "\end{table}" _newline
;
#delimit cr

file close tab1
dis _newline "  Saved: table_pf_estimates.tex"

*===============================================================================
* 4. LaTeX TABLE 2: Procurement Premium
*===============================================================================

use "$outpath/paper_markups.dta", clear

sort id year
by id: gen pp_L1 = pp_dummy[_n-1]
by id: gen pp_L2 = pp_dummy[_n-2]
gen pp_ever_3y = (pp_dummy == 1 | pp_L1 == 1 | pp_L2 == 1)
replace pp_ever_3y = . if mi(pp_dummy)

foreach spec in A B C D E OLS {
    gen lmu_`spec' = ln(markup_`spec') if markup_`spec' > 0
}

cap file close tab2
file open tab2 using "$outpath/table_premium.tex", write replace

#delimit ;
file write tab2
    "\begin{table}[htbp]\centering" _newline
    "\caption{Procurement Markup Premium Across Specifications}" _newline
    "\label{tab:premium}" _newline
    "\begin{tabular}{l*{2}{c}}" _newline
    "\toprule" _newline
    " & $pp_t$ & $pp^{3y}_t$ \\" _newline
    "\midrule" _newline
    "\multicolumn{3}{l}{\textit{Raw log markup difference}} \\" _newline
;
#delimit cr

* Raw premiums
foreach spec in A B C D E OLS {
    if "`spec'" == "A" local slabel "Base (surv.\ + pp)"
    if "`spec'" == "B" local slabel "No survival"
    if "`spec'" == "C" local slabel "No pp in Markov"
    if "`spec'" == "D" local slabel "Plain ACF"
    if "`spec'" == "E" local slabel "Translog"
    if "`spec'" == "OLS" local slabel "OLS"

    * pp_dummy
    qui ttest lmu_`spec', by(pp_dummy)
    local raw = r(mu_2) - r(mu_1)
    local se  = r(se)

    * pp_ever_3y
    qui ttest lmu_`spec', by(pp_ever_3y)
    local raw3 = r(mu_2) - r(mu_1)
    local se3  = r(se)

    file write tab2 "`slabel'"
    file write tab2 " & " %6.3f (`raw') " & " %6.3f (`raw3') " \\" _newline
    file write tab2 " & (" %5.3f (`se') ") & (" %5.3f (`se3') ") \\" _newline
}

#delimit ;
file write tab2
    "\midrule" _newline
    "\multicolumn{3}{l}{\textit{Regression (+ $k$, cogs, year$\times$nace2 FE)}} \\" _newline
;
#delimit cr

* Regression premiums
foreach spec in A B C D E OLS {
    if "`spec'" == "A" local slabel "Base (surv.\ + pp)"
    if "`spec'" == "B" local slabel "No survival"
    if "`spec'" == "C" local slabel "No pp in Markov"
    if "`spec'" == "D" local slabel "Plain ACF"
    if "`spec'" == "E" local slabel "Translog"
    if "`spec'" == "OLS" local slabel "OLS"

    xi: qui reg lmu_`spec' pp_dummy k cogs i.year*i.nace2, cluster(id)
    local reg = _b[pp_dummy]
    local se  = _se[pp_dummy]
    local N   = e(N)

    xi: qui reg lmu_`spec' pp_ever_3y k cogs i.year*i.nace2, cluster(id)
    local reg3 = _b[pp_ever_3y]
    local se3  = _se[pp_ever_3y]

    file write tab2 "`slabel'"
    file write tab2 " & " %6.3f (`reg') " & " %6.3f (`reg3') " \\" _newline
    file write tab2 " & (" %5.3f (`se') ") & (" %5.3f (`se3') ") \\" _newline
}

#delimit ;
file write tab2
    "\midrule" _newline
    "$N$ & \multicolumn{2}{c}{" %6.0fc (`N') "} \\" _newline
    "\bottomrule" _newline
    "\multicolumn{3}{p{10cm}}{\footnotesize " _newline
    "Notes: Premium = log markup difference between procurement and "
    "non-procurement firms. Regression controls for $k$, cogs, "
    "year$\times$nace2 FE; SEs clustered by firm. "
    "$pp_t$ = active supplier in year $t$; "
    "$pp^{3y}_t$ = active in $t$, $t{-}1$, or $t{-}2$.} \\" _newline
    "\end{tabular}" _newline
    "\end{table}" _newline
;
#delimit cr

file close tab2
dis _newline "  Saved: table_premium.tex"

*===============================================================================
* 5. CROSS-VALIDATION WITH PYTHON
*===============================================================================

dis _newline(2) "============================================================"
dis "  CROSS-VALIDATION: Stata vs Python"
dis "============================================================"

use "$outpath/paper_coefficients.dta", clear

dis _newline "  Stata base spec (A):"
dis "  NACE 41: b_k=" %7.4f b_k_A[1] "  b_cogs=" %7.4f b_cogs_A[1]
dis "  NACE 42: b_k=" %7.4f b_k_A[2] "  b_cogs=" %7.4f b_cogs_A[2]
dis "  NACE 43: b_k=" %7.4f b_k_A[3] "  b_cogs=" %7.4f b_cogs_A[3]

dis _newline "  Python base spec (from paper_pf_estimates.csv):"
preserve
import delimited "$outpath/paper_pf_estimates.csv", clear varnames(1)
keep if spec == "A"
dis "  NACE 41: b_k=" %7.4f b_k[1] "  b_cogs=" %7.4f b_cogs[1]
dis "  NACE 42: b_k=" %7.4f b_k[2] "  b_cogs=" %7.4f b_cogs[2]
dis "  NACE 43: b_k=" %7.4f b_k[3] "  b_cogs=" %7.4f b_cogs[3]
restore

*===============================================================================

dis _newline "============================================================"
dis "  DONE."
dis "  table_pf_estimates.tex"
dis "  table_premium.tex"
dis "============================================================"

log close
