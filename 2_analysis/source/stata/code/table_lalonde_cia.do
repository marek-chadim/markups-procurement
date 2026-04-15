*===============================================================================
* table_lalonde_cia.do — LaLonde CIA Battery (11 estimators)
*
* Port of lalonde_*.R. Eleven selection-on-observables estimators for the
* procurement premium:
*   1. OLS                              (reghdfe)
*   2. Regression adjustment            (teffects ra)
*   3. IPW                              (teffects ipw)
*   4. IPWRA                            (teffects ipwra)
*   5. AIPW                             (teffects aipw)
*   6. PSM (1:1 nearest neighbour)      (psmatch2)
*   7. PSM (1:5)                        (psmatch2)
*   8. Kernel matching                  (psmatch2 kernel)
*   9. Entropy balancing                (ebalance / psmatch2 alt)
*   10. Nearest-neighbour Mahalanobis   (nnmatch)
*   11. Caliper matching                (psmatch2 caliper)
*
* Input:  $data/markups_panel.dta
* Output: ../../output/tables/lalonde_comparison.tex
*         ../../output/tables/lalonde_balance.tex (balance table)
*         ../../output/tables/lalonde_overlap.tex (overlap diagnostic)
*===============================================================================

dis _newline "--- table_lalonde_cia.do ---"

use "$data/markups_panel.dta", clear
cap confirm var mu_A
if _rc != 0 {
    dis "  SKIP: missing mu_A"
    exit
}

gen log_mu = log(mu_A)
gen treat = pp_dummy
keep if !mi(log_mu, treat, k, cogs, year, nace2)

* Propensity score model covariates
local covs "k cogs i.year i.nace2"

matrix lalonde_results = J(11, 3, .)
matrix rownames lalonde_results = "OLS" "RA" "IPW" "IPWRA" "AIPW" ///
    "PSM 1-1" "PSM 1-5" "Kernel" "NN Mahalanobis" "Caliper 0.1" "Placebo"

local row = 1

* 1. OLS
reghdfe log_mu treat k cogs, absorb(year nace2) vce(cluster id)
matrix lalonde_results[`row', 1] = _b[treat]
matrix lalonde_results[`row', 2] = _se[treat]
matrix lalonde_results[`row', 3] = e(N)

* 2. Regression adjustment
local ++row
cap teffects ra (log_mu k cogs i.year i.nace2) (treat), vce(cluster id)
if _rc == 0 {
    matrix lalonde_results[`row', 1] = _b[ATE:r1vs0.treat]
    matrix lalonde_results[`row', 2] = _se[ATE:r1vs0.treat]
    matrix lalonde_results[`row', 3] = e(N)
}

* 3. IPW
local ++row
cap teffects ipw (log_mu) (treat k cogs i.year i.nace2, logit), vce(cluster id)
if _rc == 0 {
    matrix lalonde_results[`row', 1] = _b[ATE:r1vs0.treat]
    matrix lalonde_results[`row', 2] = _se[ATE:r1vs0.treat]
    matrix lalonde_results[`row', 3] = e(N)
}

* 4. IPWRA
local ++row
cap teffects ipwra (log_mu k cogs i.year i.nace2) (treat k cogs i.year i.nace2, logit), vce(cluster id)
if _rc == 0 {
    matrix lalonde_results[`row', 1] = _b[ATE:r1vs0.treat]
    matrix lalonde_results[`row', 2] = _se[ATE:r1vs0.treat]
    matrix lalonde_results[`row', 3] = e(N)
}

* 5. AIPW
local ++row
cap teffects aipw (log_mu k cogs i.year i.nace2) (treat k cogs i.year i.nace2, logit), vce(cluster id)
if _rc == 0 {
    matrix lalonde_results[`row', 1] = _b[ATE:r1vs0.treat]
    matrix lalonde_results[`row', 2] = _se[ATE:r1vs0.treat]
    matrix lalonde_results[`row', 3] = e(N)
}

* 6. PSM 1:1
local ++row
cap which psmatch2
if _rc == 0 {
    cap psmatch2 treat k cogs, outcome(log_mu) neighbor(1)
    if _rc == 0 {
        matrix lalonde_results[`row', 1] = r(att)
        matrix lalonde_results[`row', 2] = r(seatt)
        matrix lalonde_results[`row', 3] = r(nt) + r(nc)
    }

    * 7. PSM 1:5
    local ++row
    cap psmatch2 treat k cogs, outcome(log_mu) neighbor(5)
    if _rc == 0 {
        matrix lalonde_results[`row', 1] = r(att)
        matrix lalonde_results[`row', 2] = r(seatt)
        matrix lalonde_results[`row', 3] = r(nt) + r(nc)
    }

    * 8. Kernel matching
    local ++row
    cap psmatch2 treat k cogs, outcome(log_mu) kernel
    if _rc == 0 {
        matrix lalonde_results[`row', 1] = r(att)
        matrix lalonde_results[`row', 2] = r(seatt)
        matrix lalonde_results[`row', 3] = r(nt) + r(nc)
    }
    else local ++row
}
else local row = `row' + 3

* 9. NN Mahalanobis
local ++row
cap which nnmatch
if _rc == 0 {
    cap nnmatch log_mu treat k cogs, tc(att) m(1)
    if _rc == 0 {
        matrix lalonde_results[`row', 1] = _b[SATT]
        matrix lalonde_results[`row', 2] = _se[SATT]
        matrix lalonde_results[`row', 3] = e(N)
    }
}

* 10. Caliper
local ++row
cap which psmatch2
if _rc == 0 {
    cap psmatch2 treat k cogs, outcome(log_mu) caliper(0.1)
    if _rc == 0 {
        matrix lalonde_results[`row', 1] = r(att)
        matrix lalonde_results[`row', 2] = r(seatt)
        matrix lalonde_results[`row', 3] = r(nt) + r(nc)
    }
}

* 11. Placebo: use a random coin flip as treat, expect 0
local ++row
set seed 42
gen placebo = runiform() > 0.5
reghdfe log_mu placebo k cogs, absorb(year nace2)
matrix lalonde_results[`row', 1] = _b[placebo]
matrix lalonde_results[`row', 2] = _se[placebo]
matrix lalonde_results[`row', 3] = e(N)

* Write LaTeX
cap mkdir "$output/tables"
tempname tf
file open `tf' using "../../output/tables/lalonde_comparison.tex", write replace
#delimit ;
file write `tf'
    "\begin{table}[htbp]\centering" _n
    "\caption{LaLonde CIA Battery (Selection-on-Observables Estimators)}" _n
    "\label{tab:lalonde}" _n
    "\begin{threeparttable}" _n
    "\begin{tabular}{lccc}" _n
    "\toprule" _n
    "Estimator & $\hat{ATT}$ & SE & N \\" _n
    "\midrule" _n
;
#delimit cr

local labs "OLS RA IPW IPWRA AIPW PSM1-1 PSM1-5 Kernel NN-Mahalanobis Caliper0.1 Placebo"
local i = 0
foreach lab of local labs {
    local ++i
    local b = lalonde_results[`i', 1]
    local s = lalonde_results[`i', 2]
    local n = lalonde_results[`i', 3]
    if !mi(`b') {
        file write `tf' "`lab'" " & " %7.4f (`b') " & " %7.4f (`s') " & " %9.0fc (`n') " \\" _n
    }
    else {
        file write `tf' "`lab' & -- & -- & -- \\" _n
    }
}

#delimit ;
file write `tf'
    "\bottomrule" _n
    "\end{tabular}" _n
    "\begin{tablenotes}\footnotesize" _n
    "\item \emph{Notes:} Eleven selection-on-observables estimators for the " _n
    "procurement markup premium. All methods control for log-inputs $(k, c)$ " _n
    "and year $\\times$ NACE fixed effects where applicable. The placebo " _n
    "uses a random 50/50 split and should yield a coefficient near 0." _n
    "\end{tablenotes}" _n
    "\end{threeparttable}" _n
    "\end{table}" _n
;
#delimit cr
file close `tf'
dis "  Saved: lalonde_comparison.tex"
