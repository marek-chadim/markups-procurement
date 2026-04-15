*===============================================================================
* table_klms_benchmark.do — KLMS (2025 AER) Benchmark Analysis
*
* Port of klms_analysis.py. Computes the Kroft-Luo-Mogstad-Setzler double
* market power benchmark for Czech construction: price markup μ̂_P and
* labour wedge (1-ε) from within-firm/firm-year FE regressions.
*
* Input:  $data/markups_panel.dta
* Output: ../../output/tables/klms_benchmark.tex (if not present in final paper)
*===============================================================================

dis _newline "--- table_klms_benchmark.do ---"

use "$data/markups_panel.dta", clear
cap confirm var mu_A
if _rc != 0 {
    dis "  SKIP: missing mu_A"
    exit
}
cap confirm var empl_mid
local has_empl = !_rc

* Log markup
gen lmu = log(mu_A)

* KLMS Table 3 analog: log markup on reform dummies with firm FE
eststo clear

* 2012 single-bid ban event study
gen post2012 = year >= 2012
gen treat_2012 = pp_dummy * post2012
eststo klms_2012: reghdfe lmu treat_2012 pp_dummy, absorb(id year) vce(cluster id)
local b_2012 = _b[treat_2012]
local se_2012 = _se[treat_2012]

* 2016 MEAT criteria event study
gen post2016 = year >= 2016
gen treat_2016 = pp_dummy * post2016
eststo klms_2016: reghdfe lmu treat_2016 pp_dummy, absorb(id year) vce(cluster id)
local b_2016 = _b[treat_2016]
local se_2016 = _se[treat_2016]

* Calibrated price markup (following KLMS eq. 8): identify from firm-level variation
* μ̂_P = (1+ε) × α_COGS / α_COGS^* where α_COGS^* is the counterfactual share
* Approximation: just report sales-weighted mean of mu_A
qui summ mu_A [aw=exp(go)]
local mu_P = r(mean)

* Labour wedge (1-ε) from regression of log employment on log output
if `has_empl' {
    cap gen log_empl = log(empl_mid)
    cap reghdfe log_empl go, absorb(id year) vce(cluster id)
    if _rc == 0 {
        local eps_l = _b[go]
        local se_eps_l = _se[go]
        local markdown_l = 1 - `eps_l'
    }
}

* Write LaTeX
cap mkdir "$output/tables"
tempname tf
file open `tf' using "../../output/tables/klms_benchmark.tex", write replace
#delimit ;
file write `tf'
    "\begin{table}[htbp]\centering" _n
    "\caption{KLMS (2025) Benchmark: Czech Construction}" _n
    "\label{tab:klms_benchmark}" _n
    "\begin{threeparttable}" _n
    "\begin{tabular}{lcc}" _n
    "\toprule" _n
    "Statistic & Value & SE \\" _n
    "\midrule" _n
;
#delimit cr
file write `tf' "Sales-weighted price markup $\hat\mu_P$" " & " %7.3f (`mu_P') " & -- \\" _n
if `has_empl' {
    file write `tf' "Labour elasticity $\hat\varepsilon_L$" " & " %7.3f (`eps_l') " & " %7.3f (`se_eps_l') " \\" _n
    file write `tf' "Labour markdown $(1-\hat\varepsilon_L)$" " & " %7.3f (`markdown_l') " & -- \\" _n
}
file write `tf' "2012 single-bid-ban ATT" " & " %7.4f (`b_2012') " & " %7.4f (`se_2012') " \\" _n
file write `tf' "2016 MEAT ATT" " & " %7.4f (`b_2016') " & " %7.4f (`se_2016') " \\" _n

#delimit ;
file write `tf'
    "\bottomrule" _n
    "\end{tabular}" _n
    "\begin{tablenotes}\footnotesize" _n
    "\item \emph{Notes:} KLMS (2025) double market power analog. Price markup "
    "is the sales-weighted mean of ACF translog Spec A markups. Labour "
    "elasticity is from a firm-year FE regression of $\log(\text{empl})$ "
    "on $\log(\text{sales})$ — a proxy for the inverse wage-elasticity "
    "under labour market power. Reform ATTs are from event-study " _n
    "regressions with firm + year FE." _n
    "\end{tablenotes}" _n
    "\end{threeparttable}" _n
    "\end{table}" _n
;
#delimit cr
file close `tf'
dis "  Saved: klms_benchmark.tex"
