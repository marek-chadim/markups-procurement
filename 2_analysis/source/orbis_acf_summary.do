*===============================================================================
* orbis_acf_summary.do
*
* Generates Stata figures and LaTeX tables from Orbis ACF estimation results
* produced by orbis_acf_estimation.py (Python ACF estimator).
*
* Input:  orbis_acf_results.csv (from Python pipeline)
* Output: orbis_markup_by_industry.pdf, orbis_acf_results.tex
*
* Usage:
*   cd markups-procurement/2_analysis/source
*   stata-mp -e do orbis_acf_summary.do
*===============================================================================

clear all
set more off
set scheme s2color

local base  ".."
local tabs  "`base'/output/tables"
local figs  "`base'/output/figures"

cap log close _all
log using "`base'/output/data/orbis_acf_summary.log", text replace

dis "============================================================"
dis "  Orbis ACF Summary — Figures & Tables"
dis "============================================================"

* Try pipeline output first, then thesis output
local csvpath "`tabs'/orbis_acf_results.csv"
cap confirm file "`csvpath'"
if _rc {
    local csvpath "/Users/marek/Desktop/io/project/thesis/replication/outputs/tables/orbis_acf_results.csv"
}

import delimited "`csvpath'", clear
dis "  Loaded: `csvpath'"
dis "  Industries: " _N

*===============================================================================
* Figure: Markup by industry (bar chart)
*===============================================================================

* Sort by markup for better visualization
gsort -markup_mean

gen order = _n
labmask order, values(industry)

* Bar chart: mean markup + premium
twoway (bar markup_mean order, barwidth(0.6) fcolor(cranberry%70) ///
            lcolor(cranberry)) ///
       (rcap markup_p10 markup_p90 order, lcolor(navy) lwidth(medthick)), ///
       ytitle("Markup") xtitle("") ///
       ylabel(0(0.5)2) ///
       xlabel(1/`=_N', valuelabel angle(30) labsize(small)) ///
       title("ACF Markup Estimates by Industry (Orbis)") ///
       legend(order(1 "Mean markup" 2 "P10-P90 range") ///
              ring(0) pos(1) size(small))
graph export "`figs'/orbis_markup_by_industry.pdf", replace
dis "  Saved: orbis_markup_by_industry.pdf"

* Bar chart: procurement premium
gsort -pp_premium
replace order = _n
labmask order, values(industry)

gen zero = 0
twoway (bar pp_premium order if pp_premium != ., barwidth(0.6) ///
            fcolor(forest_green%70) lcolor(forest_green)), ///
       ytitle("Procurement premium (log markup)") xtitle("") ///
       xlabel(1/`=_N', valuelabel angle(30) labsize(small)) ///
       title("Procurement Markup Premium by Industry (Orbis)") ///
       legend(off) yline(0, lcolor(gs10))
graph export "`figs'/orbis_premium_by_industry.pdf", replace
dis "  Saved: orbis_premium_by_industry.pdf"

*===============================================================================
* LaTeX table
*===============================================================================

local fn "`tabs'/orbis_acf_results.tex"
cap file close tab
file open tab using "`fn'", write replace

file write tab "\begin{table}[htbp]\centering" _newline
file write tab "\caption{ACF Production Function Estimates by Industry (Orbis)}\label{tab:orbis_acf}" _newline
file write tab "\begin{threeparttable}" _newline
file write tab "\begin{tabular}{llccccc}" _newline
file write tab "\toprule" _newline
file write tab "Industry & NACE & $\hat{\theta}_{\text{cogs}}$ & Mean $\mu$ & P50 $\mu$ & PP Premium & $N$ \\" _newline
file write tab "\midrule" _newline

gsort nace2
local nrows = _N
forvalues i = 1/`nrows' {
    local ind = industry[`i']
    local nc  = nace2[`i']
    local th  : di %5.3f b_cogs[`i']
    local mu  : di %5.3f markup_mean[`i']
    local p50 : di %5.3f markup_p50[`i']
    local pp  : di %5.3f pp_premium[`i']
    local nn  : di %9.0fc n_obs[`i']
    file write tab "`ind' & `nc' & `th' & `mu' & `p50' & `pp' & `nn' \\" _newline
}

file write tab "\bottomrule" _newline
file write tab "\end{tabular}" _newline
file write tab "\begin{tablenotes}\footnotesize" _newline
file write tab "\item \textit{Notes:} Cobb-Douglas ACF estimates (Ackerberg, Caves \& Frazer, 2015). $\hat{\theta}_{\text{cogs}}$: output elasticity of variable input (COGS). PP Premium: coefficient on procurement dummy in $\ln\mu = \beta_0 + \beta_1 D + \text{year FE} + k + \text{cogs}$, HC1 SEs. Data: Orbis Bureau van Dijk via WRDS." _newline
file write tab "\end{tablenotes}" _newline
file write tab "\end{threeparttable}" _newline
file write tab "\end{table}" _newline
file close tab
dis "  Saved: orbis_acf_results.tex"

dis _newline "Done."
log close
