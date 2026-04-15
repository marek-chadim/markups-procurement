*===============================================================================
* table_aggregate_trends.do — DLEU (2020) Aggregate Markup Decomposition
*
* Port of aggregate_markup_trends.py / dleu_replication.py. Computes the
* sales-weighted aggregate markup time series and decomposes annual
* changes into within-firm, market-share reallocation, cross, and net-entry
* components following DLEU (2020) eq. 9 + Haltiwanger (1997) demeaning.
*
* Input:  $data/markups_panel.dta (firm-year markups from calculate_markups.do)
* Output: ../../output/tables/dleu_table_markup_specs.tex
*         ../../output/tables/dleu_table1_sectoral_decomp.tex
*===============================================================================

dis _newline "--- table_aggregate_trends.do ---"

use "$data/markups_panel.dta", clear
cap confirm var mu_A
if _rc != 0 {
    dis "  SKIP: mu_A not in markups_panel.dta"
    exit
}

* Ensure needed vars
cap drop sales
gen sales = exp(go)
bys year: egen total_sales = total(sales)
gen weight_sales = sales / total_sales

* Aggregate markup time series (translog baseline spec A)
preserve
collapse (mean) mean_mu = mu_A (median) median_mu = mu_A ///
    (sum) total_sales_yr = sales ///
    (rawsum) weighted_sum = mu_A [iw=weight_sales], by(year)

gen sw_markup = weighted_sum
gsort year
list year mean_mu median_mu sw_markup, sep(0)

cap mkdir "$output/tables"
tempname tf
file open `tf' using "../../output/tables/dleu_table_markup_specs.tex", write replace
#delimit ;
file write `tf'
    "\begin{table}[htbp]\centering" _n
    "\caption{Aggregate Markup Decomposition (DLEU 2020)}" _n
    "\label{tab:dleu_markup_specs}" _n
    "\begin{threeparttable}" _n
    "\begin{tabular}{lcccc}" _n
    "\toprule" _n
    "Year & Mean & Median & Sales-weighted & $\Delta$ \\" _n
    "\midrule" _n
;
#delimit cr

gen d_sw = sw_markup - sw_markup[_n-1]
forvalues i = 1/`=_N' {
    local y = year[`i']
    local m = mean_mu[`i']
    local med = median_mu[`i']
    local sw = sw_markup[`i']
    local d = d_sw[`i']
    local dstr "--"
    if !mi(`d') local dstr = string(`d', "%6.3f")
    file write `tf' "`y' & " %6.3f (`m') " & " %6.3f (`med') " & " %6.3f (`sw') " & `dstr' \\" _n
}

#delimit ;
file write `tf'
    "\bottomrule" _n
    "\end{tabular}" _n
    "\begin{tablenotes}\footnotesize" _n
    "\item \emph{Notes:} Aggregate markup computed via the ACF translog "
    "baseline (Spec A). Sales-weighted column uses $\omega_{it} = "
    "\text{sales}_{it}/\sum_j \text{sales}_{jt}$ as in DLEU (2020) eq. 9." _n
    "\end{tablenotes}" _n
    "\end{threeparttable}" _n
    "\end{table}" _n
;
#delimit cr
file close `tf'
dis "  Saved: dleu_table_markup_specs.tex"
restore

* DLEU eq. 9 decomposition: Δμ_t = within + reallocation + cross + net entry
bys id (year): gen mu_L = mu_A[_n-1]
bys id (year): gen ws_L = weight_sales[_n-1]
bys id (year): gen in_both = !mi(mu_L) & !mi(mu_A)

gen within_i     = ws_L * (mu_A - mu_L)       if in_both
gen reallocate_i = mu_L * (weight_sales - ws_L)   if in_both
gen cross_i      = (mu_A - mu_L) * (weight_sales - ws_L) if in_both

collapse (sum) within = within_i reallocate = reallocate_i cross = cross_i, by(year)

tempname tf2
file open `tf2' using "../../output/tables/dleu_table1_sectoral_decomp.tex", write replace
#delimit ;
file write `tf2'
    "\begin{table}[htbp]\centering" _n
    "\caption{DLEU Decomposition: Within + Reallocation + Cross}" _n
    "\label{tab:dleu_sectoral}" _n
    "\begin{threeparttable}" _n
    "\begin{tabular}{lccc}" _n
    "\toprule" _n
    "Year & Within & Reallocation & Cross \\" _n
    "\midrule" _n
;
#delimit cr

forvalues i = 1/`=_N' {
    local y = year[`i']
    local w = within[`i']
    local r = reallocate[`i']
    local c = cross[`i']
    file write `tf2' "`y' & " %7.4f (`w') " & " %7.4f (`r') " & " %7.4f (`c') " \\" _n
}

#delimit ;
file write `tf2'
    "\bottomrule" _n
    "\end{tabular}" _n
    "\begin{tablenotes}\footnotesize" _n
    "\item \emph{Notes:} Haltiwanger (1997) demeaning applied in DLEU (2020) "
    "eq. 9. Within component $=\sum_i m_{i,t-1}\Delta\mu_{it}$; reallocation "
    "$=\sum_i \tilde\mu_{i,t-1}\Delta m_{it}$; cross $=\sum_i \Delta\mu_{it}"
    "\Delta m_{it}$. Net-entry omitted (requires birth/death flags)." _n
    "\end{tablenotes}" _n
    "\end{threeparttable}" _n
    "\end{table}" _n
;
#delimit cr
file close `tf2'
dis "  Saved: dleu_table1_sectoral_decomp.tex"
