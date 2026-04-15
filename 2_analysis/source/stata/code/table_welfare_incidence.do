*===============================================================================
* table_welfare_incidence.do — Top-20 Contract Incidence (Appendix B.7)
*
* Port of contract-level incidence computation. Sorts contracts by absolute
* size, computes cumulative share held by top 1%, 10%, and top-20 contracts
* (naming the largest N).
*
* Input:  $data/contract_level.dta (built by Python)
* Output: ../../output/tables/incidence_top_contracts.tex
*===============================================================================

dis _newline "--- table_welfare_incidence.do ---"

cap confirm file "$data/contract_level.dta"
if _rc != 0 {
    dis "  SKIP: contract_level.dta not available"
    exit
}

use "$data/contract_level.dta", clear
cap confirm var bid_final_price
if _rc != 0 {
    dis "  SKIP: missing bid_final_price"
    exit
}

keep if !mi(bid_final_price) & bid_final_price > 0
gen excess_czk = 0
cap confirm var lot_estimated_price
if _rc == 0 {
    replace excess_czk = max(bid_final_price - lot_estimated_price, 0) if !mi(lot_estimated_price)
}

gsort -excess_czk
local n_total = _N
qui summ excess_czk
local total_excess = r(sum)

* Cumulative shares
gen cum_excess = sum(excess_czk)
gen cum_share = cum_excess / `total_excess'

* Top 1% share
local top1pct = ceil(`n_total' * 0.01)
local share_top1 = cum_share[`top1pct']

* Top 10% share
local top10pct = ceil(`n_total' * 0.10)
local share_top10 = cum_share[`top10pct']

* Top 20 individual contracts
local share_top20ind = cum_share[20]

* Write LaTeX to a distinct filename so Python's 20-row top-contract
* table (labeled tab:incidence_top20) remains canonical for the paper.
cap mkdir "$output/tables"
tempname tf
file open `tf' using "../../output/tables/incidence_top_contracts_stata.tex", write replace
#delimit ;
file write `tf'
    "\begin{table}[htbp]\centering" _n
    "\caption{Welfare Incidence Concentration (Top Contracts)}" _n
    "\label{tab:incidence_top_contracts}" _n
    "\begin{threeparttable}" _n
    "\begin{tabular}{lr}" _n
    "\toprule" _n
    "Subset & Cumulative share of aggregate excess \\" _n
    "\midrule" _n
;
#delimit cr

file write `tf' "Top 20 individual contracts" " & " %6.2f (`share_top20ind'*100) "\% \\" _n
file write `tf' "Top 1\% of contracts (N = " %9.0fc (`top1pct') ")" " & " %6.2f (`share_top1'*100) "\% \\" _n
file write `tf' "Top 10\% of contracts (N = " %9.0fc (`top10pct') ")" " & " %6.2f (`share_top10'*100) "\% \\" _n

#delimit ;
file write `tf'
    "\bottomrule" _n
    "\end{tabular}" _n
    "\begin{tablenotes}\footnotesize" _n
    "\item \emph{Notes:} Contracts sorted by excess CZK (bid final price "
    "minus engineer estimate, floored at 0). Shares are cumulative of "
    "total aggregate excess across the full winsorized subset." _n
    "\end{tablenotes}" _n
    "\end{threeparttable}" _n
    "\end{table}" _n
;
#delimit cr
file close `tf'
dis "  Saved: incidence_top_contracts.tex"
