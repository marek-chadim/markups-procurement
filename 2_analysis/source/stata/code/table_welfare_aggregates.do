*===============================================================================
* table_welfare_aggregates.do — Fiscal Welfare Aggregates (§7 Table 22)
*
* Port of fiscal_welfare_tenders.py. Computes Kaldor-Hicks welfare benchmarks
* for Czech public procurement 2003-2022 using the Datlab master tender
* register. Two benchmarks:
*   1. Engineer-estimate: Σ max(0, bid_final_price − lot_estimated_price)
*   2. Marginal-cost:     ((μ−1)/μ) × Σ bid_final_price  with μ = 1.14
*
* Input:  ../../1_data/input/datlab/master_tender_analytics.csv
* Output: ../../output/tables/fiscal_welfare_aggregates.tex
*===============================================================================

dis _newline "--- table_welfare_aggregates.do ---"

local tenders "../../../1_data/input/datlab/master_tender_analytics.csv"
cap confirm file "`tenders'"
if _rc != 0 {
    dis "  SKIP: tender data not found at `tenders'"
    exit
}

import delimited "`tenders'", clear varn(1) bindquote(strict) stripquote(yes) ///
    case(lower) encoding(UTF-8)
cap destring year, replace force
cap destring lot_estimated_price bid_final_price, replace force

dis "  Loaded tenders: " _N " rows"

* Drop rows with missing year/prices or zeros
drop if mi(year) | mi(bid_final_price) | bid_final_price <= 0

summ year
local ymin = r(min)
local ymax = r(max)
local nyears = `ymax' - `ymin' + 1

* FULL universe spend
qui summ bid_final_price
local total_spend_bn = r(sum) / 1e9

* Subset: both fields populated, winsorized
preserve
keep if !mi(lot_estimated_price) & lot_estimated_price > 0
gen rel_price = bid_final_price / lot_estimated_price
keep if rel_price > 0.01 & rel_price < 10
local n_sub = _N
qui summ bid_final_price
local sub_spend_bn = r(sum) / 1e9
qui summ lot_estimated_price
local sub_estimate_bn = r(sum) / 1e9

* Engineer-estimate benchmark: sum of max(final - est, 0).
* Share is computed against sub_spend (procurement spending), matching
* the paper prose claim "X% of procurement spending".
gen excess = max(bid_final_price - lot_estimated_price, 0) / 1e9
qui summ excess
local eng_bench_bn = r(sum)
local eng_bench_share = `eng_bench_bn' / `sub_spend_bn'

* Per-year average
local eng_bench_per_year = `eng_bench_bn' / `nyears'
restore

* Marginal-cost benchmark across full universe
* μ = 1.14 → deadweight share = (μ-1)/μ = 0.1228
local mu = 1.14
local mc_share = (`mu' - 1) / `mu'
local mc_bench_bn = `mc_share' * `total_spend_bn'
local mc_bench_per_year = `mc_bench_bn' / `nyears'

* Gap
local gap_pp = (`mc_share' - `eng_bench_share') * 100

* Write LaTeX table
cap mkdir "$output/tables"
tempname tf
file open `tf' using "../../output/tables/fiscal_welfare_aggregates.tex", write replace
#delimit ;
file write `tf'
    "\begin{table}[htbp]\centering" _n
    "\caption{Fiscal Welfare Aggregates (2003--2022)}" _n
    "\label{tab:fiscal_welfare}" _n
    "\begin{threeparttable}" _n
    "\begin{tabular}{lrr}" _n
    "\toprule" _n
    " & \multicolumn{1}{c}{Total (bn CZK)} & \multicolumn{1}{c}{Per year (bn CZK)} \\" _n
    "\midrule" _n
    "\multicolumn{3}{l}{\textit{Panel A: Engineer-estimate benchmark (narrow)}} \\" _n
    "Sub-universe spend "
;
#delimit cr
file write `tf' " & " %9.1f (`sub_spend_bn') " & " %9.1f (`sub_spend_bn'/`nyears') " \\" _n
file write `tf' "Engineer estimate" " & " %9.1f (`sub_estimate_bn') " & " %9.1f (`sub_estimate_bn'/`nyears') " \\" _n
file write `tf' "Excess (bid $>$ estimate)" " & " %9.1f (`eng_bench_bn') " & " %9.1f (`eng_bench_per_year') " \\" _n
file write `tf' "\quad\% of spend" " & " %9.2f (`eng_bench_share'*100) "\% & \\" _n

file write `tf' "\addlinespace" _n
file write `tf' "\multicolumn{3}{l}{\textit{Panel B: Marginal-cost benchmark (structural, $\mu = 1.14$)}} \\" _n
file write `tf' "Universe spend" " & " %9.1f (`total_spend_bn') " & " %9.1f (`total_spend_bn'/`nyears') " \\" _n
file write `tf' "Deadweight ($(\mu-1)/\mu$)" " & " %9.1f (`mc_bench_bn') " & " %9.1f (`mc_bench_per_year') " \\" _n
file write `tf' "\quad\% of spend" " & " %9.2f (`mc_share'*100) "\% & \\" _n

file write `tf' "\addlinespace" _n
file write `tf' "Gap (B-A, pp)" " & " %9.2f (`gap_pp') "\% & \\" _n

#delimit ;
file write `tf'
    "\bottomrule" _n
    "\end{tabular}" _n
    "\begin{tablenotes}\footnotesize" _n
    "\item \emph{Notes:} Panel A is computed on the subset of tenders with "
    "both engineer estimate and final bid, winsorized on the estimate/final "
    "ratio at (0.01, 10). Panel B uses the full universe of tenders with "
    "final price $>$ 0 and the paper's headline firm-level procurement "
    "markup premium $\mu = 1.14$. The gap (B-A) captures the baseline "
    "industry markup pre-priced into engineer estimates, which "
    "procurement-competition reforms alone cannot recover." _n
    "\end{tablenotes}" _n
    "\end{threeparttable}" _n
    "\end{table}" _n
;
#delimit cr
file close `tf'
dis "  Saved: fiscal_welfare_aggregates.tex"
dis "  Headlines: engineer-benchmark = " %6.2f (`eng_bench_share'*100) ///
    "% of estimate; marginal-cost = " %6.2f (`mc_share'*100) "% of spend"
