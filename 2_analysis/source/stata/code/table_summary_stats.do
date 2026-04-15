*===============================================================================
* table_summary_stats.do — Summary Stats (Panel A full sample, Panel B split)
*
* Port of summary_stats.py. Computes mean/SD/P10/P50/P90/N for 10 variables
* in Panel A (full sample) and 6 variables in Panel B (split by PP).
*
* Input:  $data/analysis_panel.dta
* Output: ../../output/tables/summary_stats.tex
*===============================================================================

dis _newline "--- table_summary_stats.do ---"

use "$data/analysis_panel.dta", clear

* Derived variables
cap drop sales_cogs_ratio
gen sales_cogs_ratio = exp(go) / exp(cogs)
label var sales_cogs_ratio "Sales/COGS ratio"

* Panel A variables — include all firm-years; procurement-only vars are
* conditional on pp_dummy == 1.
cap confirm var empl_mid
local has_empl = !_rc
cap confirm var n_contracts
local has_contracts = !_rc
cap confirm var avg_bids
local has_bids = !_rc
cap confirm var single_bid_share
local has_single = !_rc

local vars_full "go cogs k"
if `has_empl' local vars_full "`vars_full' empl_mid"
local vars_full "`vars_full' sales_cogs_ratio mktshare pp_dummy"
local vars_pp ""
if `has_contracts' local vars_pp "`vars_pp' n_contracts"
if `has_bids' local vars_pp "`vars_pp' avg_bids"
if `has_single' local vars_pp "`vars_pp' single_bid_share"

cap mkdir "$output/tables"
tempname tf
file open `tf' using "../../output/tables/summary_stats.tex", write replace
#delimit ;
file write `tf'
    "\begin{table}[htbp]\centering" _n
    "\caption{Summary Statistics --- Czech Construction Panel}" _n
    "\label{tab:summary}" _n
    "\begin{threeparttable}" _n
    "\begin{tabular}{lrrrrrr}" _n
    "\toprule" _n
    "Variable & Mean & SD & P10 & Median & P90 & N \\" _n
    "\midrule" _n
    "\multicolumn{7}{l}{\textit{Panel A: Full sample}} \\" _n
;
#delimit cr

* Panel A rows
foreach v of local vars_full {
    qui summ `v', d
    local mu = r(mean)
    local sd = r(sd)
    local p10 = r(p10)
    local p50 = r(p50)
    local p90 = r(p90)
    local n = r(N)
    local lab : variable label `v'
    if "`lab'" == "" local lab "`v'"
    * Escape underscores for LaTeX
    local lab_tex = subinstr("`lab'", "_", "\_", .)
    file write `tf' "`lab_tex'" " & " %9.3f (`mu') " & " %9.3f (`sd') ///
        " & " %9.3f (`p10') " & " %9.3f (`p50') " & " %9.3f (`p90') ///
        " & " %9.0fc (`n') " \\" _n
}

if "`vars_pp'" != "" {
    file write `tf' "\addlinespace" _n
    file write `tf' "\multicolumn{7}{l}{\textit{Panel B: Procurement firms only (pp\_dummy = 1)}} \\" _n
    foreach v of local vars_pp {
        qui summ `v' if pp_dummy == 1, d
        if r(N) == 0 continue
        local mu = r(mean)
        local sd = r(sd)
        local p10 = r(p10)
        local p50 = r(p50)
        local p90 = r(p90)
        local n = r(N)
        local lab : variable label `v'
        if "`lab'" == "" local lab "`v'"
        local lab_tex = subinstr("`lab'", "_", "\_", .)
        file write `tf' "`lab_tex'" " & " %9.3f (`mu') " & " %9.3f (`sd') ///
            " & " %9.3f (`p10') " & " %9.3f (`p50') " & " %9.3f (`p90') ///
            " & " %9.0fc (`n') " \\" _n
    }
}

#delimit ;
file write `tf'
    "\bottomrule" _n
    "\end{tabular}" _n
    "\begin{tablenotes}\footnotesize" _n
    "\item \emph{Notes:} Panel A uses the full rebuilt MagnusWeb panel "
    "(1,521 firms, 9,164 firm-years, 2005--2021). Panel B conditions on "
    "firm-years where the firm is active in public procurement "
    "(pp\_dummy = 1). Sales, COGS, and capital are deflated log values; "
    "mktshare is the firm's share of year $\times$ nace2 sales." _n
    "\end{tablenotes}" _n
    "\end{threeparttable}" _n
    "\end{table}" _n
;
#delimit cr
file close `tf'
dis "  Saved: summary_stats.tex"
