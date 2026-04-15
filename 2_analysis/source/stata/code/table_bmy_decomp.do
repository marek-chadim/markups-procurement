*===============================================================================
* table_bmy_decomp.do — BMY (Bond-Morris-Yurukoglu 2021 JME) Decomposition
*
* Port of bmy_czech_analysis.py. Decomposes aggregate Czech construction
* markup change into within-firm, reallocation, and net-entry components.
* Compares to DLEU aggregate (identity check).
*
* Input:  $data/markups_panel.dta
* Output: ../../output/tables/bmy_sample_comparison.tex
*===============================================================================

dis _newline "--- table_bmy_decomp.do ---"

use "$data/markups_panel.dta", clear

foreach v in mu_A year id go {
    cap confirm var `v'
    if _rc != 0 {
        dis "  SKIP: missing `v'"
        exit
    }
}

cap drop sales
gen sales = exp(go)

* Entry/exit flags (firm-first-year / firm-last-year)
bys id (year): gen first_year = year[1]
bys id (year): gen last_year = year[_N]
gen entrant = (year == first_year)
gen exitor  = (year == last_year)

* Sales weights within year
bys year: egen yr_sales = total(sales)
gen mshare = sales / yr_sales

* Period comparison: first vs last year of sample
summ year
local y0 = r(min)
local y1 = r(max)

* Sales-weighted aggregate markup at start and end
qui summ mu_A [aw=mshare] if year == `y0'
local mu0 = r(mean)
qui summ mu_A [aw=mshare] if year == `y1'
local mu1 = r(mean)
local d_mu = `mu1' - `mu0'

* Continuing-firm set: firms present in both `y0` and `y1`
bys id: egen in_start = max(cond(year == `y0', 1, 0))
bys id: egen in_end   = max(cond(year == `y1', 1, 0))
gen continuing = in_start & in_end

* Within-continuing component: Σ_i ms_{i,y0} × (μ_{i,y1} − μ_{i,y0})
preserve
keep if continuing
keep if year == `y0' | year == `y1'
bys id (year): gen mu_L = mu_A[_n-1]
bys id (year): gen ms_L = mshare[_n-1]
bys id (year): gen d_mu_i = mu_A - mu_L
bys id (year): gen d_ms_i = mshare - ms_L
drop if year == `y0'
qui summ d_mu_i [aw=ms_L]
local within_cont = r(mean) * r(sum_w)
qui summ d_ms_i [aw=mu_L]
local realloc_cont = r(mean) * r(sum_w)
restore

* Net entry = (mu1_entrants × ms_entrants) − (mu0_exits × ms_exits)
qui summ mu_A [aw=mshare] if year == `y1' & !continuing
local mu_ent = r(mean)
qui summ mshare if year == `y1' & !continuing
local ms_ent = r(sum)

qui summ mu_A [aw=mshare] if year == `y0' & !continuing
local mu_xit = r(mean)
qui summ mshare if year == `y0' & !continuing
local ms_xit = r(sum)

local net_entry = (`mu_ent' - `mu_xit') * ((`ms_ent' + `ms_xit') / 2)

local total_decomp = `within_cont' + `realloc_cont' + `net_entry'

* Write LaTeX
cap mkdir "$output/tables"
tempname tf
file open `tf' using "../../output/tables/bmy_sample_comparison.tex", write replace
#delimit ;
file write `tf'
    "\begin{table}[htbp]\centering" _n
    "\caption{BMY (2021) Sample Decomposition of $\Delta\bar\mu$}" _n
    "\label{tab:bmy_sample}" _n
    "\begin{threeparttable}" _n
    "\begin{tabular}{lr}" _n
    "\toprule" _n
    "Component & $\Delta\bar\mu$ contribution \\" _n
    "\midrule" _n
;
#delimit cr
file write `tf' "Start-year mean ($\bar\mu_{" "`y0'" "}$)" " & " %7.3f (`mu0') " \\" _n
file write `tf' "End-year mean ($\bar\mu_{" "`y1'" "}$)" " & " %7.3f (`mu1') " \\" _n
file write `tf' "Total change $\Delta\bar\mu$" " & " %7.3f (`d_mu') " \\" _n
file write `tf' "\addlinespace" _n
file write `tf' "Within-firm (continuing)" " & " %7.3f (`within_cont') " \\" _n
file write `tf' "Reallocation (continuing)" " & " %7.3f (`realloc_cont') " \\" _n
file write `tf' "Net entry" " & " %7.3f (`net_entry') " \\" _n
file write `tf' "Sum of components" " & " %7.3f (`total_decomp') " \\" _n
#delimit ;
file write `tf'
    "\bottomrule" _n
    "\end{tabular}" _n
    "\begin{tablenotes}\footnotesize" _n
    "\item \emph{Notes:} Bond, Morris \& Yurukoglu (2021) decomposition. "
    "Continuing firms are those present in both the start (`y0') and end "
    "(`y1') years. Net entry = (weighted mean of markups of period-" _n
    "`y1' entrants) $-$ (weighted mean of markups of pre-`y0' exitors), "
    "both scaled by their average sales shares." _n
    "\end{tablenotes}" _n
    "\end{threeparttable}" _n
    "\end{table}" _n
;
#delimit cr
file close `tf'
dis "  Saved: bmy_sample_comparison.tex"
dis "  Δμ̄ = " %6.3f (`d_mu') " | within=" %6.3f (`within_cont') ///
    " realloc=" %6.3f (`realloc_cont') " netent=" %6.3f (`net_entry')
