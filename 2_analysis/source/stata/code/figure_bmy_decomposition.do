*===============================================================================
* figure_bmy_decomposition.do — BMY (Bond-Morris-Yurukoglu 2021 JME)
*                               decomposition of Czech aggregate markup
*
* Faithful port of the decomposition in
* `references/replications/BMYReplication/Markup_Figures.do`, lines 166-218,
* adapted to the Czech construction panel (markups_panel.dta).
*
* Formulas (firm-year i,t level):
*   Lms     = L.share           (lagged sales share)
*   dmu     = D.mu               (year-over-year markup change)
*   dms     = D.share            (year-over-year share change)
*   demean  = mu - MU_agg        (markup minus contemporaneous aggregate)
*   Lmu     = L.demean           (lagged demeaned markup)
*
* Year-level aggregations (sum over firms within year):
*   Daggmu  = sum(dmu * Lms)     WITHIN-firm component
*   Daggms  = sum(dms * Lmu)     REALLOCATION (Olley-Pakes covariance)
*   Cross   = sum(dms * dmu)     CROSS term
*
* Total change and net-entry residual:
*   DMU_agg   = MU_agg_t - MU_agg_{t-1}
*   net_entry = DMU_agg - Daggmu - Daggms - Cross
*
* Cumulative series (running total from base year):
*   Dagg_sum    = cumulative sum(Daggmu)                             [WITHIN]
*   reall_inc_s = cumulative sum(Daggms + Cross)                      [REALLOC]
*   Net_sum     = cumulative sum(net_entry)                            [ENTRY]
*
* All three cumulative series are initialized at the base-year aggregate
* mu so they share a common start.
*
* Input:  $data/markups_panel.dta
* Output: $output/data/bmy_decomposition_czech.dta
*         $output/tables/bmy_decomposition_czech.csv
*         $output/tables/bmy_decomposition_summary.tex
*         4_paper/input/analysis_figures/bmy_decomposition_czech.pdf
*
* All key numbers are printed to the log and exported; prose that cites
* these values reads them from the log / CSV, not from the rendered figure.
*===============================================================================

dis _newline "--- figure_bmy_decomposition.do ---"

do "$code/graph_markups.do"

use "$data/markups_panel.dta", clear

foreach v in mu_A year id go {
    cap confirm var `v'
    if _rc != 0 {
        dis "  SKIP: missing `v'"
        exit
    }
}

drop if mi(mu_A) | mu_A <= 0 | mi(go)

cap drop sales
gen double sales = exp(go)
bysort year: egen double TOTSALES = sum(sales)
gen double msagg = sales / TOTSALES

* Aggregate sales-weighted markup by year (contemporaneous)
bysort year: egen double MU_agg = sum(msagg * mu_A)

* Set panel for lag operators
xtset id year, yearly

* Firm-year level changes and demeaned markups
gen double dmu    = D.mu_A
gen double dms    = D.msagg
gen double Lms    = L.msagg
gen double demean = mu_A - MU_agg
gen double Lmu    = L.demean

* Year-level aggregate components (sum over firms within year)
bysort year: egen double Daggmu  = sum(dmu * Lms)
bysort year: egen double Daggms  = sum(dms * Lmu)
bysort year: egen double Cross   = sum(dms * dmu)

* Reduce to one row per year
sort year
drop if year == year[_n-1]
keep year MU_agg Daggmu Daggms Cross

* Year-over-year aggregate change and net-entry residual
gen double DMU_agg   = MU_agg - MU_agg[_n-1]
gen double net_entry = DMU_agg - Daggmu - Daggms - Cross
gen double reall_inc = Daggms + Cross

* Replace missing (first year) with 0 so the cumulative sums start cleanly
replace Daggmu    = 0 if mi(Daggmu)
replace Daggms    = 0 if mi(Daggms)
replace Cross     = 0 if mi(Cross)
replace net_entry = 0 if mi(net_entry)
replace reall_inc = 0 if mi(reall_inc)

* Cumulative series (running totals)
gen double Dagg_sum    = sum(Daggmu)
gen double reall_inc_s = sum(reall_inc)
gen double Net_sum     = sum(net_entry)

* Initialize at base-year aggregate mu
qui summ year
local y_first = r(min)
local y_last  = r(max)
qui summ MU_agg if year == `y_first'
local mu_base = r(mean)
foreach x in Dagg_sum reall_inc_s Net_sum {
    replace `x' = `x' + `mu_base'
}

* Labels for readability
label var MU_agg      "Aggregate M_t (actual)"
label var Dagg_sum    "Within-firm counterfactual"
label var reall_inc_s "Reallocation counterfactual (incl. cross)"
label var Net_sum     "Net-entry counterfactual"

* ── Read key numbers from the cumulative series ──
qui summ MU_agg if year == `y_first'
local agg_first = r(mean)
qui summ MU_agg if year == `y_last'
local agg_last = r(mean)
local d_agg = `agg_last' - `agg_first'

qui summ Dagg_sum if year == `y_last'
local within_last = r(mean)
local d_within = `within_last' - `agg_first'

qui summ reall_inc_s if year == `y_last'
local realloc_last = r(mean)
local d_realloc = `realloc_last' - `agg_first'

qui summ Net_sum if year == `y_last'
local netentry_last = r(mean)
local d_netentry = `netentry_last' - `agg_first'

local sum_channels = `d_within' + `d_realloc' + `d_netentry'

dis _newline "  === BMY decomposition: `y_first' -> `y_last' ==="
dis "  Aggregate M_t:        " %6.3f `agg_first' " -> " %6.3f `agg_last' "  (change " %6.3f `d_agg' ")"
dis "  Within counterfactual:" %6.3f `agg_first' " -> " %6.3f `within_last' "  (cum within " %6.3f `d_within' ")"
dis "  Realloc counterfactual:"%6.3f `agg_first' " -> " %6.3f `realloc_last' "  (cum realloc " %6.3f `d_realloc' ")"
dis "  Net-entry counterfact:" %6.3f `agg_first' " -> " %6.3f `netentry_last' "  (cum net-entry " %6.3f `d_netentry' ")"
dis "  Sum of channels: " %6.3f `sum_channels' "  (should = aggregate change " %6.3f `d_agg' ")"

* ── Save outputs ──
cap mkdir "$output/tables"
cap mkdir "$output/data"
save "$output/data/bmy_decomposition_czech.dta", replace
export delimited year MU_agg Dagg_sum reall_inc_s Net_sum Daggmu Daggms Cross net_entry ///
    using "$output/tables/bmy_decomposition_czech.csv", replace
dis "  Saved: bmy_decomposition_czech.{dta,csv}"

* ── Summary table ──
tempname tfs
file open `tfs' using "$output/tables/bmy_decomposition_summary.tex", write replace
#delimit ;
file write `tfs'
    "\begin{table}[htbp]\centering" _n
    "\caption{BMY Decomposition of the Czech Aggregate Markup, `y_first'--`y_last'}" _n
    "\label{tab:bmy_decomposition_summary}" _n
    "\begin{threeparttable}" _n
    "\begin{tabular}{lrr}" _n
    "\toprule" _n
    "Series & `y_first' & `y_last' \\" _n
    "\midrule" _n
;
#delimit cr
file write `tfs' "Aggregate \$M_t\$ & " %6.3f (`agg_first') " & " %6.3f (`agg_last') " \\" _n
file write `tfs' "Within-firm counterfactual & " %6.3f (`agg_first') " & " %6.3f (`within_last') " \\" _n
file write `tfs' "Reallocation counterfactual & " %6.3f (`agg_first') " & " %6.3f (`realloc_last') " \\" _n
file write `tfs' "Net-entry counterfactual & " %6.3f (`agg_first') " & " %6.3f (`netentry_last') " \\" _n
file write `tfs' "\midrule" _n
file write `tfs' "Aggregate change \(\Delta M\) & \multicolumn{2}{r}{" %6.3f (`d_agg') "} \\" _n
file write `tfs' "\quad Within-firm contribution & \multicolumn{2}{r}{" %6.3f (`d_within') "} \\" _n
file write `tfs' "\quad Reallocation contribution & \multicolumn{2}{r}{" %6.3f (`d_realloc') "} \\" _n
file write `tfs' "\quad Net-entry contribution & \multicolumn{2}{r}{" %6.3f (`d_netentry') "} \\" _n
#delimit ;
file write `tfs'
    "\bottomrule" _n
    "\end{tabular}" _n
    "\begin{tablenotes}\footnotesize" _n
    "\item \emph{Notes:} Haltiwanger year-over-year decomposition following "
    "Bond, Morris, and Yurukoglu (2021, `Some Unpleasant Markup Arithmetic', "
    "JME; replication code \texttt{Markup\_Figures.do}) applied to the Czech "
    "construction panel. The within-firm counterfactual cumulates "
    "\$\\sum_i s_{i,t-1}(\\mu_{i,t}-\\mu_{i,t-1})\$ over continuing firms; the "
    "reallocation counterfactual adds the Olley-Pakes covariance term "
    "\$\\sum_i (\\mu_{i,t-1}-M_{t-1})(s_{i,t}-s_{i,t-1})\$ and the cross term; "
    "the net-entry counterfactual accumulates the residual attributable to "
    "firms entering or exiting. Each counterfactual is initialized at the "
    "base-year aggregate markup \$M_{" "`y_first'" "}\$; the sum of channel "
    "contributions equals the aggregate change by construction." _n
    "\end{tablenotes}" _n
    "\end{threeparttable}" _n
    "\end{table}" _n
;
#delimit cr
file close `tfs'
dis "  Saved: bmy_decomposition_summary.tex"

* ── Figure ──
twoway ///
    (line MU_agg year, lcolor("${markups_blue}") lwidth(thick) lpattern(solid)) ///
    (line Dagg_sum year, lcolor("${markups_red}") lpattern(dash) lwidth(medthick)) ///
    (line reall_inc_s year, lcolor("${markups_green}") lpattern(shortdash) lwidth(medthick)) ///
    (line Net_sum year, lcolor("${markups_yellow}") lpattern(dot) lwidth(medthick)) ///
    , ${markups_gropts} ///
    title("Czech aggregate markup decomposition", size(medium) position(11)) ///
    subtitle("Haltiwanger year-over-year cumulation (BMY 2021, replication port)", size(small)) ///
    ytitle("Sales-weighted aggregate markup") ///
    xtitle("Year") ///
    xline(2012, lcolor(gs10) lpattern(dot)) ///
    xline(2016, lcolor(gs10) lpattern(dot)) ///
    legend(order(1 "Aggregate (actual)" ///
                 2 "Within-firm only" ///
                 3 "Reallocation only" ///
                 4 "Net-entry only") ///
           cols(2) size(small) region(lwidth(none))) ///
    name(g_bmy, replace)

cap mkdir "../../../4_paper/input/analysis_figures"
graph export "../../../4_paper/input/analysis_figures/bmy_decomposition_czech.pdf", replace
dis "  Saved: bmy_decomposition_czech.pdf"

graph drop g_bmy
