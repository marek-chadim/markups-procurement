*===============================================================================
* figure_event_study_2012.do — Event-study plot around the 2012 single-bidding
*                              ban (Act 55/2012)
*
* Reads the event-study matrices saved by table_reforms_mechanism.do and
* produces a 2-panel figure (log markup + EBIT-like margin) with 95% CIs.
* Uses the Paul Tol style helper (graph_markups.do) for consistent palette.
*
* Input:  $output/data/reform2012_event_study.dta
*           (cols: tau, b_mu, se_mu, b_pr, se_pr — 11 rows, tau in -5..+5)
* Output: 4_paper/input/analysis_figures/event_study_2012.pdf
*===============================================================================

dis _newline "--- figure_event_study_2012.do ---"

cap confirm file "$output/data/reform2012_event_study.dta"
if _rc != 0 {
    dis "  SKIP: event-study matrix not found"
    exit
}

do "$code/graph_markups.do"

use "$output/data/reform2012_event_study.dta", clear
dis "  loaded rows: " _N

* 95% confidence intervals
gen lo_mu = b_mu - 1.96 * se_mu
gen hi_mu = b_mu + 1.96 * se_mu
gen lo_pr = b_pr - 1.96 * se_pr
gen hi_pr = b_pr + 1.96 * se_pr

* Log markup panel
twoway ///
    (rcap lo_mu hi_mu tau, lcolor("${markups_blue}") lwidth(medium)) ///
    (scatter b_mu tau, mcolor("${markups_blue}") msymbol(O) msize(medium)) ///
    (line b_mu tau, lcolor("${markups_blue}") lwidth(medthick)) ///
    , ${markups_gropts} ///
    title("(a) Log ACF markup", size(medium) position(11)) ///
    ytitle("DiD coefficient (log points)") ///
    xtitle("Years from 2012 reform") ///
    xlabel(-5(1)5) ///
    yline(0, lcolor(gs10) lpattern(dash)) ///
    xline(0, lcolor(gs10) lpattern(dash)) ///
    legend(off) ///
    name(g_mu, replace)

* EBIT-like margin panel
twoway ///
    (rcap lo_pr hi_pr tau, lcolor("${markups_pink}") lwidth(medium)) ///
    (scatter b_pr tau, mcolor("${markups_pink}") msymbol(O) msize(medium)) ///
    (line b_pr tau, lcolor("${markups_pink}") lwidth(medthick)) ///
    , ${markups_gropts} ///
    title("(b) EBIT-like margin (VA-W)/GO", size(medium) position(11)) ///
    ytitle("DiD coefficient (ratio points)") ///
    xtitle("Years from 2012 reform") ///
    xlabel(-5(1)5) ///
    yline(0, lcolor(gs10) lpattern(dash)) ///
    xline(0, lcolor(gs10) lpattern(dash)) ///
    legend(off) ///
    name(g_pr, replace)

* Combine (graph combine doesn't accept bgcolor, so pass only region options)
graph combine g_mu g_pr, col(2) ///
    graphregion(color(white)) plotregion(color(white)) ///
    ysize(3.5) xsize(9) ///
    name(g_comb, replace)

cap mkdir "../../../4_paper/input/analysis_figures"
graph export "../../../4_paper/input/analysis_figures/event_study_2012.pdf", ///
    replace
dis "  Saved: event_study_2012.pdf"

graph drop g_mu g_pr g_comb
