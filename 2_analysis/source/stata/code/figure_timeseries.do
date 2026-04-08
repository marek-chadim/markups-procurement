*===============================================================================
* figure_timeseries.do — Figure 1: Aggregate Markup Time Series
*
* Analog of DGM Figure 5:
*   DGM: cost-weighted aggregate markup over time, quantity vs revenue
*   Us:  cost-weighted aggregate markup by spec, with procurement decomposition
*
* Output: output/figure_timeseries.pdf
*===============================================================================

dis _newline "--- figure_timeseries.do ---"

use "$data/markups_panel.dta", clear

* Cost weights: firm's cogs share of total cogs in year
bys year: egen tot_cogs = total(exp(cogs))
gen w_cost = exp(cogs) / tot_cogs

*-----------------------------------------------------------------------
* Panel (a): Aggregate markup by specification
*-----------------------------------------------------------------------

preserve
collapse (mean) mu_simple_A=mu_A mu_simple_D=mu_D mu_simple_E=mu_E mu_simple_OLS=mu_OLS ///
    [aw=w_cost], by(year)

* Normalize to first year
foreach s in A D E OLS {
    qui sum mu_simple_`s' if year == 2006
    if r(N) > 0 {
        local base = r(mean)
        replace mu_simple_`s' = mu_simple_`s' / `base'
    }
}

twoway (line mu_simple_A year, lcolor(navy) lwidth(medthick)) ///
    (line mu_simple_D year, lcolor(cranberry) lpattern(dash)) ///
    (line mu_simple_E year, lcolor(forest_green) lpattern(shortdash)) ///
    (line mu_simple_OLS year, lcolor(orange) lpattern(longdash_dot)), ///
    legend(order(1 "Base" 2 "Plain" 3 "Translog" 4 "OLS") ///
        rows(1) position(6)) ///
    ytitle("Cost-weighted markup (normalized)") xtitle("Year") ///
    title("Aggregate Markup Over Time") ///
    scheme(s2color) graphregion(color(white))

graph export "$output/figure_timeseries_a.pdf", replace
restore

*-----------------------------------------------------------------------
* Panel (b): Procurement vs non-procurement firms (base spec)
*-----------------------------------------------------------------------

preserve
collapse (mean) mu_pp=mu_A [aw=w_cost], by(year pp_dummy)

reshape wide mu_pp, i(year) j(pp_dummy)

twoway (line mu_pp0 year, lcolor(navy) lwidth(medthick)) ///
    (line mu_pp1 year, lcolor(cranberry) lwidth(medthick) lpattern(dash)), ///
    legend(order(1 "Non-procurement" 2 "Procurement") ///
        rows(1) position(6)) ///
    ytitle("Cost-weighted markup (Base spec)") xtitle("Year") ///
    title("Markup by Procurement Status") ///
    scheme(s2color) graphregion(color(white))

graph export "$output/figure_timeseries_b.pdf", replace
restore

dis "  Saved: figure_timeseries_a.pdf, figure_timeseries_b.pdf"
