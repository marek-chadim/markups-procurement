*===============================================================================
* figure_binscatter.do — Figure 2: Cross-Specification Binscatters
*
* Analog of DGM Figure 4:
*   DGM: binscatter of quantity vs revenue markups (levels & FD)
*   Us:  binscatter of base vs alternative spec markups (levels & FD)
*
* Requires: ssc install binscatter (if not installed)
* Output: output/figure_binscatter_*.pdf
*===============================================================================

dis _newline "--- figure_binscatter.do ---"

* Install binscatter if needed
cap which binscatter
if _rc ssc install binscatter, replace

use "$data/markups_panel.dta", clear

* Balance on non-missing for key specs
drop if mi(l_mu_A) | mi(l_mu_D) | mi(l_mu_E) | mi(l_mu_OLS)

*-----------------------------------------------------------------------
* Panel (a): Base vs Plain (levels)
*-----------------------------------------------------------------------

binscatter l_mu_D l_mu_A, nq(50) ///
    xtitle("Log markup: Base") ytitle("Log markup: Plain") ///
    title("Base vs. Plain (levels)") ///
    lcolor(cranberry) mcolor(navy) ///
    scheme(s2color) graphregion(color(white)) ///
    absorb(nace2)

graph export "$output/figure_binscatter_base_plain.pdf", replace

*-----------------------------------------------------------------------
* Panel (b): Base vs Translog (levels)
*-----------------------------------------------------------------------

binscatter l_mu_E l_mu_A, nq(50) ///
    xtitle("Log markup: Base") ytitle("Log markup: Translog") ///
    title("Base vs. Translog (levels)") ///
    lcolor(forest_green) mcolor(navy) ///
    scheme(s2color) graphregion(color(white)) ///
    absorb(nace2)

graph export "$output/figure_binscatter_base_tl.pdf", replace

*-----------------------------------------------------------------------
* Panel (c): Base vs OLS (levels)
*-----------------------------------------------------------------------

binscatter l_mu_OLS l_mu_A, nq(50) ///
    xtitle("Log markup: Base") ytitle("Log markup: OLS") ///
    title("Base vs. OLS (levels)") ///
    lcolor(orange) mcolor(navy) ///
    scheme(s2color) graphregion(color(white)) ///
    absorb(nace2)

graph export "$output/figure_binscatter_base_ols.pdf", replace

*-----------------------------------------------------------------------
* Panel (d): First differences — Base vs Plain
*-----------------------------------------------------------------------

drop if mi(FD_l_mu_A) | mi(FD_l_mu_D)

binscatter FD_l_mu_D FD_l_mu_A, nq(50) ///
    xtitle("{&Delta} Log markup: Base") ytitle("{&Delta} Log markup: Plain") ///
    title("Base vs. Plain (first differences)") ///
    lcolor(cranberry) mcolor(navy) ///
    scheme(s2color) graphregion(color(white)) ///
    absorb(nace2)

graph export "$output/figure_binscatter_fd_base_plain.pdf", replace

dis "  Saved: figure_binscatter_*.pdf"
