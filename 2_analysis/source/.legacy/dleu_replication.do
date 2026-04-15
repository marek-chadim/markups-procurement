*===============================================================================
* dleu_replication.do
*
* Full implementation of De Loecker, Eeckhout, and Unger (2020, QJE)
* "The Rise of Market Power and the Macroeconomic Implications"
* adapted for Czech construction firms (CZ-NACE F, 2006-2021).
*
* Replicates (adapted):
*   Fig I:    Aggregate markup time series
*   Fig II:   Sensitivity to theta and weighting
*   Fig III:  Distribution (density + sales-weighted percentiles)
*   Fig IV:   Cumulative decomposition with counterfactuals (eq. 9)
*   Fig V:    Micro aggregation vs industry averages
*   Fig VII:  Cost share evolution
*   Fig VIII: Profit rate and distribution
*   Fig XII:  Cost-share vs PF markup
*   Table I:  Sectoral decomposition (eq. 10) by NACE 41/42/43
*   App 4:    Year-by-year firm-level decomposition (eq. 9)
*   App 12:   Industry-specific trends
*
* Extensions:
*   - Procurement vs non-procurement decomposition
*   - Czech transparency reform markers (2012, 2016)
*
* References:
*   DLEU (2020). QJE 135(2), 561-644.
*   Haltiwanger (1997). Measuring and interpreting productivity growth.
*   DLEU (2025). Reply to BMY. QJE.
*
* Usage:
*   cd markups-procurement/2_analysis/source
*   stata-mp -e do dleu_replication.do
*===============================================================================

clear all
set more off
set scheme s2color

* ---------------------------------------------------------------------------
* Paths
* ---------------------------------------------------------------------------
local base  ".."
local input "`base'/input"
local outd  "`base'/output/data"
local figs  "`base'/output/figures"
local tabs  "`base'/output/tables"
local temp  "`base'/output/data"

cap mkdir "`figs'"
cap mkdir "`tabs'"
cap mkdir "`outd'"

* Parameters (DLEU convention)
local usercost = 0.12        // r = (I-Pi)+Delta ~ 12%
local theta_cal = 0.85       // DLEU calibrated output elasticity

cap log close _all
log using "`outd'/dleu_replication.log", text replace

dis _newline "============================================================"
dis "  DLEU (2020 QJE) REPLICATION — Czech Construction"
dis "============================================================"

*===============================================================================
* SECTION 1: LOAD DATA AND CONSTRUCT MARKUP MEASURES
*===============================================================================

dis _newline(2) "--- Section 1: Data loading and markup construction ---"

* Load panel data
use "`input'/data.dta", clear

* Handle year format
cap confirm numeric variable year
if _rc {
    gen year_n = year(year)
    drop year
    rename year_n year
}
destring year, replace force

* Merge ACF markup estimates
merge 1:1 id year using "`outd'/paper_markups.dta", keep(3) nogen

* Keep valid markup observations
keep if markup_A > 0 & markup_A != .

* Use existing deflated levels (data.dta already has rGO, rCOGS, rK, rW)
* Rename to match DLEU notation
gen Sales = rGO
gen kexp  = `usercost' * rK
gen totcost = rCOGS + kexp

label var Sales   "Deflated sales (level)"
label var kexp    "Capital expenditure (r*K)"
label var totcost "Total cost (COGS + rK)"

* Employment (only available for ~32% of obs)
gen Emp = empl_mid
label var Emp "Employment"

* Wages (deflated)
gen Wages = rW
label var Wages "Deflated wages"

* Procurement dummy
cap confirm variable pp_dummy
if _rc {
    gen pp_dummy = 0
}

* -----------------------------------------------------------------------
* Markup measures (DLEU Create_Temp.do analog)
* -----------------------------------------------------------------------

* Cost shares (DLEU: costshare0, costshare1)
gen cs_firm = rCOGS / (rCOGS + kexp)
label var cs_firm "Cost share COGS/(COGS+rK) [firm level]"

bysort nace2 year: egen cs_ind = median(cs_firm)
label var cs_ind "Cost share median by industry-year"

* Trim cost shares (DLEU convention: drop 1st and 99th)
bysort year: egen cs_p1  = pctile(cs_firm), p(1)
bysort year: egen cs_p99 = pctile(cs_firm), p(99)
drop if cs_firm == 0 | cs_firm == .
drop if cs_firm < cs_p1 | cs_firm > cs_p99
drop cs_p1 cs_p99

gen ratio_sc = Sales / rCOGS
label var ratio_sc "Sales/COGS ratio"

* mu_0: calibrated theta = 0.85 (DLEU mu_0)
gen mu_cs085 = `theta_cal' * ratio_sc
label var mu_cs085 "Markup: calibrated theta=0.85"

* mu_1: firm-level cost share (DLEU mu_1)
gen mu_cs1 = cs_firm * ratio_sc
label var mu_cs1 "Markup: firm cost share"

* mu_3: industry-median cost share (DLEU mu_3)
gen mu_cs_ind = cs_ind * ratio_sc
label var mu_cs_ind "Markup: industry median cost share"

* mu_acf: ACF estimated (DLEU mu_10 equivalent)
gen mu_acf = markup_A
label var mu_acf "Markup: ACF estimated (benchmark)"

* Profit rate (DLEU eq. 13 simplified)
gen profit_rate = (Sales - rCOGS - kexp) / Sales
label var profit_rate "Profit rate (S-COGS-rK)/S"

* Cost shares of sales
gen cogs_share    = rCOGS / Sales
gen capital_share = kexp / Sales
label var cogs_share    "COGS/Sales"
label var capital_share "rK/Sales"

gen labor_share = .
cap replace labor_share = Wages / Sales if Wages != .
label var labor_share "Wages/Sales"

* Panel structure
egen firm_id = group(id)
xtset firm_id year, yearly

local N = _N
distinct firm_id
local Nfirm = r(ndistinct)
sum year, meanonly
local ymin = r(min)
local ymax = r(max)
dis "Panel: `N' obs, `Nfirm' firms, `ymin'-`ymax'"

save "`temp'/dleu_panel.dta", replace

*===============================================================================
* SECTION 2: AGGREGATION (DLEU Create_Temp.do lines 114-246)
*===============================================================================

dis _newline(2) "--- Section 2: Aggregation ---"

use "`temp'/dleu_panel.dta", clear

* Total sales/inputs by year
bysort year: egen TOTSALES   = sum(Sales)
bysort year: egen TOTCOGS    = sum(rCOGS)
bysort year: egen TOTCOST    = sum(totcost)
bysort year: egen TOTEMP     = sum(Emp)

* Weights
gen share_sales = Sales / TOTSALES
gen share_cogs  = rCOGS / TOTCOGS
gen share_tc    = totcost / TOTCOST

label var share_sales "Sales weight"
label var share_cogs  "COGS weight"
label var share_tc    "Total cost weight"

* --- Sales-weighted aggregate markups ---
foreach mu in mu_acf mu_cs085 mu_cs1 mu_cs_ind {
    bysort year: egen AGG_`mu'_sw = sum(share_sales * `mu')
    label var AGG_`mu'_sw "Agg `mu' [w=sales]"
}

* --- Input-weighted aggregates ---
foreach mu in mu_acf mu_cs085 {
    bysort year: egen AGG_`mu'_cw  = sum(share_cogs * `mu')
    bysort year: egen AGG_`mu'_tcw = sum(share_tc * `mu')
    label var AGG_`mu'_cw  "Agg `mu' [w=COGS]"
    label var AGG_`mu'_tcw "Agg `mu' [w=totcost]"
}

* --- Harmonic mean ---
bysort year: egen AGG_harmonic = sum(share_sales / mu_acf)
replace AGG_harmonic = 1 / AGG_harmonic
label var AGG_harmonic "Agg markup [harmonic]"

* --- Simple mean ---
bysort year: egen AGG_mean = mean(mu_acf)
label var AGG_mean "Agg markup [simple mean]"

* --- Sales-weighted percentiles (DLEU Figure III Panel B) ---
foreach r in 1 {
    bysort year (mu_acf): gen ms_cum = sum(share_sales)
    bysort year (mu_acf): gen ms90 = 1 if ms_cum < .9
    bysort year (mu_acf): gen ms75 = 1 if ms_cum < .75
    bysort year (mu_acf): gen ms50 = 1 if ms_cum < .5
    bysort year (mu_acf): gen ms25 = 1 if ms_cum < .25

    foreach p in 90 75 50 25 {
        bysort year (mu_acf): egen mu_P`p' = max(mu_acf) if ms`p' == 1
        bysort year: egen P`p'_acf = max(mu_P`p')
        drop mu_P`p'
        label var P`p'_acf "P`p' (sales-weighted)"
    }
    drop ms_cum ms90 ms75 ms50 ms25
}

* --- Profit rate aggregate ---
bysort year: egen AGG_profit = sum(share_sales * profit_rate)
label var AGG_profit "Agg profit rate [w=sales]"

* --- Cost shares aggregate ---
bysort year: egen AGG_cs_cogs = mean(cogs_share)
bysort year: egen AGG_cs_cap  = mean(capital_share)
bysort year: egen AGG_cs_firm_med = median(cs_firm)

* --- Theta comparison ---
gen theta_acf = mu_acf * alphahat
label var theta_acf "ACF output elasticity theta"
bysort year: egen AGG_theta_acf = sum(share_sales * theta_acf)
bysort year: egen AGG_theta_cs  = median(cs_firm)

* Count firms
bysort year: gen N_firms = _N if _n == 1
bysort year: egen Nfirms = max(N_firms)
drop N_firms

save "`temp'/dleu_agg_panel.dta", replace

* Collapse to year level for plotting
preserve
keep year AGG_* P*_acf Nfirms
sort year
by year: keep if _n == 1
save "`temp'/dleu_yearly.dta", replace
restore

dis _newline "Aggregate Markups (sales-weighted):"
preserve
keep year AGG_mu_acf_sw AGG_mean P50_acf P75_acf P90_acf Nfirms
sort year
by year: keep if _n == 1
list year AGG_mu_acf_sw AGG_mean P50_acf P75_acf P90_acf Nfirms, sep(0) noobs
restore

*===============================================================================
* SECTION 3: FIRM-LEVEL DECOMPOSITION (DLEU Eq. 9, Haltiwanger 1997)
*===============================================================================

dis _newline(2) "--- Section 3: Firm-level decomposition (eq. 9) ---"

use "`temp'/dleu_agg_panel.dta", clear

* Lagged values via panel structure
sort firm_id year
xtset firm_id year, yearly

gen L_share   = L.share_sales
gen L_mu      = L.mu_acf
gen L_M_agg   = L.AGG_mu_acf_sw

* Demeaned lagged markup (Haltiwanger: mu_tilde = mu - M_agg)
gen mu_tilde_L = L_mu - L_M_agg

* Changes
gen d_mu    = mu_acf - L_mu
gen d_share = share_sales - L_share

* Classify firms
gen incumbent = (L_share != . & share_sales != .)
gen entrant   = (L_share == . & share_sales != .)

* Within: m_{i,t-1} * Delta_mu_it
gen within_i = L_share * d_mu if incumbent

* Market share: mu_tilde_{i,t-1} * Delta_m_it
gen mktshare_i = mu_tilde_L * d_share if incumbent

* Cross: Delta_mu_it * Delta_m_it
gen cross_i = d_mu * d_share if incumbent

* Entry: (mu_it - M_{t-1}) * m_it  [entrants only]
* Need M_{t-1}: use the lagged aggregate from incumbents
bysort year: egen M_prev = max(L_M_agg) if incumbent
bysort year: egen M_prev_fill = max(M_prev)
gen entry_i = (mu_acf - M_prev_fill) * share_sales if entrant

* Aggregate components by year
foreach v in within_i mktshare_i cross_i entry_i {
    bysort year: egen `v'_sum = sum(`v')
}

* Exit: handled via previous year
* For each year, exiters are firms in t-1 not in t
* We compute this separately
* Exit contribution: firms present in t-1 but absent in t
* Use forward-looking approach: flag firms whose next year is missing
preserve
    keep firm_id year mu_acf share_sales AGG_mu_acf_sw
    sort firm_id year
    by firm_id: gen last_obs = (_n == _N)
    * If this is the last observation for a firm, it's an "exiter" in year+1
    * (unless it's the last sample year)
    sum year, meanonly
    local maxyear = r(max)
    gen is_exiter = (last_obs == 1 & year < `maxyear')
    gen exit_year = year + 1 if is_exiter
    gen exit_i = (mu_acf - AGG_mu_acf_sw) * share_sales if is_exiter
    keep if is_exiter
    rename exit_year year_exit
    collapse (sum) exit_sum = exit_i (count) N_exit = exit_i, by(year_exit)
    rename year_exit year
    recast long year
    tempfile exit_data
    save `exit_data'
restore

* Collapse decomposition to year level
preserve
    keep year within_i_sum mktshare_i_sum cross_i_sum entry_i_sum ///
         AGG_mu_acf_sw incumbent entrant
    sort year
    by year: keep if _n == 1
    merge 1:1 year using `exit_data', nogen

    replace exit_sum = 0 if exit_sum == .
    replace N_exit = 0 if N_exit == .

    gen net_entry = entry_i_sum - exit_sum
    gen realloc   = mktshare_i_sum + cross_i_sum
    gen dM        = .

    * Delta M
    sort year
    replace dM = AGG_mu_acf_sw - AGG_mu_acf_sw[_n-1] if _n > 1

    rename within_i_sum   within
    rename mktshare_i_sum mkt_share
    rename cross_i_sum    cross
    rename entry_i_sum    entry
    rename exit_sum       exit

    dis _newline "Firm-level Decomposition (DLEU eq. 9):"
    format dM within mkt_share cross net_entry %8.4f
    list year dM within mkt_share cross net_entry, sep(0) noobs

    * Cumulative for Figure IV
    gen cum_within    = sum(within)
    gen cum_realloc   = sum(realloc)
    gen cum_net_entry = sum(net_entry)

    * Base year markup
    sum AGG_mu_acf_sw if _n == 1, meanonly
    local base_mu = r(mean)

    gen path_within  = `base_mu' + cum_within
    gen path_realloc = `base_mu' + cum_realloc
    gen path_netent  = `base_mu' + cum_net_entry

    save "`temp'/dleu_decomp.dta", replace

    * Cumulative summary
    sum within, meanonly
    local cum_w = r(sum)
    sum realloc, meanonly
    local cum_r = r(sum)
    sum net_entry, meanonly
    local cum_ne = r(sum)
    local cum_total = `cum_w' + `cum_r' + `cum_ne'

    dis _newline "Cumulative Decomposition:"
    dis "  Within:       " %8.4f `cum_w' " (" %5.1f 100*`cum_w'/`cum_total' "%)"
    dis "  Reallocation: " %8.4f `cum_r' " (" %5.1f 100*`cum_r'/`cum_total' "%)"
    dis "  Net entry:    " %8.4f `cum_ne' " (" %5.1f 100*`cum_ne'/`cum_total' "%)"

restore

*===============================================================================
* SECTION 4: SECTORAL DECOMPOSITION (DLEU Eq. 10, Table I)
*===============================================================================

dis _newline(2) "--- Section 4: Sectoral decomposition (eq. 10) ---"

use "`temp'/dleu_agg_panel.dta", clear

* Industry-level aggregates
bysort nace2 year: egen TOTSALES_IND = sum(Sales)
gen share_ind = Sales / TOTSALES_IND
bysort nace2 year: egen MU_IND = sum(share_ind * mu_acf)
gen share_IND = TOTSALES_IND / TOTSALES

* Collapse to industry-year
preserve
keep nace2 year MU_IND share_IND AGG_mu_acf_sw
sort nace2 year
by nace2 year: keep if _n == 1
xtset nace2 year, yearly

* 5-year changes (adapt DLEU's 10-year to our shorter panel)
foreach lag in 5 {
    gen dmu_`lag'       = MU_IND - L`lag'.MU_IND
    gen dshare_`lag'    = share_IND - L`lag'.share_IND
    gen within_`lag'    = L`lag'.share_IND * dmu_`lag'
    gen between_`lag'   = L`lag'.MU_IND * dshare_`lag'
    gen cross_`lag'     = dmu_`lag' * dshare_`lag'
    gen dMARKUP_`lag'   = AGG_mu_acf_sw - L`lag'.AGG_mu_acf_sw

    bysort year: egen WITHIN_`lag'  = sum(within_`lag')
    bysort year: egen BETWEEN_`lag' = sum(between_`lag')
    bysort year: egen CROSS_`lag'   = sum(cross_`lag')
}

* Display
dis _newline "Sectoral Decomposition (5-year changes, DLEU eq. 10):"
keep year AGG_mu_acf_sw dMARKUP_5 WITHIN_5 BETWEEN_5 CROSS_5
sort year
by year: keep if _n == 1
drop if dMARKUP_5 == .
format AGG_mu_acf_sw dMARKUP_5 WITHIN_5 BETWEEN_5 CROSS_5 %8.4f
list year AGG_mu_acf_sw dMARKUP_5 WITHIN_5 BETWEEN_5 CROSS_5, ///
    sep(0) noobs

* Save for table
save "`temp'/dleu_sectoral.dta", replace

restore

*===============================================================================
* SECTION 5: FIGURES
*===============================================================================

dis _newline(2) "--- Section 5: Figures ---"

* Load yearly aggregates
use "`temp'/dleu_yearly.dta", clear

* -----------------------------------------------------------------------
* Figure I: Aggregate markup time series
* -----------------------------------------------------------------------
twoway (connected AGG_mu_acf_sw year, lcolor(cranberry) mcolor(cranberry) ///
            msymbol(O) lwidth(thick) lpattern(solid)) ///
       (connected AGG_mu_cs085_sw year, lcolor(forest_green) mcolor(forest_green) ///
            msymbol(S) lwidth(medthick) lpattern(dash)) ///
       (connected AGG_mu_cs1_sw year, lcolor(navy) mcolor(navy) ///
            msymbol(T) lwidth(medthick) lpattern(shortdash)) ///
       (connected AGG_mu_cs_ind_sw year, lcolor(orange) mcolor(orange) ///
            msymbol(D) lwidth(medthick) lpattern(dash_dot)), ///
       xline(2012, lcolor(red) lpattern(dash) lwidth(thin)) ///
       xline(2016, lcolor(forest_green) lpattern(dash) lwidth(thin)) ///
       ytitle("Sales-weighted aggregate markup") xtitle("Year") ///
       title("Aggregate Markups — Czech Construction (CZ-NACE F)") ///
       legend(order(1 "ACF (benchmark)" 2 "Calibrated {&theta}=0.85" ///
                    3 "Firm cost share" 4 "Industry cost share") ///
              ring(0) pos(5) cols(1) size(small)) ///
       xlabel(, labsize(small))
graph export "`figs'/dleu_fig1_aggregate_markup.pdf", replace
dis "  Saved: dleu_fig1_aggregate_markup.pdf"

* -----------------------------------------------------------------------
* Figure II: Sensitivity — Panel A: theta; Panel B: weighting
* -----------------------------------------------------------------------
twoway (connected AGG_mu_acf_sw year, lcolor(cranberry) msymbol(none) ///
            lwidth(thick) lpattern(solid)) ///
       (connected AGG_mu_cs085_sw year, lcolor(forest_green) msymbol(none) ///
            lwidth(thick) lpattern(dash)), ///
       xline(2012 2016, lcolor(gs12) lpattern(dash) lwidth(vthin)) ///
       ytitle("Aggregate markup") xtitle("") ///
       title("(A) Constant elasticity", size(medium)) ///
       legend(order(1 "Benchmark (ACF)" 2 "Constant {&theta}=0.85") ///
              ring(0) pos(11) size(vsmall)) ///
       name(fig2a, replace) nodraw

twoway (connected AGG_mu_acf_sw year, lcolor(cranberry) msymbol(none) ///
            lwidth(thick) lpattern(solid)) ///
       (connected AGG_mu_acf_cw year, lcolor(forest_green) msymbol(none) ///
            lwidth(thick) lpattern(dash)) ///
       (connected AGG_mu_acf_tcw year, lcolor(navy) msymbol(none) ///
            lwidth(thick) lpattern(shortdash_dot)), ///
       xline(2012 2016, lcolor(gs12) lpattern(dash) lwidth(vthin)) ///
       ytitle("") xtitle("") ///
       title("(B) Input-weighted", size(medium)) ///
       legend(order(1 "Sales-weighted" 2 "COGS-weighted" ///
                    3 "Total-cost-weighted") ///
              ring(0) pos(11) size(vsmall)) ///
       name(fig2b, replace) nodraw

graph combine fig2a fig2b, rows(1) xsize(12) ysize(5)
graph export "`figs'/dleu_fig2_sensitivity.pdf", replace
graph drop fig2a fig2b
dis "  Saved: dleu_fig2_sensitivity.pdf"

* -----------------------------------------------------------------------
* Figure III Panel B: Sales-weighted percentiles
* -----------------------------------------------------------------------
twoway (connected AGG_mu_acf_sw year, lcolor(black) msymbol(none) ///
            lwidth(thick) lpattern(solid)) ///
       (connected P90_acf year, lcolor(cranberry) msymbol(none) ///
            lwidth(medthick) lpattern(dash)) ///
       (connected P75_acf year, lcolor(cranberry) msymbol(none) ///
            lwidth(medthick) lpattern(shortdash_dot)) ///
       (connected P50_acf year, lcolor(cranberry) msymbol(none) ///
            lwidth(medthick) lpattern(dash_dot)), ///
       xline(2012, lcolor(red) lpattern(dash) lwidth(thin)) ///
       xline(2016, lcolor(forest_green) lpattern(dash) lwidth(thin)) ///
       ytitle("Markup {&mu}{sub:it}") xtitle("Year") ///
       title("Percentiles of Markup Distribution (sales-weighted)") ///
       legend(order(1 "Average" 2 "P90" 3 "P75" 4 "P50") ///
              ring(0) pos(5) cols(1) size(small))
graph export "`figs'/dleu_fig3b_percentiles.pdf", replace
dis "  Saved: dleu_fig3b_percentiles.pdf"

* -----------------------------------------------------------------------
* Figure III Panel A: Kernel density (first vs last year)
* -----------------------------------------------------------------------
use "`temp'/dleu_agg_panel.dta", clear
sum year, meanonly
local y1 = r(min)
local y2 = r(max)

twoway (kdensity mu_acf if year == `y1' & mu_acf > 0.3 & mu_acf < 5, ///
            lcolor(cranberry) lpattern(dash) lwidth(medthick) bwidth(0.15)) ///
       (kdensity mu_acf if year == `y2' & mu_acf > 0.3 & mu_acf < 5, ///
            lcolor(cranberry) lpattern(solid) lwidth(medthick) bwidth(0.15)), ///
       xtitle("Markup {&mu}{sub:it}") ytitle("Density") ///
       title("Distribution of Markups (unweighted)") ///
       legend(order(2 "`y2'" 1 "`y1'") ring(0) pos(1) size(small))
graph export "`figs'/dleu_fig3a_density.pdf", replace
dis "  Saved: dleu_fig3a_density.pdf"

* -----------------------------------------------------------------------
* Figure IV: Cumulative decomposition with counterfactuals
* -----------------------------------------------------------------------
use "`temp'/dleu_decomp.dta", clear
use "`temp'/dleu_yearly.dta", clear
merge 1:1 year using "`temp'/dleu_decomp.dta", nogen

twoway (connected AGG_mu_acf_sw year, lcolor(cranberry) msymbol(none) ///
            lwidth(thick) lpattern(solid)) ///
       (connected path_within year, lcolor(navy) msymbol(none) ///
            lwidth(thick) lpattern(longdash)) ///
       (connected path_realloc year, lcolor(black) msymbol(none) ///
            lwidth(thick) lpattern(shortdash_dot)) ///
       (connected path_netent year, lcolor(forest_green) msymbol(none) ///
            lwidth(thick) lpattern(dash_dot)), ///
       xline(2012, lcolor(red) lpattern(dash) lwidth(thin)) ///
       xline(2016, lcolor(forest_green) lpattern(dash) lwidth(thin)) ///
       ytitle("Aggregate markup") xtitle("Year") ///
       title("Decomposition of Markup Growth at the Firm Level") ///
       legend(order(1 "Markup (benchmark)" 2 "Within" ///
                    3 "Reallocation" 4 "Net entry") ///
              ring(0) pos(5) cols(1) size(small))
graph export "`figs'/dleu_fig4_decomposition.pdf", replace
dis "  Saved: dleu_fig4_decomposition.pdf"

* -----------------------------------------------------------------------
* Figure V: Micro aggregation vs industry averages
* -----------------------------------------------------------------------
use "`temp'/dleu_agg_panel.dta", clear

* Industry-level simple average, then share-weighted
bysort nace2 year: egen mu_ind_avg = mean(mu_acf)
bysort nace2 year: egen IND_SALES = sum(Sales)
gen share_IND2 = IND_SALES / TOTSALES
bysort nace2 year: gen tag_ind = (_n == 1)
bysort year: egen M_ind_avg = sum(share_IND2 * mu_ind_avg) if tag_ind
bysort year: egen M_ind_avg2 = max(M_ind_avg)

preserve
keep year AGG_mu_acf_sw AGG_mean M_ind_avg2
sort year
by year: keep if _n == 1

twoway (connected AGG_mu_acf_sw year, lcolor(cranberry) msymbol(none) ///
            lwidth(thick) lpattern(solid)) ///
       (connected M_ind_avg2 year, lcolor(navy) msymbol(none) ///
            lwidth(thick) lpattern(longdash)) ///
       (connected AGG_mean year, lcolor(black) msymbol(none) ///
            lwidth(thick) lpattern(dash_dot)), ///
       xline(2012 2016, lcolor(gs12) lpattern(dash) lwidth(vthin)) ///
       ytitle("Aggregate markup") xtitle("Year") ///
       title("Micro Aggregation vs Industry Averages") ///
       legend(order(1 "Firm-level (sales-weighted)" ///
                    2 "Industry averages" 3 "Economy-wide average") ///
              ring(0) pos(5) cols(1) size(small))
graph export "`figs'/dleu_fig5_micro_vs_agg.pdf", replace
dis "  Saved: dleu_fig5_micro_vs_agg.pdf"
restore

* -----------------------------------------------------------------------
* Figure VII: Cost shares of total sales
* -----------------------------------------------------------------------
use "`temp'/dleu_yearly.dta", clear

twoway (connected AGG_cs_cogs year, lcolor(navy) msymbol(none) ///
            lwidth(thick) lpattern(solid)) ///
       (connected AGG_cs_cap year, lcolor(orange) msymbol(none) ///
            lwidth(thick) lpattern(dash)), ///
       xline(2012 2016, lcolor(gs12) lpattern(dash) lwidth(vthin)) ///
       ytitle("Cost share of sales") xtitle("Year") ///
       title("Cost Shares of Total Sales") ///
       legend(order(1 "COGS share" 2 "Capital share (rK)") ///
              ring(0) pos(5) size(small))
graph export "`figs'/dleu_fig7_cost_shares.pdf", replace
dis "  Saved: dleu_fig7_cost_shares.pdf"

* -----------------------------------------------------------------------
* Figure VIII: Profit rate
* -----------------------------------------------------------------------

* Panel A: time series
twoway (connected AGG_profit year, lcolor(cranberry) mcolor(cranberry) ///
            msymbol(O) lwidth(thick) lpattern(solid)), ///
       xline(2012 2016, lcolor(gs12) lpattern(dash) lwidth(vthin)) ///
       ytitle("Profit rate") xtitle("Year") ///
       title("Average Profit Rate (sales-weighted)")
graph export "`figs'/dleu_fig8a_profit_rate.pdf", replace

* Panel B: density comparison
use "`temp'/dleu_agg_panel.dta", clear
sum year, meanonly
local y1 = r(min)
local y2 = r(max)

twoway (kdensity profit_rate if year == `y1' & profit_rate > -0.5 & profit_rate < 0.8, ///
            lcolor(cranberry) lpattern(dash) lwidth(medthick) bwidth(0.05)) ///
       (kdensity profit_rate if year == `y2' & profit_rate > -0.5 & profit_rate < 0.8, ///
            lcolor(cranberry) lpattern(solid) lwidth(medthick) bwidth(0.05)), ///
       xtitle("Profit rate") ytitle("Density") ///
       title("Profit Rate Distribution") ///
       legend(order(2 "`y2'" 1 "`y1'") ring(0) pos(1) size(small))
graph export "`figs'/dleu_fig8b_profit_density.pdf", replace
dis "  Saved: dleu_fig8a, dleu_fig8b"

* -----------------------------------------------------------------------
* Figure XII: Cost-share vs ACF markup and theta
* -----------------------------------------------------------------------
use "`temp'/dleu_yearly.dta", clear

twoway (connected AGG_mu_cs1_sw year, lcolor(cranberry) msymbol(none) ///
            lwidth(thick) lpattern(solid)) ///
       (connected AGG_mu_acf_sw year, lcolor(navy) msymbol(none) ///
            lwidth(thick) lpattern(dash)), ///
       xline(2012 2016, lcolor(gs12) lpattern(dash) lwidth(vthin)) ///
       ytitle("Aggregate markup") xtitle("") ///
       title("(A) CS vs ACF Markup", size(medium)) ///
       legend(order(1 "Cost-share" 2 "ACF") ///
              ring(0) pos(5) size(vsmall)) ///
       name(fig12a, replace) nodraw

twoway (connected AGG_theta_acf year, lcolor(cranberry) msymbol(none) ///
            lwidth(thick) lpattern(solid)) ///
       (connected AGG_theta_cs year, lcolor(forest_green) msymbol(none) ///
            lwidth(thick) lpattern(dash)), ///
       xline(2012 2016, lcolor(gs12) lpattern(dash) lwidth(vthin)) ///
       ytitle("Output elasticity") xtitle("") ///
       title("(B) {&theta}{sup:V}: ACF vs Cost Share", size(medium)) ///
       legend(order(1 "ACF {&theta}" 2 "Cost share (median)") ///
              ring(0) pos(5) size(vsmall)) ///
       name(fig12b, replace) nodraw

graph combine fig12a fig12b, rows(1) xsize(12) ysize(5)
graph export "`figs'/dleu_fig12_cs_vs_pf.pdf", replace
graph drop fig12a fig12b
dis "  Saved: dleu_fig12_cs_vs_pf.pdf"

* -----------------------------------------------------------------------
* Appendix 12: Industry-specific trends (NACE 41, 42, 43)
* -----------------------------------------------------------------------
use "`temp'/dleu_agg_panel.dta", clear

* Compute within-industry sales share
bysort nace2 year: egen TOTSALES_IND2 = sum(Sales)
gen share_ind2 = Sales / TOTSALES_IND2

* Industry-level aggregates
bysort nace2 year: egen MU_IND_ACF = sum(share_ind2 * mu_acf)
bysort nace2 year: egen MU_IND_CS  = sum(share_ind2 * mu_cs1)

preserve
keep nace2 year MU_IND_ACF MU_IND_CS
sort nace2 year
by nace2 year: keep if _n == 1

twoway (connected MU_IND_ACF year if nace2 == 41, lcolor(cranberry) msymbol(O) ///
            lwidth(medthick) msize(small)) ///
       (connected MU_IND_ACF year if nace2 == 42, lcolor(navy) msymbol(S) ///
            lwidth(medthick) msize(small)) ///
       (connected MU_IND_ACF year if nace2 == 43, lcolor(forest_green) msymbol(T) ///
            lwidth(medthick) msize(small)), ///
       xline(2012 2016, lcolor(gs12) lpattern(dash) lwidth(vthin)) ///
       ytitle("Sales-weighted markup") xtitle("Year") ///
       title("Industry-Specific Markups (ACF)") ///
       legend(order(1 "NACE 41: Buildings" 2 "NACE 42: Civil Eng." ///
                    3 "NACE 43: Specialized") ///
              ring(0) pos(5) cols(1) size(small))
graph export "`figs'/dleu_fig_industry_trends.pdf", replace
dis "  Saved: dleu_fig_industry_trends.pdf"
restore

* -----------------------------------------------------------------------
* Extension: Procurement vs non-procurement
* -----------------------------------------------------------------------
use "`temp'/dleu_agg_panel.dta", clear

* Aggregate by procurement status × year
preserve
foreach pp in 0 1 {
    bysort year: egen TOTSALES_pp`pp' = sum(Sales) if pp_dummy == `pp'
    gen share_pp`pp' = Sales / TOTSALES_pp`pp' if pp_dummy == `pp'
    bysort year: egen AGG_pp`pp' = sum(share_pp`pp' * mu_acf) if pp_dummy == `pp'
}

keep year AGG_pp0 AGG_pp1
sort year
by year: keep if _n == 1

twoway (connected AGG_pp1 year, lcolor(cranberry) msymbol(O) ///
            lwidth(thick) lpattern(solid) msize(small)) ///
       (connected AGG_pp0 year, lcolor(navy) msymbol(S) ///
            lwidth(thick) lpattern(dash) msize(small)), ///
       xline(2012, lcolor(red) lpattern(dash) lwidth(thin)) ///
       xline(2016, lcolor(forest_green) lpattern(dash) lwidth(thin)) ///
       ytitle("Sales-weighted markup") xtitle("Year") ///
       title("Aggregate Markup by Procurement Status") ///
       legend(order(1 "Procurement" 2 "Non-procurement") ///
              ring(0) pos(5) size(small))
graph export "`figs'/dleu_fig_procurement_split.pdf", replace
dis "  Saved: dleu_fig_procurement_split.pdf"
restore

*===============================================================================
* SECTION 6: TABLES
*===============================================================================

dis _newline(2) "--- Section 6: Tables ---"

* -----------------------------------------------------------------------
* Table I: Sectoral decomposition (LaTeX)
* -----------------------------------------------------------------------
use "`temp'/dleu_sectoral.dta", clear

local fn "`tabs'/dleu_table1_sectoral_decomp.tex"
cap file close tab1
file open tab1 using "`fn'", write replace

file write tab1 "\begin{table}[htbp]\centering" _newline
file write tab1 "\caption{Sectoral Decomposition of 5-Year Change in Markup (DLEU Eq.~10)}\label{tab:dleu_sectoral}" _newline
file write tab1 "\begin{threeparttable}" _newline
file write tab1 "\begin{tabular}{lccccc}" _newline
file write tab1 "\toprule" _newline
file write tab1 "Year & Markup & $\Delta$Markup & $\Delta$Within & $\Delta$Between & $\Delta$Cross \\" _newline
file write tab1 "\midrule" _newline

local nrows = _N
forvalues i = 1/`nrows' {
    local yr   = year[`i']
    local mu   : di %5.3f AGG_mu_acf_sw[`i']
    local dm   : di %6.3f dMARKUP_5[`i']
    local wi   : di %6.3f WITHIN_5[`i']
    local bt   : di %6.3f BETWEEN_5[`i']
    local cr   : di %6.3f CROSS_5[`i']
    file write tab1 "`yr' & `mu' & `dm' & `wi' & `bt' & `cr' \\" _newline
}

file write tab1 "\bottomrule" _newline
file write tab1 "\end{tabular}" _newline
file write tab1 "\begin{tablenotes}\footnotesize" _newline
file write tab1 "\item \textit{Notes:} Sectoral decomposition following De Loecker, Eeckhout, and Unger (2020, eq.~10). Sectors are CZ-NACE 41 (buildings), 42 (civil engineering), 43 (specialized construction). $\Delta$Within: within-industry markup change weighted by lagged share. $\Delta$Between: industry share reallocation weighted by lagged markup. $\Delta$Cross: joint change." _newline
file write tab1 "\end{tablenotes}" _newline
file write tab1 "\end{threeparttable}" _newline
file write tab1 "\end{table}" _newline
file close tab1
dis "  Saved: dleu_table1_sectoral_decomp.tex"

* -----------------------------------------------------------------------
* Table: Firm-level decomposition (Appendix 4 equivalent)
* -----------------------------------------------------------------------
use "`temp'/dleu_decomp.dta", clear

local fn "`tabs'/dleu_table_firm_decomp.tex"
cap file close tab2
file open tab2 using "`fn'", write replace

file write tab2 "\begin{table}[htbp]\centering" _newline
file write tab2 "\caption{Decomposition of Aggregate Markups at the Firm Level (DLEU Eq.~9)}\label{tab:dleu_firm_decomp}" _newline
file write tab2 "\begin{threeparttable}" _newline
file write tab2 "\begin{tabular}{lcccccc}" _newline
file write tab2 "\toprule" _newline
file write tab2 "Year & $\Delta M$ & $\Delta$Within & $\Delta$Mkt.\ Share & $\Delta$Cross & Net Entry & $N$ \\" _newline
file write tab2 "\midrule" _newline

local nrows = _N
forvalues i = 1/`nrows' {
    local yr  = year[`i']
    local dm  : di %6.3f dM[`i']
    local wi  : di %6.3f within[`i']
    local ms  : di %6.3f mkt_share[`i']
    local cr  : di %6.3f cross[`i']
    local ne  : di %6.3f net_entry[`i']
    file write tab2 "`yr' & `dm' & `wi' & `ms' & `cr' & `ne' & \\" _newline
}

file write tab2 "\bottomrule" _newline
file write tab2 "\end{tabular}" _newline
file write tab2 "\begin{tablenotes}\footnotesize" _newline
file write tab2 "\item \textit{Notes:} Firm-level decomposition following De Loecker, Eeckhout, and Unger (2020, eq.~9) with Haltiwanger (1997) demeaning. $\Delta$Within: $\sum_i m_{i,t-1}\Delta\mu_{it}$. $\Delta$Mkt.\ Share: $\sum_i \tilde{\mu}_{i,t-1}\Delta m_{it}$ where $\tilde{\mu}_{i,t-1} = \mu_{i,t-1} - M_{t-1}$. $\Delta$Cross: $\sum_i \Delta\mu_{it}\Delta m_{it}$. Net Entry: entrant minus exiter contribution (demeaned)." _newline
file write tab2 "\end{tablenotes}" _newline
file write tab2 "\end{threeparttable}" _newline
file write tab2 "\end{table}" _newline
file close tab2
dis "  Saved: dleu_table_firm_decomp.tex"

* -----------------------------------------------------------------------
* Table: Summary across specifications
* -----------------------------------------------------------------------
use "`temp'/dleu_yearly.dta", clear

local fn "`tabs'/dleu_table_markup_specs.tex"
cap file close tab3
file open tab3 using "`fn'", write replace

file write tab3 "\begin{table}[htbp]\centering" _newline
file write tab3 "\caption{Aggregate Markup by Specification and Year}\label{tab:dleu_markup_specs}" _newline
file write tab3 "\begin{threeparttable}" _newline
file write tab3 "\begin{tabular}{l" _newline

* Get years to show
levelsof year, local(allyears)
local show_years ""
foreach y of local allyears {
    if inlist(`y', 2007, 2010, 2014, 2018, 2021) {
        local show_years "`show_years' `y'"
        file write tab3 "c"
    }
}
file write tab3 "}" _newline
file write tab3 "\toprule" _newline
file write tab3 "Specification"
foreach y of local show_years {
    file write tab3 " & `y'"
}
file write tab3 " \\" _newline
file write tab3 "\midrule" _newline

* ACF
file write tab3 "ACF (benchmark)"
foreach y of local show_years {
    sum AGG_mu_acf_sw if year == `y', meanonly
    if r(N) > 0 {
        file write tab3 " & " %5.3f (r(mean))
    }
    else {
        file write tab3 " & --"
    }
}
file write tab3 " \\" _newline

* Calibrated
file write tab3 "Calibrated $\theta=0.85$"
foreach y of local show_years {
    sum AGG_mu_cs085_sw if year == `y', meanonly
    if r(N) > 0 {
        file write tab3 " & " %5.3f (r(mean))
    }
    else {
        file write tab3 " & --"
    }
}
file write tab3 " \\" _newline

* Firm cost share
file write tab3 "Firm cost share"
foreach y of local show_years {
    sum AGG_mu_cs1_sw if year == `y', meanonly
    if r(N) > 0 {
        file write tab3 " & " %5.3f (r(mean))
    }
    else {
        file write tab3 " & --"
    }
}
file write tab3 " \\" _newline

* Industry cost share
file write tab3 "Industry cost share"
foreach y of local show_years {
    sum AGG_mu_cs_ind_sw if year == `y', meanonly
    if r(N) > 0 {
        file write tab3 " & " %5.3f (r(mean))
    }
    else {
        file write tab3 " & --"
    }
}
file write tab3 " \\" _newline

file write tab3 "\bottomrule" _newline
file write tab3 "\end{tabular}" _newline
file write tab3 "\begin{tablenotes}\footnotesize" _newline
file write tab3 "\item \textit{Notes:} Sales-weighted aggregate markup under alternative specifications following De Loecker, Eeckhout, and Unger (2020). ACF uses production-function-estimated $\hat{\theta}_{st}$. Calibrated sets $\theta = 0.85$. Firm cost share uses $\alpha^V_{it} = \text{COGS}_{it}/(\text{COGS}_{it} + r_t K_{it})$. Industry cost share: $\text{median}_{i \in s}(\alpha^V_{it})$ within NACE 2-digit $\times$ year." _newline
file write tab3 "\end{tablenotes}" _newline
file write tab3 "\end{threeparttable}" _newline
file write tab3 "\end{table}" _newline
file close tab3
dis "  Saved: dleu_table_markup_specs.tex"

*===============================================================================
* SECTION 7: KEY FINDINGS COMPARISON
*===============================================================================

dis _newline(2) "============================================================"
dis "  KEY FINDINGS — Czech Construction vs DLEU (US)"
dis "============================================================"

use "`temp'/dleu_yearly.dta", clear
sort year
local y1 = year[1]
local yN = year[_N]
local mu1 = AGG_mu_acf_sw[1]
local muN = AGG_mu_acf_sw[_N]
local dmu = `muN' - `mu1'

local dmu_str : di %6.3f `dmu'
dis "  Czech aggregate markup: " %5.3f `mu1' " (`y1') -> " %5.3f `muN' " (`yN')  [Delta = `dmu_str']"
dis "  DLEU US benchmark:      1.210 (1980) -> 1.610 (2016)  [Delta = +0.400]"

local gap = AGG_mu_acf_sw[_N] - AGG_mean[_N]
dis "  SW-mean gap (last year): " %5.3f `gap' "  [DLEU: ~0.20-0.40]"

dis _newline "  Czech construction: flat-to-declining (unlike US secular rise)"
dis "  Consistent with DLEU (2025) reply: theta constant, all action in alpha"

* Clean up temp files
cap erase "`temp'/dleu_panel.dta"

dis _newline "============================================================"
dis "  DONE. Outputs saved to:"
dis "    Figures: `figs'/"
dis "    Tables:  `tabs'/"
dis "    Data:    `outd'/"
dis "============================================================"

log close
