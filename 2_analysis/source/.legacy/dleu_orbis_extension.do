*===============================================================================
* dleu_orbis_extension.do
*
* Cross-industry DLEU aggregate markup analysis using Orbis all-industry panel.
* Tests whether the declining Czech construction markup pattern generalizes
* across industries.
*
* Extends dleu_replication.do (MagnusWeb construction) to:
*   - Top-5 procurement industries (NACE 41, 42, 43, 46, 71)
*   - All-industry aggregate
*   - Procurement vs non-procurement split by industry
*   - Sectoral decomposition (DLEU eq. 10) across broad sectors
*
* Data: orbis_panel.dta (1.3M obs, 275K firms, 77 NACE-2d, 2006-2023)
* Markups: cost-share approach (no PF estimation needed)
*   mu = [COGS/(COGS+rK)] * (Sales/COGS)
*
* Usage:
*   cd markups-procurement/2_analysis/source
*   stata-mp -e do dleu_orbis_extension.do
*===============================================================================

clear all
set more off
set scheme s2color

local base  ".."
local input "`base'/input"
local outd  "`base'/output/data"
local figs  "`base'/output/figures"
local tabs  "`base'/output/tables"

local usercost = 0.12

cap log close _all
log using "`outd'/dleu_orbis_extension.log", text replace

dis _newline "============================================================"
dis "  DLEU Cross-Industry Extension — Orbis All-Industry Panel"
dis "============================================================"

*===============================================================================
* SECTION 1: LOAD AND PREPARE
*===============================================================================

use "`input'/orbis_panel.dta", clear

* Firm ID
rename ico id

* Deflated levels already in data: rGO, rCOGS, rK
gen Sales = rGO
gen kexp  = `usercost' * rK
gen totcost = rCOGS + kexp

* Cost-share markup (DLEU mu_1 equivalent)
gen cs_firm = rCOGS / (rCOGS + kexp)
gen mu_cs = cs_firm * (Sales / rCOGS)
label var mu_cs "Cost-share markup"

* Trim extreme markups
drop if mu_cs <= 0 | mu_cs == .
egen mu_p1 = pctile(mu_cs), p(1)
egen mu_p99 = pctile(mu_cs), p(99)
drop if mu_cs < mu_p1 | mu_cs > mu_p99
drop mu_p1 mu_p99

* Panel structure
egen firm_id = group(id)

* Industry labels
gen ind_label = ""
replace ind_label = "Buildings (41)" if nace2 == 41
replace ind_label = "Civil Eng. (42)" if nace2 == 42
replace ind_label = "Spec. Constr. (43)" if nace2 == 43
replace ind_label = "Wholesale (46)" if nace2 == 46
replace ind_label = "Architecture (71)" if nace2 == 71

* Flag top-5 procurement industries
gen top5 = inlist(nace2, 41, 42, 43, 46, 71)

dis "Panel: " _N " obs, " ///
    string(firm_id[_N], "%12.0fc") " firms"

save "`outd'/dleu_orbis_panel.dta", replace

*===============================================================================
* SECTION 2: AGGREGATE MARKUPS BY INDUSTRY
*===============================================================================

dis _newline "--- Aggregate markups by industry ---"

* Within-industry sales weights
bysort nace2 year: egen TOTSALES_IND = sum(Sales)
gen share_ind = Sales / TOTSALES_IND

* Industry-level sales-weighted markup
bysort nace2 year: egen MU_IND = sum(share_ind * mu_cs)
label var MU_IND "Sales-weighted industry markup"

* Economy-wide
bysort year: egen TOTSALES = sum(Sales)
gen share_econ = Sales / TOTSALES

bysort year: egen MU_ALL = sum(share_econ * mu_cs)
label var MU_ALL "Economy-wide aggregate markup"

* Collapse to industry-year
preserve
keep nace2 year MU_IND MU_ALL TOTSALES_IND ind_label top5
sort nace2 year
by nace2 year: keep if _n == 1

* Count firms per industry-year
save "`outd'/dleu_orbis_indyr.dta", replace

* Display top-5 summary
dis _newline "Top-5 Procurement Industries — Aggregate Markup:"
list nace2 year MU_IND if top5 & inlist(year, 2008, 2012, 2016, 2020), ///
    sep(0) noobs

restore

*===============================================================================
* SECTION 3: CROSS-INDUSTRY FIGURE
*===============================================================================

dis _newline "--- Figures ---"

use "`outd'/dleu_orbis_indyr.dta", clear

* Figure: Multi-panel markup trends for top-5 industries
twoway (connected MU_IND year if nace2 == 41, lcolor(cranberry) msymbol(O) ///
            lwidth(medthick) msize(small)) ///
       (connected MU_IND year if nace2 == 42, lcolor(navy) msymbol(S) ///
            lwidth(medthick) msize(small)) ///
       (connected MU_IND year if nace2 == 43, lcolor(forest_green) msymbol(T) ///
            lwidth(medthick) msize(small)) ///
       (connected MU_IND year if nace2 == 46, lcolor(orange) msymbol(D) ///
            lwidth(medthick) msize(small)) ///
       (connected MU_IND year if nace2 == 71, lcolor(purple) msymbol(+) ///
            lwidth(medthick) msize(medlarge)) ///
       (connected MU_ALL year, lcolor(black) msymbol(none) ///
            lwidth(thick) lpattern(dash)), ///
       xline(2012, lcolor(red) lpattern(dash) lwidth(thin)) ///
       xline(2016, lcolor(forest_green) lpattern(dash) lwidth(thin)) ///
       ytitle("Sales-weighted markup (cost-share)") xtitle("Year") ///
       title("Aggregate Markups by Industry — Czech Orbis Panel") ///
       legend(order(1 "Buildings (41)" 2 "Civil Eng. (42)" ///
                    3 "Spec. Constr. (43)" 4 "Wholesale (46)" ///
                    5 "Architecture (71)" 6 "All-industry") ///
              ring(0) pos(1) cols(2) size(vsmall)) ///
       xlabel(2008(2)2022, labsize(small))
graph export "`figs'/orbis_dleu_cross_industry.pdf", replace
dis "  Saved: orbis_dleu_cross_industry.pdf"

* Separate panels per industry
foreach n in 41 42 43 46 71 {
    local lab = ""
    if `n' == 41 local lab "Buildings"
    if `n' == 42 local lab "Civil Engineering"
    if `n' == 43 local lab "Specialized Construction"
    if `n' == 46 local lab "Wholesale Trade"
    if `n' == 71 local lab "Architecture & Engineering"

    twoway (connected MU_IND year if nace2 == `n', lcolor(cranberry) ///
                msymbol(O) lwidth(medthick) msize(small)), ///
           xline(2012 2016, lcolor(gs12) lpattern(dash) lwidth(vthin)) ///
           ytitle("Markup") xtitle("") ///
           title("NACE `n': `lab'", size(medium)) ///
           name(ind`n', replace) nodraw
}
graph combine ind41 ind42 ind43 ind46 ind71, ///
    rows(2) xsize(14) ysize(8) ///
    title("Industry-Specific Markup Trends (Orbis)")
graph export "`figs'/orbis_dleu_industry_panels.pdf", replace
graph drop ind41 ind42 ind43 ind46 ind71
dis "  Saved: orbis_dleu_industry_panels.pdf"

*===============================================================================
* SECTION 4: PROCUREMENT SPLIT BY INDUSTRY
*===============================================================================

dis _newline "--- Procurement split by industry ---"

use "`outd'/dleu_orbis_panel.dta", clear

* Aggregate by pp_dummy × nace2 × year
preserve
foreach pp in 0 1 {
    bysort nace2 year: egen TOTSALES_pp`pp' = sum(Sales) if pp_dummy == `pp'
    gen share_pp`pp' = Sales / TOTSALES_pp`pp' if pp_dummy == `pp'
    bysort nace2 year: egen MU_pp`pp' = sum(share_pp`pp' * mu_cs) if pp_dummy == `pp'
}

keep nace2 year MU_pp0 MU_pp1
sort nace2 year
by nace2 year: keep if _n == 1

* Construction subsector comparison
twoway (connected MU_pp1 year if nace2 == 41, lcolor(cranberry) msymbol(O) ///
            lwidth(medthick) msize(small)) ///
       (connected MU_pp0 year if nace2 == 41, lcolor(navy) msymbol(S) ///
            lwidth(medthick) lpattern(dash) msize(small)) ///
       (connected MU_pp1 year if nace2 == 43, lcolor(cranberry) msymbol(T) ///
            lwidth(thin) msize(small)) ///
       (connected MU_pp0 year if nace2 == 43, lcolor(navy) msymbol(D) ///
            lwidth(thin) lpattern(dash) msize(small)), ///
       xline(2012, lcolor(red) lpattern(dash) lwidth(thin)) ///
       xline(2016, lcolor(forest_green) lpattern(dash) lwidth(thin)) ///
       ytitle("Sales-weighted markup") xtitle("Year") ///
       title("Markup by Procurement Status (Orbis)") ///
       legend(order(1 "NACE 41 PP" 2 "NACE 41 non-PP" ///
                    3 "NACE 43 PP" 4 "NACE 43 non-PP") ///
              ring(0) pos(1) cols(2) size(vsmall))
graph export "`figs'/orbis_procurement_split_by_nace.pdf", replace
dis "  Saved: orbis_procurement_split_by_nace.pdf"

restore

*===============================================================================
* SECTION 5: SECTORAL DECOMPOSITION (DLEU eq. 10)
*===============================================================================

dis _newline "--- Sectoral decomposition ---"

use "`outd'/dleu_orbis_indyr.dta", clear

* Sector share of total economy
bysort year: egen TOTSALES_YR = sum(TOTSALES_IND)
gen share_s = TOTSALES_IND / TOTSALES_YR

xtset nace2 year, yearly

* 5-year changes
gen dmu    = MU_IND - L5.MU_IND
gen dshare = share_s - L5.share_s
gen within_5  = L5.share_s * dmu
gen between_5 = L5.MU_IND * dshare
gen cross_5   = dmu * dshare

* Aggregate decomposition
bysort year: egen WITHIN  = sum(within_5) if top5
bysort year: egen BETWEEN = sum(between_5) if top5
bysort year: egen CROSS   = sum(cross_5) if top5

preserve
keep if top5
keep year MU_ALL WITHIN BETWEEN CROSS
sort year
by year: keep if _n == 1
drop if WITHIN == .

dis _newline "Sectoral Decomposition (top-5, 5-year changes):"
format MU_ALL WITHIN BETWEEN CROSS %8.4f
list year MU_ALL WITHIN BETWEEN CROSS, sep(0) noobs

* Save
save "`outd'/dleu_orbis_sectoral.dta", replace
restore

*===============================================================================
* SECTION 6: LATEX TABLE
*===============================================================================

dis _newline "--- Tables ---"

use "`outd'/dleu_orbis_indyr.dta", clear

* Table: markup by industry and selected years
local fn "`tabs'/orbis_dleu_cross_industry.tex"
cap file close tab
file open tab using "`fn'", write replace

file write tab "\begin{table}[htbp]\centering" _newline
file write tab "\caption{Aggregate Markup by Industry (Orbis Panel)}\label{tab:orbis_industry}" _newline
file write tab "\begin{threeparttable}" _newline
file write tab "\begin{tabular}{llcccc}" _newline
file write tab "\toprule" _newline
file write tab "NACE & Industry & 2008 & 2012 & 2016 & 2020 \\" _newline
file write tab "\midrule" _newline

foreach n in 41 42 43 46 71 {
    local lab = ""
    if `n' == 41 local lab "Buildings"
    if `n' == 42 local lab "Civil Engineering"
    if `n' == 43 local lab "Specialized Constr."
    if `n' == 46 local lab "Wholesale Trade"
    if `n' == 71 local lab "Architecture \& Eng."

    file write tab "`n' & `lab'"
    foreach y in 2008 2012 2016 2020 {
        sum MU_IND if nace2 == `n' & year == `y', meanonly
        if r(N) > 0 {
            local val : di %5.3f r(mean)
            file write tab " & `val'"
        }
        else {
            file write tab " & ---"
        }
    }
    file write tab " \\" _newline
}

file write tab "\midrule" _newline
file write tab " & All industries"
foreach y in 2008 2012 2016 2020 {
    sum MU_ALL if year == `y', meanonly
    if r(N) > 0 {
        local val : di %5.3f r(mean)
        file write tab " & `val'"
    }
}
file write tab " \\" _newline

file write tab "\bottomrule" _newline
file write tab "\end{tabular}" _newline
file write tab "\begin{tablenotes}\footnotesize" _newline
file write tab `"\item \textit{Notes:} Sales-weighted aggregate markup by NACE 2-digit industry using cost-share output elasticity \$\alpha^V_{it} = \text{COGS}_{it}/(\text{COGS}_{it} + r_t K_{it})\$. Data: Orbis Bureau van Dijk via WRDS, 275,000 Czech firms, 2006--2023. User cost of capital \$r = 0.12\$."' _newline
file write tab "\end{tablenotes}" _newline
file write tab "\end{threeparttable}" _newline
file write tab "\end{table}" _newline
file close tab
dis "  Saved: orbis_dleu_cross_industry.tex"

* Also save CSV for Python companion
use "`outd'/dleu_orbis_indyr.dta", clear
keep if top5
keep nace2 year MU_IND MU_ALL ind_label
sort nace2 year
export delimited "`outd'/orbis_dleu_cross_industry.csv", replace
dis "  Saved: orbis_dleu_cross_industry.csv"

*===============================================================================
* DONE
*===============================================================================

dis _newline "============================================================"
dis "  DONE. Outputs:"
dis "    Figures: orbis_dleu_cross_industry.pdf"
dis "            orbis_dleu_industry_panels.pdf"
dis "            orbis_procurement_split_by_nace.pdf"
dis "    Tables:  orbis_dleu_cross_industry.tex"
dis "    Data:    orbis_dleu_cross_industry.csv"
dis "============================================================"

log close
