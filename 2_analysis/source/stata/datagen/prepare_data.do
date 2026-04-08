*===============================================================================
* prepare_data.do — Prepare analysis datasets
*
* Analog of DGM prepare_datasets.do:
*   DGM: merge FICUS-FARE with EAP price data, generate s/v/k/m/o
*   Us:  load rebuilt data, generate analysis variables, winsorize
*
* Input:  data_rebuilt.dta (from rebuild_data.py)
* Output: analysisdata/analysis_panel.dta
*===============================================================================

dis _newline "--- prepare_data.do ---"

use "$srcdata", clear
xtset id year, yearly

dis "Source data: " _N " obs, " year[1] "-" year[_N]
tab nace2

*-----------------------------------------------------------------------
* Variable construction (analog of DGM variable definitions)
*-----------------------------------------------------------------------

* DGM variables → our analogs:
*   s = ln(catotalR)  →  go  (already in data, = ln deflated sales)
*   v = ln(salR)      →  w   (ln wages, if available) or le (ln employees)
*   k = ln(immocorR)  →  k   (already in data, = ln deflated fixed assets)
*   m = ln(acha4R)    →  cogs (already in data, = ln deflated COGS)
*   o = ln(autachaR)  →  o   (ln overhead, if available)
*   ratio = catotal/acha4 → ratio = exp(go)/exp(cogs)

* Expenditure ratio (for markup calculation: μ = θ/α, α = exp(cogs)/exp(go))
cap drop ratio
gen ratio = exp(go) / exp(cogs)
label var ratio "Sales/COGS ratio (= 1/alpha)"

* Market share by year × nace2 (analog of DGM ms5d)
cap drop total_sales
cap drop mktshare
bys year nace2: egen total_sales = total(exp(go))
gen mktshare = exp(go) / total_sales
label var mktshare "Market share (year × nace2)"
drop total_sales

* Cost share (for cost-weighted aggregation, DGM Figure 5)
cap drop total_cogs
cap drop costshare
bys year nace2: egen total_cogs = total(exp(cogs))
gen costshare = exp(cogs) / total_cogs
label var costshare "Cost share (year × nace2)"
drop total_cogs

* Age (years since first observed in panel)
cap drop entry_year
cap drop age
bys id: egen entry_year = min(year)
gen age = year - entry_year
label var age "Firm age (years since panel entry)"

* Firm size categories
cap drop size_cat
xtile size_cat = go, nq(4)
label var size_cat "Size quartile (by go)"

*-----------------------------------------------------------------------
* Winsorization (DGM: 2% tails by sector)
*-----------------------------------------------------------------------

* DGM winsorize at 2% by naf2d for s, m, v, k, o, p, q
* We winsorize at 2% by nace2 for go, k, cogs (already done in rebuild_data.py)
* Just verify and add ratio winsorization

foreach var in go k cogs {
    * Check for extreme values
    qui sum `var', d
    dis "  `var': p1=" %6.2f r(p1) " p99=" %6.2f r(p99) " range=" %6.2f r(max)-r(min)
}

* Winsorize ratio at 2% by nace2 (DGM equivalent of their winsor command)
cap drop ratio_raw
gen ratio_raw = ratio
levelsof nace2, local(naces)
foreach n of local naces {
    qui sum ratio if nace2 == `n', d
    local lo = r(p1)
    local hi = r(p99)
    replace ratio = max(min(ratio, `hi'), `lo') if nace2 == `n'
}

*-----------------------------------------------------------------------
* Treatment variables
*-----------------------------------------------------------------------

* Ensure procurement indicators exist
cap confirm var pp_dummy
if _rc {
    dis "ERROR: pp_dummy not found in data"
    error 198
}

* Lagged treatment
cap drop pp_lag
sort id year
by id: gen pp_lag = pp_dummy[_n-1]
label var pp_lag "L.pp_dummy"

* 3-year treatment window
cap drop pp_L1
cap drop pp_L2
cap drop pp_ever_3y
by id: gen pp_L1 = pp_dummy[_n-1]
by id: gen pp_L2 = pp_dummy[_n-2]
gen pp_ever_3y = (pp_dummy == 1 | pp_L1 == 1 | pp_L2 == 1)
replace pp_ever_3y = . if mi(pp_dummy)
label var pp_ever_3y "PP active in t, t-1, or t-2"
drop pp_L1 pp_L2

*-----------------------------------------------------------------------
* Survival probit (CWDL 2015)
*-----------------------------------------------------------------------

cap drop next_yr
cap drop survival_1
cap drop phat_survival
by id: gen next_yr = year[_n+1]
gen survival_1 = (next_yr - year == 1) if next_yr != .
replace survival_1 = 0 if survival_1 == .

xi: logit survival_1 k cogs pp_dummy i.year*i.nace2, iterate(100)
predict phat_survival, pr
label var phat_survival "Predicted survival probability"

dis "  Survival phat: mean=" %6.4f r(mean)

drop next_yr survival_1 _I*

*-----------------------------------------------------------------------
* Save
*-----------------------------------------------------------------------

order id year nace2 go k cogs ratio mktshare costshare ///
    pp_dummy pp_lag pp_ever_3y phat_survival age size_cat

compress
save "$data/analysis_panel.dta", replace
dis "  Saved: analysis_panel.dta (" _N " obs)"
