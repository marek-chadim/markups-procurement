*===============================================================================
* build_contract_level.do — Build contract-level panel for welfare analysis
*
* Port of contract_level_welfare.py merge step. Joins Datlab master tender
* analytics to firm-year markups from markups_panel.dta on (id, year).
*
* Input:
*   ../../../1_data/input/datlab/master_tender_analytics.csv
*   $data/markups_panel.dta (must have mu_A, nace2, pp_dummy)
*
* Output:
*   $data/contract_level.dta  (contract-year panel with firm_markup merged)
*===============================================================================

dis _newline "--- build_contract_level.do ---"

local tenders "../../../1_data/input/datlab/master_tender_analytics.csv"
cap confirm file "`tenders'"
if _rc != 0 {
    dis "  SKIP: tender data not found at `tenders'"
    exit
}

* ---- Load tenders ---------------------------------------------------------
import delimited "`tenders'", clear varn(1) bindquote(strict) stripquote(yes) ///
    case(lower) encoding(UTF-8)
cap destring year, replace force
cap destring lot_estimated_price bid_final_price, replace force

dis "  Loaded " _N " raw contract rows"

* Keep rows with both estimate and final price
keep if !mi(lot_estimated_price) & !mi(bid_final_price) ///
    & lot_estimated_price > 0 & bid_final_price > 0 & !mi(year)
dis "  After estimate+final filter: " _N

gen rel_price = bid_final_price / lot_estimated_price
keep if rel_price >= 0.2 & rel_price <= 5.0
dis "  After winsorization [0.2, 5.0]: " _N

* Cast bidder_id to integer IČO (Czech firm identifier)
cap tostring bidder_id, replace force
gen double id_num = real(bidder_id)
drop if mi(id_num) | id_num <= 0
gen long id = id_num
drop id_num bidder_id
cap destring year, replace force
gen int year_int = year
drop year
rename year_int year

tempfile contracts_tmp
save `contracts_tmp'
dis "  Saved `contracts_tmp' with " _N " contracts"

* ---- Load markups ---------------------------------------------------------
* Prefer Python's paper_markups.dta (which has markup_A — the canonical
* headline ACF translog markup used in the paper prose). Fall back to
* Stata's markups_panel.dta (mu_A) if paper_markups.dta is unavailable.
local pm "../../output/data/paper_markups.dta"
cap confirm file "`pm'"
if _rc == 0 {
    use "`pm'", clear
    cap confirm var markup_A
    if _rc == 0 {
        keep id year markup_A pp_dummy nace2
        rename markup_A firm_markup
    }
    else {
        dis "  paper_markups.dta has no markup_A — falling back to markups_panel.dta"
        use "$data/markups_panel.dta", clear
        keep id year mu_A pp_dummy nace2
        rename mu_A firm_markup
    }
}
else {
    use "$data/markups_panel.dta", clear
    cap confirm var mu_A
    if _rc != 0 {
        dis "  SKIP: no markup source available"
        exit
    }
    keep id year mu_A pp_dummy nace2
    rename mu_A firm_markup
}
cap gen long firm_id = id
duplicates drop id year, force

tempfile markups_tmp
save `markups_tmp'

* ---- Merge ----------------------------------------------------------------
use `contracts_tmp', clear
merge m:1 id year using `markups_tmp', keep(master matched) keepusing(firm_markup pp_dummy nace2 firm_id)
keep if _merge == 3
drop _merge

* NACE F construction only
keep if inlist(nace2, 41, 42, 43)
dis "  After NACE filter (41/42/43): " _N

* Drop exact duplicates (Datlab has some near-identical contract rows)
duplicates drop id year bid_final_price lot_estimated_price, force
dis "  After dedup: " _N

dis "  Contracts merged: " _N "  firms: "
qui distinct id
dis r(ndistinct)

* ---- Save -----------------------------------------------------------------
order id year firm_id firm_markup pp_dummy nace2 rel_price ///
    bid_final_price lot_estimated_price
save "$data/contract_level.dta", replace
dis "  Saved: contract_level.dta (" _N " rows)"
