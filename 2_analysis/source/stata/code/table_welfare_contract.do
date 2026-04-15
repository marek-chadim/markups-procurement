*===============================================================================
* table_welfare_contract.do — Contract-Level Pass-Through (§7 Table 23)
*
* Port of contract_level_welfare.py. Regresses contract-level Rel_Price on
* firm-level markup with firm FE, year FE, and within-firm demeaned markup.
* Tests whether firm markup premium translates to contract-level overbidding.
*
* Input:  $data/contract_level.dta (built by Python pipeline)
*         OR $data/analysis_panel.dta merged with tender register
* Output: ../../output/tables/contract_level_welfare.tex
*===============================================================================

dis _newline "--- table_welfare_contract.do ---"

cap confirm file "$data/contract_level.dta"
if _rc != 0 {
    dis "  SKIP: contract_level.dta not available"
    dis "  (Requires Python contract_level_welfare.py to build the firm-contract panel)"
    exit
}

use "$data/contract_level.dta", clear
dis "  Contracts: " _N

* Check variables exist
foreach v in rel_price firm_markup firm_id year nace2 {
    cap confirm var `v'
    if _rc != 0 {
        dis "  SKIP: missing required variable `v'"
        exit
    }
}

* Match Python specification exactly: log(rel_price) regressed on log(firm_markup)
keep if rel_price > 0 & firm_markup > 0
gen log_rel_price = log(rel_price)
gen log_firm_markup = log(firm_markup)

* Cross-firm OLS (year+NACE FE, no firm FE)
eststo clear
eststo crossfirm: reghdfe log_rel_price log_firm_markup, absorb(year nace2) vce(cluster firm_id)
local beta_cross = _b[log_firm_markup]
local se_cross = _se[log_firm_markup]

* Within-firm OLS (firm+year FE) — the key specification
eststo withinfirm: reghdfe log_rel_price log_firm_markup, absorb(firm_id year) vce(cluster firm_id)
local beta_within = _b[log_firm_markup]
local se_within = _se[log_firm_markup]
local n_within = e(N)
local n_firms = e(N_clust)

* NACE heterogeneity
foreach n in 41 42 43 {
    cap eststo nace`n': reghdfe log_rel_price log_firm_markup if nace2 == `n', ///
        absorb(firm_id year) vce(cluster firm_id)
    if _rc == 0 {
        local beta_`n' = _b[log_firm_markup]
        local se_`n' = _se[log_firm_markup]
    }
}

* Pre/post 2012 reform
foreach era in "pre" "post" {
    if "`era'" == "pre" local cond "year < 2012"
    if "`era'" == "post" local cond "year >= 2012"
    cap eststo era_`era': reghdfe log_rel_price log_firm_markup if `cond', ///
        absorb(firm_id year) vce(cluster firm_id)
    if _rc == 0 {
        local beta_`era' = _b[log_firm_markup]
        local se_`era' = _se[log_firm_markup]
    }
}

* Write LaTeX table — distinct filename so Python's canonical
* contract_level_welfare.tex (used by the paper prose claims) is not
* overwritten. Stata's version is a cross-check written to
* contract_level_welfare_stata.tex.
cap mkdir "$output/tables"
tempname tf
file open `tf' using "../../output/tables/contract_level_welfare_stata.tex", write replace
#delimit ;
file write `tf'
    "\begin{table}[htbp]\centering" _n
    "\caption{Contract-Level Pass-Through of Firm Markup Premium}" _n
    "\label{tab:contract_welfare}" _n
    "\begin{threeparttable}" _n
    "\begin{tabular}{lcc}" _n
    "\toprule" _n
    " & (1) Cross-firm & (2) Within-firm \\" _n
    " & year + NACE FE & firm + year FE \\" _n
    "\midrule" _n
;
#delimit cr
file write `tf' "Firm markup $\mu_{it}$"
file write `tf' " & " %9.3f (`beta_cross') " & " %9.3f (`beta_within') " \\" _n
file write `tf' " "
file write `tf' " & (" %6.3f (`se_cross') ") & (" %6.3f (`se_within') ") \\" _n
file write `tf' "N" " & " %9.0fc (`n_within') " & " %9.0fc (`n_within') " \\" _n

#delimit ;
file write `tf'
    "\bottomrule" _n
    "\end{tabular}" _n
    "\begin{tablenotes}\footnotesize" _n
    "\item \emph{Notes:} Dependent variable is the contract-level ratio "
    "$\text{bid\_final\_price} / \text{engineer\_estimate}$. "
    "Column (1) absorbs year and NACE fixed effects. Column (2) absorbs "
    "firm and year fixed effects, isolating the within-firm correlation "
    "between the firm's ACF markup and contract-level Rel\_Price. "
    "Standard errors clustered by firm." _n
    "\end{tablenotes}" _n
    "\end{threeparttable}" _n
    "\end{table}" _n
;
#delimit cr
file close `tf'
dis "  Saved: contract_level_welfare.tex"
dis "  Within-firm β = " %6.3f (`beta_within') " (SE " %6.3f (`se_within') ")"
