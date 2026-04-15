*===============================================================================
* table_favoritism.do — Favoritism Decomposition (§6 Table, 6 panels)
*
* Port of favoritism_decomposition.py. Six panels testing whether procurement
* premium is driven by favoritism rather than genuine productivity:
*   A. Baseline premium on pp_dummy
*   B. Controlling for firm fixed effects (within firm)
*   C. Controlling for firm × year FE (identification via crosswalk)
*   D. Single-bidder interaction (if favoritism, effect larger under single bid)
*   E. Pre/post 2012 single-bidding ban (should shrink premium)
*   F. Direct test: Rel_Price control (if Rel_Price fully explains, β→0)
*
* Input:  $data/markups_panel.dta
* Output: ../../output/tables/favoritism_decomposition.tex
*===============================================================================

dis _newline "--- table_favoritism.do ---"

use "$data/markups_panel.dta", clear

cap confirm var mu_A
if _rc != 0 {
    dis "  SKIP: missing mu_A"
    exit
}

gen lmu = log(mu_A)

eststo clear

* Panel A: baseline
eststo panA: reghdfe lmu pp_dummy, absorb(year nace2) vce(cluster id)
local bA = _b[pp_dummy]
local seA = _se[pp_dummy]
local nA = e(N)

* Panel B: firm FE
eststo panB: reghdfe lmu pp_dummy, absorb(id year) vce(cluster id)
local bB = _b[pp_dummy]
local seB = _se[pp_dummy]
local nB = e(N)

* Panel C: firm × year FE (can't identify main effect of pp_dummy, so drop)
* Alternative: firm FE + year × nace2 FE
eststo panC: reghdfe lmu pp_dummy, absorb(id year##nace2) vce(cluster id)
local bC = _b[pp_dummy]
local seC = _se[pp_dummy]
local nC = e(N)

* Panel D: interaction with single-bid share (if available)
cap confirm var single_bid_share
if _rc == 0 {
    gen pp_single = pp_dummy * single_bid_share
    eststo panD: reghdfe lmu pp_dummy pp_single, absorb(id year) vce(cluster id)
    local bD = _b[pp_dummy]
    local seD = _se[pp_dummy]
    local bD_int = _b[pp_single]
    local seD_int = _se[pp_single]
    local nD = e(N)
}

* Panel E: pre/post 2012 single-bid ban
gen post2012 = (year >= 2012)
gen pp_post = pp_dummy * post2012
eststo panE: reghdfe lmu pp_dummy pp_post, absorb(id year) vce(cluster id)
local bE_pre = _b[pp_dummy]
local bE_diff = _b[pp_post]
local seE_pre = _se[pp_dummy]
local seE_diff = _se[pp_post]
local nE = e(N)

* Panel F: direct Rel_Price control (if contract_level data merges in)
cap confirm var rel_price_firm_avg
if _rc == 0 {
    eststo panF: reghdfe lmu pp_dummy rel_price_firm_avg, absorb(id year) vce(cluster id)
    local bF = _b[pp_dummy]
    local seF = _se[pp_dummy]
    local nF = e(N)
}

* Write LaTeX
cap mkdir "$output/tables"
tempname tf
file open `tf' using "../../output/tables/favoritism_decomposition.tex", write replace
#delimit ;
file write `tf'
    "\begin{table}[htbp]\centering" _n
    "\caption{Favoritism Decomposition: 6-Panel Test of Markup Premium}" _n
    "\label{tab:favoritism}" _n
    "\begin{threeparttable}" _n
    "\begin{tabular}{lccc}" _n
    "\toprule" _n
    "Panel & $\hat\beta_{pp}$ & SE & N \\" _n
    "\midrule" _n
;
#delimit cr

file write `tf' "A: Baseline (year + NACE FE)" " & " %7.4f (`bA') " & " %7.4f (`seA') " & " %9.0fc (`nA') " \\" _n
file write `tf' "B: Firm + year FE" " & " %7.4f (`bB') " & " %7.4f (`seB') " & " %9.0fc (`nB') " \\" _n
file write `tf' "C: Firm + year\(\times\)NACE FE" " & " %7.4f (`bC') " & " %7.4f (`seC') " & " %9.0fc (`nC') " \\" _n

cap confirm var single_bid_share
if _rc == 0 {
    file write `tf' "D: + single-bid interaction" " & " %7.4f (`bD') " & " %7.4f (`seD') " & " %9.0fc (`nD') " \\" _n
    file write `tf' "\quad \(pp \times\) single-bid share" " & " %7.4f (`bD_int') " & " %7.4f (`seD_int') " & \\" _n
}

file write `tf' "E: Pre-2012 baseline" " & " %7.4f (`bE_pre') " & " %7.4f (`seE_pre') " & " %9.0fc (`nE') " \\" _n
file write `tf' "\quad $\times$ Post-2012" " & " %7.4f (`bE_diff') " & " %7.4f (`seE_diff') " & \\" _n

cap confirm var rel_price_firm_avg
if _rc == 0 {
    file write `tf' "F: + Rel\_Price firm avg" " & " %7.4f (`bF') " & " %7.4f (`seF') " & " %9.0fc (`nF') " \\" _n
}

#delimit ;
file write `tf'
    "\bottomrule" _n
    "\end{tabular}" _n
    "\begin{tablenotes}\footnotesize" _n
    "\item \emph{Notes:} Dependent variable is $\log\mu_{it}$ from ACF "
    "translog Spec A. Panels B--C absorb firm FE so the coefficient "
    "identifies from within-firm switches into/out of procurement. Panel D "
    "interacts with single-bid share — if favoritism drives the premium, "
    "the interaction should be positive and large. Panel E tests whether "
    "the 2012 Czech single-bidding ban shrank the premium. Panel F directly "
    "controls for Rel\_Price (ratio of bid to engineer estimate); if the "
    "premium shrinks to zero, it was just reflecting tender-level overbid." _n
    "\end{tablenotes}" _n
    "\end{threeparttable}" _n
    "\end{table}" _n
;
#delimit cr
file close `tf'
dis "  Saved: favoritism_decomposition.tex"
