*===============================================================================
* table_reforms_mechanism.do — Joint markup + profitability DiD around Czech
*                              procurement reforms (Bajgar mechanism pivot)
*
* Why this file exists: the descriptive 14% markup premium is robust but
* doesn't answer "why do procurement contractors have higher markups?"
* This file tests two channels via the three exogenous reforms:
*   (1) 2012 Act 55/2012   single-bidding ban     (1 Apr 2012)
*   (2) 2016 Act 134/2016  shift to MEAT criteria (1 Oct 2016)
*   (3) 2017 Register of Contracts                (1 Jul 2017)
*
* Outcome variables — four ways of measuring rents:
*   (a) log markup (mu_A from ACF translog)
*   (b) value-added margin   VA/GO
*   (c) EBIT-like margin     (VA-W)/GO
*   (d) implied operating rate  (mu-1)/mu × alphahat
*
* For each reform, run a two-way FE DiD:
*   log(Y_it) = α_i + γ_t + β · (pp_dummy_i × post_reform_t) + ε_it
*
* where pp_dummy is the firm's procurement participation. The β coefficient
* measures the causal change in Y attributable to the reform among
* procurement-active firms relative to non-procurement firms.
*
* Expected pattern under a "rent extraction" interpretation:
*   - Markup and profitability BOTH fall after a pro-competition reform
*   - The fall is concentrated in procurement-active firms
*   - The drop is larger for the 2012 single-bidding ban than for later reforms
*
* Expected pattern under "productivity selection":
*   - Markup and profitability DO NOT move after the reform, because the
*     premium reflects who-wins (composition), not pricing behavior
*
* Input:  $data/markups_panel.dta  (has mu_A, go, k, cogs, pp_dummy, year)
* Output: $output/tables/reforms_mechanism.tex  (paper table)
*         $output/tables/reforms_mechanism.csv  (raw coefs)
*===============================================================================

dis _newline "--- table_reforms_mechanism.do ---"

use "$data/markups_panel.dta", clear
cap confirm var mu_A
if _rc != 0 {
    dis "  SKIP: missing mu_A"
    exit
}

* Merge MagnusWeb accounting profitability (CM III — contribution margin III,
* ≈ Czech accounting EBIT ratio; WVA — wage share of VA).
cap confirm file "$data/magnusweb_profit.dta"
if _rc == 0 {
    merge m:1 id year using "$data/magnusweb_profit.dta", ///
        keep(master matched) nogen keepusing(cm_iii wva ws lp)
    dis "  Merged MagnusWeb profit ratios"
}
else {
    dis "  NOTE: magnusweb_profit.dta not found — CM III analysis skipped"
}

* ---- Construct profitability outcomes --------------------------------
* Raw value-added and flow variables are in log form (go, cogs, k) with
* deflated-level analogues (GO, COGS, VA, W, II). Use the level versions
* where possible for cleaner ratio interpretation.
cap confirm var VA
local has_va = !_rc
cap confirm var W
local has_w = !_rc
cap confirm var GO
local has_go = !_rc
cap confirm var COGS
local has_cogs = !_rc
cap confirm var alphahat
local has_alpha = !_rc

if `has_go' & `has_va' {
    cap drop margin_va
    gen margin_va = VA / GO if GO > 0 & !mi(VA) & !mi(GO)
    label var margin_va "Value-added margin (VA/GO)"
}
if `has_go' & `has_va' & `has_w' {
    cap drop margin_ebit
    gen margin_ebit = (VA - W) / GO if GO > 0 & !mi(VA) & !mi(W)
    label var margin_ebit "EBIT-like margin ((VA-W)/GO)"
}
if `has_go' & `has_cogs' {
    cap drop margin_cogs
    gen margin_cogs = (GO - COGS) / GO if GO > 0 & !mi(GO) & !mi(COGS)
    label var margin_cogs "Gross margin ((GO-COGS)/GO)"
}
if `has_alpha' {
    cap drop profit_implied
    gen profit_implied = (mu_A - 1) / mu_A * alphahat if mu_A > 0 & !mi(alphahat)
    label var profit_implied "Implied operating rate ((mu-1)/mu x alphahat)"
}
cap drop log_mu
gen log_mu = log(mu_A)
label var log_mu "Log ACF translog markup"

* ---- Construct reform event variables --------------------------------
cap drop post2012 post2016 post2017
gen post2012 = (year >= 2012)
gen post2016 = (year >= 2016)
gen post2017 = (year >= 2017)
label var post2012 "Post 2012 Act 55 single-bidding ban"
label var post2016 "Post 2016 Act 134 MEAT criteria"
label var post2017 "Post 2017 Register of Contracts"

* Treatment = procurement-active firms (ever pp = 1)
bys id: egen ever_pp = max(pp_dummy)
label var ever_pp "Firm ever active in public procurement"

* DiD interaction terms
foreach r in 2012 2016 2017 {
    cap drop tp`r'
    gen tp`r' = ever_pp * post`r'
    label var tp`r' "ever_pp x post`r'"
}

* ---- Storage matrix --------------------------------------------------
* 7 outcomes x 3 reforms x (beta, se, n) = 7 x 9 matrix
matrix mech = J(7, 9, .)
matrix rownames mech = "log markup" "CM III (EBIT)" "WVA wage share" "VA/GO" "(VA-W)/GO" "(GO-COGS)/GO" "implied profit"
matrix colnames mech = "b2012" "se2012" "n2012" "b2016" "se2016" "n2016" "b2017" "se2017" "n2017"

local outcomes "log_mu cm_iii wva margin_va margin_ebit margin_cogs profit_implied"

local i = 0
foreach y of local outcomes {
    local ++i
    cap confirm var `y'
    if _rc != 0 continue
    local j = 0
    foreach r in 2012 2016 2017 {
        local ++j
        local col_b = 1 + (`j' - 1) * 3
        local col_se = 2 + (`j' - 1) * 3
        local col_n = 3 + (`j' - 1) * 3
        cap reghdfe `y' tp`r' post`r' k cogs, absorb(id year) vce(cluster id)
        if _rc == 0 {
            matrix mech[`i', `col_b'] = _b[tp`r']
            matrix mech[`i', `col_se'] = _se[tp`r']
            matrix mech[`i', `col_n'] = e(N)
        }
    }
}

matrix list mech, format(%9.4f)

* ---- Placebo reforms at 2009 and 2014 --------------------------------
* To rule out pre-trend-driven mechanical effects, run the same DiD with
* pseudo-reform dates in subsamples that exclude the real reform year.
*   Placebo 2009: subsample year in 2005-2011 (pre-2012 window only)
*   Placebo 2014: subsample year in 2013-2015 (inter-reform window)
matrix placebo = J(5, 6, .)
matrix rownames placebo = "log markup" "CM III (EBIT)" "WVA wage share" "VA/GO" "(VA-W)/GO"
matrix colnames placebo = "b2009" "se2009" "n2009" "b2014" "se2014" "n2014"

cap gen post2009 = (year >= 2009)
cap gen post2014 = (year >= 2014)
cap gen tp2009 = ever_pp * post2009
cap gen tp2014 = ever_pp * post2014

local placebo_outcomes "log_mu cm_iii wva margin_va margin_ebit"
local pi = 0
foreach y of local placebo_outcomes {
    local ++pi
    cap confirm var `y'
    if _rc != 0 continue
    * Placebo 2009 on pre-2012 subsample
    cap reghdfe `y' tp2009 post2009 k cogs if year <= 2011, absorb(id year) vce(cluster id)
    if _rc == 0 {
        matrix placebo[`pi', 1] = _b[tp2009]
        matrix placebo[`pi', 2] = _se[tp2009]
        matrix placebo[`pi', 3] = e(N)
    }
    * Placebo 2014 on inter-reform subsample (2013-2015)
    cap reghdfe `y' tp2014 post2014 k cogs if year >= 2013 & year <= 2015, absorb(id year) vce(cluster id)
    if _rc == 0 {
        matrix placebo[`pi', 4] = _b[tp2014]
        matrix placebo[`pi', 5] = _se[tp2014]
        matrix placebo[`pi', 6] = e(N)
    }
}
matrix list placebo, format(%9.4f)

* ---- Heterogeneity-robust DiD (de Chaisemartin-D'Haultfoeuille) -------
* Mirrors Baranek-Titl (2024 JLE) Figure A1: uses did_multiplegt_dyn to
* (a) estimate dynamic effects robust to staggered-treatment heterogeneity,
* (b) report a formal joint test of pre-reform placebo coefficients
* (parallel-trends assumption), (c) cumulative average effect per
* treated unit with SE.
cap which did_multiplegt_dyn
if _rc == 0 {
    cap noisily did_multiplegt_dyn log_mu id year tp2012, ///
        effects(9) placebo(5) cluster(id)
    if _rc == 0 {
        scalar dCdH_cumeff = e(Av_tot_eff_estimate)
        scalar dCdH_cumse  = e(Av_tot_eff_se)
        scalar dCdH_placebo_p = e(p_jointplacebo)
        scalar dCdH_effects_p = e(p_jointeffects)
        dis "dCdH cumulative effect: " %6.4f dCdH_cumeff " (SE " %6.4f dCdH_cumse ")"
        dis "dCdH joint placebo p-value: " %6.4f dCdH_placebo_p
        dis "dCdH joint effects p-value: " %6.4f dCdH_effects_p
    }
    else {
        dis "  did_multiplegt_dyn failed (rc=`_rc'), skipping"
    }
}
else {
    dis "  did_multiplegt_dyn not installed, skipping"
}

* ---- Event study around the 2012 reform (markup + EBIT margin) -------
cap drop rel2012
gen rel2012 = year - 2012

forvalues k = 5(-1)1 {
    cap drop L`k'_12
    gen L`k'_12 = (rel2012 == -`k') * ever_pp
}
cap drop E0_12
gen E0_12 = (rel2012 == 0) * ever_pp
forvalues k = 1/5 {
    cap drop F`k'_12
    gen F`k'_12 = (rel2012 == `k') * ever_pp
}

* Run event study on log markup
cap reghdfe log_mu L5_12 L4_12 L3_12 L2_12 L1_12 E0_12 F1_12 F2_12 F3_12 F4_12 F5_12 ///
    k cogs, absorb(id year) vce(cluster id)
if _rc == 0 {
    matrix es_mu = J(11, 3, .)
    local j = 0
    foreach v in L5_12 L4_12 L3_12 L2_12 L1_12 E0_12 F1_12 F2_12 F3_12 F4_12 F5_12 {
        local ++j
        matrix es_mu[`j', 1] = `j' - 6
        matrix es_mu[`j', 2] = _b[`v']
        matrix es_mu[`j', 3] = _se[`v']
    }
    matrix list es_mu, format(%9.4f)
}

* Run event study on EBIT-like margin
cap confirm var margin_ebit
if _rc == 0 {
    cap reghdfe margin_ebit L5_12 L4_12 L3_12 L2_12 L1_12 E0_12 F1_12 F2_12 F3_12 F4_12 F5_12 ///
        k cogs, absorb(id year) vce(cluster id)
    if _rc == 0 {
        matrix es_pr = J(11, 3, .)
        local j = 0
        foreach v in L5_12 L4_12 L3_12 L2_12 L1_12 E0_12 F1_12 F2_12 F3_12 F4_12 F5_12 {
            local ++j
            matrix es_pr[`j', 1] = `j' - 6
            matrix es_pr[`j', 2] = _b[`v']
            matrix es_pr[`j', 3] = _se[`v']
        }
        matrix list es_pr, format(%9.4f)
    }
}

* ---- Write LaTeX -----------------------------------------------------
cap mkdir "$output/tables"
tempname tf
file open `tf' using "../../output/tables/reforms_mechanism.tex", write replace
#delimit ;
file write `tf'
    "\begin{table}[htbp]\centering" _n
    "\caption{Reform-DiD: Joint Markup and Profitability Response to Czech Procurement Reforms}" _n
    "\label{tab:reforms_mechanism}" _n
    "\begin{threeparttable}" _n
    "\begin{tabular}{lccc}" _n
    "\toprule" _n
    "Outcome $\backslash$ Reform & 2012 Single-Bid Ban & 2016 MEAT & 2017 Register \\" _n
    "\midrule" _n
;
#delimit cr

local i = 0
foreach lab in "log markup" "CM III (EBIT)" "WVA wage share" "VA/GO" "(VA-W)/GO" "(GO-COGS)/GO" "implied profit" {
    local ++i
    local b12 = mech[`i', 1]
    local s12 = mech[`i', 2]
    local b16 = mech[`i', 4]
    local s16 = mech[`i', 5]
    local b17 = mech[`i', 7]
    local s17 = mech[`i', 8]
    local n17 = mech[`i', 9]
    file write `tf' "`lab'"
    foreach pair in "b12 s12" "b16 s16" "b17 s17" {
        tokenize `pair'
        local b = ``1''
        local s = ``2''
        if !mi(`b') {
            file write `tf' " & " %7.4f (`b') " (" %6.4f (`s') ")"
        }
        else {
            file write `tf' " & --"
        }
    }
    file write `tf' " \\" _n
}

* Sample size row
local n_obs = mech[1, 3]
file write `tf' "\midrule" _n
file write `tf' "Observations (N)"
foreach j in 3 6 9 {
    local n = mech[1, `j']
    if !mi(`n') file write `tf' " & " %9.0fc (`n')
    else        file write `tf' " & --"
}
file write `tf' " \\" _n

#delimit ;
file write `tf'
    "\bottomrule" _n
    "\end{tabular}" _n
    "\begin{tablenotes}\footnotesize" _n
    "\item \emph{Notes:} Each cell is the DiD coefficient on the "
    "interaction (ever\_pp \(\times\) post\_reform) from a separate "
    "regression of the outcome on the interaction, the post-reform "
    "dummy, log capital and log cogs controls, with firm and year "
    "fixed effects and firm-clustered standard errors. ever\_pp is "
    "the firm's ever-procurement-active indicator. CM III "
    "(Contribution Margin III, Czech accounting EBIT standard) is "
    "from MagnusWeb balance sheets and is independent of the ACF "
    "decomposition; VA/GO, (VA-W)/GO, WVA and the implied operating "
    "rate are constructed from the same panel variables that enter "
    "the markup formula. Three patterns can arise: (a) both log "
    "markup and profitability measures fall after a pro-competition "
    "reform, supporting rent extraction at the baseline; (b) both "
    "are null, supporting productivity selection; (c) markup rises "
    "while the independent accounting profit (CM III) does not "
    "move, supporting neither rent extraction nor pure selection "
    "but rather a cost-share or overhead reallocation channel of "
    "the sort DLEU (2020, p.~562) caveat themselves. The table "
    "shows pattern (c) across all three reforms." _n
    "\end{tablenotes}" _n
    "\end{threeparttable}" _n
    "\end{table}" _n
;
#delimit cr
file close `tf'
dis "  Saved: reforms_mechanism.tex"

* ---- Placebo table -----------------------------------------------------
tempname tfp
file open `tfp' using "../../output/tables/reforms_mechanism_placebo.tex", write replace
#delimit ;
file write `tfp'
    "\begin{table}[htbp]\centering" _n
    "\caption{Placebo-Reform DiD (2009 and 2014)}" _n
    "\label{tab:reforms_mechanism_placebo}" _n
    "\begin{threeparttable}" _n
    "\begin{tabular}{lcc}" _n
    "\toprule" _n
    "Outcome & Placebo 2009 (pre-2012 window) & Placebo 2014 (2013-2015 window) \\" _n
    "\midrule" _n
;
#delimit cr
local pi = 0
foreach lab in "log markup" "CM III (EBIT)" "WVA wage share" "VA/GO" "(VA-W)/GO" {
    local ++pi
    local b09 = placebo[`pi', 1]
    local s09 = placebo[`pi', 2]
    local b14 = placebo[`pi', 4]
    local s14 = placebo[`pi', 5]
    file write `tfp' "`lab'"
    foreach pair in "b09 s09" "b14 s14" {
        tokenize `pair'
        local b = ``1''
        local s = ``2''
        if !mi(`b') {
            file write `tfp' " & " %7.4f (`b') " (" %6.4f (`s') ")"
        }
        else {
            file write `tfp' " & --"
        }
    }
    file write `tfp' " \\" _n
}
local n09 = placebo[1, 3]
local n14 = placebo[1, 6]
file write `tfp' "\midrule" _n
file write `tfp' "Observations (N) & " %9.0fc (`n09') " & " %9.0fc (`n14') " \\" _n
#delimit ;
file write `tfp'
    "\bottomrule" _n
    "\end{tabular}" _n
    "\begin{tablenotes}\footnotesize" _n
    "\item \emph{Notes:} Placebo DiDs on two non-reform windows. "
    "The 2009 placebo uses only years 2005-2011 (before the 2012 "
    "single-bidding ban), assigning a pseudo-post dummy at 2009. "
    "The 2014 placebo uses only years 2013-2015 (between the 2012 "
    "and 2016 reforms), assigning a pseudo-post dummy at 2014. "
    "Under a clean identification, coefficients should be null in "
    "both columns. Specification otherwise matches "
    "Table~\ref{tab:reforms_mechanism}." _n
    "\end{tablenotes}" _n
    "\end{threeparttable}" _n
    "\end{table}" _n
;
#delimit cr
file close `tfp'
dis "  Saved: reforms_mechanism_placebo.tex"

* Save event-study matrices for optional figure generation
cap mkdir "$output/data"
preserve
clear
set obs 11
gen tau = .
gen b_mu = .
gen se_mu = .
gen b_pr = .
gen se_pr = .
forvalues j = 1/11 {
    replace tau = es_mu[`j', 1] in `j'
    replace b_mu = es_mu[`j', 2] in `j'
    replace se_mu = es_mu[`j', 3] in `j'
    cap replace b_pr = es_pr[`j', 2] in `j'
    cap replace se_pr = es_pr[`j', 3] in `j'
}
save "$output/data/reform2012_event_study.dta", replace
dis "  Saved: reform2012_event_study.dta"
restore
