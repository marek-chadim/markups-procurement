*===============================================================================
* table_dlw_eval.do — DLW (2012) Treatment Evaluation Framework
*
* Replicates De Loecker & Warzynski (2012, AER) treatment evaluation
* with procurement status replacing export status:
*   1. Cross-sectional premium (DLW Table 3, eq. 20)
*   2. Entry/exit dynamics (DLW Table 4, eq. 22)
*   3. Controlling for productivity (DLW eq. 21)
*   4. Intensive margin (procurement share)
*   5. Productivity-markup correlation (DLW Section V.C)
*   6. Sub-industry heterogeneity (DLW destination analog)
*   7. Event study (dynamic entry effects)
*
* Parallels Python dlw_treatment_eval.py
*
* Input:  $data/markups_panel.dta (from calculate_markups.do)
* Output: $output/dlw_treatment_eval_stata.tex
*         $output/dlw_event_study_stata.csv
*===============================================================================

dis _newline "--- table_dlw_eval.do ---"

*-----------------------------------------------------------------------
* Setup: check for reghdfe/estout, fall back gracefully
*-----------------------------------------------------------------------

cap which reghdfe
local has_reghdfe = (_rc == 0)
if !`has_reghdfe' {
    dis as text "  Note: reghdfe not installed; using areg with group FE"
}

cap which esttab
local has_estout = (_rc == 0)
if !`has_estout' {
    dis as text "  Note: estout not installed; using manual LaTeX output"
}

*-----------------------------------------------------------------------
* Load data
*-----------------------------------------------------------------------

use "$data/markups_panel.dta", clear

* Verify key variables exist
foreach v in l_mu_A omega_A pp_dummy k cogs nace2 year id {
    cap confirm var `v'
    if _rc {
        dis as error "ERROR: variable `v' not found in markups_panel.dta"
        error 198
    }
}

* Generate year × nace2 FE group
cap drop yr_nace
egen yr_nace = group(year nace2)

* Log markup shorthand
local lmu "l_mu_A"

dis "  Loaded: " _N " obs, " r(N) " firms"
distinct id
local N_firms = r(ndistinct)
dis "  Firms: `N_firms'"

*-----------------------------------------------------------------------
* Construct entry/exit/always variables
*-----------------------------------------------------------------------

sort id year

* First procurement year per firm
bys id (year): gen _first_pp = year if pp_dummy == 1 & ///
    (pp_dummy[_n-1] == 0 | _n == 1)
bys id: egen first_pp_year = min(_first_pp)
drop _first_pp

* Last procurement year per firm
bys id: egen last_pp_year = max(cond(pp_dummy == 1, year, .))

* Last observed year per firm
bys id: egen last_year = max(year)

* Always-procurer: pp_dummy = 1 in every year
bys id: egen _min_pp = min(pp_dummy)
gen always_pp = (_min_pp == 1)
drop _min_pp

* Entry: non-always firm, post first procurement year
gen entry_pp = 0
replace entry_pp = 1 if always_pp == 0 & first_pp_year != . & ///
    year >= first_pp_year

* Exit: non-always firm, post last procurement year (firm stopped procuring)
gen exit_pp = 0
replace exit_pp = 1 if always_pp == 0 & last_pp_year != . & ///
    last_pp_year < last_year & year > last_pp_year

* Count switches for multi-switcher filter (DLW fn 52)
bys id (year): gen _switch = abs(pp_dummy - pp_dummy[_n-1]) if _n > 1
bys id: egen n_switches = total(_switch)
drop _switch

* Procurement share (if available; fill 0 otherwise)
cap confirm var pp_share
if _rc {
    gen pp_share = pp_dummy
    label var pp_share "Procurement share (proxy = pp_dummy)"
}
replace pp_share = 0 if mi(pp_share)

* Interactions for analysis 4
gen pp_x_share = pp_dummy * pp_share
label var pp_x_share "pp_dummy × pp_share"

dis _newline "  Entry/exit construction:"
dis "    Always-procurers: " %6.0f _N " obs where always_pp==1"
count if always_pp == 1
dis "    Entry obs:        " r(N)
count if entry_pp == 1
dis "    Exit obs:         " r(N)

* Preserve full sample before dropping multi-switchers
tempfile fulldata
save `fulldata'


*=======================================================================
*  ANALYSIS 1: Cross-sectional premium (DLW Table 3, eq. 20)
*=======================================================================

dis _newline(2) "========================================"
dis "  1. CROSS-SECTIONAL PREMIUM (DLW Table 3)"
dis "========================================"

if `has_reghdfe' {
    reghdfe `lmu' pp_dummy k cogs, absorb(yr_nace) cluster(id)
}
else {
    areg `lmu' pp_dummy k cogs, absorb(yr_nace) cluster(id)
}
est store dlw1_xsec

local b1 = _b[pp_dummy]
local se1 = _se[pp_dummy]
local N1 = e(N)
dis _newline "  Premium (pp_dummy): " %7.4f `b1' " (SE " %6.4f `se1' ")"
dis "  N = " %6.0fc `N1'


*=======================================================================
*  ANALYSIS 2: Entry/exit dynamics (DLW Table 4, eq. 22)
*=======================================================================

dis _newline(2) "========================================"
dis "  2. ENTRY/EXIT DYNAMICS (DLW Table 4)"
dis "========================================"

* Drop multi-switchers (DLW fn 52)
preserve
drop if n_switches > 2
local N_dropped = _N
use `fulldata', clear
local N_full = _N
restore

preserve
drop if n_switches > 2
dis "  Dropped multi-switchers (>2 switches), " _N " obs remain"

if `has_reghdfe' {
    reghdfe `lmu' entry_pp exit_pp always_pp k cogs, absorb(yr_nace) cluster(id)
}
else {
    areg `lmu' entry_pp exit_pp always_pp k cogs, absorb(yr_nace) cluster(id)
}
est store dlw2_entry

local b2_entry  = _b[entry_pp]
local se2_entry = _se[entry_pp]
local b2_exit   = _b[exit_pp]
local se2_exit  = _se[exit_pp]
local b2_always = _b[always_pp]
local se2_always = _se[always_pp]
local N2 = e(N)

dis _newline "  Entry (pp):  " %7.4f `b2_entry'  " (SE " %6.4f `se2_entry'  ")"
dis "  Exit (pp):   " %7.4f `b2_exit'   " (SE " %6.4f `se2_exit'   ")"
dis "  Always (pp): " %7.4f `b2_always' " (SE " %6.4f `se2_always' ")"
dis "  N = " %6.0fc `N2'

restore


*=======================================================================
*  ANALYSIS 3: Controlling for productivity (DLW eq. 21)
*=======================================================================

dis _newline(2) "========================================"
dis "  3. CONTROLLING FOR PRODUCTIVITY (DLW eq. 21)"
dis "========================================"

* Without omega
if `has_reghdfe' {
    reghdfe `lmu' pp_dummy k cogs, absorb(yr_nace) cluster(id)
}
else {
    areg `lmu' pp_dummy k cogs, absorb(yr_nace) cluster(id)
}
est store dlw3_noomega

local b3a = _b[pp_dummy]
local se3a = _se[pp_dummy]

* With omega
if `has_reghdfe' {
    reghdfe `lmu' pp_dummy omega_A k cogs, absorb(yr_nace) cluster(id)
}
else {
    areg `lmu' pp_dummy omega_A k cogs, absorb(yr_nace) cluster(id)
}
est store dlw3_omega

local b3b = _b[pp_dummy]
local se3b = _se[pp_dummy]
local b3_omega = _b[omega_A]
local se3_omega = _se[omega_A]
local N3 = e(N)

local pct_explained = (1 - `b3b'/`b3a') * 100

dis _newline "  Without omega: pp = " %7.4f `b3a' " (SE " %6.4f `se3a' ")"
dis "  With omega:    pp = " %7.4f `b3b' " (SE " %6.4f `se3b' ")"
dis "  omega_A coef:       " %7.4f `b3_omega' " (SE " %6.4f `se3_omega' ")"
dis "  Productivity explains " %5.1f `pct_explained' "% of premium"
dis "  (DLW: ~70% for Slovenian exporters)"
dis "  N = " %6.0fc `N3'


*=======================================================================
*  ANALYSIS 4: Intensive margin (procurement share)
*=======================================================================

dis _newline(2) "========================================"
dis "  4. INTENSIVE MARGIN (procurement share)"
dis "========================================"

if `has_reghdfe' {
    reghdfe `lmu' pp_dummy pp_x_share k cogs, absorb(yr_nace) cluster(id)
}
else {
    areg `lmu' pp_dummy pp_x_share k cogs, absorb(yr_nace) cluster(id)
}
est store dlw4_intensive

local b4_ext  = _b[pp_dummy]
local se4_ext = _se[pp_dummy]
local b4_int  = _b[pp_x_share]
local se4_int = _se[pp_x_share]
local N4 = e(N)

dis _newline "  Extensive (pp_dummy):   " %7.4f `b4_ext'  " (SE " %6.4f `se4_ext'  ")"
dis "  Intensive (pp×share):   " %7.4f `b4_int'  " (SE " %6.4f `se4_int'  ")"
dis "  N = " %6.0fc `N4'


*=======================================================================
*  ANALYSIS 5: Productivity-markup correlation (DLW Section V.C)
*=======================================================================

dis _newline(2) "========================================"
dis "  5. PRODUCTIVITY-MARKUP CORRELATION (DLW V.C)"
dis "========================================"

if `has_reghdfe' {
    reghdfe `lmu' omega_A k cogs, absorb(yr_nace) cluster(id)
}
else {
    areg `lmu' omega_A k cogs, absorb(yr_nace) cluster(id)
}
est store dlw5_omega

local b5 = _b[omega_A]
local se5 = _se[omega_A]
local N5 = e(N)

dis _newline "  omega_A coef: " %7.4f `b5' " (SE " %6.4f `se5' ")"
dis "  (DLW find beta = 0.3 for Slovenian manufacturing)"
dis "  N = " %6.0fc `N5'


*=======================================================================
*  ANALYSIS 6: Sub-industry heterogeneity (DLW destination analog)
*=======================================================================

dis _newline(2) "========================================"
dis "  6. SUB-INDUSTRY HETEROGENEITY"
dis "========================================"

gen pp_x_n41 = pp_dummy * (nace2 == 41)
gen pp_x_n42 = pp_dummy * (nace2 == 42)
gen pp_x_n43 = pp_dummy * (nace2 == 43)

label var pp_x_n41 "pp × NACE 41 (Buildings)"
label var pp_x_n42 "pp × NACE 42 (Civil eng.)"
label var pp_x_n43 "pp × NACE 43 (Specialized)"

if `has_reghdfe' {
    reghdfe `lmu' pp_x_n41 pp_x_n42 pp_x_n43 k cogs, ///
        absorb(yr_nace) cluster(id)
}
else {
    areg `lmu' pp_x_n41 pp_x_n42 pp_x_n43 k cogs, ///
        absorb(yr_nace) cluster(id)
}
est store dlw6_nace

local b6_41  = _b[pp_x_n41]
local se6_41 = _se[pp_x_n41]
local b6_42  = _b[pp_x_n42]
local se6_42 = _se[pp_x_n42]
local b6_43  = _b[pp_x_n43]
local se6_43 = _se[pp_x_n43]
local N6 = e(N)

dis _newline "  NACE 41 (Buildings):     " %7.4f `b6_41' " (SE " %6.4f `se6_41' ")"
dis "  NACE 42 (Civil eng.):    " %7.4f `b6_42' " (SE " %6.4f `se6_42' ")"
dis "  NACE 43 (Specialized):   " %7.4f `b6_43' " (SE " %6.4f `se6_43' ")"
dis "  N = " %6.0fc `N6'


*=======================================================================
*  ANALYSIS 7: Event study (dynamic entry effects)
*=======================================================================

dis _newline(2) "========================================"
dis "  7. EVENT STUDY (dynamic entry effects)"
dis "========================================"

* Restrict to firms with observed entry
preserve
keep if first_pp_year != .
gen rel_year = year - first_pp_year

* Event window: [-3, +5], omit τ = -1
keep if inrange(rel_year, -3, 5)

dis "  Event study sample: " _N " obs"
distinct id
dis "  Event study firms:  " r(ndistinct)

* Create relative-year dummies (Stata-safe names: m = minus, p = plus)
gen tau_m3 = (rel_year == -3)
gen tau_m2 = (rel_year == -2)
* tau = -1 omitted (reference period)
gen tau_p0 = (rel_year == 0)
gen tau_p1 = (rel_year == 1)
gen tau_p2 = (rel_year == 2)
gen tau_p3 = (rel_year == 3)
gen tau_p4 = (rel_year == 4)
gen tau_p5 = (rel_year == 5)

if `has_reghdfe' {
    reghdfe `lmu' tau_m3 tau_m2 tau_p0 tau_p1 tau_p2 tau_p3 tau_p4 tau_p5 ///
        k cogs, absorb(yr_nace) cluster(id)
}
else {
    areg `lmu' tau_m3 tau_m2 tau_p0 tau_p1 tau_p2 tau_p3 tau_p4 tau_p5 ///
        k cogs, absorb(yr_nace) cluster(id)
}
est store dlw7_event
local N7 = e(N)

* Collect event study coefficients
dis _newline "  tau   Coef       SE         95% CI"
dis "  " _dup(55) "-"

* Reference period
dis "   -1  (reference)"

matrix C = e(b)
matrix V = e(V)

* Save event study coefficients to CSV
cap file close evf
file open evf using "$output/dlw_event_study_stata.csv", write replace
file write evf "tau,coef,se,ci_lo,ci_hi" _n

* Pre-treatment
foreach tau_name in tau_m3 tau_m2 {
    local tau_num = cond("`tau_name'" == "tau_m3", -3, -2)
    local b = _b[`tau_name']
    local se = _se[`tau_name']
    local ci_lo = `b' - 1.96 * `se'
    local ci_hi = `b' + 1.96 * `se'
    dis "  " %3.0f `tau_num' "   " %7.4f `b' "   " %7.4f `se' ///
        "   [" %7.4f `ci_lo' ", " %7.4f `ci_hi' "]"
    file write evf "`tau_num'," %12.8f (`b') "," %12.8f (`se') ///
        "," %12.8f (`ci_lo') "," %12.8f (`ci_hi') _n
}

* Reference period (tau = -1)
file write evf "-1,0,0,0,0" _n

* Post-treatment
foreach tau_name in tau_p0 tau_p1 tau_p2 tau_p3 tau_p4 tau_p5 {
    local tau_num = cond("`tau_name'" == "tau_p0", 0, ///
                   cond("`tau_name'" == "tau_p1", 1, ///
                   cond("`tau_name'" == "tau_p2", 2, ///
                   cond("`tau_name'" == "tau_p3", 3, ///
                   cond("`tau_name'" == "tau_p4", 4, 5)))))
    local b = _b[`tau_name']
    local se = _se[`tau_name']
    local ci_lo = `b' - 1.96 * `se'
    local ci_hi = `b' + 1.96 * `se'
    dis "  " %3.0f `tau_num' "   " %7.4f `b' "   " %7.4f `se' ///
        "   [" %7.4f `ci_lo' ", " %7.4f `ci_hi' "]"
    file write evf "`tau_num'," %12.8f (`b') "," %12.8f (`se') ///
        "," %12.8f (`ci_lo') "," %12.8f (`ci_hi') _n
}

file close evf
dis _newline "  Saved: $output/dlw_event_study_stata.csv"

restore


*=======================================================================
*  OUTPUT: LaTeX table
*=======================================================================

dis _newline(2) "========================================"
dis "  Writing LaTeX table"
dis "========================================"

cap file close tf
file open tf using "$output/dlw_treatment_eval_stata.tex", write replace

#delimit ;
file write tf
    "\begin{table}[htbp]\centering" _n
    "\caption{DLW (2012) Treatment Evaluation: Procurement Premium}"
    "\label{tab:dlw_eval}" _n
    "\begin{tabular}{lcccc}" _n
    "\hline\hline" _n
    " & (1) & (2) & (3) & (4) \\" _n
    " & Cross-section & Entry/Exit & With "
    _char(36) "\omega" _char(36) " & Intensive \\" _n
    "\hline" _n
;
#delimit cr

* Row: pp_dummy
file write tf _char(36) "pp_t" _char(36)
file write tf " & " %6.4f (`b1')
file write tf " & "
file write tf " & " %6.4f (`b3b')
file write tf " & " %6.4f (`b4_ext') " \\" _n
file write tf " & (" %5.4f (`se1') ")"
file write tf " & "
file write tf " & (" %5.4f (`se3b') ")"
file write tf " & (" %5.4f (`se4_ext') ") \\" _n

* Row: entry_pp
file write tf "Entry"
file write tf " & "
file write tf " & " %6.4f (`b2_entry')
file write tf " & "
file write tf " & \\" _n
file write tf " & "
file write tf " & (" %5.4f (`se2_entry') ")"
file write tf " & "
file write tf " & \\" _n

* Row: exit_pp
file write tf "Exit"
file write tf " & "
file write tf " & " %6.4f (`b2_exit')
file write tf " & "
file write tf " & \\" _n
file write tf " & "
file write tf " & (" %5.4f (`se2_exit') ")"
file write tf " & "
file write tf " & \\" _n

* Row: always_pp
file write tf "Always"
file write tf " & "
file write tf " & " %6.4f (`b2_always')
file write tf " & "
file write tf " & \\" _n
file write tf " & "
file write tf " & (" %5.4f (`se2_always') ")"
file write tf " & "
file write tf " & \\" _n

* Row: omega_A
file write tf _char(36) "\omega_A" _char(36)
file write tf " & "
file write tf " & "
file write tf " & " %6.4f (`b3_omega')
file write tf " & \\" _n
file write tf " & "
file write tf " & "
file write tf " & (" %5.4f (`se3_omega') ")"
file write tf " & \\" _n

* Row: pp_x_share
file write tf _char(36) "pp_t \times \text{share}" _char(36)
file write tf " & "
file write tf " & "
file write tf " & "
file write tf " & " %6.4f (`b4_int') " \\" _n
file write tf " & "
file write tf " & "
file write tf " & "
file write tf " & (" %5.4f (`se4_int') ") \\" _n

* Separator and sub-industry panel
#delimit ;
file write tf
    "\hline" _n
    "\multicolumn{5}{l}{\textit{Sub-industry heterogeneity (Analysis 6)}} \\" _n
;
#delimit cr

file write tf "NACE 41 (Buildings)"
file write tf " & \multicolumn{4}{c}{" %6.4f (`b6_41') " (" %5.4f (`se6_41') ")} \\" _n
file write tf "NACE 42 (Civil eng.)"
file write tf " & \multicolumn{4}{c}{" %6.4f (`b6_42') " (" %5.4f (`se6_42') ")} \\" _n
file write tf "NACE 43 (Specialized)"
file write tf " & \multicolumn{4}{c}{" %6.4f (`b6_43') " (" %5.4f (`se6_43') ")} \\" _n

* Footer
#delimit ;
file write tf
    "\hline" _n
    "Year " _char(36) "\times" _char(36) " NACE FE & Yes & Yes & Yes & Yes \\" _n
    "Controls (" _char(36) "k, \text{cogs}" _char(36) ") & Yes & Yes & Yes & Yes \\" _n
    _char(36) "N" _char(36) " & " %6.0fc (`N1')
        " & " %6.0fc (`N2')
        " & " %6.0fc (`N3')
        " & " %6.0fc (`N4') " \\" _n
    "\hline\hline" _n
    "\multicolumn{5}{p{12cm}}{\footnotesize" _n
    "\emph{Notes:} All regressions include year " _char(36) "\times" _char(36)
    " NACE 2-digit FE and controls for " _char(36) "k_{it}" _char(36)
    " and " _char(36) "\text{cogs}_{it}" _char(36) ". "
    "Firm-clustered standard errors in parentheses. "
    "Following De Loecker and Warzynski (2012), "
    "Entry" _char(36) "_{it}" _char(36) " = 1 post-first-procurement, "
    "Exit" _char(36) "_{it}" _char(36) " = 1 post-last-procurement, "
    "Always" _char(36) "_i" _char(36) " = 1 for permanent procurement firms. "
    "Multi-switchers ($>$2 status changes) excluded from col.\ (2). "
    "Analysis 5 (productivity-markup correlation): "
    _char(36) "\hat{\beta}_\omega" _char(36)
    " = " %6.4f (`b5') " (SE " %5.4f (`se5') "). "
    "Event study coefficients in separate CSV.} \\" _n
    "\end{tabular}" _n
    "\end{table}" _n
;
#delimit cr

file close tf
dis "  Saved: $output/dlw_treatment_eval_stata.tex"


*=======================================================================
*  Summary
*=======================================================================

dis _newline(2) "========================================"
dis "  DLW TREATMENT EVALUATION SUMMARY"
dis "========================================"
dis "  1. Cross-section premium:   " %7.4f `b1'  " (SE " %6.4f `se1'  ")"
dis "  2. Entry effect:            " %7.4f `b2_entry'  " (SE " %6.4f `se2_entry'  ")"
dis "     Exit effect:             " %7.4f `b2_exit'   " (SE " %6.4f `se2_exit'   ")"
dis "     Always-procurer:         " %7.4f `b2_always' " (SE " %6.4f `se2_always' ")"
dis "  3. Premium with omega:      " %7.4f `b3b' " (omega explains " %5.1f `pct_explained' "%)"
dis "  4. Intensive margin:        " %7.4f `b4_int'  " (SE " %6.4f `se4_int'  ")"
dis "  5. Omega-markup corr:       " %7.4f `b5'  " (SE " %6.4f `se5'  ")"
dis "  6. NACE 41/42/43:           " %5.3f `b6_41' " / " %5.3f `b6_42' " / " %5.3f `b6_43'
dis "  7. Event study:             see dlw_event_study_stata.csv"
dis "========================================"

* Clean up generated variables
cap drop yr_nace first_pp_year last_pp_year last_year
cap drop always_pp entry_pp exit_pp n_switches
cap drop pp_x_share pp_x_n41 pp_x_n42 pp_x_n43
