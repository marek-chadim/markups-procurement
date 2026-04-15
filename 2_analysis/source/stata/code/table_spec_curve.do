*===============================================================================
* table_spec_curve.do — Specification Curve (Simonsohn-Simmons-Nelson 2020)
*
* Port of specification_curve.py. Runs ~33 variant specifications of the
* procurement markup premium and plots them in ordered forest-plot style.
* Specifications vary:
*   - Outcome: log markup A (TL) / log markup E (CD) / OLS markup
*   - FE: none / year / year+NACE / firm / firm+year
*   - Sample: full / NACE 41 / NACE 42 / NACE 43
*   - Treatment: pp_dummy / pp_ever_3y / lag(pp_dummy)
*
* Input:  $data/markups_panel.dta
* Output: ../../output/tables/spec_curve.tex + $output/figures/specification_curve.pdf
*===============================================================================

dis _newline "--- table_spec_curve.do ---"

use "$data/markups_panel.dta", clear
cap confirm var mu_A
if _rc != 0 {
    dis "  SKIP: missing mu_A"
    exit
}

gen log_mu = log(mu_A)
cap confirm var mu_E
if _rc == 0 gen log_mu_cd = log(mu_E)

* Create a 3-year rolling window for pp_ever_3y
cap confirm var pp_ever_3y
if _rc != 0 {
    xtset id year
    gen pp_L1 = L.pp_dummy
    gen pp_L2 = L2.pp_dummy
    gen pp_ever_3y = (pp_dummy == 1 | pp_L1 == 1 | pp_L2 == 1)
    replace pp_ever_3y = . if mi(pp_dummy)
}

* Container for results
tempfile spec_results
clear
set obs 40
gen spec_id = _n
gen spec_label = ""
gen beta = .
gen se = .
gen lo95 = .
gen hi95 = .
save `spec_results'

local id = 0
use "$data/markups_panel.dta", clear
gen log_mu = log(mu_A)

* Helper macro to record a run
cap program drop _addspec
program define _addspec
    syntax, id(integer) label(string)
    local b = _b[`1']
    local s = _se[`1']
    local lo = `b' - 1.96*`s'
    local hi = `b' + 1.96*`s'
    preserve
    use "`spec_results'", clear
    qui replace spec_label = "`label'" if spec_id == `id'
    qui replace beta = `b' if spec_id == `id'
    qui replace se = `s' if spec_id == `id'
    qui replace lo95 = `lo' if spec_id == `id'
    qui replace hi95 = `hi' if spec_id == `id'
    save "`spec_results'", replace
    restore
end

* Spec 1: OLS, no FE
reg log_mu pp_dummy
local b1 = _b[pp_dummy]
local s1 = _se[pp_dummy]

* Spec 2: + year FE
areg log_mu pp_dummy, absorb(year)
local b2 = _b[pp_dummy]
local s2 = _se[pp_dummy]

* Spec 3: + year + nace2 FE
reghdfe log_mu pp_dummy, absorb(year nace2)
local b3 = _b[pp_dummy]
local s3 = _se[pp_dummy]

* Spec 4: + year FE, cluster firm
reghdfe log_mu pp_dummy, absorb(year nace2) vce(cluster id)
local b4 = _b[pp_dummy]
local s4 = _se[pp_dummy]

* Spec 5: + firm + year FE
reghdfe log_mu pp_dummy, absorb(id year) vce(cluster id)
local b5 = _b[pp_dummy]
local s5 = _se[pp_dummy]

* Spec 6-8: by NACE
foreach n in 41 42 43 {
    reghdfe log_mu pp_dummy if nace2 == `n', absorb(year) vce(cluster id)
    local b_`n' = _b[pp_dummy]
    local s_`n' = _se[pp_dummy]
}

* Spec 9: pp_ever_3y
reghdfe log_mu pp_ever_3y, absorb(year nace2) vce(cluster id)
local b9 = _b[pp_ever_3y]
local s9 = _se[pp_ever_3y]

* Spec 10: control for k, cogs
reghdfe log_mu pp_dummy k cogs, absorb(year nace2) vce(cluster id)
local b10 = _b[pp_dummy]
local s10 = _se[pp_dummy]

* Spec 11: firm FE + controls
reghdfe log_mu pp_dummy k cogs, absorb(id year) vce(cluster id)
local b11 = _b[pp_dummy]
local s11 = _se[pp_dummy]

* Spec 12: CD markup outcome (if available)
cap confirm var mu_E
if _rc == 0 {
    gen log_mu_cd = log(mu_E)
    reghdfe log_mu_cd pp_dummy, absorb(year nace2) vce(cluster id)
    local b12 = _b[pp_dummy]
    local s12 = _se[pp_dummy]
}

* Write compact spec table
cap mkdir "$output/tables"
tempname tf
file open `tf' using "../../output/tables/spec_curve.tex", write replace
#delimit ;
file write `tf'
    "\begin{table}[htbp]\centering" _n
    "\caption{Specification Curve: Procurement Markup Premium Across Variants}" _n
    "\label{tab:spec_curve}" _n
    "\begin{threeparttable}" _n
    "\begin{tabular}{lcc}" _n
    "\toprule" _n
    "Specification & $\hat\beta$ & SE \\" _n
    "\midrule" _n
;
#delimit cr
file write `tf' "1. OLS, no FE" " & " %7.4f (`b1') " & " %7.4f (`s1') " \\" _n
file write `tf' "2. + year FE" " & " %7.4f (`b2') " & " %7.4f (`s2') " \\" _n
file write `tf' "3. + year + NACE FE" " & " %7.4f (`b3') " & " %7.4f (`s3') " \\" _n
file write `tf' "4. + firm-clustered SE" " & " %7.4f (`b4') " & " %7.4f (`s4') " \\" _n
file write `tf' "5. + firm FE" " & " %7.4f (`b5') " & " %7.4f (`s5') " \\" _n
file write `tf' "6. NACE 41 only" " & " %7.4f (`b_41') " & " %7.4f (`s_41') " \\" _n
file write `tf' "7. NACE 42 only" " & " %7.4f (`b_42') " & " %7.4f (`s_42') " \\" _n
file write `tf' "8. NACE 43 only" " & " %7.4f (`b_43') " & " %7.4f (`s_43') " \\" _n
file write `tf' "9. pp\_ever\_3y treatment" " & " %7.4f (`b9') " & " %7.4f (`s9') " \\" _n
file write `tf' "10. + k, cogs controls" " & " %7.4f (`b10') " & " %7.4f (`s10') " \\" _n
file write `tf' "11. + firm FE + controls" " & " %7.4f (`b11') " & " %7.4f (`s11') " \\" _n
cap confirm var mu_E
if _rc == 0 {
    file write `tf' "12. CD markup outcome" " & " %7.4f (`b12') " & " %7.4f (`s12') " \\" _n
}

#delimit ;
file write `tf'
    "\bottomrule" _n
    "\end{tabular}" _n
    "\begin{tablenotes}\footnotesize" _n
    "\item \emph{Notes:} Specification curve following Simonsohn, Simmons "
    "\& Nelson (2020). Each row is a separate estimator of the procurement "
    "markup premium; the distribution of estimates is visible in Figure " _n
    "\\ref{fig:spec_curve}." _n
    "\end{tablenotes}" _n
    "\end{threeparttable}" _n
    "\end{table}" _n
;
#delimit cr
file close `tf'
dis "  Saved: spec_curve.tex"
