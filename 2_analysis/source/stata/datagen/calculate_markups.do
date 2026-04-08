*===============================================================================
* calculate_markups.do — Compute firm-level markups from PF coefficients
*
* Analog of DGM markup_calculation.do:
*   DGM: elasticity × (sales/materials), 4 specs × 10 treatments
*   Us:  θ/α for CD and TL, 6 specs (A-E + OLS)
*
* Markup formula (DLW 2012):
*   CD:  μ = β_cogs / α,  where α = exp(cogs)/exp(y_c)
*   TL:  μ = (β_cogs + 2·β_cogs2·cogs + β_kcogs·k) / α
*
* Input:  analysisdata/analysis_panel.dta + coefficients_byind.dta
* Output: analysisdata/markups_panel.dta
*===============================================================================

dis _newline "--- calculate_markups.do ---"

use "$data/analysis_panel.dta", clear
xtset id year, yearly

* Load coefficients
preserve
use "$data/coefficients_byind.dta", clear
tempfile coeff
save `coeff'
restore

* Merge coefficients by nace2
merge m:1 nace2 using `coeff', assert(3) nogen

*-----------------------------------------------------------------------
* First stage (needed for corrected α)
*-----------------------------------------------------------------------

* We need phi and epsilon. Re-run first stage by industry.
gen phi = .
gen epsilon = .

levelsof nace2, local(naces)
foreach n of local naces {
    xi: qui reg go c.k*#pp_dummy c.cogs*#pp_dummy i.year if nace2 == `n'
    predict phi_`n' if nace2 == `n'
    predict eps_`n' if nace2 == `n', res
    replace phi = phi_`n' if nace2 == `n'
    replace epsilon = eps_`n' if nace2 == `n'
    drop phi_`n' eps_`n' _I*
}

* Corrected expenditure share (DGM: ratio = catotal/acha4; we use α = exp(cogs)/exp(y_c))
gen y_c = go - epsilon
gen Y_c = exp(y_c)
gen alphahat = exp(cogs) / Y_c
label var alphahat "Corrected expenditure share (cogs/output)"

*-----------------------------------------------------------------------
* Markup computation (analog of DGM's 4-spec × 10-treatment loop)
*-----------------------------------------------------------------------

* Spec A: Base CD
gen mu_A = b_cogs_A / alphahat
label var mu_A "Markup: Base (surv + pp Markov)"

* Spec B: No survival CD
gen mu_B = b_cogs_B / alphahat
label var mu_B "Markup: No survival"

* Spec C: No pp Markov CD
gen mu_C = b_cogs_C / alphahat
label var mu_C "Markup: No pp in Markov"

* Spec D: Plain CD
gen mu_D = b_cogs_D / alphahat
label var mu_D "Markup: Plain ACF"

* Spec E: Translog
gen theta_E = b_cogs_E + 2*b_cogs2_E*cogs + b_kcogs_E*k
gen mu_E = theta_E / alphahat
label var mu_E "Markup: Translog (surv + pp)"
label var theta_E "Output elasticity (TL)"

* OLS
gen mu_OLS = b_cogs_OLS / alphahat
label var mu_OLS "Markup: OLS"

*-----------------------------------------------------------------------
* Log markups (DGM: l_mu_*)
*-----------------------------------------------------------------------

foreach s in A B C D E OLS {
    gen l_mu_`s' = ln(mu_`s') if mu_`s' > 0
    label var l_mu_`s' "Log markup: spec `s'"
}

*-----------------------------------------------------------------------
* Productivity (omega, from base spec A)
*-----------------------------------------------------------------------

gen omega_A = phi - b_k_A*k - b_cogs_A*cogs
* No constant subtraction needed — it's absorbed in phi

label var omega_A "Productivity (omega, base spec)"

*-----------------------------------------------------------------------
* First differences (DGM Figure 4)
*-----------------------------------------------------------------------

sort id year
foreach s in A B C D E OLS {
    by id: gen FD_l_mu_`s' = l_mu_`s' - l_mu_`s'[_n-1]
    label var FD_l_mu_`s' "FD log markup: spec `s'"
}

*-----------------------------------------------------------------------
* Summary statistics
*-----------------------------------------------------------------------

dis _newline "=== Markup Distributions ==="
foreach s in A B C D E OLS {
    dis _newline "  Spec `s':"
    tabstat mu_`s', by(nace2) stat(mean sd p10 p50 p90 N) format(%9.3f)
}

dis _newline "=== By Procurement Status (Base Spec A) ==="
tabstat mu_A l_mu_A omega_A, by(pp_dummy) stat(mean sd p50 N) format(%9.4f)

*-----------------------------------------------------------------------
* Save
*-----------------------------------------------------------------------

drop b_k_* b_cogs_* b_k2_* b_cogs2_* b_kcogs_* se_* crit_* N_obs r2_first
drop phi epsilon y_c Y_c theta_E

order id year nace2 pp_dummy pp_lag pp_ever_3y ///
    mu_A mu_B mu_C mu_D mu_E mu_OLS ///
    l_mu_A l_mu_B l_mu_C l_mu_D l_mu_E l_mu_OLS ///
    omega_A alphahat ratio

compress
save "$data/markups_panel.dta", replace
dis _newline "  Saved: markups_panel.dta (" _N " obs)"
