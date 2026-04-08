*===============================================================================
* table_bh_ri.do — Borusyak-Hull (2021) Randomization Inference
*
* Uses ri_pvalue.ado and ri_ci.ado (copied from BH replication package)
* to test whether the procurement markup premium survives randomization
* inference under counterfactual treatment assignments.
*===============================================================================

dis _newline "--- table_bh_ri.do ---"

* Add ado path
adopath + "$root/ado"

use "$data/analysis_panel.dta", clear
xtset id year, yearly

* Ensure log markup exists
cap confirm variable l_mu_A
if _rc {
    cap confirm variable markup_A
    if !_rc gen l_mu_A = ln(markup_A) if markup_A > 0
}

drop if missing(l_mu_A, pp_dummy, k, cogs)

* ── Generate counterfactual treatment assignments ──
* Strategy: permute pp_dummy within (year, nace2) strata
* K=500 permutations (Stata is slower than Python's 1000)

set seed 42
local K = 500

* Create strata
egen strata = group(year nace2)

* Sort-based permutation: sorting by random number within strata
* effectively shuffles pp_dummy within each stratum
forvalues s = 1/`K' {
    tempvar rand
    gen `rand' = runiform()
    bys strata (`rand'): gen pp_sim_`s' = pp_dummy[_n]
    drop `rand'
    if mod(`s', 100) == 0 dis "  Permutation `s'/`K' done"
}

dis "  Generated `K' permuted treatment assignments"

* ── Run BH randomization inference ──

* RI p-value (null: premium = 0)
ri_pvalue l_mu_A pp_dummy, zsim(pp_sim_*) beta(0) controls(k cogs) absorb(yr_nace)
local ri_p = r(pvalue)
local ri_stat = r(stat)
local ri_nsims = r(nsims)
dis _newline "  RI p-value: " `ri_p' " (stat=" `ri_stat' ", nsims=" `ri_nsims' ")"

* RI confidence interval (search range: -0.1 to 0.4)
ri_ci l_mu_A pp_dummy, zsim(pp_sim_*) range(-0.1 0.4) controls(k cogs) absorb(yr_nace)
local ci_lo = r(ci_left)
local ci_hi = r(ci_right)
dis "  RI 95% CI: [" `ci_lo' ", " `ci_hi' "]"

* Standard OLS for comparison
reghdfe l_mu_A pp_dummy k cogs, absorb(yr_nace) cluster(id)
local ols_b = _b[pp_dummy]
local ols_se = _se[pp_dummy]
local ols_lo = `ols_b' - 1.96 * `ols_se'
local ols_hi = `ols_b' + 1.96 * `ols_se'

dis _newline "  === Comparison ==="
dis "  OLS:     " %6.4f `ols_b' " [" %6.4f `ols_lo' ", " %6.4f `ols_hi' "]"
dis "  BH RI:   " %6.4f `ols_b' " [" %6.4f `ci_lo' ", " %6.4f `ci_hi' "] p=" %6.4f `ri_p'

* Save results to CSV
tempname fh
file open `fh' using "$output/bh_ri_stata.csv", write replace
file write `fh' "method,estimate,ci_lo,ci_hi,pvalue" _n
file write `fh' "OLS clustered," %6.4f (`ols_b') "," %6.4f (`ols_lo') "," %6.4f (`ols_hi') ",." _n
file write `fh' "BH RI (within yr x nace)," %6.4f (`ols_b') "," %6.4f (`ci_lo') "," %6.4f (`ci_hi') "," %6.4f (`ri_p') _n
file close `fh'

dis _newline "--- table_bh_ri.do DONE ---"
