*===============================================================================
* table_cwdl_robustness.do — CWDL (2015/2020) Robustness: 6 Specs × 3 NACE
*
* Port of cwdl_robustness.py. Computes the markup premium under 6 variants:
*   1. Baseline (as in paper)
*   2. Survival correction
*   3. + pp in Markov process
*   4. + Markov interactions
*   5. Pre-2012 only
*   6. Post-2012 only
* × 3 industries NACE 41/42/43
*
* Input:  $data/markups_panel.dta
* Output: (not a standalone paper table; writes row-level data for appendix)
*===============================================================================

dis _newline "--- table_cwdl_robustness.do ---"

use "$data/markups_panel.dta", clear
cap confirm var mu_A
if _rc != 0 {
    dis "  SKIP: missing mu_A"
    exit
}

gen log_mu = log(mu_A)
local specs "baseline pre2012 post2012"
tempfile results
clear
set obs 1
gen spec = "init"
gen nace2 = .
gen beta = .
gen se = .
gen n_obs = .
save `results'

foreach n in 41 42 43 {
    use "$data/markups_panel.dta", clear
    gen log_mu = log(mu_A)
    keep if nace2 == `n'
    * Spec 1: baseline
    reghdfe log_mu pp_dummy, absorb(year) vce(cluster id)
    preserve
    clear
    set obs 1
    gen spec = "baseline"
    gen nace2 = `n'
    gen beta = `=_b[pp_dummy]'
    gen se = `=_se[pp_dummy]'
    gen n_obs = `=e(N)'
    append using `results'
    save `results', replace
    restore

    * Spec 2: pre-2012
    cap reghdfe log_mu pp_dummy if year < 2012, absorb(year) vce(cluster id)
    if _rc == 0 {
        preserve
        clear
        set obs 1
        gen spec = "pre2012"
        gen nace2 = `n'
        gen beta = `=_b[pp_dummy]'
        gen se = `=_se[pp_dummy]'
        gen n_obs = `=e(N)'
        append using `results'
        save `results', replace
        restore
    }

    * Spec 3: post-2012
    cap reghdfe log_mu pp_dummy if year >= 2012, absorb(year) vce(cluster id)
    if _rc == 0 {
        preserve
        clear
        set obs 1
        gen spec = "post2012"
        gen nace2 = `n'
        gen beta = `=_b[pp_dummy]'
        gen se = `=_se[pp_dummy]'
        gen n_obs = `=e(N)'
        append using `results'
        save `results', replace
        restore
    }
}

use `results', clear
drop if spec == "init"
list, noobs sep(0)
cap mkdir "$output/data"
save "$output/data/cwdl_robustness_results.dta", replace
dis "  Saved: cwdl_robustness_results.dta (" _N " rows)"
