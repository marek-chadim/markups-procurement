*===============================================================================
* table_cwdl_acf.do — CWDL (2020) ACF cross-validation
*
* Independent Stata ACF estimation using Collard-Wexler & De Loecker (2020)
* Mata GMM procedure. Cross-validates Python acf_estimator.py results.
*===============================================================================

dis _newline "--- table_cwdl_acf.do ---"

use "$data/analysis_panel.dta", clear
xtset id year, yearly

*-----------------------------------------------------------------------
* Regenerate ACF first-stage variables locally.
*
* estimate_pf.do creates phi/Lphi/Lk/Lcogs inside a `preserve ... restore`
* block, so they do not persist to analysis_panel.dta. Rather than modify
* the upstream pipeline (which would ripple through multiple scripts), we
* reconstruct the variables here using the same first-stage regression
* specification as estimate_pf.do (xi: reg go c.k*#pp_dummy c.cogs*#pp_dummy
* i.year, by nace2).
*-----------------------------------------------------------------------
cap drop phi Lphi Lk Lcogs Lpp yr_nace const
gen phi = .
levelsof nace2, local(nace_list_init)
foreach n of local nace_list_init {
    xi: qui reg go c.k*#pp_dummy c.cogs*#pp_dummy i.year if nace2 == `n'
    cap drop phi_tmp_`n'
    predict phi_tmp_`n' if nace2 == `n'
    replace phi = phi_tmp_`n' if nace2 == `n'
    drop phi_tmp_`n' _I*
}

* Lags via the panel structure (xtset above)
gen Lphi  = L.phi
gen Lk    = L.k
gen Lcogs = L.cogs
gen Lpp   = L.pp_dummy

* Year × NACE absorbing interaction (for reghdfe / areg)
egen yr_nace = group(year nace2)

gen const = 1

* Sanity check (should never trigger after the block above)
cap confirm variable phi Lphi k cogs Lk Lcogs Lpp pp_dummy yr_nace
if _rc {
    dis as error "Required variables missing after regeneration — aborting."
    exit 198
}

* CWDL Mata GMM programs
clear mata
mata: mata set matastrict off
mata:

// Plain ACF (no pp in Markov) — CWDL original
// Note: X and X_lag carry only (k, cogs), no const. The intercept is soaked
// up by the second-stage Markov regression's constant column C, so putting a
// const in the first-stage X creates identification redundancy (and mismatches
// the 2-element optimizer init below, producing a Mata conformability error).
void CritACF_plain(todo, b, crit, g, H)
{
    PHI=st_data(.,("phi")); PHI_LAG=st_data(.,("Lphi"))
    Z=st_data(.,("const","k","Lcogs"))
    X=st_data(.,("k","cogs"))
    X_lag=st_data(.,("Lk","Lcogs"))
    W=invsym(Z'Z)
    C=st_data(.,("const"))
    OMEGA=PHI-X*b'; OL=PHI_LAG-X_lag*b'
    P=(C,OL); gb=invsym(P'P)*P'OMEGA; XI=OMEGA-P*gb
    crit=(Z'XI)'*W*(Z'XI)
}

// Base ACF (pp in Markov) — De Loecker (2013) extension
void CritACF_base(todo, b, crit, g, H)
{
    PHI=st_data(.,("phi")); PHI_LAG=st_data(.,("Lphi"))
    Z=st_data(.,("const","k","Lcogs"))
    X=st_data(.,("k","cogs"))
    X_lag=st_data(.,("Lk","Lcogs"))
    W=invsym(Z'Z)
    C=st_data(.,("const")); PP=st_data(.,("Lpp"))
    OMEGA=PHI-X*b'; OL=PHI_LAG-X_lag*b'
    P=(C,OL,PP); gb=invsym(P'P)*P'OMEGA; XI=OMEGA-P*gb
    crit=(Z'XI)'*W*(Z'XI)
}

// Run optimizer
void run_cwdl(pointer(function) scalar fn, string scalar label)
{
    S = optimize_init()
    optimize_init_params(S, (0.05, 0.95))
    optimize_init_evaluator(S, fn)
    optimize_init_which(S, "min")
    optimize_init_conv_warning(S, "off")
    optimize_init_technique(S, "nm")
    optimize_init_nmsimplexdeltas(S, 0.1)
    optimize_init_conv_maxiter(S, 5000)
    p = optimize(S)
    // Refine
    optimize_init_nmsimplexdeltas(S, 1e-5)
    optimize_init_params(S, p)
    p = optimize(S)
    optimize_init_params(S, p)
    p = optimize(S)

    cb = optimize_result_value(S)
    printf("\n  %s: b_k=%9.4f, b_cogs=%9.4f, crit=%12.6f\n", label, p[1], p[2], cb)
    st_matrix("cwdl_" + label, p)
}

end

* Run for each NACE
levelsof nace2, local(nace_list)
foreach n of local nace_list {
    preserve
    keep if nace2 == `n'
    dis _newline "=== NACE `n' (N = " _N ") ==="

    * Ensure estimation sample (need lags)
    drop if missing(phi, Lphi, k, cogs, Lk, Lcogs)

    mata: run_cwdl(&CritACF_plain(), "plain_`n'")
    mata: run_cwdl(&CritACF_base(), "base_`n'")

    * Compute markups from CWDL estimates
    mata: b_base = st_matrix("cwdl_base_`n'")
    mata: st_numscalar("bk", b_base[1])
    mata: st_numscalar("bc", b_base[2])

    gen markup_cwdl = scalar(bc) / (exp(cogs) / exp(go))
    gen lmu_cwdl = ln(markup_cwdl) if markup_cwdl > 0

    * Premium
    reghdfe lmu_cwdl pp_dummy k cogs, absorb(yr_nace) cluster(id)
    dis "  CWDL premium (NACE `n'): " _b[pp_dummy] " (SE " _se[pp_dummy] ")"

    drop markup_cwdl lmu_cwdl
    restore
}

dis _newline "--- table_cwdl_acf.do DONE ---"
