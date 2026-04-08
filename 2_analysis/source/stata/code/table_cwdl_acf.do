*===============================================================================
* table_cwdl_acf.do — CWDL (2020) ACF cross-validation
*
* Independent Stata ACF estimation using Collard-Wexler & De Loecker (2020)
* Mata GMM procedure. Cross-validates Python acf_estimator.py results.
*===============================================================================

dis _newline "--- table_cwdl_acf.do ---"

use "$data/analysis_panel.dta", clear
xtset id year, yearly

* Ensure variables exist
cap confirm variable phi Lphi k cogs Lk Lcogs pp_dummy
if _rc {
    dis as error "Required variables missing. Run estimate_pf.do first."
    exit 198
}

gen const = 1

* CWDL Mata GMM programs
clear mata
mata:

// Plain ACF (no pp in Markov) — CWDL original
void CritACF_plain(todo, b, crit, g, H)
{
    PHI=st_data(.,("phi")); PHI_LAG=st_data(.,("Lphi"))
    Z=st_data(.,("const","k","Lcogs")); X=st_data(.,("const","k","cogs"))
    X_lag=st_data(.,("const","Lk","Lcogs")); W=invsym(Z'Z)
    C=st_data(.,("const"))
    OMEGA=PHI-X*b'; OL=PHI_LAG-X_lag*b'
    P=(C,OL); gb=invsym(P'P)*P'OMEGA; XI=OMEGA-P*gb
    crit=(Z'XI)'*W*(Z'XI)
}

// Base ACF (pp in Markov) — De Loecker (2013) extension
void CritACF_base(todo, b, crit, g, H)
{
    PHI=st_data(.,("phi")); PHI_LAG=st_data(.,("Lphi"))
    Z=st_data(.,("const","k","Lcogs")); X=st_data(.,("const","k","cogs"))
    X_lag=st_data(.,("const","Lk","Lcogs")); W=invsym(Z'Z)
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
