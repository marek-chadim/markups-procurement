*===============================================================================
* estimate_pf.do — ACF Production Function Estimation
*
* Analog of DGM markup_estimation_dgmacf.do:
*   DGM: 4-input TL/CD by industry via Mata GMM (21 industries)
*   Us:  2-input CD/TL by nace2 via Mata GMM (3 industries)
*        with CWDL (2015) extensions (survival, pp in Markov)
*
* Specifications estimated:
*   A. Base: CD, linear Markov with pp_lag + survival
*   B. No survival: CD, linear Markov with pp_lag
*   C. No pp Markov: CD, linear Markov with survival only
*   D. Plain: CD, linear Markov only
*   E. Translog: TL, linear Markov with pp_lag + survival
*
* Input:  analysisdata/analysis_panel.dta
* Output: analysisdata/coefficients_byind.dta (one row per nace2 × spec)
*===============================================================================

dis _newline "--- estimate_pf.do ---"

use "$data/analysis_panel.dta", clear
xtset id year, yearly

levelsof nace2, local(nace2_list)

*-----------------------------------------------------------------------
* Mata GMM programs (4 Markov variants + translog)
*-----------------------------------------------------------------------

clear mata
mata:

// --- CD: base (pp + survival in Markov) ---
void GMM_CD_BASE(todo, b, PHI, PHI_LAG, PP_lag, PHAT_lag, Z, X, X_lag, W, crit, g, H)
{
    PHI=st_data(.,("phi")); PHI_LAG=st_data(.,("Lphi"))
    PP_lag=st_data(.,("Lpp")); PHAT_lag=st_data(.,("Lphat"))
    Z=st_data(.,("const","k","Lcogs")); X=st_data(.,("const","k","cogs"))
    X_lag=st_data(.,("const","Lk","Lcogs")); W=invsym(Z'Z)
    C=st_data(.,("const"))
    OMEGA=PHI-X*b'; OL=PHI_LAG-X_lag*b'
    P=(C,OL,PP_lag,PHAT_lag); gb=invsym(P'P)*P'OMEGA; XI=OMEGA-P*gb
    crit=(Z'XI)'*W*(Z'XI)
}

// --- CD: pp only ---
void GMM_CD_PP(todo, b, PHI, PHI_LAG, PP_lag, PHAT_lag, Z, X, X_lag, W, crit, g, H)
{
    PHI=st_data(.,("phi")); PHI_LAG=st_data(.,("Lphi"))
    PP_lag=st_data(.,("Lpp"))
    Z=st_data(.,("const","k","Lcogs")); X=st_data(.,("const","k","cogs"))
    X_lag=st_data(.,("const","Lk","Lcogs")); W=invsym(Z'Z)
    C=st_data(.,("const"))
    OMEGA=PHI-X*b'; OL=PHI_LAG-X_lag*b'
    P=(C,OL,PP_lag); gb=invsym(P'P)*P'OMEGA; XI=OMEGA-P*gb
    crit=(Z'XI)'*W*(Z'XI)
}

// --- CD: survival only ---
void GMM_CD_SURV(todo, b, PHI, PHI_LAG, PP_lag, PHAT_lag, Z, X, X_lag, W, crit, g, H)
{
    PHI=st_data(.,("phi")); PHI_LAG=st_data(.,("Lphi"))
    PHAT_lag=st_data(.,("Lphat"))
    Z=st_data(.,("const","k","Lcogs")); X=st_data(.,("const","k","cogs"))
    X_lag=st_data(.,("const","Lk","Lcogs")); W=invsym(Z'Z)
    C=st_data(.,("const"))
    OMEGA=PHI-X*b'; OL=PHI_LAG-X_lag*b'
    P=(C,OL,PHAT_lag); gb=invsym(P'P)*P'OMEGA; XI=OMEGA-P*gb
    crit=(Z'XI)'*W*(Z'XI)
}

// --- CD: plain ---
void GMM_CD_PLAIN(todo, b, PHI, PHI_LAG, PP_lag, PHAT_lag, Z, X, X_lag, W, crit, g, H)
{
    PHI=st_data(.,("phi")); PHI_LAG=st_data(.,("Lphi"))
    Z=st_data(.,("const","k","Lcogs")); X=st_data(.,("const","k","cogs"))
    X_lag=st_data(.,("const","Lk","Lcogs")); W=invsym(Z'Z)
    C=st_data(.,("const"))
    OMEGA=PHI-X*b'; OL=PHI_LAG-X_lag*b'
    P=(C,OL); gb=invsym(P'P)*P'OMEGA; XI=OMEGA-P*gb
    crit=(Z'XI)'*W*(Z'XI)
}

// --- TL: base ---
void GMM_TL_BASE(todo, b, PHI, PHI_LAG, PP_lag, PHAT_lag, Z, X, X_lag, W, crit, g, H)
{
    PHI=st_data(.,("phi")); PHI_LAG=st_data(.,("Lphi"))
    PP_lag=st_data(.,("Lpp")); PHAT_lag=st_data(.,("Lphat"))
    Z=st_data(.,("const","k","Lcogs","k2","Lcogs2","kLcogs"))
    X=st_data(.,("const","k","cogs","k2","cogs2","kcogs"))
    X_lag=st_data(.,("const","Lk","Lcogs","Lk2","Lcogs2","LkLcogs"))
    W=invsym(Z'Z); C=st_data(.,("const"))
    OMEGA=PHI-X*b'; OL=PHI_LAG-X_lag*b'
    P=(C,OL,PP_lag,PHAT_lag); gb=invsym(P'P)*P'OMEGA; XI=OMEGA-P*gb
    crit=(Z'XI)'*W*(Z'XI)
}

// --- Generic optimizer (NM only, robust) ---
void run_opt(pointer(function) scalar fn, real rowvector sv,
    string scalar bn, string scalar cn)
{
    S = optimize_init()
    for (i=1; i<=8; i++) optimize_init_argument(S, i, .)
    optimize_init_params(S, sv)
    optimize_init_evaluator(S, fn)
    optimize_init_which(S, "min")
    optimize_init_conv_warning(S, "off")
    optimize_init_conv_nrtol(S, 1e-10)
    optimize_init_conv_maxiter(S, 10000)
    optimize_init_technique(S, "nm")
    optimize_init_tracelevel(S, "none")

    // Coarse NM
    optimize_init_nmsimplexdeltas(S, 0.1)
    p = optimize(S)

    // Refine NM (3 rounds with tight deltas)
    optimize_init_nmsimplexdeltas(S, 1e-5)
    optimize_init_params(S, p)
    p = optimize(S)
    optimize_init_params(S, p)
    p = optimize(S)
    optimize_init_params(S, p)
    p = optimize(S)

    cb = optimize_result_value(S)
    printf("  Criterion = %12.6f\n", cb)
    printf("  Params: ")
    for (j=1; j<=cols(p); j++) printf("%9.4f ", p[j])
    printf("\n")
    st_matrix(bn, p)
    st_numscalar(cn, cb)
    printf("  Matrix %s saved.\n", bn)
}

// --- Analytical SEs (ACH 2012, clustered) ---
void compute_se(string scalar bn, string scalar mt,
    string rowvector zv, string rowvector xv,
    string rowvector xlv)
{
    PHI = st_data(.,("phi"))
    PHIL = st_data(.,("Lphi"))
    C = st_data(.,("const"))
    Z = st_data(.,zv)
    X = st_data(.,xv)
    XL = st_data(.,xlv)
    cl = st_data(.,("id"))
    b = st_matrix(bn)
    N = rows(X)
    Kz = cols(Z)

    if (mt=="base") {
        MV = (C, st_data(.,("Lpp")), st_data(.,("Lphat")))
    }
    else if (mt=="pp") {
        MV = (C, st_data(.,("Lpp")))
    }
    else if (mt=="surv") {
        MV = (C, st_data(.,("Lphat")))
    }
    else {
        MV = C
    }

    // Omega and XI
    OM = PHI - X*b'
    OL = PHIL - XL*b'
    P = (MV, OL)
    gb = invsym(P'P) * P'OM
    XI = OM - P*gb

    // Jacobian (numerical)
    h = 1e-7
    K = cols(b)
    G = J(Kz, K, 0)
    for (j=1; j<=K; j++) {
        bp = b
        bm = b
        bp[j] = bp[j] + h
        bm[j] = bm[j] - h
        OMp = PHI - X*bp'
        OLp = PHIL - XL*bp'
        Pp = (MV, OLp)
        gbp = invsym(Pp'Pp) * Pp'OMp
        XIp = OMp - Pp*gbp
        OMm = PHI - X*bm'
        OLm = PHIL - XL*bm'
        Pm = (MV, OLm)
        gbm = invsym(Pm'Pm) * Pm'OMm
        XIm = OMm - Pm*gbm
        G[.,j] = (Z'XIp/N - Z'XIm/N) / (2*h)
    }

    // Clustered meat (ACH 2012)
    W = invsym(Z'Z/N)
    ZXI = Z :* XI
    cids = uniqrows(cl)
    Nc = rows(cids)
    S = J(Kz, Kz, 0)
    for (c=1; c<=Nc; c++) {
        sel = selectindex(cl :== cids[c])
        mc = colsum(ZXI[sel,.])
        S = S + mc'*mc
    }
    S = S/N * (Nc/(Nc-1))
    GWGi = invsym(G'*W*G)
    V = GWGi * G'*W*S*W*G * GWGi / N
    st_matrix("V_an", V)
    st_matrix("se_an", sqrt(diagonal(V))')
}

end

*-----------------------------------------------------------------------
* Industry loop
*-----------------------------------------------------------------------

local row = 0

foreach n of local nace2_list {

    dis _newline(2) "=========================================="
    dis "  NACE2 = `n'"
    dis "=========================================="

    preserve
    keep if nace2 == `n'

    * Polynomial terms
    gen k2 = k^2
    gen cogs2 = cogs^2
    gen kcogs = k * cogs

    * First stage OLS (with pp_dummy interactions + year FE)
    xi: reg go c.k*#pp_dummy c.cogs*#pp_dummy i.year
    predict phi
    predict epsilon, res
    local r2_1st = e(r2)

    * Corrected expenditure share
    gen y_c = go - epsilon
    gen Y_c = exp(y_c)
    gen alphahat = exp(cogs) / Y_c

    * Lags
    xtset id year
    gen Lphi  = L.phi
    gen Lpp   = L.pp_dummy
    gen Lphat = L.phat_survival
    gen Lk    = L.k
    gen Lcogs = L.cogs
    gen Lk2   = Lk^2
    gen Lcogs2 = Lcogs^2
    gen LkLcogs = Lk * Lcogs
    gen kLcogs = k * Lcogs
    gen const  = 1

    drop if mi(k, cogs, Lk, Lcogs, phi, Lphi, Lpp)
    qui sum Lphat
    replace Lphat = r(mean) if mi(Lphat)

    local nobs = _N
    dis "  N = `nobs'"

    * OLS starting values
    qui reg go k cogs
    local sv_c = _b[_cons]
    local sv_k = _b[k]
    local sv_m = _b[cogs]

    * --- Spec A: Base ---
    sort id year
    mata: run_opt(&GMM_CD_BASE(), (`sv_c',`sv_k',`sv_m'), "bA", "cA")
    local b0A = bA[1,1]
    replace phi = phi - `b0A'
    replace Lphi = Lphi - `b0A'
    matrix bAnc = bA[1,2..3]
    mata: compute_se("bAnc","base",("k","Lcogs"),("k","cogs"),("Lk","Lcogs"))
    replace phi = phi + `b0A'
    replace Lphi = Lphi + `b0A'
    matrix seA = se_an

    * --- Spec B: No survival ---
    sort id year
    mata: run_opt(&GMM_CD_PP(), (`sv_c',`sv_k',`sv_m'), "bB", "cB")
    local b0B = bB[1,1]
    replace phi = phi - `b0B'
    replace Lphi = Lphi - `b0B'
    matrix bBnc = bB[1,2..3]
    mata: compute_se("bBnc","pp",("k","Lcogs"),("k","cogs"),("Lk","Lcogs"))
    replace phi = phi + `b0B'
    replace Lphi = Lphi + `b0B'
    matrix seB = se_an

    * --- Spec C: No pp ---
    sort id year
    mata: run_opt(&GMM_CD_SURV(), (`sv_c',`sv_k',`sv_m'), "bC", "cC")
    local b0C = bC[1,1]
    replace phi = phi - `b0C'
    replace Lphi = Lphi - `b0C'
    matrix bCnc = bC[1,2..3]
    mata: compute_se("bCnc","surv",("k","Lcogs"),("k","cogs"),("Lk","Lcogs"))
    replace phi = phi + `b0C'
    replace Lphi = Lphi + `b0C'
    matrix seC = se_an

    * --- Spec D: Plain ---
    sort id year
    mata: run_opt(&GMM_CD_PLAIN(), (`sv_c',`sv_k',`sv_m'), "bD", "cD")
    local b0D = bD[1,1]
    replace phi = phi - `b0D'
    replace Lphi = Lphi - `b0D'
    matrix bDnc = bD[1,2..3]
    mata: compute_se("bDnc","plain",("k","Lcogs"),("k","cogs"),("Lk","Lcogs"))
    replace phi = phi + `b0D'
    replace Lphi = Lphi + `b0D'
    matrix seD = se_an

    * --- Spec E: Translog ---
    qui reg go k cogs k2 cogs2 kcogs
    matrix sv_tl = (_b[_cons], _b[k], _b[cogs], _b[k2], _b[cogs2], _b[kcogs])
    sort id year
    mata: run_opt(&GMM_TL_BASE(), st_matrix("sv_tl"), "bE", "cE")
    local b0E = bE[1,1]
    replace phi = phi - `b0E'
    replace Lphi = Lphi - `b0E'
    matrix bEnc = bE[1,2..6]
    mata: compute_se("bEnc","base", ///
        ("k","Lcogs","k2","Lcogs2","kLcogs"), ///
        ("k","cogs","k2","cogs2","kcogs"), ///
        ("Lk","Lcogs","Lk2","Lcogs2","LkLcogs"))
    replace phi = phi + `b0E'
    replace Lphi = Lphi + `b0E'
    matrix seE = se_an

    * --- OLS ---
    qui reg go k cogs
    local ols_k = _b[k]
    local ols_m = _b[cogs]

    * --- Save coefficients ---
    clear
    set obs 1
    gen nace2 = `n'
    gen N_obs = `nobs'
    gen r2_first = `r2_1st'

    foreach s in A B C D {
        gen b_k_`s' = b`s'[1,2]
        gen b_cogs_`s' = b`s'[1,3]
        gen se_k_`s' = se`s'[1,1]
        gen se_cogs_`s' = se`s'[1,2]
        gen crit_`s' = c`s'
    }
    gen b_k_E = bE[1,2]
    gen b_cogs_E = bE[1,3]
    gen b_k2_E = bE[1,4]
    gen b_cogs2_E = bE[1,5]
    gen b_kcogs_E = bE[1,6]
    gen se_k_E = seE[1,1]
    gen se_cogs_E = seE[1,2]
    gen se_k2_E = seE[1,3]
    gen se_cogs2_E = seE[1,4]
    gen se_kcogs_E = seE[1,5]
    gen crit_E = cE
    gen b_k_OLS = `ols_k'
    gen b_cogs_OLS = `ols_m'

    save "$temp/coeff_`n'.dta", replace

    restore
}

* Combine
use "$temp/coeff_41.dta", clear
foreach n in 42 43 {
    append using "$temp/coeff_`n'.dta"
}
save "$data/coefficients_byind.dta", replace
dis _newline "  Saved: coefficients_byind.dta"
list nace2 b_k_A b_cogs_A se_k_A se_cogs_A N_obs, noobs
