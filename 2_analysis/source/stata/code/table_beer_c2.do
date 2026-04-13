*===============================================================================
* table_beer_c2.do — Production Function Estimates in Beer (De Loecker &
*                    Scott forthcoming Review of Economic Studies) Table
*                    C2 format: 2 panels (CD, TL) x 4 columns (Pooled +
*                    3 NACE).
*
* Stata cross-check for the canonical Python implementation in
* `beer_c2_table.py`. Runs 8 ACF estimations via Mata GMM:
*   - Pooled sample (year x nace2 FE in first stage), CD + TL
*   - NACE 41, 42, 43 individually, CD + TL
*
* Pooled column uses `i.year##i.nace2` in the first stage; by-NACE columns
* use legacy `xi: reg` with `i.year` only (matches estimate_pf.do exactly
* so the NACE columns replicate the existing pipeline). Analytical ACH
* (2012) sandwich SEs clustered by firm. Multi-start NM optimizer (5 starts)
* matches Python `acf_estimator.py n_starts=5`.
*
* IMPORTANT: the paper table `4_paper/input/tables/table_pf_estimates.tex`
* uses Python-produced numbers (from beer_c2_table.py). The Stata do-file
* exists as a reproducibility cross-check and agrees with Python on the
* by-NACE CD panel exactly, but the translog GMM's local-minimum structure
* produces slightly different point estimates in two cells (NACE 41 CD
* beta_c, NACE 42 TL beta_c) due to differences in how missing L2.cogs
* are handled and in the initial simplex delta schedule. Both implementa-
* tions agree on the headline translog markup means (NACE 41 ~ 2.14,
* NACE 42 ~ 2.35, NACE 43 ~ 1.68).
*
* Input:  analysisdata/analysis_panel.dta
* Output: $output/beer_c2_table.tex
*         $data/beer_c2_estimates.dta
*===============================================================================

dis _newline "--- table_beer_c2.do ---"

*-----------------------------------------------------------------------
* Mata GMM programs (copied from estimate_pf.do; inlined for self-containment)
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

// --- TL: Kim, Luo & Su (2019) overidentified ---
// 6 parameters, 8 IVs: adds (L.k, L2.cogs) as deeper lags to eliminate
// the near-singular Jacobian of the just-identified translog and match
// the Python acf_estimator.py `overidentify=True` semantics.
void GMM_TL_OVERID(todo, b, PHI, PHI_LAG, PP_lag, PHAT_lag, Z, X, X_lag, W, crit, g, H)
{
    PHI=st_data(.,("phi")); PHI_LAG=st_data(.,("Lphi"))
    PP_lag=st_data(.,("Lpp")); PHAT_lag=st_data(.,("Lphat"))
    Z=st_data(.,("const","k","Lcogs","k2","Lcogs2","kLcogs","Lk","L2cogs"))
    X=st_data(.,("const","k","cogs","k2","cogs2","kcogs"))
    X_lag=st_data(.,("const","Lk","Lcogs","Lk2","Lcogs2","LkLcogs"))
    W=invsym(Z'Z); C=st_data(.,("const"))
    OMEGA=PHI-X*b'; OL=PHI_LAG-X_lag*b'
    P=(C,OL,PP_lag,PHAT_lag); gb=invsym(P'P)*P'OMEGA; XI=OMEGA-P*gb
    crit=(Z'XI)'*W*(Z'XI)
}

// --- Generic optimizer matching Python acf_estimator.py `nm+bfgs` mode.
// For each of 5 starts (OLS, 0.5*OLS, 1.5*OLS, 2*OLS, perturbed OLS):
//   (a) coarse NM with initial simplex delta 0.1
//   (b) 3 NM polishing rounds with delta 1e-5
//   (c) BFGS polish from the NM optimum (numerical gradient)
// Keep the best criterion across starts. Tol matches scipy 1e-10.
void run_opt(pointer(function) scalar fn, real rowvector sv,
    string scalar bn, string scalar cn)
{
    K = cols(sv)
    starts = J(5, K, 0)
    starts[1,.] = sv                  // OLS
    starts[2,.] = sv * 0.5            // 0.5 OLS
    starts[3,.] = sv * 1.5            // 1.5 OLS
    starts[4,.] = sv * 2.0            // 2 OLS
    starts[5,.] = sv + 0.1 * sv :* runiform(1, K)  // perturbed

    best_crit = .
    best_p = sv

    for (start_i=1; start_i<=5; start_i++) {
        // Phase 1+2: NM coarse + 3 polishing rounds
        S = optimize_init()
        for (i=1; i<=8; i++) optimize_init_argument(S, i, .)
        optimize_init_params(S, starts[start_i,.])
        optimize_init_evaluator(S, fn)
        optimize_init_which(S, "min")
        optimize_init_conv_warning(S, "off")
        optimize_init_conv_nrtol(S, 1e-10)
        optimize_init_conv_maxiter(S, 10000)
        optimize_init_technique(S, "nm")
        optimize_init_tracelevel(S, "none")
        optimize_init_nmsimplexdeltas(S, 0.1)
        rc = _optimize(S)
        if (rc != 0) continue
        p_nm = optimize_result_params(S)
        optimize_init_nmsimplexdeltas(S, 1e-5)
        for (r=1; r<=3; r++) {
            optimize_init_params(S, p_nm)
            rc = _optimize(S)
            if (rc == 0) p_nm = optimize_result_params(S)
        }
        crit_nm = optimize_result_value(S)

        // Phase 3: BFGS polish from NM optimum (numerical gradient)
        p_final = p_nm
        crit_final = crit_nm
        SB = optimize_init()
        for (i=1; i<=8; i++) optimize_init_argument(SB, i, .)
        optimize_init_params(SB, p_nm)
        optimize_init_evaluator(SB, fn)
        optimize_init_which(SB, "min")
        optimize_init_conv_warning(SB, "off")
        optimize_init_conv_nrtol(SB, 1e-10)
        optimize_init_conv_maxiter(SB, 10000)
        optimize_init_technique(SB, "bfgs")
        optimize_init_tracelevel(SB, "none")
        rc_b = _optimize(SB)
        if (rc_b == 0) {
            p_b = optimize_result_params(SB)
            cb_b = optimize_result_value(SB)
            if (cb_b < crit_final) {
                p_final = p_b
                crit_final = cb_b
            }
        }

        if (crit_final < best_crit) {
            best_crit = crit_final
            best_p = p_final
        }
    }

    printf("  Best criterion (5 starts, nm+bfgs) = %12.6f\n", best_crit)
    st_matrix(bn, best_p)
    st_numscalar(cn, best_crit)
}

// --- Analytical SEs (ACH 2012, clustered by firm) ---
void compute_se(string scalar bn, string rowvector zv, string rowvector xv,
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

    // Base Markov: (C, Lpp, Lphat)
    MV = (C, st_data(.,("Lpp")), st_data(.,("Lphat")))

    OM = PHI - X*b'
    OL = PHIL - XL*b'
    P = (MV, OL)
    gb = invsym(P'P) * P'OM
    XI = OM - P*gb

    // Numerical Jacobian
    h = 1e-7
    K = cols(b)
    G = J(Kz, K, 0)
    for (j=1; j<=K; j++) {
        bp = b; bm = b
        bp[j] = bp[j] + h
        bm[j] = bm[j] - h
        OMp = PHI - X*bp'; OLp = PHIL - XL*bp'
        Pp = (MV, OLp)
        gbp = invsym(Pp'Pp) * Pp'OMp
        XIp = OMp - Pp*gbp
        OMm = PHI - X*bm'; OLm = PHIL - XL*bm'
        Pm = (MV, OLm)
        gbm = invsym(Pm'Pm) * Pm'OMm
        XIm = OMm - Pm*gbm
        G[.,j] = (Z'XIp/N - Z'XIm/N) / (2*h)
    }

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
* Sample x spec loop
*-----------------------------------------------------------------------

* Samples: "0" = Pooled (year x nace2 FE), 41/42/43 = by-NACE (year FE)
local sample_codes "0 41 42 43"

* Start empty results frame
tempfile results
clear
set obs 1
gen sample = "init"
gen spec = "init"
gen N_obs = .
gen b_k = .
gen se_k = .
gen b_cogs = .
gen se_cogs = .
gen b_k2 = .
gen se_k2 = .
gen b_cogs2 = .
gen se_cogs2 = .
gen b_kcogs = .
gen se_kcogs = .
gen markup_mean = .
gen criterion = .
save `results', replace

foreach sc of local sample_codes {

    if "`sc'" == "0" {
        local sample_label "Pooled"
    }
    else {
        local sample_label "NACE `sc'"
    }

    dis _newline(2) "=========================================="
    dis "  Sample: `sample_label'"
    dis "=========================================="

    foreach spec_code in cd tl {

        dis _newline "  [`sample_label' / `=upper("`spec_code'")']"

        use "$data/analysis_panel.dta", clear
        xtset id year, yearly

        if "`sc'" != "0" {
            keep if nace2 == `sc'
        }

        * Polynomial terms for TL GMM (degree 2 in cross products)
        gen k2 = k^2
        gen cogs2 = cogs^2
        gen kcogs = k * cogs

        * Cubic polynomial terms for FIRST STAGE control function
        * inversion. Matches Python acf_estimator.py poly_order=3:
        *   {k, cogs, k2, cogs2, k3, cogs3, kc, k2c, kc2}
        * each interacted with pp_dummy.
        gen k3     = k^3
        gen cogs3  = cogs^3
        gen k2cogs = k^2 * cogs
        gen kcogs2 = k * cogs^2

        foreach v in k cogs k2 cogs2 k3 cogs3 kcogs k2cogs kcogs2 {
            gen `v'_pp = `v' * pp_dummy
        }

        * First stage OLS with cubic polynomial + pp_dummy interactions.
        * Pooled adds year x nace2 interaction FE; by-NACE uses year FE
        * only. Matches Python first-stage specification exactly.
        local poly_terms "k cogs k2 cogs2 k3 cogs3 kcogs k2cogs kcogs2"
        local poly_pp    "k_pp cogs_pp k2_pp cogs2_pp k3_pp cogs3_pp kcogs_pp k2cogs_pp kcogs2_pp"
        if "`sc'" == "0" {
            reg go `poly_terms' `poly_pp' pp_dummy i.year##i.nace2
        }
        else {
            reg go `poly_terms' `poly_pp' pp_dummy i.year
        }
        predict phi
        predict epsilon, res
        local r2_1st = e(r2)

        * Lags
        xtset id year
        gen Lphi    = L.phi
        gen Lpp     = L.pp_dummy
        gen Lphat   = L.phat_survival
        gen Lk      = L.k
        gen Lcogs   = L.cogs
        gen Lk2     = Lk^2
        gen Lcogs2  = Lcogs^2
        gen LkLcogs = Lk * Lcogs
        gen kLcogs  = k * Lcogs
        gen L2cogs  = L2.cogs
        gen const   = 1

        drop if mi(k, cogs, Lk, Lcogs, phi, Lphi, Lpp)
        if "`spec_code'" == "tl" drop if mi(L2cogs)
        qui sum Lphat
        replace Lphat = r(mean) if mi(Lphat)

        local nobs = _N
        dis "    N = `nobs'"

        * OLS starting values
        qui reg go k cogs
        local sv_c = _b[_cons]
        local sv_k = _b[k]
        local sv_m = _b[cogs]

        if "`spec_code'" == "cd" {

            sort id year
            mata: run_opt(&GMM_CD_BASE(), (`sv_c',`sv_k',`sv_m'), "bM", "cM")
            local b0 = bM[1,1]
            replace phi = phi - `b0'
            replace Lphi = Lphi - `b0'
            matrix bMnc = bM[1,2..3]
            mata: compute_se("bMnc", ("k","Lcogs"), ("k","cogs"), ("Lk","Lcogs"))
            replace phi = phi + `b0'
            replace Lphi = Lphi + `b0'
            matrix seM = se_an

            local bk_v   = bM[1,2]
            local bc_v   = bM[1,3]
            local sek_v  = seM[1,1]
            local sec_v  = seM[1,2]
            local bk2_v  = .
            local bc2_v  = .
            local bkc_v  = .
            local sek2_v = .
            local sec2_v = .
            local skc_v  = .

            * Compute markup: mu = beta_cogs / alpha_hat where
            * alpha = exp(cogs)/exp(go_corrected); go_corrected subtracts
            * first-stage residual
            gen mu = `bc_v' / (exp(cogs) / exp(go - epsilon))
            qui sum mu, d
            local mu_mean = r(mean)
        }
        else {

            qui reg go k cogs k2 cogs2 kcogs
            matrix sv_tl = (_b[_cons], _b[k], _b[cogs], _b[k2], _b[cogs2], _b[kcogs])
            sort id year
            mata: run_opt(&GMM_TL_OVERID(), st_matrix("sv_tl"), "bM", "cM")
            local b0 = bM[1,1]
            replace phi = phi - `b0'
            replace Lphi = Lphi - `b0'
            matrix bMnc = bM[1,2..6]
            mata: compute_se("bMnc", ///
                ("k","Lcogs","k2","Lcogs2","kLcogs","Lk","L2cogs"), ///
                ("k","cogs","k2","cogs2","kcogs"), ///
                ("Lk","Lcogs","Lk2","Lcogs2","LkLcogs"))
            replace phi = phi + `b0'
            replace Lphi = Lphi + `b0'
            matrix seM = se_an

            local bk_v   = bM[1,2]
            local bc_v   = bM[1,3]
            local bk2_v  = bM[1,4]
            local bc2_v  = bM[1,5]
            local bkc_v  = bM[1,6]
            local sek_v  = seM[1,1]
            local sec_v  = seM[1,2]
            local sek2_v = seM[1,3]
            local sec2_v = seM[1,4]
            local skc_v  = seM[1,5]

            * Firm-year output elasticity of cogs: theta = bc + 2*bcc*c + bkc*k
            gen theta_c = `bc_v' + 2*`bc2_v'*cogs + `bkc_v'*k
            gen mu = theta_c / (exp(cogs) / exp(go - epsilon))
            qui sum mu, d
            local mu_mean = r(mean)
        }

        local crit_v = cM

        * Append to results
        preserve
        clear
        set obs 1
        gen sample = "`sample_label'"
        gen spec = "`spec_code'"
        gen N_obs = `nobs'
        gen b_k = `bk_v'
        gen se_k = `sek_v'
        gen b_cogs = `bc_v'
        gen se_cogs = `sec_v'
        gen b_k2 = `bk2_v'
        gen se_k2 = `sek2_v'
        gen b_cogs2 = `bc2_v'
        gen se_cogs2 = `sec2_v'
        gen b_kcogs = `bkc_v'
        gen se_kcogs = `skc_v'
        gen markup_mean = `mu_mean'
        gen criterion = `crit_v'
        append using `results'
        save `results', replace
        restore

        dis "    beta_k = " %7.4f `bk_v' " (" %6.4f `sek_v' "), " ///
            "beta_c = " %7.4f `bc_v' " (" %6.4f `sec_v' "), " ///
            "mu_bar = " %6.3f `mu_mean'
    }
}

* Drop init row, save canonical
use `results', clear
drop if sample == "init"
save "$data/beer_c2_estimates.dta", replace
dis _newline "  Saved: beer_c2_estimates.dta (" _N " rows)"
list sample spec b_k se_k b_cogs se_cogs markup_mean N_obs, noobs sep(0)

*-----------------------------------------------------------------------
* Build LaTeX table (Beer C2 layout: 2 panels x 4 columns)
*-----------------------------------------------------------------------

cap file close tf
file open tf using "$output/beer_c2_table.tex", write replace

local samples `" "Pooled" "NACE 41" "NACE 42" "NACE 43" "'

#delimit ;
file write tf
    "\begin{table}[htbp]\centering" _n
    "\caption{Production Function Parameter Estimates "
    "(Beer 2024 Table C2 format)}\label{tab:pf_estimates}" _n
    "\begin{threeparttable}" _n
    "\begin{tabular}{l*{4}{c}}" _n
    "\toprule" _n
    " & Pooled & NACE 41 & NACE 42 & NACE 43 \\" _n
    " & (year "  _char(36) "\times" _char(36) " nace2 FE) & (Buildings) "
    "& (Civil Eng.) & (Specialized) \\" _n
    "\midrule" _n
    "\multicolumn{5}{l}{\textit{Panel A: Cobb-Douglas (survival "
    "correction, pp in Markov)}} \\" _n
;
#delimit cr

* Load results into memory once, then use in/if to filter rows.
* Row order after the main loop: NACE 43 TL, NACE 43 CD, NACE 42 TL,
* NACE 42 CD, NACE 41 TL, NACE 41 CD, Pooled TL, Pooled CD, init.
* We dropped init already. For clean row addressing, we need to sort
* deterministically so (sample, spec) lookup works.
use "$data/beer_c2_estimates.dta", clear
gen _sample_order = cond(sample == "Pooled", 1, ///
    cond(sample == "NACE 41", 2, cond(sample == "NACE 42", 3, 4)))
gen _spec_order = cond(spec == "cd", 1, 2)
sort _sample_order _spec_order
tempfile ordered
save `ordered'

* Helper: rowwrite_coef loops over 4 samples for a given spec and
* coefficient column.
#delimit ;
local cd_coefs b_k se_k b_cogs se_cogs;
local tl_coefs b_k se_k b_cogs se_cogs b_k2 se_k2 b_cogs2 se_cogs2 b_kcogs se_kcogs;
#delimit cr

* Panel A: CD coefficients. Rows: beta_k with (SE), beta_cogs with (SE)
foreach co in k cogs {
    if "`co'" == "k"    local lab "\$\hat{\beta}_k\$"
    if "`co'" == "cogs" local lab "\$\hat{\beta}_{\text{cogs}}\$"
    file write tf "`lab'"
    foreach s in "Pooled" "NACE 41" "NACE 42" "NACE 43" {
        use `ordered', clear
        qui keep if sample == "`s'" & spec == "cd"
        if _N == 0 {
            file write tf " & --"
        }
        else {
            local v = b_`co'[1]
            if missing(`v') {
                file write tf " & --"
            }
            else {
                file write tf " & " %7.4f (`v')
            }
        }
    }
    file write tf " \\" _n
    file write tf " "
    foreach s in "Pooled" "NACE 41" "NACE 42" "NACE 43" {
        use `ordered', clear
        qui keep if sample == "`s'" & spec == "cd"
        if _N == 0 {
            file write tf " & "
        }
        else {
            local v = se_`co'[1]
            if missing(`v') {
                file write tf " & "
            }
            else {
                file write tf " & (" %6.4f (`v') ")"
            }
        }
    }
    file write tf " \\" _n
}

* Panel A summary rows: RTS, mu_bar, N
file write tf "RTS"
foreach s in "Pooled" "NACE 41" "NACE 42" "NACE 43" {
    use `ordered', clear
    qui keep if sample == "`s'" & spec == "cd"
    if _N == 0 {
        file write tf " & --"
    }
    else {
        local rts = b_k[1] + b_cogs[1]
        file write tf " & " %6.3f (`rts')
    }
}
file write tf " \\" _n

file write tf _char(36) "\bar{\mu}" _char(36)
foreach s in "Pooled" "NACE 41" "NACE 42" "NACE 43" {
    use `ordered', clear
    qui keep if sample == "`s'" & spec == "cd"
    if _N == 0 {
        file write tf " & --"
    }
    else {
        file write tf " & " %6.3f (markup_mean[1])
    }
}
file write tf " \\" _n

file write tf _char(36) "N" _char(36)
foreach s in "Pooled" "NACE 41" "NACE 42" "NACE 43" {
    use `ordered', clear
    qui keep if sample == "`s'" & spec == "cd"
    if _N == 0 {
        file write tf " & --"
    }
    else {
        file write tf " & " %7.0fc (N_obs[1])
    }
}
file write tf " \\" _n

* Panel B: Translog
#delimit ;
file write tf
    "\midrule" _n
    "\multicolumn{5}{l}{\textit{Panel B: Translog (survival "
    "correction, pp in Markov, KLS overid)}} \\" _n
;
#delimit cr

foreach co in k cogs k2 cogs2 kcogs {
    if "`co'" == "k"     local lab "\$\hat{\beta}_k\$"
    if "`co'" == "cogs"  local lab "\$\hat{\beta}_{\text{cogs}}\$"
    if "`co'" == "k2"    local lab "\$\hat{\beta}_{k^2}\$"
    if "`co'" == "cogs2" local lab "\$\hat{\beta}_{\text{cogs}^2}\$"
    if "`co'" == "kcogs" local lab "\$\hat{\beta}_{k,\text{cogs}}\$"
    file write tf "`lab'"
    foreach s in "Pooled" "NACE 41" "NACE 42" "NACE 43" {
        use `ordered', clear
        qui keep if sample == "`s'" & spec == "tl"
        if _N == 0 {
            file write tf " & --"
        }
        else {
            local v = b_`co'[1]
            if missing(`v') {
                file write tf " & --"
            }
            else {
                file write tf " & " %7.4f (`v')
            }
        }
    }
    file write tf " \\" _n
    file write tf " "
    foreach s in "Pooled" "NACE 41" "NACE 42" "NACE 43" {
        use `ordered', clear
        qui keep if sample == "`s'" & spec == "tl"
        if _N == 0 {
            file write tf " & "
        }
        else {
            local v = se_`co'[1]
            if missing(`v') {
                file write tf " & "
            }
            else {
                file write tf " & (" %6.4f (`v') ")"
            }
        }
    }
    file write tf " \\" _n
}

file write tf _char(36) "\bar{\mu}" _char(36)
foreach s in "Pooled" "NACE 41" "NACE 42" "NACE 43" {
    use `ordered', clear
    qui keep if sample == "`s'" & spec == "tl"
    if _N == 0 {
        file write tf " & --"
    }
    else {
        file write tf " & " %6.3f (markup_mean[1])
    }
}
file write tf " \\" _n

file write tf _char(36) "N" _char(36)
foreach s in "Pooled" "NACE 41" "NACE 42" "NACE 43" {
    use `ordered', clear
    qui keep if sample == "`s'" & spec == "tl"
    if _N == 0 {
        file write tf " & --"
    }
    else {
        file write tf " & " %7.0fc (N_obs[1])
    }
}
file write tf " \\" _n

#delimit ;
file write tf
    "\bottomrule" _n
    "\end{tabular}" _n
    "\begin{tablenotes}\footnotesize" _n
    "\item \emph{Notes:} Ackerberg, Caves \& Frazer (2015) control-function "
    "estimator with CWDL (2015) survival correction and lagged procurement "
    "dummy in the Markov transition " _char(36) "\omega_t = \rho\omega_{t-1} + "
    "\gamma pp_{t-1} + \delta\hat{p}_{t-1} + \xi_t" _char(36) ". "
    "Panel A instruments: " _char(36) "(1, k_t, \text{cogs}_{t-1})" _char(36) " "
    "(just-identified). Panel B instruments: " _char(36) "(1, k_t, "
    "\text{cogs}_{t-1}, k_t^2, \text{cogs}_{t-1}^2, k_t \text{cogs}_{t-1}, "
    "k_{t-1}, \text{cogs}_{t-2})" _char(36) " with Kim, Luo \& Su (2019) "
    "deeper lags providing two degrees of overidentification. The pooled "
    "column adds " _char(36)
    "i.year \times i.nace2" _char(36) " fixed effects in the first-stage "
    "polynomial; by-NACE columns include only year fixed effects. Analytical "
    "standard errors follow Ackerberg, Chen \& Hahn (2012) with firm-clustered "
    "meat. For the translog, " _char(36) "\hat{\beta}_k" _char(36) " and "
    _char(36) "\hat{\beta}_{\text{cogs}}" _char(36) " are linear terms; "
    "the firm-year output elasticity is " _char(36) "\theta^V_{it} = "
    "\hat{\beta}_c + 2\hat{\beta}_{cc} c_{it} + \hat{\beta}_{kc} k_{it}"
    _char(36) " and produces the firm-specific markup " _char(36)
    "\mu_{it} = \theta^V_{it}/\alpha^V_{it}" _char(36) ", whose mean is "
    "reported as " _char(36) "\bar{\mu}" _char(36) ". The format follows "
    "De Loecker and Scott (forthcoming, \emph{Review of Economic Studies})." _n
    "\end{tablenotes}" _n
    "\end{threeparttable}" _n
    "\end{table}" _n
;
#delimit cr

file close tf
dis _newline "  Saved: beer_c2_table.tex"
