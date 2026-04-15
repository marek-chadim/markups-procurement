*===============================================================================
* table_ags_twostep.do — AGS (2017) Numerical Λ Sensitivity Matrix
*
* Port of ags_twostep_identification.py. Computes the Andrews-Gentzkow-Shapiro
* (2017) scaled sensitivity matrix
*     Λ̂ = −(Ĝ' Ŵ Ĝ)^{-1} Ĝ' Ŵ diag(σ̂)
* at the translog ACF GMM optimum for NACE 41. Identifies load-bearing
* moments by magnitude of Λ̂_{k,m}.
*
* Input:  $data/coefficients_byind.dta (from estimate_pf.do) + analysis_panel
* Output: ../../output/tables/ags_twostep_sensitivity.tex
*         ../../output/tables/ags_twostep_ci.tex (identification-robust CI via S-stat)
*===============================================================================

dis _newline "--- table_ags_twostep.do ---"

cap confirm file "$data/coefficients_byind.dta"
if _rc != 0 {
    dis "  SKIP: coefficients_byind.dta not found"
    exit
}

* Load coefficients for NACE 41 TL baseline
use "$data/coefficients_byind.dta", clear
keep if nace2 == 41
if _N == 0 {
    dis "  SKIP: no NACE 41 row in coefficients_byind.dta"
    exit
}

local b_k    = b_k_E[1]
local b_cogs = b_cogs_E[1]
local b_k2   = b_k2_E[1]
local b_c2   = b_cogs2_E[1]
local b_kc   = b_kcogs_E[1]

dis "  Loaded NACE 41 TL: β_k=" %6.4f `b_k' ", β_c=" %6.4f `b_cogs'

* Re-open panel and build instruments / Jacobian numerically
use "$data/analysis_panel.dta", clear
keep if nace2 == 41
xtset id year
gen k2 = k^2
gen cogs2 = cogs^2
gen kcogs = k * cogs
gen Lk    = L.k
gen Lcogs = L.cogs
gen L2cogs = L2.cogs
gen Lk2 = Lk^2
gen Lcogs2 = Lcogs^2
gen LkLcogs = Lk * Lcogs
gen kLcogs = k * Lcogs
gen Lpp = L.pp_dummy
keep if !mi(k, cogs, Lk, Lcogs, L2cogs, Lpp)

* Rebuild PHI via first-stage polynomial
reg go k cogs k2 cogs2 kcogs pp_dummy c.k#pp_dummy c.cogs#pp_dummy i.year
predict phi
predict eps, residuals
gen Lphi = L.phi
gen const = 1

drop if mi(phi, Lphi)

* Mata: numerical Jacobian of moment function at the stored coefficients
clear mata
mata:
b = st_matrix("r(params)")
// Use Stata locals via st_local
bk    = strtoreal(st_local("b_k"))
bc    = strtoreal(st_local("b_cogs"))
bk2   = strtoreal(st_local("b_k2"))
bc2   = strtoreal(st_local("b_c2"))
bkc   = strtoreal(st_local("b_kc"))
b_hat = (0, bk, bc, bk2, bc2, bkc)

Z = st_data(., ("const", "k", "Lcogs", "k2", "Lcogs2", "kLcogs", "Lk", "L2cogs"))
X = st_data(., ("const", "k", "cogs", "k2", "cogs2", "kcogs"))
X_lag = st_data(., ("const", "Lk", "Lcogs", "Lk2", "Lcogs2", "LkLcogs"))
PHI = st_data(., ("phi"))
PHIL = st_data(., ("Lphi"))
Lpp = st_data(., ("Lpp"))
C = st_data(., ("const"))
cl = st_data(., ("id"))
N = rows(X)
K = cols(b_hat)
Kz = cols(Z)

real matrix compute_xi(real rowvector b)
{
    external PHI, PHIL, X, X_lag, Lpp, C, Z, N
    OMEGA = PHI - X*b'
    OL = PHIL - X_lag*b'
    P = (C, OL, Lpp)
    gb = invsym(P'P) * P'OMEGA
    return(OMEGA - P*gb)
}

// Moment vector g(β) = Z'ξ/N
real matrix compute_g(real rowvector b)
{
    external Z, N
    xi = compute_xi(b)
    return(Z'xi / N)
}

g0 = compute_g(b_hat)

// Numerical Jacobian G (Kz × K)
h = 1e-6
G = J(Kz, K, 0)
for (j=1; j<=K; j++) {
    bp = b_hat
    bm = b_hat
    bp[j] = bp[j] + h
    bm[j] = bm[j] - h
    G[.,j] = (compute_g(bp) - compute_g(bm)) / (2*h)
}

// Weight matrix W = (Z'Z/N)^{-1}
W = invsym(Z'Z / N)

// AGS Λ = -(G'WG)^{-1} G'W diag(σ)
// where σ = firm-clustered SD of moments
xi_hat = compute_xi(b_hat)
ZXI = Z :* xi_hat
cids = uniqrows(cl)
Nc = rows(cids)
S = J(Kz, Kz, 0)
for (c=1; c<=Nc; c++) {
    sel = selectindex(cl :== cids[c])
    mc = colsum(ZXI[sel,.])
    S = S + mc'*mc
}
S = S / N * (Nc/(Nc-1))
sigma_m = sqrt(diagonal(S))

Lambda = -invsym(G'*W*G) * G'*W * diag(sigma_m)
st_matrix("AGS_lambda", Lambda)
printf("  Λ matrix computed (%gx%g)\n", rows(Lambda), cols(Lambda))
end

* Dump to table
cap mkdir "$output/tables"
tempname tf
file open `tf' using "../../output/tables/ags_twostep_sensitivity.tex", write replace
#delimit ;
file write `tf'
    "\begin{table}[htbp]\centering" _n
    "\caption{AGS (2017) Scaled Sensitivity \$\hat\Lambda\$ --- NACE 41 TL}" _n
    "\label{tab:ags_twostep_sensitivity}" _n
    "\begin{threeparttable}" _n
    "\begin{tabular}{lcccccc}" _n
    "\toprule" _n
    "Moment \$\downarrow\$ / Parameter \$\rightarrow\$ & const & \$k\$ & cogs & \$k^2\$ & cogs\$^2\$ & \$k\cdot\$cogs \\" _n
    "\midrule" _n
;
#delimit cr

local mnames "const k Lcogs k2 Lcogs2 kLcogs Lk L2cogs"
forvalues i = 1/8 {
    local m : word `i' of `mnames'
    local mname = subinstr("`m'", "_", "\_", .)
    file write `tf' "\texttt{`mname'}"
    forvalues j = 1/6 {
        local v = AGS_lambda[`i',`j']
        file write `tf' " & " %10.2f (`v')
    }
    file write `tf' " \\" _n
}

#delimit ;
file write `tf'
    "\bottomrule" _n
    "\end{tabular}" _n
    "\begin{tablenotes}\footnotesize" _n
    "\item \emph{Notes:} $\hat\Lambda = -(\hat G'\hat W \hat G)^{-1}\hat G' "
    "\hat W \cdot \text{diag}(\hat\sigma)$ evaluated at the translog ACF " _n
    "GMM estimates for NACE 41. Rows are moments (instruments $\times$ " _n
    "residual); columns are translog parameters. Large $|\Lambda|$ entries "
    "identify load-bearing moments per Andrews, Gentzkow \& Shapiro (2017)." _n
    "\end{tablenotes}" _n
    "\end{threeparttable}" _n
    "\end{table}" _n
;
#delimit cr
file close `tf'
dis "  Saved: ags_twostep_sensitivity.tex"
