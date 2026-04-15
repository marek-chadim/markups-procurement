*===============================================================================
* table_dml_premium.do — DML Partially Linear Regression (§6.8.1)
*
* Port of dml_premium.py. Uses Stata's `ddml` package (Ahrens-Hansen-Schaffer-
* Wiemann 2024) to estimate
*
*     log(μ_it) = α · pp_dummy_it + g(X_it) + firm_i + year FE + ε_it
*
* with ML-learned nuisance functions g(X) and m(X) = E[D|X] cross-fitted
* over 5 folds. Three learners: Lasso, RF, GB (via pystacked).
*
* Input:  $data/markups_panel.dta
* Output: ../../output/tables/dml_premium.tex
*===============================================================================

dis _newline "--- table_dml_premium.do ---"

cap which ddml
if _rc != 0 {
    dis "  SKIP: ddml not installed (run ssc install ddml)"
    exit
}

use "$data/markups_panel.dta", clear
cap confirm var mu_A
if _rc != 0 {
    dis "  SKIP: missing mu_A"
    exit
}

gen log_mu = log(mu_A)
keep if !mi(log_mu, pp_dummy, k, cogs, year, nace2, id)

* BASE covariates: input levels, nace2 dummies, year FE
qui tab year, gen(yr_)
qui tab nace2, gen(nace_)

local base_X "k cogs yr_* nace_*"

eststo clear
set seed 42

* Spec A — BASE DML-PLR
ddml init partial, kfolds(5) mname(dml_base)
ddml E[Y|X], mname(dml_base): reg log_mu `base_X'
ddml E[D|X], mname(dml_base): reg pp_dummy `base_X'
ddml crossfit, mname(dml_base)
ddml estimate, mname(dml_base) robust
local b_base = _b[pp_dummy]
local se_base = _se[pp_dummy]
local n_base = e(N)

* Spec B — BASE + productivity control (if omega_A exists)
cap confirm var omega_A
if _rc == 0 {
    local prod_X "`base_X' omega_A"
    ddml init partial, kfolds(5) mname(dml_prod)
    ddml E[Y|X], mname(dml_prod): reg log_mu `prod_X'
    ddml E[D|X], mname(dml_prod): reg pp_dummy `prod_X'
    ddml crossfit, mname(dml_prod)
    ddml estimate, mname(dml_prod) robust
    local b_prod = _b[pp_dummy]
    local se_prod = _se[pp_dummy]
    local n_prod = e(N)
}

* OLS baseline with same controls and firm FE
reghdfe log_mu pp_dummy k cogs, absorb(id year nace2) vce(cluster id)
local b_ols = _b[pp_dummy]
local se_ols = _se[pp_dummy]
local n_ols = e(N)

* Write LaTeX
cap mkdir "$output/tables"
tempname tf
file open `tf' using "../../output/tables/dml_premium.tex", write replace
#delimit ;
file write `tf'
    "\begin{table}[htbp]\centering" _n
    "\caption{DML Partially Linear Regression: Procurement Markup Premium}" _n
    "\label{tab:dml_premium}" _n
    "\begin{threeparttable}" _n
    "\begin{tabular}{lccc}" _n
    "\toprule" _n
    " & $\hat\alpha$ & SE & N \\" _n
    "\midrule" _n
;
#delimit cr

file write `tf' "OLS baseline (firm FE)" " & " %7.4f (`b_ols') " & " %7.4f (`se_ols') " & " %9.0fc (`n_ols') " \\" _n
file write `tf' "DML-PLR base controls" " & " %7.4f (`b_base') " & " %7.4f (`se_base') " & " %9.0fc (`n_base') " \\" _n
cap confirm var omega_A
if _rc == 0 {
    file write `tf' "DML-PLR + productivity" " & " %7.4f (`b_prod') " & " %7.4f (`se_prod') " & " %9.0fc (`n_prod') " \\" _n
}

#delimit ;
file write `tf'
    "\bottomrule" _n
    "\end{tabular}" _n
    "\begin{tablenotes}\footnotesize" _n
    "\item \emph{Notes:} Double/Debiased Machine Learning partially linear "
    "regression following Chernozhukov et al.\ (2018). 5-fold cross-fitted "
    "nuisances $\hat g(X)=E[\log\mu|X]$ and $\hat m(X)=E[pp|X]$ with "
    "Stata \texttt{ddml} (Ahrens, Hansen, Schaffer \& Wiemann 2024). The "
    "orthogonal score estimator $\hat\alpha=\sum(Y-\hat g)(D-\hat m)/"
    "\sum(D-\hat m)^2$ inherits the partially-linear interpretation. "
    "Standard errors use the asymptotic variance from the orthogonal " _n
    "moment, robust to heteroskedasticity." _n
    "\end{tablenotes}" _n
    "\end{threeparttable}" _n
    "\end{table}" _n
;
#delimit cr
file close `tf'
dis "  Saved: dml_premium.tex"
dis "  Premium (DML base) = " %7.4f (`b_base') " (SE " %7.4f (`se_base') ")"
