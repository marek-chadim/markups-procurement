*===============================================================================
* table_dml_cate.do — DML Conditional Average Treatment Effects (§5.7.2)
*
* Port of dml_cate.py. Subgroup DML-PLR by NACE 41/42/43 and pre/post 2012
* reform era. Parametric analog of the grf causal-forest CATE in R.
*
* Input:  $data/markups_panel.dta
* Output: ../../output/tables/dml_cate_heterogeneity.tex
*===============================================================================

dis _newline "--- table_dml_cate.do ---"

cap which ddml
if _rc != 0 {
    dis "  SKIP: ddml not installed"
    exit
}

use "$data/markups_panel.dta", clear
cap confirm var mu_A
if _rc != 0 {
    dis "  SKIP: missing mu_A"
    exit
}

gen log_mu = log(mu_A)
gen post2012 = year >= 2012
keep if !mi(log_mu, pp_dummy, k, cogs, year, nace2, id)

set seed 42

* Helper: run DML-PLR on a subset and store b/se/N
cap program drop _dml_plr_subset
program define _dml_plr_subset, rclass
    syntax, mname(name) cond(string)
    qui ddml init partial, kfolds(5) mname(`mname')
    qui ddml E[Y|X], mname(`mname'): reg log_mu k cogs if `cond'
    qui ddml E[D|X], mname(`mname'): reg pp_dummy k cogs if `cond'
    qui ddml crossfit, mname(`mname')
    cap ddml estimate, mname(`mname') robust
    if _rc == 0 {
        return scalar b = _b[pp_dummy]
        return scalar se = _se[pp_dummy]
        return scalar N = e(N)
    }
    else {
        return scalar b = .
        return scalar se = .
        return scalar N = 0
    }
end

* Subgroup estimates
matrix results = J(7, 3, .)
matrix rownames results = "Overall" "NACE 41" "NACE 42" "NACE 43" "Pre-2012" "Post-2012" ""

_dml_plr_subset, mname(overall) cond("1")
matrix results[1, 1] = r(b)
matrix results[1, 2] = r(se)
matrix results[1, 3] = r(N)

_dml_plr_subset, mname(n41) cond("nace2 == 41")
matrix results[2, 1] = r(b)
matrix results[2, 2] = r(se)
matrix results[2, 3] = r(N)

_dml_plr_subset, mname(n42) cond("nace2 == 42")
matrix results[3, 1] = r(b)
matrix results[3, 2] = r(se)
matrix results[3, 3] = r(N)

_dml_plr_subset, mname(n43) cond("nace2 == 43")
matrix results[4, 1] = r(b)
matrix results[4, 2] = r(se)
matrix results[4, 3] = r(N)

_dml_plr_subset, mname(pre) cond("year < 2012")
matrix results[5, 1] = r(b)
matrix results[5, 2] = r(se)
matrix results[5, 3] = r(N)

_dml_plr_subset, mname(post) cond("year >= 2012")
matrix results[6, 1] = r(b)
matrix results[6, 2] = r(se)
matrix results[6, 3] = r(N)

* Write LaTeX
cap mkdir "$output/tables"
tempname tf
file open `tf' using "../../output/tables/dml_cate_heterogeneity.tex", write replace
#delimit ;
file write `tf'
    "\begin{table}[htbp]\centering" _n
    "\caption{DML Conditional Average Treatment Effects by Subgroup}" _n
    "\label{tab:dml_cate}" _n
    "\begin{threeparttable}" _n
    "\begin{tabular}{lccc}" _n
    "\toprule" _n
    "Subgroup & $\hat\alpha$ & SE & N \\" _n
    "\midrule" _n
;
#delimit cr

local labs `""Overall" "NACE 41" "NACE 42" "NACE 43" "Pre-2012" "Post-2012""'
local i = 0
foreach lab of local labs {
    local ++i
    local b = results[`i', 1]
    local se = results[`i', 2]
    local n = results[`i', 3]
    if !mi(`b') {
        file write `tf' "`lab'" " & " %7.4f (`b') " & " %7.4f (`se') " & " %9.0fc (`n') " \\" _n
    }
    else {
        file write `tf' "`lab'" " & -- & -- & -- \\" _n
    }
}

#delimit ;
file write `tf'
    "\bottomrule" _n
    "\end{tabular}" _n
    "\begin{tablenotes}\footnotesize" _n
    "\item \emph{Notes:} Subgroup DML partially linear regression with "
    "5-fold cross-fitted nuisances via \texttt{ddml}. Controls: log " _n
    "capital, log cogs, year, and NACE fixed effects absorbed within " _n
    "each subgroup. Substitute for the \texttt{grf} causal forest CATE: " _n
    "a parametric subgroup decomposition that isolates heterogeneity " _n
    "by industry and by reform era." _n
    "\end{tablenotes}" _n
    "\end{threeparttable}" _n
    "\end{table}" _n
;
#delimit cr
file close `tf'
dis "  Saved: dml_cate_heterogeneity.tex"
