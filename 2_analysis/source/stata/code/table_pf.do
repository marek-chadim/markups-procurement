*===============================================================================
* table_pf.do — Table 1: Production Function Estimates
*
* No direct DGM analog (DGM report coefficients inline).
* We create a standalone coefficient table with robustness panel.
*
* Output: output/table_pf.tex
*===============================================================================

dis _newline "--- table_pf.do ---"

use "$data/coefficients_byind.dta", clear

* LaTeX output
cap file close tf
file open tf using "$output/table_pf.tex", write replace

#delimit ;
file write tf
    "\begin{table}[htbp]\centering" _n
    "\caption{Production Function Estimates by Industry}" _n
    "\label{tab:pf_estimates}" _n
    "\begin{tabular}{l*{3}{c}}" _n
    "\hline\hline" _n
    " & NACE 41 & NACE 42 & NACE 43 \\" _n
    " & (Buildings) & (Civil Eng.) & (Specialized) \\" _n
    "\hline" _n
    "\multicolumn{4}{l}{\textit{Panel A: Baseline (survival + pp in Markov)}} \\" _n
;
#delimit cr

* Panel A: beta_k
file write tf _char(36) "\hat{\beta}_k" _char(36)
forvalues i = 1/3 {
    local b = b_k_A[`i']
    file write tf " & " %6.4f (`b')
}
file write tf " \\" _n " "
forvalues i = 1/3 {
    local s = se_k_A[`i']
    file write tf " & (" %6.4f (`s') ")"
}
file write tf " \\" _n

* Panel A: beta_cogs
file write tf _char(36) "\hat{\beta}_{\text{cogs}}" _char(36)
forvalues i = 1/3 {
    local b = b_cogs_A[`i']
    file write tf " & " %6.4f (`b')
}
file write tf " \\" _n " "
forvalues i = 1/3 {
    local s = se_cogs_A[`i']
    file write tf " & (" %6.4f (`s') ")"
}
file write tf " \\" _n

* RTS and N
file write tf "RTS"
forvalues i = 1/3 {
    local rts = b_k_A[`i'] + b_cogs_A[`i']
    file write tf " & " %5.3f (`rts')
}
file write tf " \\" _n
file write tf _char(36) "N" _char(36)
forvalues i = 1/3 {
    local nn = N_obs[`i']
    file write tf " & " %6.0fc (`nn')
}
file write tf " \\" _n

* Panel B
#delimit ;
file write tf
    "\hline" _n
    "\multicolumn{4}{l}{\textit{Panel B: Robustness --- "
    _char(36) "\hat{\beta}_{\text{cogs}}" _char(36) " (SE)}} \\" _n
;
#delimit cr

foreach spec in B C D OLS {
    if "`spec'"=="B" local lab "No survival"
    if "`spec'"=="C" local lab "No pp in Markov"
    if "`spec'"=="D" local lab "Plain ACF"
    if "`spec'"=="OLS" local lab "OLS"

    file write tf "`lab'"
    forvalues i = 1/3 {
        local b = b_cogs_`spec'[`i']
        file write tf " & " %6.4f (`b')
    }
    file write tf " \\" _n
    if "`spec'" != "OLS" {
        file write tf " "
        forvalues i = 1/3 {
            local s = se_cogs_`spec'[`i']
            file write tf " & (" %6.4f (`s') ")"
        }
        file write tf " \\" _n
    }
}

* Translog
file write tf "Translog"
forvalues i = 1/3 {
    local b = b_cogs_E[`i']
    file write tf " & " %6.4f (`b')
}
file write tf " \\" _n " "
forvalues i = 1/3 {
    local s = se_cogs_E[`i']
    file write tf " & (" %6.4f (`s') ")"
}
file write tf " \\" _n

#delimit ;
file write tf
    "\hline\hline" _n
    "\multicolumn{4}{p{10cm}}{\footnotesize" _n
    "Notes: Two-step ACF estimator (Ackerberg et al.\ 2015). "
    "Baseline Markov: " _char(36) "\omega_t = \rho\omega_{t-1} + \gamma pp_{t-1}"
    " + \delta\hat{p}_{t-1} + \xi_t" _char(36) ". "
    "Instruments: " _char(36) "(1, k_t, \text{cogs}_{t-1})" _char(36) ". "
    "Analytical SEs (ACH 2012), clustered by firm.} \\" _n
    "\end{tabular}" _n
    "\end{table}" _n
;
#delimit cr

file close tf
dis "  Saved: table_pf.tex"
