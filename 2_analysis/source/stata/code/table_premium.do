*===============================================================================
* table_premium.do — Table 2: Procurement Premium
*
* No DGM analog. Novel to our application.
* Reports raw and regression premiums across all specifications.
*
* Output: output/table_premium.tex
*===============================================================================

dis _newline "--- table_premium.do ---"

use "$data/markups_panel.dta", clear

cap file close tf
file open tf using "$output/table_premium.tex", write replace

#delimit ;
file write tf
    "\begin{table}[htbp]\centering" _n
    "\caption{Procurement Markup Premium Across Specifications}" _n
    "\label{tab:premium}" _n
    "\begin{tabular}{l*{2}{c}}" _n
    "\hline\hline" _n
    " & " _char(36) "pp_t" _char(36) " & " _char(36) "pp^{3y}_t" _char(36) " \\" _n
    "\hline" _n
    "\multicolumn{3}{l}{\textit{Raw log markup difference}} \\" _n
;
#delimit cr

foreach spec in A B C D E OLS {
    if "`spec'"=="A" local lab "Base (surv.\ + pp)"
    if "`spec'"=="B" local lab "No survival"
    if "`spec'"=="C" local lab "No pp in Markov"
    if "`spec'"=="D" local lab "Plain ACF"
    if "`spec'"=="E" local lab "Translog"
    if "`spec'"=="OLS" local lab "OLS"

    * pp_dummy raw
    qui ttest l_mu_`spec', by(pp_dummy)
    local raw = r(mu_2) - r(mu_1)
    local se_raw = r(se)
    * pp_ever_3y raw
    qui ttest l_mu_`spec', by(pp_ever_3y)
    local raw3 = r(mu_2) - r(mu_1)
    local se3 = r(se)

    file write tf "`lab'"
    file write tf " & " %6.3f (`raw') " & " %6.3f (`raw3') " \\" _n
    file write tf " & (" %5.3f (`se_raw') ") & (" %5.3f (`se3') ") \\" _n
}

#delimit ;
file write tf
    "\hline" _n
    "\multicolumn{3}{l}{\textit{Regression (" _char(36) "k" _char(36) ", cogs, year" _char(36) "\times" _char(36) "nace2 FE)}} \\" _n
;
#delimit cr

foreach spec in A B C D E OLS {
    if "`spec'"=="A" local lab "Base (surv.\ + pp)"
    if "`spec'"=="B" local lab "No survival"
    if "`spec'"=="C" local lab "No pp in Markov"
    if "`spec'"=="D" local lab "Plain ACF"
    if "`spec'"=="E" local lab "Translog"
    if "`spec'"=="OLS" local lab "OLS"

    xi: qui reg l_mu_`spec' pp_dummy k cogs i.year*i.nace2, cluster(id)
    local reg = _b[pp_dummy]
    local se = _se[pp_dummy]
    local N = e(N)

    xi: qui reg l_mu_`spec' pp_ever_3y k cogs i.year*i.nace2, cluster(id)
    local reg3 = _b[pp_ever_3y]
    local se3 = _se[pp_ever_3y]

    file write tf "`lab'"
    file write tf " & " %6.3f (`reg') " & " %6.3f (`reg3') " \\" _n
    file write tf " & (" %5.3f (`se') ") & (" %5.3f (`se3') ") \\" _n
}

#delimit ;
file write tf
    "\hline" _n
    _char(36) "N" _char(36) " & \multicolumn{2}{c}{" %6.0fc (`N') "} \\" _n
    "\hline\hline" _n
    "\multicolumn{3}{p{10cm}}{\footnotesize" _n
    "Notes: Premium is the log markup difference between procurement"
    " and non-procurement firms. Regression controls for " _char(36) "k" _char(36) ", cogs,"
    " year" _char(36) "\times" _char(36) "nace2 FE; SEs clustered by firm."
    " " _char(36) "pp_t" _char(36) " = active supplier in year " _char(36) "t" _char(36) ";"
    " " _char(36) "pp^{3y}_t" _char(36) " = active in " _char(36) "t" _char(36) ", " _char(36) "t{-}1" _char(36) ", or " _char(36) "t{-}2" _char(36) ".} \\" _n
    "\end{tabular}" _n
    "\end{table}" _n
;
#delimit cr

file close tf
dis "  Saved: table_premium.tex"
