*===============================================================================
* table_moments.do — Table 3: Markup Distributional Moments
*
* Analog of DGM Table 4:
*   DGM: mean/sd/p25/p50/p75 by industry for FSQ/BFSR/NoFSQ/NoFSR
*   Us:  mean/sd/p10/p50/p90 by nace2 for Base/NoSurv/NoPP/Plain/TL/OLS
*
* Output: output/table_moments.tex
*===============================================================================

dis _newline "--- table_moments.do ---"


use "$data/markups_panel.dta", clear

cap file close tf
file open tf using "$output/table_moments.tex", write replace

#delimit ;
file write tf
    "\begin{table}[htbp]\centering" _n
    "\caption{Markup Distributional Moments by Industry and Specification}" _n
    "\label{tab:moments}" _n
    "\small" _n
    "\begin{tabular}{ll*{5}{c}}" _n
    "\hline\hline" _n
    " & Spec & Mean & SD & p10 & p50 & p90 \\" _n
    "\hline" _n
;
#delimit cr

levelsof nace2, local(naces)
foreach n of local naces {
    if `n'==41 local nlab "NACE 41 (Buildings)"
    if `n'==42 local nlab "NACE 42 (Civil Eng.)"
    if `n'==43 local nlab "NACE 43 (Specialized)"

    file write tf "\multicolumn{7}{l}{\textit{`nlab'}} \\" _n

    foreach s in A D E OLS {
        if "`s'"=="A" local slab "Base"
        if "`s'"=="D" local slab "Plain"
        if "`s'"=="E" local slab "Translog"
        if "`s'"=="OLS" local slab "OLS"

        qui sum mu_`s' if nace2==`n', d
        local mn = r(mean)
        local sd = r(sd)
        local p10 = r(p10)
        local p50 = r(p50)
        local p90 = r(p90)

        file write tf " & `slab'"
        file write tf " & " %5.2f (`mn') " & " %5.2f (`sd')
        file write tf " & " %5.2f (`p10') " & " %5.2f (`p50') " & " %5.2f (`p90')
        file write tf " \\" _n
    }
    file write tf "\hline" _n
}

* Pooled
file write tf "\multicolumn{7}{l}{\textit{All industries}} \\" _n
foreach s in A D E OLS {
    if "`s'"=="A" local slab "Base"
    if "`s'"=="D" local slab "Plain"
    if "`s'"=="E" local slab "Translog"
    if "`s'"=="OLS" local slab "OLS"

    qui sum mu_`s', d
    file write tf " & `slab'"
    file write tf " & " %5.2f (r(mean)) " & " %5.2f (r(sd))
    file write tf " & " %5.2f (r(p10)) " & " %5.2f (r(p50)) " & " %5.2f (r(p90))
    file write tf " \\" _n
}

#delimit ;
file write tf
    "\hline\hline" _n
    "\multicolumn{7}{p{12cm}}{\footnotesize" _n
    "Notes: Markup " _char(36) "\mu = \hat{\theta}^V / \hat{\alpha}^V" _char(36)
    " where " _char(36) "\hat{\theta}^V" _char(36) " is the estimated output elasticity"
    " and " _char(36) "\hat{\alpha}^V" _char(36) " the corrected expenditure share."
    " Base = CD with survival + pp in Markov;"
    " Plain = CD without Markov controls;"
    " Translog = TL with survival + pp.} \\" _n
    "\end{tabular}" _n
    "\end{table}" _n
;
#delimit cr

file close tf
dis "  Saved: table_moments.tex"
