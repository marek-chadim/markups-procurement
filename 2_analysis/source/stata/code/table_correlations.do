*===============================================================================
* table_correlations.do — Table 4: Cross-Specification Correlations
*
* Analog of DGM Table 5:
*   DGM: Pearson/Spearman correlations between quantity & revenue markups
*   Us:  correlations between Base/Plain/TL/OLS markups (levels & FD)
*
* Output: output/table_correlations.tex
*===============================================================================

dis _newline "--- table_correlations.do ---"


use "$data/markups_panel.dta", clear

* Balance panel on non-missing across specs
drop if mi(l_mu_A) | mi(l_mu_D) | mi(l_mu_E) | mi(l_mu_OLS)

cap file close tf
file open tf using "$output/table_correlations.tex", write replace

#delimit ;
file write tf
    "\begin{table}[htbp]\centering" _n
    "\caption{Cross-Specification Markup Correlations}" _n
    "\label{tab:correlations}" _n
    "\begin{tabular}{l*{4}{c}}" _n
    "\hline\hline" _n
    " & Base--Plain & Base--TL & Base--OLS & Plain--OLS \\" _n
    "\hline" _n
    "\multicolumn{5}{l}{\textit{Levels (log markups)}} \\" _n
;
#delimit cr

* Pearson correlations (levels)
local pairs "A-D A-E A-OLS D-OLS"
local labs  `""Base--Plain" "Base--TL" "Base--OLS" "Plain--OLS""'

* Pooled Pearson
file write tf "Pearson (pooled)"
foreach pair in A_D A_E A_OLS D_OLS {
    local s1 = substr("`pair'", 1, strpos("`pair'","_")-1)
    local s2 = substr("`pair'", strpos("`pair'","_")+1, .)
    qui corr l_mu_`s1' l_mu_`s2'
    file write tf " & " %5.3f (r(rho))
}
file write tf " \\" _n

* Spearman
file write tf "Spearman (pooled)"
foreach pair in A_D A_E A_OLS D_OLS {
    local s1 = substr("`pair'", 1, strpos("`pair'","_")-1)
    local s2 = substr("`pair'", strpos("`pair'","_")+1, .)
    qui spearman l_mu_`s1' l_mu_`s2'
    file write tf " & " %5.3f (r(rho))
}
file write tf " \\" _n

* By industry
levelsof nace2, local(naces)
foreach n of local naces {
    file write tf "Pearson (NACE `n')"
    foreach pair in A_D A_E A_OLS D_OLS {
        local s1 = substr("`pair'", 1, strpos("`pair'","_")-1)
        local s2 = substr("`pair'", strpos("`pair'","_")+1, .)
        qui corr l_mu_`s1' l_mu_`s2' if nace2==`n'
        file write tf " & " %5.3f (r(rho))
    }
    file write tf " \\" _n
}

* First differences
#delimit ;
file write tf
    "\hline" _n
    "\multicolumn{5}{l}{\textit{First differences (" _char(36) "\Delta" _char(36) " log markups)}} \\" _n
;
#delimit cr

* Pooled FD Pearson
file write tf "Pearson (pooled)"
foreach pair in A_D A_E A_OLS D_OLS {
    local s1 = substr("`pair'", 1, strpos("`pair'","_")-1)
    local s2 = substr("`pair'", strpos("`pair'","_")+1, .)
    qui corr FD_l_mu_`s1' FD_l_mu_`s2'
    file write tf " & " %5.3f (r(rho))
}
file write tf " \\" _n

* FD Spearman
file write tf "Spearman (pooled)"
foreach pair in A_D A_E A_OLS D_OLS {
    local s1 = substr("`pair'", 1, strpos("`pair'","_")-1)
    local s2 = substr("`pair'", strpos("`pair'","_")+1, .)
    qui spearman FD_l_mu_`s1' FD_l_mu_`s2'
    file write tf " & " %5.3f (r(rho))
}
file write tf " \\" _n

local N = _N
#delimit ;
file write tf
    "\hline" _n
    _char(36) "N" _char(36) " & \multicolumn{4}{c}{" %6.0fc (`N') "} \\" _n
    "\hline\hline" _n
    "\multicolumn{5}{p{11cm}}{\footnotesize" _n
    "Notes: Correlations between log markups across specifications."
    " Base = CD with survival + pp in Markov;"
    " Plain = CD without Markov controls;"
    " TL = translog; OLS = first-stage only."
    " First differences: " _char(36) "\Delta\ln\mu_{it} = \ln\mu_{it} - \ln\mu_{it-1}" _char(36) ".} \\" _n
    "\end{tabular}" _n
    "\end{table}" _n
;
#delimit cr

file close tf
dis "  Saved: table_correlations.tex"
