cap confirm file "../../output/stata/magnusweb_profit.dta"
if _rc != 0 {
    dis "NO FILE"
    exit
}
use "../../output/stata/magnusweb_profit.dta", clear
dis "Rows: " _N
ds
summ cm_iii wva ws lp, d
count if !mi(cm_iii)
dis "Non-missing cm_iii: " r(N)
