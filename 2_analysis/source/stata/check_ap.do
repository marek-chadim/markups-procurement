clear all
use "../../output/stata/analysis_panel.dta", clear
dis "analysis_panel.dta: " _N " obs"
count if !mi(k) & !mi(cogs) & !mi(go)
dis "  non-mi go,k,cogs: " r(N)
xtset id year
gen Lcogs_ = L.cogs
gen Lk_ = L.k
gen Lpp_ = L.pp_dummy
count if !mi(k) & !mi(cogs) & !mi(go) & !mi(Lk_) & !mi(Lcogs_) & !mi(Lpp_)
dis "  + first lags: " r(N)
