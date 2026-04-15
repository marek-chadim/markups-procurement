clear all
use "../../input/data_rebuilt.dta", clear
dis "data_rebuilt.dta: " _N " obs"
count if !mi(k) & !mi(cogs) & !mi(go)
dis "  non-mi go,k,cogs: " r(N)

xtset id year
gen Lcogs = L.cogs
gen Lk = L.k
gen L2cogs = L2.cogs
gen Lpp = L.pp_dummy
count if !mi(k) & !mi(cogs) & !mi(Lk) & !mi(Lcogs) & !mi(Lpp) & !mi(go)
dis "  + first lags + Lpp: " r(N)
count if !mi(k) & !mi(cogs) & !mi(Lk) & !mi(Lcogs) & !mi(L2cogs) & !mi(Lpp) & !mi(go)
dis "  + L2cogs: " r(N)

tab nace2 if !mi(k) & !mi(cogs) & !mi(Lk) & !mi(Lcogs) & !mi(Lpp) & !mi(go)
