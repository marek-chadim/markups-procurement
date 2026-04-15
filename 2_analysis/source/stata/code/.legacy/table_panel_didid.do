*===============================================================================
* table_panel_didid.do — Panel DiD Battery (csdid + bacondecomp)
*
* Port of panel_treatment.py / panel_treatment_effects.R / sunab_event_study.R
* Uses three panel-DiD estimators robust to TWFE bias in staggered designs:
*   1. csdid        — Callaway-Sant'Anna (2021) doubly robust
*   2. bacondecomp  — Goodman-Bacon (2021) decomposition
*   3. reghdfe event-study with pre/post leads+lags (Sun-Abraham analog)
*
* Input:  $data/markups_panel.dta
* Output: ../../output/tables/panel_bacon.tex
*         ../../output/tables/panel_sdid.tex
*         ../../output/tables/panel_treatment_effects.tex
*===============================================================================

dis _newline "--- table_panel_didid.do ---"

use "$data/markups_panel.dta", clear
cap confirm var mu_A
if _rc != 0 {
    dis "  SKIP: missing mu_A"
    exit
}

gen log_mu = log(mu_A)
xtset id year

* Construct first-treatment year (for CS and SA)
bys id (year): gen first_pp = year if pp_dummy == 1 & (pp_dummy[_n-1] == 0 | _n == 1)
bys id: egen gvar = min(first_pp)
replace gvar = 0 if mi(gvar)

* ---- csdid: Callaway-Sant'Anna ----
cap which csdid
if _rc == 0 {
    cap csdid log_mu k cogs, ivar(id) time(year) gvar(gvar) ///
        notyet long2 agg(simple)
    if _rc == 0 {
        local att_cs = e(b)[1,1]
        local se_cs = sqrt(e(V)[1,1])
        matrix b_cs = e(b)
    }
}

* ---- bacondecomp: Goodman-Bacon ----
cap which bacondecomp
if _rc == 0 {
    cap bacondecomp log_mu pp_dummy, ddetail
    if _rc == 0 {
        matrix bacon_weights = e(sumbacon)
        local bacon_att = e(avgdd)
    }
}

* ---- reghdfe event study ----
* Build relative-time dummies around first_pp
bys id (year): egen first_treat = min(cond(pp_dummy == 1, year, .))
gen rel_time = year - first_treat if !mi(first_treat)
forvalues k = 4(-1)1 {
    gen L`k'_ev = (rel_time == -`k')
}
gen E0_ev = (rel_time == 0)
forvalues k = 1/4 {
    gen F`k'_ev = (rel_time == `k')
}
gen F5p_ev = (rel_time >= 5 & !mi(rel_time))
gen L5p_ev = (rel_time <= -5 & !mi(rel_time))

cap reghdfe log_mu L5p_ev L4_ev L3_ev L2_ev L1_ev E0_ev F1_ev F2_ev F3_ev F4_ev F5p_ev, ///
    absorb(id year) vce(cluster id)

* Write Bacon decomposition table
cap mkdir "$output/tables"
tempname tf
file open `tf' using "../../output/tables/panel_bacon.tex", write replace
#delimit ;
file write `tf'
    "\begin{table}[htbp]\centering" _n
    "\caption{Goodman-Bacon (2021) Decomposition of TWFE DiD}" _n
    "\label{tab:panel_bacon}" _n
    "\begin{threeparttable}" _n
    "\begin{tabular}{lc}" _n
    "\toprule" _n
    "Comparison & Weight \\" _n
    "\midrule" _n
;
#delimit cr
cap confirm matrix bacon_weights
if _rc == 0 {
    forvalues i = 1/`=rowsof(bacon_weights)' {
        local rn : word `i' of `: rownames bacon_weights'
        local v = bacon_weights[`i', 1]
        file write `tf' "`rn' & " %7.4f (`v') " \\" _n
    }
}
else {
    file write `tf' "(bacondecomp unavailable) & -- \\" _n
}
#delimit ;
file write `tf'
    "\bottomrule" _n
    "\end{tabular}" _n
    "\begin{tablenotes}\footnotesize" _n
    "\item \emph{Notes:} Goodman-Bacon (2021) decomposition of the "
    "standard TWFE DiD coefficient into its component 2×2 comparisons " _n
    "and the associated weights." _n
    "\end{tablenotes}" _n
    "\end{threeparttable}" _n
    "\end{table}" _n
;
#delimit cr
file close `tf'

* Write CS table
tempname tf2
file open `tf2' using "../../output/tables/panel_sdid.tex", write replace
#delimit ;
file write `tf2'
    "\begin{table}[htbp]\centering" _n
    "\caption{Callaway-Sant'Anna (2021) Doubly Robust ATT}" _n
    "\label{tab:panel_sdid}\label{tab:panel_treatment}" _n
    "\begin{threeparttable}" _n
    "\begin{tabular}{lcc}" _n
    "\toprule" _n
    "Estimator & $\hat{ATT}$ & SE \\" _n
    "\midrule" _n
;
#delimit cr
cap confirm scalar att_cs
if _rc == 0 {
    file write `tf2' "CS Doubly Robust (simple agg.)" " & " %7.4f (`att_cs') " & " %7.4f (`se_cs') " \\" _n
}
else {
    file write `tf2' "(csdid not run) & -- & -- \\" _n
}
#delimit ;
file write `tf2'
    "\bottomrule" _n
    "\end{tabular}" _n
    "\begin{tablenotes}\footnotesize" _n
    "\item \emph{Notes:} Callaway \& Sant'Anna (2021) doubly robust " _n
    "heterogeneous-cohort ATT. Aggregation `simple` averages all " _n
    "year-since-treatment cells." _n
    "\end{tablenotes}" _n
    "\end{threeparttable}" _n
    "\end{table}" _n
;
#delimit cr
file close `tf2'

* Write event-study table
tempname tf3
file open `tf3' using "../../output/tables/panel_treatment_effects.tex", write replace
#delimit ;
file write `tf3'
    "\begin{table}[htbp]\centering" _n
    "\caption{Event-Study Estimates around First-Procurement Year}" _n
    "\label{tab:panel_treatment_effects}" _n
    "\begin{threeparttable}" _n
    "\begin{tabular}{lcc}" _n
    "\toprule" _n
    "Relative year & $\hat\beta$ & SE \\" _n
    "\midrule" _n
;
#delimit cr
foreach v in L4_ev L3_ev L2_ev L1_ev E0_ev F1_ev F2_ev F3_ev F4_ev F5p_ev {
    cap local b = _b[`v']
    cap local s = _se[`v']
    if !mi("`b'") {
        local vlabel = subinstr("`v'", "_", "\_", .)
        file write `tf3' "`vlabel' & " %7.4f (`b') " & " %7.4f (`s') " \\" _n
    }
}
#delimit ;
file write `tf3'
    "\bottomrule" _n
    "\end{tabular}" _n
    "\end{threeparttable}" _n
    "\end{table}" _n
;
#delimit cr
file close `tf3'

dis "  Saved: panel_bacon.tex, panel_sdid.tex, panel_treatment_effects.tex"
