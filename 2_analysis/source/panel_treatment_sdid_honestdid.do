/*
 * panel_treatment_sdid_honestdid.do
 *
 * Runs Synthetic Difference-in-Differences (Clarke, Pailanir, Athey & Imbens
 * 2023 Stata Journal, implementing Arkhangelsky et al. 2021 AER) and
 * HonestDiD sensitivity analysis (Caceres & Rambachan 2023, wrapping
 * Rambachan & Roth 2023 RES) on the Czech construction panel.
 *
 * Packages:
 *   sdid      (github.com/Daniel-Pailanir/sdid)
 *   honestdid (github.com/mcaceresb/stata-honestdid)
 *   reghdfe, ftools
 *
 * Absorbing-treatment reformulation: first_treat = earliest year with
 * pp_dummy = 1; treat_post = 1 if year >= first_treat; never-treated firms
 * have first_treat = 0 and treat_post = 0 throughout.
 *
 * Outputs:
 *   ../output/panel_sdid_honestdid_results.csv
 *   ../output/tables/panel_sdid.tex
 *   ../output/tables/panel_honestdid.tex
 *   ../output/panel_treatment_sdid_honestdid.log
 */

clear all
set more off
cap log close
log using "../output/panel_treatment_sdid_honestdid.log", replace

local input_dir  "../input"
local output_dir "../output"
local tables_dir "`output_dir'/tables"
cap mkdir "`tables_dir'"

* ======================================================================
*  Load & prepare panel
* ======================================================================
use "`input_dir'/data.dta", clear
keep id year pp_dummy
merge 1:1 id year using "`output_dir'/data/paper_markups.dta", ///
    keepusing(markup_A) keep(match) nogen

drop if missing(markup_A) | markup_A <= 0
gen log_mu = log(markup_A)

* Consecutive unit ID for sdid (which prefers small integer IDs)
egen id_num = group(id)

xtset id_num year

* ---------- Absorbing treatment reformulation ----------
bysort id_num (year): egen first_treat = min(cond(pp_dummy == 1, year, .))
replace first_treat = 0 if missing(first_treat)
gen ever_treat = (first_treat > 0)
gen treat_post = ever_treat & year >= first_treat

* Event time relative to first_treat (only defined for ever-treated firms)
gen event_time = year - first_treat if ever_treat == 1

* Sample diagnostics
summarize log_mu pp_dummy treat_post ever_treat
bysort id_num: gen first_row = (_n == 1)
count if first_row & ever_treat == 1
local n_ever = r(N)
count if first_row & ever_treat == 0
local n_never = r(N)
drop first_row
display "ever-treated firms = `n_ever'"
display "never-treated firms = `n_never'"

* ======================================================================
*  SDID on a forced-balance sub-panel
*  sdid requires each unit to appear in every period. We enforce balance
*  by keeping firms observed in every year of the panel window.
* ======================================================================
preserve

bysort id_num: egen n_years = count(year)
summarize n_years, meanonly
local max_years = r(max)
display "Max observed years per firm = `max_years'"

keep if n_years == `max_years'

* Re-index ids for balance
drop id_num
egen id_num = group(id)
xtset id_num year

summarize id_num
local n_balanced = r(max)
count
local n_obs_balanced = r(N)
display "Balanced sub-panel: `n_balanced' firms, `n_obs_balanced' obs"

* sdid requires that no firm be treated in the very first period of the
* panel. Drop firms whose first_treat equals the minimum year (they are
* effectively "always treated" from the sdid perspective).
summarize year, meanonly
local y_min = r(min)
count if first_treat == `y_min'
if (r(N) > 0) {
    display "Dropping firms with first_treat == `y_min' (always treated)"
    drop if first_treat == `y_min'
    drop id_num
    egen id_num = group(id)
    xtset id_num year
}

* Run synthetic DiD with bootstrap inference. The sdid package stores
* results in scalars e(ATT), e(se), e(ATT_l), e(ATT_r), e(N), e(N_clust).
sdid log_mu id_num year treat_post, ///
    vce(bootstrap) reps(50) seed(42)
local sdid_att = e(ATT)
local sdid_se  = e(se)
local sdid_lb  = e(ATT_l)
local sdid_ub  = e(ATT_r)
local sdid_tstat = `sdid_att' / `sdid_se'
local sdid_pval  = 2 * (1 - normal(abs(`sdid_tstat')))
local sdid_n = e(N)
local sdid_firms = e(N_clust)

count if first_treat > 0
local sdid_ever = r(N)
bysort id_num: gen first_row_b = (_n == 1)
count if first_row_b & first_treat > 0
local sdid_ever = r(N)

display "SDID ATT = `sdid_att' SE = `sdid_se' 95% CI [`sdid_lb', `sdid_ub']"
display "  N_obs = `sdid_n'  N_firms = `sdid_firms'  N_ever_treated = `sdid_ever'"

restore

* ======================================================================
*  Goodman-Bacon (2021) decomposition of the absorbing-treatment TWFE
*  Decomposes the weighted average of 2x2 DiD comparisons that compose
*  the naive TWFE estimator under staggered adoption. Implemented via
*  the bacondecomp Stata package (Goodrich, Goldring, Lovenheim,
*  Schmidt-Catran 2019, Stata Journal).
*
*  bacondecomp requires a balanced panel; we reuse the 45-firm balanced
*  sub-panel (same as SDID *before* dropping always-treated firms, which
*  the Bacon decomposition handles natively as a separate category).
*
*  The package prints a summary table with four comparison classes:
*    - Always-treated vs timing groups
*    - Never-treated vs timing groups
*    - Earlier-treated vs later-treated ("forbidden" under heterogeneous
*      effects)
*    - Later-treated vs earlier-treated ("forbidden" under heterogeneous
*      effects)
*  The display output is the canonical source (no return matrix with the
*  4-category aggregation); a Python post-processor parses the log file.
* ======================================================================
preserve
bysort id_num: egen n_years_bacon = count(year)
qui summarize n_years_bacon, meanonly
keep if n_years_bacon == r(max)
drop id_num
egen id_num = group(id)
xtset id_num year

count
local bacon_n = r(N)
bysort id_num: gen _first = (_n == 1)
count if _first
local bacon_firms = r(N)
drop _first

display "BEGIN BACONDECOMP TABLE"
bacondecomp log_mu treat_post, ddetail nograph
display "END BACONDECOMP TABLE"

* Store total TWFE from bacondecomp
local bacon_twfe = _b[treat_post]
local bacon_se   = _se[treat_post]
display "Bacon TWFE = `bacon_twfe' SE = `bacon_se' N = `bacon_n' firms = `bacon_firms'"
restore

* ======================================================================
*  Event study for HonestDiD
*  Use the full sample (balanced or unbalanced) and fit a standard TWFE
*  event-study with relative-year dummies. Reference = event_time == -1.
* ======================================================================

* Construct event-time dummies. Window [-5, +5] with pre1 (t=-1) omitted.
gen byte es_m5 = ever_treat == 1 & event_time == -5
gen byte es_m4 = ever_treat == 1 & event_time == -4
gen byte es_m3 = ever_treat == 1 & event_time == -3
gen byte es_m2 = ever_treat == 1 & event_time == -2
gen byte es_p0 = ever_treat == 1 & event_time ==  0
gen byte es_p1 = ever_treat == 1 & event_time ==  1
gen byte es_p2 = ever_treat == 1 & event_time ==  2
gen byte es_p3 = ever_treat == 1 & event_time ==  3
gen byte es_p4 = ever_treat == 1 & event_time ==  4
gen byte es_p5 = ever_treat == 1 & event_time >=  5

reghdfe log_mu es_m5 es_m4 es_m3 es_m2 es_p0 es_p1 es_p2 es_p3 es_p4 es_p5, ///
    absorb(id_num year) cluster(id_num)

matrix b_es = e(b)
matrix V_es = e(V)
matlist b_es
display "Event-study pre-trend F-test:"
testparm es_m5 es_m4 es_m3 es_m2

* ======================================================================
*  HonestDiD: sensitivity to parallel-trends violations
*  The reghdfe coefficient vector has pre-treatment at positions 1-4
*  (es_m5, es_m4, es_m3, es_m2) and post-treatment at positions 5-10
*  (es_p0, es_p1, es_p2, es_p3, es_p4, es_p5). Reference is event_time=-1.
*
*  Note: the honestdid Stata package does NOT store the robust CI matrix
*  in r(). The printed table is the canonical source. We print it here;
*  a shell post-processor parses the log to produce the LaTeX table.
* ======================================================================
display "BEGIN HONESTDID TABLE"
honestdid, pre(1/4) post(5/10) ///
    mvec(0 0.5 1 1.5 2) delta(rm)
display "END HONESTDID TABLE"
display "Event-study pre-trend F-test result: F(4, e(df_r)) Prob>F already reported above."

* ======================================================================
*  Save outputs
* ======================================================================

* SDID CSV
file open sdid_csv using "`output_dir'/panel_sdid_honestdid_results.csv", ///
    write replace
file write sdid_csv "estimator,att,se,ci_lb,ci_ub,n_obs,n_firms,n_ever_treated" _n
file write sdid_csv "SDID (Arkhangelsky et al. 2021),"
file write sdid_csv %6.4f (`sdid_att') ","
file write sdid_csv %6.4f (`sdid_se') ","
file write sdid_csv %6.4f (`sdid_lb') ","
file write sdid_csv %6.4f (`sdid_ub') ","
file write sdid_csv "`sdid_n',`sdid_firms',`sdid_ever'" _n
file close sdid_csv

* SDID LaTeX table
file open sdid_tex using "`tables_dir'/panel_sdid.tex", write replace
file write sdid_tex "\begin{table}[htbp]\centering" _n
file write sdid_tex "\caption{Synthetic Difference-in-Differences Estimate}\label{tab:panel_sdid}" _n
file write sdid_tex "\begin{threeparttable}" _n
file write sdid_tex "\begin{tabular}{lccc}" _n
file write sdid_tex "\toprule" _n
file write sdid_tex "Estimator & ATT & SE & 95\% CI \\" _n
file write sdid_tex "\midrule" _n
file write sdid_tex "Synthetic DiD (Arkhangelsky et al.\ 2021) & "
file write sdid_tex %6.4f (`sdid_att') " & (" %6.4f (`sdid_se') ") & ["
file write sdid_tex %6.4f (`sdid_lb') ", " %6.4f (`sdid_ub') "] \\" _n
file write sdid_tex "\midrule" _n
file write sdid_tex "\$N\$ firm-years & \multicolumn{3}{c}{`sdid_n'} \\" _n
file write sdid_tex "Balanced-panel firms & \multicolumn{3}{c}{`sdid_firms'} \\" _n
file write sdid_tex "Ever-treated firms in sub-panel & \multicolumn{3}{c}{`sdid_ever'} \\" _n
file write sdid_tex "\bottomrule" _n
file write sdid_tex "\end{tabular}" _n
file write sdid_tex "\begin{tablenotes}\footnotesize" _n
file write sdid_tex "\item \textit{Notes:} Synthetic Difference-in-Differences (SDID) following Arkhangelsky, Athey, Hirshberg, Imbens, and Wager (2021), estimated using the \texttt{sdid} Stata package of Clarke, Pailanir, Athey, and Imbens (2023). SDID is applied to a balanced sub-panel of Czech construction firms observed in every panel year under the absorbing-treatment reformulation. Standard errors from 50 bootstrap replications. The balanced-panel constraint drops firms with incomplete coverage; the small ever-treated count reflects the restriction that firms must be observed in all 15 panel years." _n
file write sdid_tex "\end{tablenotes}" _n
file write sdid_tex "\end{threeparttable}" _n
file write sdid_tex "\end{table}" _n
file close sdid_tex

display "Saved: panel_sdid_honestdid_results.csv"
display "Saved: tables/panel_sdid.tex"
display "HonestDiD LaTeX table will be produced by shell post-processor from the log"

log close
exit
