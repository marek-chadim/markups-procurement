*===============================================================================
* figure_ri_reform_2012.do — B&T-style randomization inference on the
*                            2012 reform DiD coefficient
*
* Baranek-Titl (2024 JLE) placebo.do randomly permutes the "connection"
* indicator within election cycles and re-runs the main DiD, building a
* null distribution of 200 permutation coefficients and comparing the
* observed coefficient to it via kdensity. This do-file applies the same
* logic to my 2012 single-bidding-ban DiD on log markup:
*
*   - Permute ever_pp within NACE 2-digit strata (preserving empirical
*     treatment rate by sub-industry)
*   - Re-run reghdfe log_mu tp2012_perm post2012 k cogs, absorb(id year)
*     cluster(id), storing the coefficient on tp2012_perm
*   - 500 iterations (Stata is slow; B&T uses 200)
*   - Plot histogram of permuted coefficients with observed coefficient as
*     a vertical reference line; report one-sided RI p-value = Pr(beta_perm
*     >= beta_real).
*
* This is COMPLEMENTARY to table_bh_ri.do, which tests the cross-sectional
* 14% premium. Here the object is the reform-DiD coefficient (+0.029), a
* different estimand testing a different null.
*
* Input:  $data/markups_panel.dta
* Output: 4_paper/input/analysis_figures/ri_reform_2012.pdf
*===============================================================================

dis _newline "--- figure_ri_reform_2012.do ---"

do "$code/graph_markups.do"

use "$data/markups_panel.dta", clear

bys id: egen ever_pp_real = max(pp_dummy)
gen log_mu = log(mu_A)
gen post2012 = (year >= 2012)
gen tp_real = ever_pp_real * post2012

* Observed DiD coefficient (same spec as table_reforms_mechanism.do)
cap confirm var nace2
if _rc != 0 {
    dis "  ERROR: nace2 not found — skipping"
    exit
}

qui reghdfe log_mu tp_real post2012 k cogs, absorb(id year) vce(cluster id)
local b_real = _b[tp_real]
local se_real = _se[tp_real]
dis "  Observed DiD: " %7.4f `b_real' " (SE " %7.4f `se_real' ")"

* ── Permute ever_pp within NACE 2-digit strata ──
* Collapse to one row per firm, permute within NACE via self-merge on
* (nace2, rank_key). The natural rank (pos_orig) is deterministic by id;
* the shuffled rank (pos_shuf) comes from sorting on a random uniform key.
* Lookup table: (nace2, pos_shuf → ever_pp_real); joined back on
* (nace2, pos_orig = pos_shuf) gives each firm a shuffled value.
preserve
keep id nace2 ever_pp_real
bys id: keep if _n == 1
count
local n_firms = r(N)
dis "  Firm-level obs: " `n_firms'
tempfile firm_level
save `firm_level'
restore

set seed 42
local K = 500
local count_ge = 0
matrix ri_betas = J(`K', 1, .)

* Pre-build once: the firm's deterministic pos_orig (doesn't change across iter)
preserve
use `firm_level', clear
bys nace2 (id): gen pos_orig = _n
tempfile firm_base
save `firm_base'
restore

forvalues k = 1/`K' {
    preserve
    use `firm_base', clear
    * Random position per firm within nace2
    gen u = runiform()
    sort nace2 u
    by nace2: gen pos_shuf = _n
    * Build lookup table in a separate tempfile (no nested preserve)
    tempfile lookup
    * First write the (nace2, pos_shuf, ever_pp_real) triples
    qui save `lookup', replace
    * Rename on disk via a sub-use
    use `lookup', clear
    keep nace2 pos_shuf ever_pp_real
    rename pos_shuf pos_key
    rename ever_pp_real ever_pp_perm
    qui save `lookup', replace
    * Now reload firm_base and merge: pos_orig as the join key
    use `firm_base', clear
    rename pos_orig pos_key
    qui merge 1:1 nace2 pos_key using `lookup', nogen keep(master matched)
    keep id ever_pp_perm
    tempfile firm_perm
    qui save `firm_perm'
    restore

    preserve
    qui merge m:1 id using `firm_perm', nogen keepusing(ever_pp_perm)
    gen tp_perm = ever_pp_perm * post2012
    qui reghdfe log_mu tp_perm post2012 k cogs, absorb(id year) vce(cluster id)
    matrix ri_betas[`k', 1] = _b[tp_perm]
    if abs(_b[tp_perm]) >= abs(`b_real') local ++count_ge
    restore

    if mod(`k', 50) == 0 dis "    iteration " `k' " / " `K'
}

local ri_p = `count_ge' / `K'
dis "  RI p-value (two-sided, |perm| >= |real|): " %6.4f `ri_p'
dis "  (`count_ge' out of `K' permutations at least as extreme)"

* ── Build histogram ──
preserve
clear
svmat double ri_betas, name(b_perm)
rename b_perm1 b_perm
summ b_perm, d
local mean_perm = r(mean)
local sd_perm = r(sd)
local p2 = r(p1)
local p98 = r(p99)
dis "  Null distribution: mean=" %6.4f `mean_perm' " sd=" %6.4f `sd_perm'

twoway ///
    (hist b_perm, frequency bin(40) ///
        fcolor("${markups_blue}%40") lcolor("${markups_blue}") lwidth(vthin)) ///
    , ${markups_gropts} ///
    xtitle("Permuted DiD coefficient on ever_pp x post2012") ///
    ytitle("Frequency across 500 permutations") ///
    xline(`b_real', lcolor("${markups_red}") lpattern(dash) lwidth(medthick)) ///
    xline(0, lcolor(gs10) lpattern(dot)) ///
    xscale(range(-0.04 0.04)) ///
    xlabel(-0.04(0.01)0.04) ///
    title("Randomization inference: 2012 reform DiD on log markup", size(medium) position(11)) ///
    subtitle("Observed b = `=string(`b_real', "%5.4f")', RI p = `=string(`ri_p', "%5.4f")'", size(small)) ///
    note("Null: permute ever_pp within NACE 2-digit strata (500 iter); rerun reghdfe. Dashed red = observed.", size(vsmall)) ///
    legend(off) ///
    name(g_ri, replace)

cap mkdir "../../../4_paper/input/analysis_figures"
graph export "../../../4_paper/input/analysis_figures/ri_reform_2012.pdf", ///
    replace
dis "  Saved: ri_reform_2012.pdf"

* Also save the null distribution + observed for audit
gen observed = `b_real'
gen ri_p_val = `ri_p'
save "$output/data/ri_reform_2012_null.dta", replace
dis "  Saved: ri_reform_2012_null.dta"

graph drop g_ri
restore
