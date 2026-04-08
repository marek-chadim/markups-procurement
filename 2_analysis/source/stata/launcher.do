*===============================================================================
* launcher.do — Master launcher for empirical analysis
*
* Adapted from De Ridder, Grassi & Morzenti (2026, Econometrica)
* "Hitchhiker's Guide to Markup Estimation" Replication Package
*
* Context: Czech construction procurement (Chadim 2026)
*   - 3 sub-industries: NACE 41 (buildings), 42 (civil eng), 43 (specialized)
*   - Production function inputs: go (output), k (capital), cogs (materials)
*   - Treatment: pp_dummy (public procurement participation)
*
* DGM adaptations:
*   - Their 4 inputs (v,k,m,o) → our 2 inputs (k, cogs)
*   - Their quantity vs revenue comparison → our specification robustness
*   - Their 21 industries → our 3 NACE sub-industries
*   - Added: procurement premium analysis (not in DGM)
*===============================================================================

clear all
set more off
set seed 42

* Set paths — Gentzkow Lab Template structure
* Run from: markups-procurement/2_analysis/source/stata/
global root    "."
global datagen "$root/datagen"
global code    "$root/code"
global data    "../../output/stata"
global output  "../../output/stata"
global temp    "../../output/stata/temp"

* Source data (produced by rebuild_data.py in module 1_data)
global srcdata "../../input/data.dta"

* Create output directories
cap mkdir "$data"
cap mkdir "$output"
cap mkdir "$temp"

cap log close
log using "$output/launcher.log", text replace

timer on 1

*-----------------------------------------------------------------------
* STEP 1: DATA GENERATION
*-----------------------------------------------------------------------
dis _newline(2) "========================================"
dis "  STEP 1: Data Generation"
dis "========================================"

* 1a. Prepare analysis datasets
do "$datagen/prepare_data.do"

* 1b. Estimate production functions (ACF, multiple specs)
do "$datagen/estimate_pf.do"

* 1c. Calculate firm-level markups
do "$datagen/calculate_markups.do"

*-----------------------------------------------------------------------
* STEP 2: ANALYSIS (TABLES AND FIGURES)
*-----------------------------------------------------------------------
dis _newline(2) "========================================"
dis "  STEP 2: Analysis"
dis "========================================"

do "$code/master_code.do"

*-----------------------------------------------------------------------
timer off 1
timer list

dis _newline "========================================"
dis "  DONE. Output in: $output"
dis "========================================"

log close
