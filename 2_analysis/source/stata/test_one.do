clear all
set more off
global data    "../../output/stata"
global output  "../../output/stata"
global code    "./code"

cap log close
log using "../../output/stata/test_one.log", text replace

do "./code/table_summary_stats.do"

log close
