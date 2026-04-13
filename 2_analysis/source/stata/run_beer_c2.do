* Wrapper to run table_beer_c2.do with launcher globals set up
clear all
set more off
set seed 42

global root    "."
global datagen "$root/datagen"
global code    "$root/code"
global data    "../../output/stata"
global output  "../../output/stata"
global temp    "../../output/stata/temp"

cap log close
log using "$output/table_beer_c2.log", text replace

do "$code/table_beer_c2.do"

log close
