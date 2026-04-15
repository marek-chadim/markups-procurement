use "../../output/stata/markups_panel.dta", clear
dis _N " obs"
cap confirm var markup_A
if _rc == 0 dis "markup_A: OK"
else dis "markup_A: MISSING"
ds markup*
ds go k cogs pp_dummy nace2 year id
