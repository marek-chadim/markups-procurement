* Install oster and sensemakr from GitHub
set checksum off

cap net install oster, from("https://raw.githubusercontent.com/emilyoster/oster_stata/master") replace
cap which psacalc
dis "oster/psacalc: " _rc

cap ssc install psacalc, replace
cap which psacalc
dis "psacalc via ssc: " _rc

cap net install sensemakr, from("https://raw.githubusercontent.com/chadhazlett/sensemakr/master/stata") replace
cap which sensemakr
dis "sensemakr: " _rc

* bacondecomp — verify it's loadable
cap which bacondecomp
dis "bacondecomp: " _rc

* csdid
cap which csdid
dis "csdid: " _rc

* teffects (built in)
cap which teffects
dis "teffects: " _rc
