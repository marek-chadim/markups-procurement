use "../../output/stata/beer_c2_estimates.dta", clear
format b_* %9.4f
format markup_mean criterion %9.6f
list sample spec b_k b_cogs b_k2 b_cogs2 b_kcogs markup_mean criterion, noobs sep(0)
