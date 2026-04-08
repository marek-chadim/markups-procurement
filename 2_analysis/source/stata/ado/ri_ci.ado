*! ri_ci: Computes randomzation confidence intervals for OLS & IV regressions. 
*! Reference: Borusyak and Hull (2023) "Non-Random Exposure to Exogenous Shocks"
*! Version: September 20, 2020
*! Author: Kirill Borusyak
/*
Main syntax: ri_ci y x [z], zsim(varlist) range(b0 b1) [optional parameters]
- Here y and x define the structural equation.
- If specified, z defines the instrument; otherwise z=x is assumed
- The zsim varlist captures simulated versions of z (or x, for OLS) computed from counterfactual shocks; typical usage is zsim(z*) or zsim(z1-z999).
- Since the method uses binary search, you need to specify the interval of betas where the search happens [b0 b1]. If you get a corresponding error message, extend the interval.
- Note: binary search only works well if the randomization pvalue is monotone over the beta parameter; contact the authors for a slow grid-search version if this is violated in your case (but some warning messages about non-monotonicities are fine).

Additional options:
- weights, if, and in clauses are allowed
- controls and absorb: additional controls and FE. The controls have to be predetermined; this function does not work with multiple treatments!
- noconstant: a regression without constant (cannot be combined with absorb)
- alpha: significance level (5% by default)
- tol, maxit: precision in the beta space and the max number of iterations for the binary search

Returns: r(ci_left) and r(ci_right) define the ends of the two-sided confidence interval. Please note that sometimes they can be infinite on one or both sides, or even empty (where the latter, and only it, suggests misspecification)
*/

cap program drop ri_ci
program ri_ci, rclass
syntax varlist(min=2 max=3 numeric) [aw iw fw pw] [if] [in], zsim(varlist) Range(numlist min=2 max=2 sort) [alpha(real 0.05) extreme(numlist min=2 max=2 sort) tol(real 0.0001) maxit(integer 50) CONTRols(string) Absorb(string) noCONstant]
qui {	
   // Prepare inputs
	tokenize `varlist'
	local y `1'
	local x `2'
	if (`"`3'"'=="") local z `x'
		else local z `3'
	marksample touse
	markout `touse' `varlist' `zsim'
	if ("`controls'"!="" | "`absorb'"!="") markout `touse' `controls' `absorb'
	if ("`absorb'"!="" & "`constant'"!="") {
	    di as error "Options absorb and noconstant cannot be combined"
		error 184
	}
	if ("`extreme'"=="") {
		local extreme_left = -10^7
		local extreme_right = 10^7
	}
	else {
		tokenize `extreme'
		local extreme_left = `1'
		local extreme_right = `2'
	}
	tokenize `range'
	local range_left = `1'
	local range_right = `2'
	
	// Residualize y and x
	if ("`controls'"!="" | "`absorb'"!="" | "`constant'"=="") {
		foreach v in y x {
			local vvar ``v''
			tempvar `v'perp
			if (`"`absorb'"'=="") {
				reg `vvar' `controls' [`weight'`exp'] if `touse', `constant'
				predict ``v'perp' if `touse', resid
			}
			else reghdfe `vvar' `controls' if `touse' , absorb(`absorb') resid(``v'perp')
			local `v' ``v'perp'
		}
	}
	
	local alpha1 = `alpha'/2
	local alpha2 = 1-`alpha1'
	local extrapar = "zsim(`zsim') nocon" // these parameters are passed to ri_pvalue (always nocon since already residualized)
	local jlist = ""
	
	// Check extreme values first, decide if CI is infinite on either side, potentially empty [unless non-monotonic rejrate], or good
	ri_pvalue `y' `x' `z', beta(`extreme_left') `extrapar'
		local pextrleft = r(pvalue)
		local rejextrleft = r(rejrate)
	ri_pvalue `y' `x' `z', beta(`extreme_right') `extrapar'
		local pextrright = r(pvalue)
		local rejextrright = r(rejrate)
	
	if !inrange(`alpha1',min(`rejextrleft',`rejextrright')-10^-6,max(`rejextrleft',`rejextrright')+10^-6) & ///
				!inrange(`alpha2',min(`rejextrleft',`rejextrright')-10^-6,max(`rejextrleft',`rejextrright')+10^-6) {
		if (`pextrleft'>=2*`alpha1' & `pextrright'>=2*`alpha1') { // can't reject anything ([alpha1,alpha2] \superset [rejextrleft,rejextrright])
			di as error "Warning: CI is infinite on both sides"
			return scalar ci_left = `extreme_left'
			return scalar ci_right = `extreme_right'
		}
		else if (`pextrleft'<2*`alpha1' & `pextrright'<2*`alpha1') { // reject everything ([alpha1,alpha2] doesn't overlap with [rejextrleft,rejextrright])
			di as error "Warning: CI is empty"
			return scalar ci_left = .
			return scalar ci_right = .
		}
		else di as error "Something strange happened #1"
		exit
	}
	else if !inrange(`alpha1',min(`rejextrleft',`rejextrright')-10^-6,max(`rejextrleft',`rejextrright')+10^-6) {
		local jlist = "2" // only search for intersection with alpha2
		di as error "Warning: CI is infinite on one side"
		if (`pextrleft'>=`alpha') local CI1 = `extrleft' // decide which extreme value to include in the CI	
			else if (`pextrright'>=`alpha') local CI1 = `extrright'
			else di as error "Something strange happened #2"
	}
	else if !inrange(`alpha2',min(`rejextrleft',`rejextrright')-10^-6,max(`rejextrleft',`rejextrright')+10^-6) {
		local jlist = "1" // only search for intersection with alpha1
		di as error "Warning: CI is infinite on one side"
		if (`pextrleft'>=`alpha') local CI2 = `extrleft' // decide which extreme value to include in the CI
			else if (`pextrright'>=`alpha') local CI2 = `extrright'
			else di as error "Something strange happened #3"
	}
	else local jlist = "1 2" 
	
	// Now run within the range, making sure the ends of the CI are not outside the range (but within the extreme values)
	local left = `range_left' // specify the interval on which to do binary search. The program will complain if the test is not rejected at the boundaries 
	local right = `range_right'
	ri_pvalue `y' `x' `z', beta(`left') `extrapar'
		local rejleft = r(rejrate)
	ri_pvalue `y' `x' `z', beta(`right') `extrapar'
		local rejright = r(rejrate)

	if (!inrange(`alpha1',min(`rejleft',`rejright')-10^-6,max(`rejleft',`rejright')+10^-6) & "`jlist'"!="2") | ///
					(!inrange(`alpha2',min(`rejleft',`rejright')-10^-6,max(`rejleft',`rejright')+10^-6) & "`jlist'"!="1") {
		di as error "Boundary hit. Please extend the search interval in the range option. Left: (`left', `rejleft') Right: (`right', `rejright')"
		exit 498
	}

	foreach j in `jlist' { // correspond to two ends of CI. Which one is left depend on whether the rejrate is increasing or falling in beta
		local l = `left'
		local r = `right'
		local rejl = `rejleft'
		local rejr = `rejright'
		local it`j' = 0
		while (`r'-`l'>2*`tol' & `it`j''<`maxit') { // at most `maxit' iterations
			ri_pvalue `y' `x' `z', beta(`=(`l'+`r')/2') `extrapar' // test the middle, then decide which direction to go (or detect non-monotonicity)
			local current = r(rejrate)
			if !inrange(`current',min(`rejl',`rejr')-10^-6,max(`rejl',`rejr')+10^-6) {
				di as error "Warning: non-monotonicity for `alpha`j''. Left: (`l', `rejl') Right: (`r', `rejr') Middle: `current'"
			}
			if inrange(`alpha`j'',min(`rejl',`current'),max(`rejl',`current')) {
				local r = (`l'+`r')/2
				local rejr = r(rejrate)
			}
			else {
				local l = (`l'+`r')/2
				local rejl = r(rejrate)
			}
			local ++it`j'
		}
		local CI`j' = (`l'+`r')/2
	}
	
	if (`CI1'>`CI2') {
		local CI = `CI1'
		local CI1 = `CI2'
		local CI2 = `CI'
	}
//	noi di "Converged in `it1'+`it2' iterations; CI: [`CI1',`CI2'] for `alpha0'" // could check that actually converged
	return scalar ci_left = `CI1' 
	return scalar ci_right = `CI2'
}
end
