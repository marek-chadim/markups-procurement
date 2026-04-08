*! ri_pvalue: Computes randomzation pvalue for OLS & IV regressions.
*! Reference: Borusyak and Hull (2023) "Non-Random Exposure to Exogenous Shocks"
*! Version: September 20, 2020
*! Author: Kirill Borusyak
/*
Main syntax: ri_pvalue y x [z], zsim(varlist) [beta(real) [optional parameters]
- Here y and x define the structural equation. 
- If specified, z defines the instrument; otherwise z=x is assumed
- The zsim varlist captures simulated versions of z (or x, for OLS) computed from counterfactual shocks; typical usage is zsim(z*) or zsim(z1-z999).
- beta defines the null hypothesis, 0 by default. 
- Note: When beta=0, this is a traditional randomization test

Additional options:
- weights, if, and in clauses are allowed
- controls and absorb: additional controls and FE. The controls have to be predetermined; this function does not work with multiple treatments!
- noconstant: a regression without constant (cannot be combined with absorb).
- Note: If using many times for different betas (like in the confidence interval construction), it's faster to residualize y and x on controls in advance and use the noconstant option
- onesided: if specified, a one-sided pvalue is returned; two-sided is the recommended default

Returns: 
- r(pvalue): the pvalue of the two-sided (unless onesided specified) test
- r(nsims): the number of simulations used, as defined by zsim()
- r(rejrate): the one-sided pvalue, for programming use only
- r(stat): the value of the RI statistic
*/

cap program drop ri_pvalue
program def ri_pvalue, rclass
syntax varlist(min=2 max=3 numeric) [aw iw fw pw] [if] [in] , zsim(varlist) [Beta(real 0) CONTrols(string) Absorb(string) noCONstant ONEsided]
	// To-do: adjustments for (k+1)/(n+1) (i.e. include the actual realization as one of the counterfactuals); breaking the ties; ritest functions (save & plot the simulations, produce permutations in other ways)
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
	
	// Prepare residual
	tempvar eps eps_resid
	gen `eps' = `y'-`beta'*`x'
	if ("`controls'"=="" & "`absorb'"=="" & "`constant'"!="") gen `eps_resid' = `eps'
	else { // residualize unless there are no controls, not even a constant
	    if ("`absorb'"=="") {
		    reg `eps' `controls' [`weight'`exp'] if `touse', `constant'
			predict `eps_resid' if `touse', resid
		}
		else reghdfe `eps' `controls' [`weight'`exp'] if `touse', a(`absorb') resid(`eps_resid')
	}
	
	tempvar zeps
	gen `zeps' = `z' * `eps_resid' if `touse'
	sum `zeps' [`weight'`exp'] if `touse'
	local stat = r(sum)
	
	local nlarger = 0 
	local ntotal = 0
	foreach v of varlist `zsim' { 
		replace `zeps' = `v' * `eps_resid' if `touse'
		sum `zeps' [`weight'`exp'] if `touse'
		local nlarger = `nlarger' + (`stat'< r(sum))
		local ++ntotal
	}
	local rejrate = `nlarger'/`ntotal' // fraction of simulations with the simulated statistic larger than the actual one
	return scalar rejrate = `rejrate'
	if ("`onesided'"!="") return scalar pvalue = `rejrate'
		else return scalar pvalue = 2*min(`rejrate',1-`rejrate')
	return scalar nsims = `ntotal'
	return scalar stat = `stat'
	drop `eps' `eps_resid'
}
end
