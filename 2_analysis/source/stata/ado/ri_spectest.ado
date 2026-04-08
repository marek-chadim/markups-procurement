*! ri_spectest: Conducts randomization-based specification tests by regressing the recentered instrument on predetermined controls
*! Reference: Borusyak and Hull (2023) "Non-Random Exposure to Exogenous Shocks"
*! Version: September 23, 2020
*! Author: Kirill Borusyak
/*
Syntax: ri_spectest z_rc controls, sim(varlist) [optional parameters]
- z_rc is the recentered(!) instrument
- the list of predetermined controls can be empty => then just a constant is used
- sim: varlist that captures simulated versions of z_rc computed from counterfactual shocks; typical usage is sim(z_rc1-z_rc999)

Additional options:
- weights, if, and in clauses are allowed
- absorb: predetermined FE
- INDividual: besides the joint test pvalue, report pvalues for individual coefficients
- When individual is specified, two additional options are available:
	- inddemean: demeans all controls. Useful to make the constant and its pvalue informative
	- indonesided: report one-sided test pvalues for individual coefs (two-sided is the recommended default)

Returns: 
- r(jointp): the pvalue of joint specification test, that the mean-squared fitted value is not in the tails of the simulation distribution
- r(p`control'): if `individual' is specified, pvalues for individual-coefficient tests are reported, for each control
- r(p_cons): if `individual' is specified and there are no `absorb' variables, the pvalue for the constant-based test is reported (but useful only with `inddemean'!)
- r(nsims): the number of simulations used, as defined by sim()
*/

cap program drop ri_spectest
program def ri_spectest, rclass
syntax varlist [aw iw fw pw] [if] [in] , sim(varlist) [Absorb(string) INDividual inddemean indonesided]
	/* To-do: 
		- adjustments for (k+1)/(n+1) or equal statistics (esp. important for the constant when the sims are demeaned!)
		- check for conflicts of parameters for inddemean and indonesided
		- report coef values (e.g. for the constant after demeaning) and mean-squared fitted values
		- report critical values of the RI test
	*/
qui {
   // Prepare inputs
	tokenize `varlist'
	local z `1'
	macro shift
	if ("`*'"!="") unab controls : `*'
	
	marksample touse
	markout `touse' `varlist' `absorb'
	if ("`absorb'"!="" & "`constant'"!="") {
	    di as error "Options absorb and noconstant cannot be combined"
		error 184
	}
	
	if ("`inddemean'"!="") {
		local newcontrols = ""
		foreach c in `controls' {
			tempvar c`c'
			sum `c' [`weight'`exp'] if `touse'
			gen `c`c'' = `c'-r(mean)
			local newcontrols `newcontrols' `c`c''
		}
		local oldcontrols `controls'
		local controls `newcontrols'
	}
	
	tempvar fitted
	if ("`absorb'"=="") {
		reg `z' `controls' [`weight'`exp'] if `touse'
		predict `fitted' if `touse', xb
		local includeconst "_cons"
	}
	else {
		reghdfe `z' `controls' [`weight'`exp'] if `touse', a(`absorb') resid(`fitted')
		replace `fitted' = `z'-`fitted' if `touse'
		local includeconst "" // with FE don't have the constant
	}
	replace `fitted' = `fitted'^2
	sum `fitted' [`weight'`exp'] if `touse'
	local stat = r(mean)
	drop `fitted'
	if ("`individual'"!="") {
		foreach c in `includeconst' `controls' {
			local s_`c' = _b[`c']
			local nl_`c' = 0 // nlarger for this var
		}
	}

	local nlarger = 0 
	local ntotal = 0
	foreach v of varlist `sim' { 
		if ("`absorb'"=="") {
			reg `v' `controls' [`weight'`exp'] if `touse'
			predict `fitted' if `touse', xb
		}
		else {
			reghdfe `v' `controls' [`weight'`exp'] if `touse', a(`absorb') resid(`fitted')
			replace `fitted' = `v'-`fitted' if `touse'
		}
		replace `fitted' = `fitted'^2
		sum `fitted' [`weight'`exp'] if `touse'
		local nlarger = `nlarger' + (`stat'< r(mean))
		drop `fitted'
		if ("`individual'"!="") {
			foreach c in `includeconst' `controls' {
				local nl_`c' = `nl_`c'' + (`s_`c''<_b[`c'])
			}
		}
		local ++ntotal
	}
		
	return scalar jointp = `nlarger'/`ntotal' // fraction of simulations with the simulated statistic larger than the actual one
	if ("`individual'"!="") {
		if ("`inddemean'"=="") {
			foreach c in `includeconst' `controls' {
				local rejrate = `nl_`c''/`ntotal'
				if ("`indonesided'"!="") return scalar p`c' = `rejrate'
					else return scalar p`c' = 2*min(`rejrate',1-`rejrate')
			}
		} 
		else { // renamed the control vars so need to report correct names
			foreach c in `includeconst' {
				local rejrate = `nl_`c''/`ntotal'
				if ("`indonesided'"!="") return scalar p`c' = `rejrate'
					else return scalar p`c' = 2*min(`rejrate',1-`rejrate')
			}
			foreach c in `oldcontrols' {
				local rejrate = `nl_`c`c'''/`ntotal'
				if ("`indonesided'"!="") return scalar p`c' = `rejrate'
					else return scalar p`c' = 2*min(`rejrate',1-`rejrate')
			}
		}
	}
	return scalar nsims = `ntotal'
}
end
