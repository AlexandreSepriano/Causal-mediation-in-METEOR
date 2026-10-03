* msmbalance Version 1
* Baseline / balance tables with standardised mean differences using Jackson formula
* Reference: Jackson JW. Diagnostics for confounding of time-varying and other joint exposures. Epidemiology. 2016;27(6):859-69.
* Save as msmbalanceV1.ado in ado directory to run without calling

capture program drop msmbalanceV1
program define msmbalanceV1, rclass
	version 15.0

	syntax varlist(numeric) ,	treatvar(varname)	///
					weight(varname)		///
					[			///
					wraw(varname)		///
					touse(string)		///
					label(string)		///
					minarm(integer 10)	///
					cutoff(real 0.10)	///
					folder(string)		///
					filename(string)	///
					save(string)		///
					]

	local savetable = lower(trim("`save'"))
	if "`savetable'"=="" local savetable "no"

	if "`touse'"=="" local touse "1"

	tempvar use
	qui gen byte `use' = (`touse') & !missing(`treatvar') & !missing(`weight')

	qui count if `use' & `treatvar'==1
	local n1 = r(N)
	qui count if `use' & `treatvar'==0
	local n0 = r(N)

	di ""
	di as text "{hline 78}"
	di as text "`label'"
	di as text "{hline 78}"
	di as result "  original population: treated `n1', untreated `n0', total `=`n1'+`n0''"

	if `n1' < 2 | `n0' < 2 {
		di as error "  NOT ESTIMABLE: fewer than two patients in one arm."
		exit
	}

	qui sum `weight' if `use' & `treatvar'==1
	local w1 = r(sum)
	qui sum `weight' if `use' & `treatvar'==0
	local w0 = r(sum)
	qui sum `weight' if `use'
	local wt = r(sum)
	local mw = r(mean)

	di as result "  pseudo-population:   treated `=round(`w1')', untreated `=round(`w0')', total `=round(`wt')' (mean weight " %5.3f `mw' ")"

	* approximate sampling SE of an SMD, evaluated at d = 0
	local se = sqrt((`n1'+`n0')/(`n1'*`n0'))
	di as result "  approximate SE of an SMD here: " %4.2f `se' "  (smaller arm n = `=min(`n1',`n0')')"

	if min(`n1',`n0') < `minarm' {
		di as error "  UNINFORMATIVE: the smaller arm is too small for the SMD to carry evidence."
		di as error "  Report this panel descriptively; do not apply the `cutoff' threshold to it."
	}
	else if `se' >= 0.10 {
		di as error "  CAUTION: the SE reaches the 0.10 convention, so individual rows near the"
		di as error "  threshold cannot be distinguished from sampling variation."
	}

	local nv : word count `varlist'
	local ncol = 6
	if "`wraw'" != "" local ncol = 7

	tempname M
	matrix `M' = J(`nv',`ncol',.)

	local i = 0
	local rn ""
	foreach x of local varlist {
		local ++i
		local rn "`rn' `x'"

		* is the covariate binary in this panel?
		local isbin = 1
		qui levelsof `x' if `use', local(lv)
		foreach l of local lv {
			if !inlist(`l',0,1) local isbin = 0
		}

		* unweighted means, and the unweighted pooled SD that standardises
		* every column of this row (Jackson, Table 1 footnote a)
		qui sum `x' if `use' & `treatvar'==1
		local m1 = r(mean)
		local s1 = r(Var)
		qui sum `x' if `use' & `treatvar'==0
		local m0 = r(mean)
		local s0 = r(Var)
		if `isbin'==1 {
			local s1 = `m1'*(1-`m1')
			local s0 = `m0'*(1-`m0')
		}
		local den = sqrt((`s1'+`s0')/2)

		if `den' > 0 & `den' < . {
			local usmd = (`m1'-`m0')/`den'
		}
		else {
			local usmd = .
		}

		* weighted means, same denominator
		qui sum `x' if `use' & `treatvar'==1 [aw=`weight']
		local wm1 = r(mean)
		qui sum `x' if `use' & `treatvar'==0 [aw=`weight']
		local wm0 = r(mean)
		if `den' > 0 & `den' < . {
			local wsmd = (`wm1'-`wm0')/`den'
		}
		else {
			local wsmd = .
		}

		matrix `M'[`i',1] = `m1'
		matrix `M'[`i',2] = `m0'
		matrix `M'[`i',3] = `usmd'
		matrix `M'[`i',4] = `wm1'
		matrix `M'[`i',5] = `wm0'
		matrix `M'[`i',6] = `wsmd'

		if "`wraw'" != "" {
			qui sum `x' if `use' & `treatvar'==1 [aw=`wraw']
			local rm1 = r(mean)
			qui sum `x' if `use' & `treatvar'==0 [aw=`wraw']
			local rm0 = r(mean)
			if `den' > 0 & `den' < . {
				matrix `M'[`i',7] = (`rm1'-`rm0')/`den'
			}
		}
	}

	if "`wraw'"=="" {
		matrix colnames `M' = "Mean treated" "Mean untreated" "Unweighted SMD" ///
		                      "Mean treated" "Mean untreated" "Weighted SMD"
	}
	else {
		matrix colnames `M' = "Mean treated" "Mean untreated" "Unweighted SMD" ///
		                      "Mean treated" "Mean untreated" "Weighted SMD" "Weighted SMD untrimmed"
	}
	matrix rownames `M' = `rn'

	matlist `M', format(%9.2f) rowtitle("Variable") ///
		title({bf} `label')

	* flag the rows above the cut-off in the weighted column
	local unbal ""
	local i = 0
	foreach x of local varlist {
		local ++i
		if abs(`M'[`i',6]) >= `cutoff' & abs(`M'[`i',6]) < . {
			local unbal "`unbal' `x'"
		}
	}
	if "`unbal'" != "" {
		di as text "  |weighted SMD| >= `cutoff': " as error "`unbal'"
	}
	else {
		di as result "  no |weighted SMD| >= `cutoff'"
	}

	if "`savetable'"=="yes" {
		putexcel clear
		putexcel set "`folder'/`filename'.xlsx", sheet("`filename'") replace
		putexcel A1 = matrix(`M'), names
		di as text "  saved: `folder'/`filename'.xlsx"
	}

	return scalar n1 = `n1'
	return scalar n0 = `n0'
	return scalar se = `se'
	return matrix balance = `M'
end
