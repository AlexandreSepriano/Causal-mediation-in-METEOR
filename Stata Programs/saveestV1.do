******************************************************************************************************************************************************
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
**** saveestV1: collect the Stata estimates and diagnostics for the supplementary tables
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
******************************************************************************************************************************************************

*===============================================================================
* What this program does
*===============================================================================
* Supporting program for Analysis.do. It appends every Stata estimate and every
* diagnostic needed to fill Supplementary Tables S2 to S7 and S10 (and their footnotes)
* to two csv files in the Tables folder:
*
*   estimates_stata.csv    table,estimator,effect,estimate,se,lci,uci,evalue,reps,n
*   diagnostics_stata.csv  table,estimator,diagnostic,value
*
* The two csv files are the accumulator: each run of Analysis.do is a separate
* execution, so the rows have to be appended to a plain text file as they are
* produced. Run 8 of Analysis.do then calls savetabs, which compiles them into
*
*   estimates_stata.xlsx   one sheet per supplementary table
*
* Decimals shown in the workbook are set by $xldec and $xldecd just below.
*
* Rows accumulate across runs: execute run 1 through run 6.1 in turn and the two
* files fill up. Erase both files to start a clean set. The R estimates and the
* R diagnostics are added to the tables separately.
*
* Programs defined here
*   savest      one estimate row (generic: the other programs call it)
*   savediag    one diagnostic row
*   savebs      every effect after: bootstrap ... then estat bootstrap, percentile
*               (also used after gformula, which posts e(b) and e(V))
*   savemedeff  the four effects, the total and the two averages after: medeff
*   savete      the ATE after: teffects ..., ate
*
* Globals expected from Analysis.do
*   $tables   output folder
*   $anaN     patients in the analysis
*   $anaR     bootstrap replications requested
*
* Labels passed to these programs must not contain commas: the two files are plain csv.
*===============================================================================

global estfile  "${tables}estimates_stata.csv"
global diagfile "${tables}diagnostics_stata.csv"
global xlsfile  "${tables}estimates_stata.xlsx"

* Decimals shown in the workbook. These are display formats only: the values
* themselves are written at full precision, so widening them costs nothing.
global xldec  = 2 // estimates, standard errors, confidence limits and e-values
global xldecd = 3 // diagnostics and p-values


*===============================================================================
* savest: write one estimate row
*===============================================================================

capture program drop savest
program define savest
	args tbl est eff b se lo hi ev reps n

	* Some estimators report an interval but no standard error (medeff, and
	* gformula when it posts no variance matrix). It is backed out of the
	* percentile interval so the SE column is never empty:
	*   se = (uci - lci) / (2 x 1.96)
	if "`se'" == "" local se "."
	if "`lo'" == "" local lo "."
	if "`hi'" == "" local hi "."
	if `se' >= . & `lo' < . & `hi' < . local se = (`hi' - `lo')/(2*1.96)

	capture file close sv // in case a previous run left the handle open
	capture confirm file "$estfile"
	if _rc {
		file open sv using "$estfile", write text replace
		file write sv "table,estimator,effect,estimate,se,lci,uci,evalue,reps,n" _n
	}
	else file open sv using "$estfile", write text append
	file write sv "`tbl',`est',`eff',`b',`se',`lo',`hi',`ev',`reps',`n'" _n
	file close sv
end


*===============================================================================
* savediag: write one diagnostic row
*===============================================================================

capture program drop savediag
program define savediag
	args tbl est diag val
	capture file close dg // in case a previous run left the handle open
	capture confirm file "$diagfile"
	if _rc {
		file open dg using "$diagfile", write text replace
		file write dg "table,estimator,diagnostic,value" _n
	}
	else file open dg using "$diagfile", write text append
	file write dg "`tbl',`est',`diag',`val'" _n
	file close dg
end


*===============================================================================
* savebs: write every effect held in e(b)
*===============================================================================
* Used after "bootstrap ... " followed by "estat bootstrap, percentile", where the
* percentile bounds sit in e(ci_percentile) and the bootstrap standard errors in
* e(se). Also used after gformula, which posts e(b) and e(V) but neither of those:
* the standard error is then taken from e(V) and the interval is normal-based.
* The third argument is optional and carries the e-value of the first effect.
* The fourth is optional too: a list of effect names to keep, so one estimation
* can fill two tables - the four main paths in one, the whole set in another.

capture program drop savebs
program define savebs
	args tbl est ev keep
	if "`ev'" == "" local ev "."
	tempname B S V C
	matrix `B' = e(b)
	capture matrix `S' = e(se)
	local hasse = (_rc == 0)
	capture matrix `V' = e(V)
	local hasv = (_rc == 0)
	capture matrix `C' = e(ci_percentile)
	local hasci = (_rc == 0)
	local reps = e(N_reps)
	if `reps' >= . local reps = $anaR
	local nm : colnames `B'
	forvalues j = 1/`=colsof(`B')' {
		local eff : word `j' of `nm'
		if "`keep'" != "" & strpos(" `keep' ", " `eff' ") == 0 continue
		local b = `B'[1,`j']
		local s = .
		if `hasse' local s = `S'[1,`j']
		if `hasse' == 0 & `hasv' local s = sqrt(`V'[`j',`j'])
		local lo = .
		local hi = .
		if `hasci' {
			local lo = `C'[1,`j']
			local hi = `C'[2,`j']
		}
		if `hasci' == 0 & `s' < . {
			local lo = `b' - 1.96*`s'
			local hi = `b' + 1.96*`s'
		}
		local e1 "."
		if `j' == 1 local e1 "`ev'"
		savest "`tbl'" "`est'" "`eff'" `b' `s' `lo' `hi' `e1' `reps' $anaN
	}
end


*===============================================================================
* savemedeff: write the effects returned by medeff
*===============================================================================
* medeff is r-class, so every value is read in the five commands below before any
* other command can clear r(). The suffix of the bootstrap bounds differs between
* medeff versions, so both spellings are read and the one that is not missing is
* used; if lci and uci still come out empty in the csv, neither spelling applied
* and the bounds have to be taken from the medeff output by hand.

capture program drop savemedeff
program define savemedeff
	args tbl est
	local vpt "`=r(delta0)'   `=r(delta1)'   `=r(zeta0)'   `=r(zeta1)'   `=r(tau)'"
	local vlo "`=r(delta0_l)' `=r(delta1_l)' `=r(zeta0_l)' `=r(zeta1_l)' `=r(tau_l)'"
	local vhi "`=r(delta0_u)' `=r(delta1_u)' `=r(zeta0_u)' `=r(zeta1_u)' `=r(tau_u)'"
	local wlo "`=r(delta0lo)' `=r(delta1lo)' `=r(zeta0lo)' `=r(zeta1lo)' `=r(taulo)'"
	local whi "`=r(delta0hi)' `=r(delta1hi)' `=r(zeta0hi)' `=r(zeta1hi)' `=r(tauhi)'"

	local names "ACME_control ACME_treated ADE_control ADE_treated Total"

	forvalues j = 1/5 {
		local B`j' : word `j' of `vpt'
		local L`j' : word `j' of `vlo'
		local U`j' : word `j' of `vhi'
		if `L`j'' >= . local L`j' : word `j' of `wlo'
		if `U`j'' >= . local U`j' : word `j' of `whi'
		local eff : word `j' of `names'
		savest "`tbl'" "`est'" "`eff'" `B`j'' . `L`j'' `U`j'' . $anaR $anaN
	}

	* Average indirect and average direct effect. Without a treatment-mediator
	* interaction in the model the two ACMEs are the same number and so are their
	* intervals, so averaging them returns exactly that effect and that interval.
	* The standard error is then backed out of the interval inside savest.
	savest "`tbl'" "`est'" "AIE" `=(`B1'+`B2')/2' . `=(`L1'+`L2')/2' `=(`U1'+`U2')/2' . $anaR $anaN
	savest "`tbl'" "`est'" "ADE" `=(`B3'+`B4')/2' . `=(`L3'+`L4')/2' `=(`U3'+`U4')/2' . $anaR $anaN
end


*===============================================================================
* savegf: write the effects left behind by gformula
*===============================================================================
* gformula is rclass and runs ereturn clear before it finishes, so nothing is left
* in e(). What it does leave, in the ordinary matrix namespace, are the matrices it
* built from its own bootstrap: b (point estimates), se (bootstrap standard errors)
* and ci_percentile (the percentile interval it prints). Their columns come in the
* order gformula prints them:
*
*   1 Total  total causal effect
*   2 NDE    natural direct effect    (the average direct effect when there is no interaction)
*   3 NIE    natural indirect effect  (the average indirect effect when there is no interaction)
*   4 PM     proportion mediated
*   5 CDE    controlled direct effect (only present when control() is specified)
*
* Those matrices are not temporary, so they are dropped before each gformula call:
* otherwise a failed run would leave the previous one's numbers to be saved as new.

capture program drop savegf
program define savegf
	args tbl est
	capture confirm matrix b
	if _rc {
		noi di as error "savegf: matrix b not found - gformula saved nothing for Table `tbl'"
		exit
	}
	local names "Total NDE NIE PM CDE"
	local nb = colsof(b)
	if `nb' > 5 local nb = 5
	forvalues j = 1/`nb' {
		local eff : word `j' of `names'
		local point = b[1,`j']
		local s = se[1,`j']
		local lo = .
		local hi = .
		capture local lo = ci_percentile[1,`j']
		capture local hi = ci_percentile[2,`j']
		savest "`tbl'" "`est'" "`eff'" `point' `s' `lo' `hi' . $anaR $anaN
	}
end


*===============================================================================
* savete: write the ATE after teffects
*===============================================================================

capture program drop savete
program define savete
	args tbl est eff
	tempname B V
	matrix `B' = e(b)
	matrix `V' = e(V)
	local b = `B'[1,1]
	local s = sqrt(`V'[1,1])
	savest "`tbl'" "`est'" "`eff'" `b' `s' `=`b'-1.96*`s'' `=`b'+1.96*`s'' . . $anaN
end


*===============================================================================
* savesheet: write one sheet of the workbook
*===============================================================================
* Called at the end of a run, with the table that run feeds. It rereads the two
* csv files, keeps that table's rows and rewrites the sheet, so the workbook is
* up to date after every single estimator: run 5, open Excel, sheet S6 is there.
*
* When the same estimator is run twice (say with more replications the second
* time) both sets of rows sit in the csv file, which stays an append-only log of
* everything that was ever run. Only the last one reaches the sheet.
*
* Each sheet carries the estimates from the first row and, two rows below them,
* the diagnostics. The sheet is rewritten whole, so no stale rows are left behind.
* Nothing typed into the workbook by hand survives: the R estimates go into the
* manuscript tables, not into this file. Close the workbook in Excel before running.

capture program drop savesheet
program define savesheet
	args tbl
	preserve

	* Table S9 is produced by run 7 and has a dataset of its own
	if "`tbl'" == "S9" {
		capture confirm file "${tables}Supplementary_Table_S9.dta"
		if _rc == 0 {
			use "${tables}Supplementary_Table_S9.dta", clear
			local f "%12.${xldec}f"
			gen str40 ATE_95CI = strtrim(string(ATE, "`f'")) + " (" + strtrim(string(ATE_LCL, "`f'")) + "; " + strtrim(string(ATE_UCL, "`f'")) + ")"
			gen str40 ADE_95CI = strtrim(string(ADE, "`f'")) + " (" + strtrim(string(ADE_LCL, "`f'")) + "; " + strtrim(string(ADE_UCL, "`f'")) + ")"
			gen str40 AIE_95CI = strtrim(string(AIE, "`f'")) + " (" + strtrim(string(AIE_LCL, "`f'")) + "; " + strtrim(string(AIE_UCL, "`f'")) + ")"
			order subgroup level N REPS ATE_95CI ADE_95CI AIE_95CI
			format ATE ATE_LCL ATE_UCL ADE ADE_LCL ADE_UCL AIE AIE_LCL AIE_UCL %12.${xldec}f
			gen str8 P_TXT = strtrim(string(P_INTER, "%12.${xldecd}f")) // as text: Excel keeps the decimals
			drop P_INTER
			rename P_TXT P_INTER
			format level N REPS %12.0f
			capture export excel using "$xlsfile", sheet("S9", replace) firstrow(variables)
			if _rc noi di as error "could not write $xlsfile - is it open in Excel? The estimates are safe in the csv files: close it and type savetabs"
			qui count
			noi di as text "sheet S9: " as result r(N) as text " subgroups"
		}
	}

	* Table S9.1 is produced by run 7 (two mediators within strata) and has a dataset of its own
	else if "`tbl'" == "S9_1" {
		capture confirm file "${tables}Supplementary_Table_S9_1.dta"
		if _rc == 0 {
			use "${tables}Supplementary_Table_S9_1.dta", clear
			local f "%12.${xldec}f"
			gen str40 TOTAL_95CI   = strtrim(string(TOTAL, "`f'"))   + " (" + strtrim(string(TOTAL_LCL, "`f'"))   + "; " + strtrim(string(TOTAL_UCL, "`f'"))   + ")"
			gen str40 DIRECT_95CI  = strtrim(string(DIRECT, "`f'"))  + " (" + strtrim(string(DIRECT_LCL, "`f'"))  + "; " + strtrim(string(DIRECT_UCL, "`f'"))  + ")"
			gen str40 PSE_CRP_95CI = strtrim(string(PSE_CRP, "`f'")) + " (" + strtrim(string(PSE_CRP_LCL, "`f'")) + "; " + strtrim(string(PSE_CRP_UCL, "`f'")) + ")"
			gen str40 PSE_PRO_95CI = strtrim(string(PSE_PRO, "`f'")) + " (" + strtrim(string(PSE_PRO_LCL, "`f'")) + "; " + strtrim(string(PSE_PRO_UCL, "`f'")) + ")"
			gen str40 ATE_CRP_95CI = strtrim(string(ATE_CRP, "`f'")) + " (" + strtrim(string(ATE_CRP_LCL, "`f'")) + "; " + strtrim(string(ATE_CRP_UCL, "`f'")) + ")"
			gen str40 ATE_PRO_95CI = strtrim(string(ATE_PRO, "`f'")) + " (" + strtrim(string(ATE_PRO_LCL, "`f'")) + "; " + strtrim(string(ATE_PRO_UCL, "`f'")) + ")"
			order subgroup level N REPS TOTAL_95CI DIRECT_95CI PSE_CRP_95CI PSE_PRO_95CI ATE_CRP_95CI ATE_PRO_95CI
			format TOTAL TOTAL_LCL TOTAL_UCL DIRECT DIRECT_LCL DIRECT_UCL PSE_CRP PSE_CRP_LCL PSE_CRP_UCL PSE_PRO PSE_PRO_LCL PSE_PRO_UCL ATE_CRP ATE_CRP_LCL ATE_CRP_UCL ATE_PRO ATE_PRO_LCL ATE_PRO_UCL %12.${xldec}f
			format level N REPS %12.0f
			capture export excel using "$xlsfile", sheet("S9_1", replace) firstrow(variables)
			if _rc noi di as error "could not write $xlsfile - is it open in Excel? Close it and type savetabs"
			qui count
			noi di as text "sheet S9_1: " as result r(N) as text " subgroups"
		}
	}

	* Tables S2 to S7 and S10 are built from the two csv files
	else {

		local nest = 0
		local ndiag = 0

		* Estimates, from the first row of the sheet
		capture confirm file "$estfile"
		if _rc == 0 {
			import delimited "$estfile", clear varnames(1) case(preserve) stringcols(1 2 3)
			gen long ord = _n
			bysort table estimator effect (ord): gen byte lastrow = (_n == _N) // the last run of an estimator wins
			keep if lastrow
			sort ord
			drop ord lastrow
			keep if table == "`tbl'"
			qui count
			local nest = r(N)
			if `nest' > 0 {
				drop table
				local f "%12.${xldec}f"
				gen str40 est_95CI = strtrim(string(estimate, "`f'"))
				replace est_95CI = est_95CI + " (" + strtrim(string(lci, "`f'")) + "; " + strtrim(string(uci, "`f'")) + ")" if lci < . & uci < .
				order estimator effect est_95CI se estimate lci uci evalue reps n
				format estimate se lci uci evalue %12.${xldec}f
				format reps n %12.0f
				capture export excel using "$xlsfile", sheet("`tbl'", replace) firstrow(variables)
				if _rc noi di as error "could not write $xlsfile - is it open in Excel? The estimates are safe in the csv files: close it and type savetabs"
			}
		}

		* Diagnostics, two rows below the estimates
		capture confirm file "$diagfile"
		if _rc == 0 {
			import delimited "$diagfile", clear varnames(1) case(preserve) stringcols(1 2 3)
			gen long ord = _n
			bysort table estimator diagnostic (ord): gen byte lastrow = (_n == _N) // the last run of an estimator wins
			keep if lastrow
			sort ord
			drop ord lastrow
			keep if table == "`tbl'"
			qui count
			local ndiag = r(N)
			if `ndiag' > 0 {
				drop table
				* The diagnostics sit below the estimates, so they are written with
				* sheetmodify, and sheetmodify does not carry a Stata display format into
				* the Excel cell format: the column would come out with no decimals at all.
				* The value is therefore written as text, already rounded to $xldecd.
				gen str20 v = strtrim(string(value, "%12.${xldecd}f"))
				drop value
				rename v value
				local start = `nest' + 3
				if `nest' == 0 {
					capture export excel using "$xlsfile", sheet("`tbl'", replace) firstrow(variables)
					if _rc noi di as error "could not write $xlsfile - is it open in Excel? The estimates are safe in the csv files: close it and type savetabs"
				}
				else capture export excel using "$xlsfile", sheet("`tbl'") sheetmodify cell(A`start') firstrow(variables)
				if _rc noi di as error "could not write $xlsfile - is it open in Excel? The estimates are safe in the csv files: close it and type savetabs"
			}
		}

		noi di as text "sheet `tbl': " as result `nest' as text " estimates, " as result `ndiag' as text " diagnostics"
	}

	restore
end


*===============================================================================
* savetabs: rewrite every sheet of the workbook
*===============================================================================
* Not called by Analysis.do: each run already refreshes its own sheet. Type it in
* the Command window when the whole workbook has to be rebuilt from the csv files
* without re-estimating anything, for instance after deleting it or after changing
* the number of decimals in $xldec.

capture program drop savetabs
program define savetabs
	foreach t in S2 S3 S4 S5 S6 S7 S9 S9_1 S10 T2 {
		savesheet "`t'"
	}
	noi di as text "workbook written: $xlsfile"
end
