qui { 

******************************************************************************************************************************************************
******************************************************************************************************************************************************
**************************************** CAUSAL MEDIATION ANALYSIS IN METEOR (bDMARDs -> ASDAS -> BASFI at 6 months) *********************************
******************************************************************************************************************************************************
******************************************************************************************************************************************************

cls              // Cleans result window
clear all        // Drops loaded data
macro drop _all  // Drops all macros

*===============================================================================
* Working folders
*===============================================================================

global path     "C:\Users\alexa\OneDrive\Work\Projects\Causal_axSpA\METEOR\Data\Mediation manuscript\" // code
global datasets "${path}Datasets\"  // the source csv file
global programs "${path}Programs\"  // supporting programs
global tables   "${path}Tables\"    // tables and exported results
global figures  "${path}Figures\"   // figures

*===============================================================================
* Select seed and replications for the G-formula estimators
*===============================================================================

local seed=1234 // seed for bootstraping (one for all estimators)
local sims=1000  // number of bootstraping simulations
local mcsims=10000 // Monte Carlo sample size inside gformula (runs 1.2 and 4.1)

* Note: large mcsims increases computational time substantially

*===============================================================================
* Select estimator / analysis
*===============================================================================

local run = 8.4 // Chose estimator/analysis

* 0    = Descriptive: Treated vs untreated (Table 1) and included vs excluded (Table S1)
* 1    = Box S1   Total effect, time-fixed g-formula          (manual)          -> Table S2
* 1.1  = Box S1   Total effect                                (Stata medeff)    -> Table S2
* 1.2  = Box S1   Total effect and mediation effects          (Stata gformula)  -> Tables S2 and S3
* 2    = Box S2   Mediation ignoring MOC                      (manual)          -> Table S3
* 2.1  = Box S2   Mediation ignoring MOC                      (Stata medeff)    -> Table S3
* 3    = Box S2.1 Two mediators, main paths, no interaction   (manual)          -> Table S4
* 3.1  = Box S2.1 Two mediators, all paths, no interaction    (manual)          -> Tables S4 and S5
* 3.2  = Box S2.1 Two mediators, all paths, with interaction  (manual)          -> Tables S4 and S5
* 4    = Box S6   Mediation accounting for MOC                (manual)          -> Table S10
* 4.1  = Box S6   Mediation accounting for MOC                (Stata gformula)  -> Table S10
* 5    = Box S3   MSM with IPTW                               (manual)          -> Tables S6 and S8
* 5.1  = Box S3   IPTW                                        (Stata teffects)  -> Table S6
* 6    = Box S4   TMLE, total effect                          (manual)          -> Table S7
* 6.1  = Box S4   AIPW, total effect                          (Stata teffects)  -> Table S7
* 7    = Mediation within levels of baseline characteristics  (manual)          -> Tables S9 and S9.1
* 8    = Complete vs missing 6-month data                                       -> Table S11 and sheet T2
* 8.1  = Sensitivity: >=1M of bDMARD, one & two mediators     (manual)          -> Table 2 (sheet T2)
* 8.2  = Sensitivity: >=3M of bDMARD, one & two mediators     (manual)          -> Table 2 (sheet T2)
* 8.3  = Box S7 Sensitivity: total effect, MSM (sIPTW x CW)   (manual)          -> Figure 3B (sheet T2)
* 8.4  = Box S7 Sensitivity: total effect,TMLE with censoring (manual)          -> Figure 3B (sheet T2)

* MOC = mediator-outcome confounders; PSE = path-specific effect
* Table S8 (balance) comes from run 5 via msmbalanceV1


*===============================================================================
* Collect the Stata estimates for the supplementary tables
*===============================================================================
* Each estimate is appended to ${tables}estimates_stata.csv and each diagnostic to
* ${tables}diagnostics_stata.csv, and each run then rewrites its own sheet of
* ${tables}estimates_stata.xlsx, one sheet per supplementary table. So the Excel file
* is up to date after every single estimator, and running an estimator again
* replaces its rows in the sheet. 
* Excel must be closed before running.

qui do "${programs}saveestV1.do" // Run supporting program


*===============================================================================
* Select the cohort
*===============================================================================

local cohort = 800 // Chose cohort

* 800  = baseline and 6 months,        n = 419 (paper)
* 761  = baseline, 6 and 12 months,    n = 352 (poster)


******************************************************************************************************************************************************
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
**** Prepare the dataset
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
******************************************************************************************************************************************************

*===============================================================================
* Dataset:   ${datasets}originalfulllong800.csv  (cohort 800)
*            ${datasets}originalfulllong761.csv  (cohort 761)
*
* Cohort:    population3 == 1 & allcompleters == 1
*            cohort 800: n = 419 over two visits; cohort 761: n = 352 over three.
*            High disease activity at baseline is already applied in both.
*            complete ASDAS and BASFI at every visit and complete baseline
*            covariates: the main analysis population of the manuscript.
*
* Treatment: bionew, taken as it is in the source file. It records that a
*            bDMARD was used in the interval ending at the visit, so on the
*            6-month row it is the bDMARD covering baseline to 6 months. 
*
* Visits:    t1 (baseline), t2 (6 months); cohort 761 also has t3 (12 months)
*===============================================================================

if `cohort' == 800 {
	import delimited "${datasets}originalfulllong800.csv", clear varnames(1) case(preserve)
	clonevar population3 = population2 // same role as population3, over two visits
	local n = 419
}

if `cohort' == 761 {
	import delimited "${datasets}originalfulllong761.csv", clear varnames(1) case(preserve)
	local n = 352
}

*-------------------------------------------------------------------------------
* Analysis cohort, one row per patient
*-------------------------------------------------------------------------------

keep if population3 == 1 & allcompleters == 1 // main analysis population
keep if t == 6 // one row per patient

qui count
assert r(N) == `n' // stops here if the source file or the selection changes
noi di as text "patients in the analysis: " as result r(N)
global anaN = r(N)   // patients, used by the estimate collector
global anaR = `sims' // replications, used by the estimate collector

*-------------------------------------------------------------------------------
* ASDAS without CRP, used by the two-mediator analyses (run 3, 3.1 and 3.2)
*-------------------------------------------------------------------------------

capture drop asdastotalt1_pro
gen asdastotalt1_pro = asdastotalt1 - 0.579*ln(crpt1 + 1)
capture drop asdastotalt2_pro
gen asdastotalt2_pro = asdastotalt2 - 0.579*ln(crpt2 + 1)

*-------------------------------------------------------------------------------
* Variables
*-------------------------------------------------------------------------------

**** Pre-treatment baseline confounders 
global W="age sex comorbbin mny asasmri hla pertvt1 ibdbl emmtvt1 comedtvt1 asdastotalt1 basfitotalt1"

**** Treatment (binary)
global A="bionew"

**** Mediator (continuous)
global M="asdastotalt2"

**** Mediator-outcome confounders (MOC) (all binary)
global MOCa ="pertvt2"
global MOCb ="emmtvt2"
global MOCc ="comedtvt2"

**** Outcome (continuous)
global Y="basfitotalt2"

*-------------------------------------------------------------------------------
* Keep the analysis variables only
*-------------------------------------------------------------------------------

keep id t $W $A $M $Y $MOCa $MOCb $MOCc ///
	crpt1 crpt2 crpelevatedt1 basdaitotalt1 basdaitotalt2 asdastotalt1_pro asdastotalt2_pro ///
	asdascatt1 asdascatt2 ///
	population3 allcompleters bionew1m bionew3m /// bionew1m bionew3m: treatment definitions of runs 8.1 and 8.2
	symptomsduration smokeever bmi ///
	arthtvt1 enthtvt1 psotvt1 aautvt1 nsaidnewt1 csdmardnewt1 gcnewt1

order id $W $A $M $Y 


******************************************************************************************************************************************************
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
**** 0 - DESCRIPTIVE: the analysis cohort by treatment (Table 1)
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
******************************************************************************************************************************************************

if `run' == 0 {

*-------------------------------------------------------------------------------
* Table 1 
*-------------------------------------------------------------------------------
qui do "${programs}smdtableV2.do" // Run supporting program

global descriptive="age sex symptomsduration comorbbin mny asasmri hla crpelevatedt1 arthtvt1 enthtvt1 psotvt1 ibdbl aautvt1 asdastotalt1 i.asdascatt1 basfitotalt1 nsaidnewt1 csdmardnewt1 gcnewt1" 
global groupvars "population3==1 & allcompleters==1"  
noi smdtableV2, selector($groupvars) visit(6) analysis(table1) variable($descriptive) treatvar(bionew) roundvar(0) roundsmd(2) difference(3) completers(yes) save(yes) ///
    folder($tables) filename(Table1)

*-------------------------------------------------------------------------------
* Supplementary Table S1: included vs excluded among the 800 eligible patients 
*-------------------------------------------------------------------------------
* The data in memory are the 419 analysed patients: they are set aside, the
* 800 eligible patients are read again and the analysis data are restored at the end. 

if `cohort' == 800 {

tempfile anadata0
save `anadata0'

import delimited "${datasets}originalfulllong800.csv", clear varnames(1) case(preserve)
clonevar population3 = population2
keep if t == 6 // one row per patient; baseline values are in the t1 variables

qui count
assert r(N) == 800 // all eligible patients

gen byte included = 0
replace included = 1 if population3 == 1 & allcompleters == 1

qui count if included == 1
assert r(N) == 419 // the analysis cohort

capture label drop includedlbl
label define includedlbl 1 "Included" 0 "Excluded"
label values included includedlbl

global descriptive="bionew age sex symptomsduration comorbbin mny asasmri hla crpelevatedt1 arthtvt1 enthtvt1 psotvt1 ibdbl aautvt1 asdastotalt1 i.asdascatt1 basfitotalt1 nsaidnewt1 csdmardnewt1 gcnewt1"
global groupvars "included==0 | included==1"
noi smdtableV2, selector($groupvars) visit(6) analysis(table1) variable($descriptive) treatvar(included) roundvar(0) roundsmd(2) difference(0.1) save(yes) ///
    folder($tables) filename(Supplementary_Table_S1)

use `anadata0', clear // back to the 419 analysed patients

} // close cohort 800


} // close 0


******************************************************************************************************************************************************
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
**** 1 - SUPPLEMENTARY BOX S1: total effect of bDMARDs on BASFI at 6 months, time-fixed g-formula (manual)
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
******************************************************************************************************************************************************

if `run' == 1 {

*===============================================================================
* REFERENCE:    Supplementary Box S1 (Total Effect of bDMARDs on BASFI at 6M)
* METHOD:       Parametric Time-Fixed G-Formula
* ESTIMAND:     Average Total Effect (ATE) of bDMARDs on BASFI at 6 Months
*===============================================================================

*===============================================================================
* PROGRAM:       Manual calculation
*===============================================================================


capture program drop run_gformula
program define run_gformula, rclass // Program for confidence interval (ignore if only interested in point estimate)

///////////////// Step 1 — Model the Observed Data

**** Outcome model (mandatory: include mediator only if later simulated)
regress $Y $M $A $W
** Save coeficients from the outcome model in a vector
estimates store outcome

**** Mediator model (optional: simulation can run without simulating the mediator, mandatory only if mediator included in outcome model)
regress $M $A $W
** Save coeficients from the Mediator model in a vector
estimates store mediator

**** Treatment model (optional: only for diagnostics)
logit $A $W
** Save coeficients from the Mediator model in a vector
estimates store treatment

///////////////// Step 2 — Monte Carlo Simulation (adjust for confounding)

/////////// Step 2.1. Simulate counterfactual mediator (ASDAS) at 6 months (Optional, only if included in the outcome model)

* Restore the mediator model coefficients
estimates restore mediator

* Generate M1 (Counterfactual ASDAS at 6 months if everyone were treated, bionew = 1)
gen M1 = _b[_cons] + ///
         _b[bionew]*1 + ///
         _b[age]*age + ///
         _b[sex]*sex + ///
         _b[comorbbin]*comorbbin + ///
         _b[mny]*mny + ///
         _b[asasmri]*asasmri + ///
         _b[hla]*hla + ///
         _b[pertvt1]*pertvt1 + ///
         _b[ibdbl]*ibdbl + ///
         _b[emmtvt1]*emmtvt1 + ///
         _b[comedtvt1]*comedtvt1 + ///
         _b[asdastotalt1]*asdastotalt1 + ///
         _b[basfitotalt1]*basfitotalt1

* Generate M0 (Counterfactual ASDAS at 6 months if everyone were untreated, bionew = 0)
gen M0 = _b[_cons] + ///
         _b[bionew]*0 + ///
         _b[age]*age + ///
         _b[sex]*sex + ///
         _b[comorbbin]*comorbbin + ///
         _b[mny]*mny + ///
         _b[asasmri]*asasmri + ///
         _b[hla]*hla + ///
         _b[pertvt1]*pertvt1 + ///
         _b[ibdbl]*ibdbl + ///
         _b[emmtvt1]*emmtvt1 + ///
         _b[comedtvt1]*comedtvt1 + ///
         _b[asdastotalt1]*asdastotalt1 + ///
         _b[basfitotalt1]*basfitotalt1
		 
/////////// Step 2.2. Simulate counterfactual outcome (BASFI) at 6 months

* Restore the outcome model coefficients
estimates restore outcome

* Generate Y1M1 (Counterfactual BASFI if everyone were treated, with ASDAS at its treated value M1)
gen Y1M1 = _b[_cons] + ///
           _b[asdastotalt2]*M1 + ///
           _b[bionew]*1 + ///
           _b[age]*age + ///
           _b[sex]*sex + ///
           _b[comorbbin]*comorbbin + ///
           _b[mny]*mny + ///
           _b[asasmri]*asasmri + ///
           _b[hla]*hla + ///
           _b[pertvt1]*pertvt1 + ///
           _b[ibdbl]*ibdbl + ///
           _b[emmtvt1]*emmtvt1 + ///
           _b[comedtvt1]*comedtvt1 + ///
           _b[asdastotalt1]*asdastotalt1 + ///
           _b[basfitotalt1]*basfitotalt1

* Generate Y0M0 (Counterfactual BASFI if everyone were untreated, with ASDAS at its untreated value M0)
gen Y0M0 = _b[_cons] + ///
           _b[asdastotalt2]*M0 + ///
           _b[bionew]*0 + ///
           _b[age]*age + ///
           _b[sex]*sex + ///
           _b[comorbbin]*comorbbin + ///
           _b[mny]*mny + ///
           _b[asasmri]*asasmri + ///
           _b[hla]*hla + ///
           _b[pertvt1]*pertvt1 + ///
           _b[ibdbl]*ibdbl + ///
           _b[emmtvt1]*emmtvt1 + ///
           _b[comedtvt1]*comedtvt1 + ///
           _b[asdastotalt1]*asdastotalt1 + ///
           _b[basfitotalt1]*basfitotalt1

/////////// Step 3 — Total Treatment Effect

* Calculate mean PO accross all patients and save in scalars for bootstraping

qui sum Y1M1
scalar mean_Y1M1 = r(mean)
    
qui sum Y0M0
scalar mean_Y0M0 = r(mean)
    
* Calculate the marginal total treatment effect (ATE) and save in scalar for bootstraping

scalar total_effect = mean_Y1M1 - mean_Y0M0
    
/////////// Return values for the bootstraping
return scalar ate = total_effect
return scalar mean11 = mean_Y1M1
return scalar mean00 = mean_Y0M0


/////////// Optional step. Predict treatment, mediator and outcome under observed data (only for diagnostics)

* Observed treatment
qui sum $A, meanonly 
return scalar mean_$A = r(mean)

* Generate A (predicted treatment under observed data)
estimates restore treatment

gen Ap = _b[_cons] + ///
         _b[age]*age + ///
         _b[sex]*sex + ///
         _b[comorbbin]*comorbbin + ///
         _b[mny]*mny + ///
         _b[asasmri]*asasmri + ///
         _b[hla]*hla + ///
         _b[pertvt1]*pertvt1 + ///
         _b[ibdbl]*ibdbl + ///
         _b[emmtvt1]*emmtvt1 + ///
         _b[comedtvt1]*comedtvt1 + ///
         _b[asdastotalt1]*asdastotalt1 + ///
         _b[basfitotalt1]*basfitotalt1
replace Ap= invlogit(Ap)
replace Ap= rbinomial(1, Ap) // random error added

qui sum Ap, meanonly 
return scalar mean_Ap = r(mean)

* Observed mediator
qui sum $M, meanonly 
return scalar mean_$M = r(mean)

* Restore the mediator model coefficients
estimates restore mediator

* Generate M (Predicted ASDAS at 6 months under observed treatment)
gen Mp = _b[_cons] + ///
         _b[bionew]*Ap + ///
         _b[age]*age + ///
         _b[sex]*sex + ///
         _b[comorbbin]*comorbbin + ///
         _b[mny]*mny + ///
         _b[asasmri]*asasmri + ///
         _b[hla]*hla + ///
         _b[pertvt1]*pertvt1 + ///
         _b[ibdbl]*ibdbl + ///
         _b[emmtvt1]*emmtvt1 + ///
         _b[comedtvt1]*comedtvt1 + ///
         _b[asdastotalt1]*asdastotalt1 + ///
         _b[basfitotalt1]*basfitotalt1

qui sum Mp, meanonly 
return scalar mean_Mp = r(mean)

* Observed outcome
qui sum $Y, meanonly 
return scalar mean_$Y = r(mean)

* Restore the outcome model coeficients
estimates restore outcome

gen Yp   = _b[_cons] + ///
           _b[asdastotalt2]*Mp + ///
           _b[bionew]*Ap + ///
           _b[age]*age + ///
           _b[sex]*sex + ///
           _b[comorbbin]*comorbbin + ///
           _b[mny]*mny + ///
           _b[asasmri]*asasmri + ///
           _b[hla]*hla + ///
           _b[pertvt1]*pertvt1 + ///
           _b[ibdbl]*ibdbl + ///
           _b[emmtvt1]*emmtvt1 + ///
           _b[comedtvt1]*comedtvt1 + ///
           _b[asdastotalt1]*asdastotalt1 + ///
           _b[basfitotalt1]*basfitotalt1

qui sum Yp, meanonly 
return scalar mean_Yp = r(mean)



/////////// Clean up temporary variables before the next replication

drop M1 M0 Y1M1 Y0M0 Ap Mp Yp _est_mediator _est_outcome _est_treatment

end 

* Run the bootstrap across 1000 replications (second line of bootstrap comman contains only diagnostics and is optional
set seed `seed'
bootstrap ATE=r(ate) Mean_Y1M1=r(mean11) Mean_Y0M0=r(mean00) ///
          Mean_$A=r(mean_$A) Mean_Ap=r(mean_Ap) Mean_$M=r(mean_$M) Mean_Mp=r(mean_Mp) Mean_$Y=r(mean_$Y) Mean_Yp=r(mean_Yp) ///
		  , reps(`sims') nodots: run_gformula

* Display the percentile-based bootstrap 95% Confidence Intervals
noi estat bootstrap, percentile


///////////////// e-value

* Pull ATE point estimate from the bootstrap results
matrix b = e(b)
scalar ATE_gform = b[1,1]

* SD of the observed outcome
qui sum $Y
scalar SD_Y = r(sd)

* Compute the e-value
local rr = exp(0.91*abs(scalar(ATE_gform)/scalar(SD_Y)))
scalar evalue_point = `rr' + sqrt(`rr'*(`rr'-1))
noi di "E-value (point estimate): " %6.3f scalar(evalue_point)

savebs "S2" "G-formula (manual)" `=scalar(evalue_point)'
savesheet "S2"



} // close 1


******************************************************************************************************************************************************
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
**** 1.1 - SUPPLEMENTARY BOX S1: total effect (Stata medeff)
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
******************************************************************************************************************************************************

if `run' == 1.1 {

*===============================================================================
* PROGRAM:       Stata medeff
*===============================================================================

set seed `seed' // medeff bootstraps
noi medeff (regress $M $A $W) (regress $Y $A $M $W), mediate($M) treat($A) vce(bootstrap, reps(`sims'))

savemedeff "S2" "medeff"
savesheet "S2"


} // close 1.1


******************************************************************************************************************************************************
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
**** 1.2 - SUPPLEMENTARY BOX S1: total effect and mediation effects (Stata gformula, mediation syntax)
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
******************************************************************************************************************************************************

if `run' == 1.2 {

*===============================================================================
* PROGRAM:       Stata gformula (mediation syntax)
*===============================================================================


* Stata gformula simulates the mediator, the MOC and the outcome instead of using their
* predicted means, so its point estimates carry Monte Carlo error on top of the
* sampling error. Two options bring it close to the manual code:
*   minsim:                uses the predicted mean of the outcome instead of drawing it
*   moreMC simulations():  lets the Monte Carlo sample exceed the analysis sample, so
*                          the remaining noise from the mediator and MOC draws shrinks
*                          with the square root of the Monte Carlo sample size


capture matrix drop b se ci_normal ci_percentile ci_bc ci_bca // essential to avoid old values being used

noi gformula $Y $M $A $W, ///
mediation ex($A) mediator($M) out($Y) ///
eq($M: $A $W, $Y: $A $M $W) /// 
com($M:regress, $Y:regress)   ///
obe base_confs($W) ///
minsim moreMC simulations(`mcsims') ///
seed(`seed') samples(`sims') all

* The same call gives the total effect (Table S2) and the natural direct and
* indirect effects (Table S3): the rows are written once under each table.
savegf "S2" "gformula"
savegf "S3" "gformula"
savesheet "S2"
savesheet "S3"



} // close 1.2


******************************************************************************************************************************************************
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
**** 2 - SUPPLEMENTARY BOX S2: mediation ignoring mediator-outcome confounding (manual)
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
******************************************************************************************************************************************************

if `run' == 2 {

*===============================================================================
* REFERENCE:    Supplementary Box S2 (ignoring  mediator-outcome confounding)
* METHOD:       Parametric Time-Fixed G-Formula Algorithm for causal mediation
* ESTIMAND:     Causal mediation effects of bDMARDs on BASFI at 6 Months
*===============================================================================


*===============================================================================
* METHOD:       Manual calculation
*===============================================================================

///////////////// Step 0 - Test for a treatment-mediator interaction

* The A x M product term is fitted once on the observed data. Its p-value decides
* whether the interaction is carried into the outcome models below. 

local pcut = 0.20 // significance level for keeping the interaction
gen inter = $M*$A // built once here, not rebuilt at every bootstrap replication
regress $Y $M $A $W inter

qui test inter
local pinter = r(p)

local interaction = 0
if `pinter' < `pcut' {
	local interaction = 1
}

global interX   ""
global interY11 ""
global interY10 ""
global interY01 ""
global interY00 ""

if `interaction' == 1 {
	global interX   "inter"
	global interY11 "+ _b[inter]*(M1*1)"
	global interY10 "+ _b[inter]*(M0*1)"
	global interY01 "+ _b[inter]*(M1*0)"
	global interY00 "+ _b[inter]*(M0*0)"
}

noi di ""
noi di as text "Treatment-mediator interaction ($A x $M): p = " as result %6.4f `pinter'
if `interaction' == 1 {
	noi di as text "  -> interaction kept in the outcome model and in the counterfactual outcomes"
}
if `interaction' == 0 {
	noi di as text "  -> no interaction (p >= `pcut'): the outcome models have no A x M term"
}
noi di ""

savediag "S3" "G-formula ignoring MOC (manual)" "p interaction AxM" `pinter'


capture program drop run_gformula
program define run_gformula, rclass // Program for confidence interval (ignore if only interested in point estimate)

///////////////// Step 1 — Model the Observed Data

**** Outcome model (mandatory: include mediator only if later simulated)
regress $Y $M $A $W $interX // A x M term, empty if Step 0 found no interaction
** Save coeficients from the outcome model in a vector
estimates store outcome

**** Mediator model (optional: simulation can run without simulating the mediator, mandatory only if mediator included in outcome model)
regress $M $A $W
** Save coeficients from the Mediator model in a vector
estimates store mediator

**** Treatment model (optional: only for diagnostics)
logit $A $W
** Save coeficients from the Mediator model in a vector
estimates store treatment

///////////////// Step 2 — Monte Carlo Simulation (adjust for confounding)

/////////// Step 2.1. Simulate counterfactual mediator (ASDAS) at 6 months (Optional, only if included in the outcome model)

* Restore the mediator model coefficients
estimates restore mediator

* Generate M1 (Counterfactual ASDAS at 6 months if everyone were treated, bionew = 1)
gen M1 = _b[_cons] + ///
         _b[bionew]*1 + ///
         _b[age]*age + ///
         _b[sex]*sex + ///
         _b[comorbbin]*comorbbin + ///
         _b[mny]*mny + ///
         _b[asasmri]*asasmri + ///
         _b[hla]*hla + ///
         _b[pertvt1]*pertvt1 + ///
         _b[ibdbl]*ibdbl + ///
         _b[emmtvt1]*emmtvt1 + ///
         _b[comedtvt1]*comedtvt1 + ///
         _b[asdastotalt1]*asdastotalt1 + ///
         _b[basfitotalt1]*basfitotalt1 

* Generate M0 (Counterfactual ASDAS at 6 months if everyone were untreated, bionew = 0)
gen M0 = _b[_cons] + ///
         _b[bionew]*0 + ///
         _b[age]*age + ///
         _b[sex]*sex + ///
         _b[comorbbin]*comorbbin + ///
         _b[mny]*mny + ///
         _b[asasmri]*asasmri + ///
         _b[hla]*hla + ///
         _b[pertvt1]*pertvt1 + ///
         _b[ibdbl]*ibdbl + ///
         _b[emmtvt1]*emmtvt1 + ///
         _b[comedtvt1]*comedtvt1 + ///
         _b[asdastotalt1]*asdastotalt1 + ///
         _b[basfitotalt1]*basfitotalt1 
		 
/////////// Step 2.2. Simulate counterfactual outcome (BASFI) at 6 months (differs from BOX S1 only from this point onward)

* Restore the outcome model coefficients
estimates restore outcome

* Generate Y1M1 (Counterfactual BASFI if everyone were treated, with ASDAS at its treated value M1)
gen Y1M1 = _b[_cons] + ///
           _b[asdastotalt2]*M1 + ///
           _b[bionew]*1 + ///
           _b[age]*age + ///
           _b[sex]*sex + ///
           _b[comorbbin]*comorbbin + ///
           _b[mny]*mny + ///
           _b[asasmri]*asasmri + ///
           _b[hla]*hla + ///
           _b[pertvt1]*pertvt1 + ///
           _b[ibdbl]*ibdbl + ///
           _b[emmtvt1]*emmtvt1 + ///
           _b[comedtvt1]*comedtvt1 + ///
           _b[asdastotalt1]*asdastotalt1 + ///
           _b[basfitotalt1]*basfitotalt1 $interY11 // A x M term, empty if Step 0 found no interaction

* Generate Y1M0 (Counterfactual BASFI if everyone were treated, with ASDAS at its untreated value M0)
gen Y1M0 = _b[_cons] + ///
           _b[asdastotalt2]*M0 + ///
           _b[bionew]*1 + ///
           _b[age]*age + ///
           _b[sex]*sex + ///
           _b[comorbbin]*comorbbin + ///
           _b[mny]*mny + ///
           _b[asasmri]*asasmri + ///
           _b[hla]*hla + ///
           _b[pertvt1]*pertvt1 + ///
           _b[ibdbl]*ibdbl + ///
           _b[emmtvt1]*emmtvt1 + ///
           _b[comedtvt1]*comedtvt1 + ///
           _b[asdastotalt1]*asdastotalt1 + ///
           _b[basfitotalt1]*basfitotalt1 $interY10 // A x M term, empty if Step 0 found no interaction
		
* Generate Y0M1 (Counterfactual BASFI if everyone were untreated, with ASDAS at its treated value M1)
gen Y0M1 = _b[_cons] + ///
           _b[asdastotalt2]*M1 + ///
           _b[bionew]*0 + ///
           _b[age]*age + ///
           _b[sex]*sex + ///
           _b[comorbbin]*comorbbin + ///
           _b[mny]*mny + ///
           _b[asasmri]*asasmri + ///
           _b[hla]*hla + ///
           _b[pertvt1]*pertvt1 + ///
           _b[ibdbl]*ibdbl + ///
           _b[emmtvt1]*emmtvt1 + ///
           _b[comedtvt1]*comedtvt1 + ///
           _b[asdastotalt1]*asdastotalt1 + ///
           _b[basfitotalt1]*basfitotalt1 $interY01 // A x M term, empty if Step 0 found no interaction
		   
* Generate Y0M0 (Counterfactual BASFI if everyone were untreated, with ASDAS at its untreated value M0)
gen Y0M0 = _b[_cons] + ///
           _b[asdastotalt2]*M0 + ///
           _b[bionew]*0 + ///
           _b[age]*age + ///
           _b[sex]*sex + ///
           _b[comorbbin]*comorbbin + ///
           _b[mny]*mny + ///
           _b[asasmri]*asasmri + ///
           _b[hla]*hla + ///
           _b[pertvt1]*pertvt1 + ///
           _b[ibdbl]*ibdbl + ///
           _b[emmtvt1]*emmtvt1 + ///
           _b[comedtvt1]*comedtvt1 + ///
           _b[asdastotalt1]*asdastotalt1 + ///
           _b[basfitotalt1]*basfitotalt1 $interY00 // A x M term, empty if Step 0 found no interaction

/////////// Step 3 — Mediation effects

* Calculate mean PO accross all patients and save in scalars for bootstraping

qui sum Y1M1
scalar mean_Y1M1 = r(mean)
 
qui sum Y1M0
scalar mean_Y1M0 = r(mean)

qui sum Y0M1
scalar mean_Y0M1 = r(mean)

qui sum Y0M0
scalar mean_Y0M0 = r(mean)


* Calculate mediation effects accross all patients and save in scalars for bootstraping

gen NIE = Y0M1-Y0M0
qui sum NIE
return scalar mean_NIE = r(mean)
    
gen TIE = Y1M1-Y1M0
qui sum TIE
return scalar mean_TIE = r(mean)

gen AIE = (NIE+TIE)/2
qui sum AIE
return scalar mean_AIE = r(mean)
	
gen NDE = Y1M0-Y0M0
qui sum NDE
return scalar mean_NDE = r(mean)	
 
gen TDE = Y1M1-Y0M1
qui sum TDE
return scalar mean_TDE = r(mean)	

gen ADE = (NDE+TDE)/2
qui sum ADE
return scalar mean_ADE = r(mean)


* Calculate the marginal total treatment effect (ATE) and save in scalar for bootstraping

scalar total_effect = mean_Y1M1 - mean_Y0M0
    
/////////// Return values for the bootstraping
return scalar ate = total_effect
return scalar mean11 = mean_Y1M1
return scalar mean10 = mean_Y1M0
return scalar mean01 = mean_Y0M1
return scalar mean00 = mean_Y0M0


/////////// Optional step. Predict treatment, mediator and outcome under observed data (only for diagnostics)

* Observed treatment
qui sum $A, meanonly 
return scalar mean_$A = r(mean)

* Generate A (predicted treatment under observed data)
estimates restore treatment

gen Ap = _b[_cons] + ///
         _b[age]*age + ///
         _b[sex]*sex + ///
         _b[comorbbin]*comorbbin + ///
         _b[mny]*mny + ///
         _b[asasmri]*asasmri + ///
         _b[hla]*hla + ///
         _b[pertvt1]*pertvt1 + ///
         _b[ibdbl]*ibdbl + ///
         _b[emmtvt1]*emmtvt1 + ///
         _b[comedtvt1]*comedtvt1 + ///
         _b[asdastotalt1]*asdastotalt1 + ///
         _b[basfitotalt1]*basfitotalt1
replace Ap= invlogit(Ap)
replace Ap= rbinomial(1, Ap) // random error added

qui sum Ap, meanonly 
return scalar mean_Ap = r(mean)

* Observed mediator
qui sum $M, meanonly 
return scalar mean_$M = r(mean)

* Restore the mediator model coefficients
estimates restore mediator

* Generate M (Predicted ASDAS at 6 months under observed treatment)
gen Mp = _b[_cons] + ///
         _b[bionew]*Ap + ///
         _b[age]*age + ///
         _b[sex]*sex + ///
         _b[comorbbin]*comorbbin + ///
         _b[mny]*mny + ///
         _b[asasmri]*asasmri + ///
         _b[hla]*hla + ///
         _b[pertvt1]*pertvt1 + ///
         _b[ibdbl]*ibdbl + ///
         _b[emmtvt1]*emmtvt1 + ///
         _b[comedtvt1]*comedtvt1 + ///
         _b[asdastotalt1]*asdastotalt1 + ///
         _b[basfitotalt1]*basfitotalt1

qui sum Mp, meanonly 
return scalar mean_Mp = r(mean)

* Observed outcome
qui sum $Y, meanonly 
return scalar mean_$Y = r(mean)


* Restore the outcome model coeficients
estimates restore outcome

gen Yp   = _b[_cons] + ///
           _b[asdastotalt2]*Mp + ///
           _b[bionew]*Ap + ///
           _b[age]*age + ///
           _b[sex]*sex + ///
           _b[comorbbin]*comorbbin + ///
           _b[mny]*mny + ///
           _b[asasmri]*asasmri + ///
           _b[hla]*hla + ///
           _b[pertvt1]*pertvt1 + ///
           _b[ibdbl]*ibdbl + ///
           _b[emmtvt1]*emmtvt1 + ///
           _b[comedtvt1]*comedtvt1 + ///
           _b[asdastotalt1]*asdastotalt1 + ///
           _b[basfitotalt1]*basfitotalt1 //+ _b[inter]*(Mp*Ap) // Delete last term if no interaction is intended.

qui sum Yp, meanonly 
return scalar mean_Yp = r(mean)



/////////// Clean up temporary variables before the next replication

capture drop M1 M0 Y1M1 Y1M0 Y0M1 Y0M0 NIE TIE AIE NDE TDE ADE Ap Mp Yp _est_mediator _est_outcome _est_treatment // inter is kept: it is built once in Step 0

end

* Run the bootstrap across 1000 replications (second line of bootstrap comman contains only diagnostics and is optional
set seed `seed'
bootstrap ATE=r(ate) Mean_Y1M1=r(mean11) Mean_Y1M0=r(mean10) Mean_Y0M1=r(mean01) Mean_Y0M0=r(mean00) ///
		  NIE=r(mean_NIE) TIE=r(mean_TIE) AIE=r(mean_AIE) NDE=r(mean_NDE) TDE=r(mean_TDE) ADE=r(mean_ADE) ///
          Mean_$A=r(mean_$A) Mean_Ap=r(mean_Ap) Mean_$M=r(mean_$M) Mean_Mp=r(mean_Mp) Mean_$Y=r(mean_$Y) Mean_Yp=r(mean_Yp) ///
		  , reps(`sims') nodots: run_gformula

* Display the percentile-based bootstrap 95% Confidence Intervals
noi estat bootstrap, percentile

savebs "S3" "G-formula ignoring MOC (manual)"
savesheet "S3"


} // close 2


******************************************************************************************************************************************************
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
**** 2.1 - SUPPLEMENTARY BOX S2: mediation ignoring mediator-outcome confounding (Stata medeff)
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
******************************************************************************************************************************************************

if `run' == 2.1 {

*===============================================================================
* PROGRAM:       Stata medeff
*===============================================================================
local pcut = 0.20
gen inter=$M*$A
regress $Y $M $A $W inter

qui test inter
local pinter = r(p)

local interX ""    // term in the outcome model
local interopt ""  // medeff option
if `pinter' < `pcut' {
	local interX "inter"
	local interopt "interact(inter)"
}

if "`interopt'" != "" {
	noi di as text "  -> interaction kept in the outcome modelin medeff"
}
if "`interopt'" == "" {
	noi di as text "  -> no interaction (p >= `pcut'): medeff runs without the A x M term"
}

set seed `seed'
noi medeff (regress $M $A $W) (regress $Y $A $M $W `interX'), mediate($M) treat($A) vce(bootstrap, reps(`sims')) `interopt'

savemedeff "S3" "medeff"
savesheet "S3"


*===============================================================================
* PROGRAM:       Stata gformula (mediation syntax) 
*===============================================================================

* Not supported (can only include post-treatment mediator-outcome confounders)




} // close 2.1


******************************************************************************************************************************************************
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
**** 3 - SUPPLEMENTARY BOX S2.1: two causally ordered mediators, paths type 1 decomposition (manual)
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
******************************************************************************************************************************************************

if `run' == 3 {

*===============================================================================
* REFERENCE:    Supplementary Box S2.1 (ignoring  mediator-outcome confounding)
* METHOD:       Parametric Time-Fixed G-Formula Algorithm for causal mediation
* ESTIMAND:     Causal mediation effects of bDMARDs on BASFI at 6 Months
*               with two causally ordered mediators
*===============================================================================


*===============================================================================
* METHOD:       Manual calculation of the R paths package (type 1 decomposition)
*===============================================================================

*===============================================================================
*
* METHOD:       Parametric edge g-formula for path-specific effects (PSEs)
*               with two causally ordered mediators: M1 -> M2 -> Y ("Paths")
*
* ESTIMAND:     Direct effect (A->Y), path-specific effect via M1 (A->M1~>Y,
*               which may also travel through M2), and path-specific effect
*               via M2 (A->M2->Y, direct arrow because M2 is last in the chain)
*
* REFERENCE:    Generalizes the single-mediator logic in Box S2 to 2 
*               causally ordered mediators, following the same decomposition
*               logic implemented in the R package "paths" (Zhou & Yamamoto
*               2020, "Tracing Causal Paths from Experimental and
*               Observational Data") — Type I decomposition.
*
* IDENTIFYING ASSUMPTIONS (stronger than the single- mediator case):
*   1. No unmeasured confounding of A-Y, A-M1, A-M2, M1-Y, M2-Y (as before)
*   2. No unmeasured confounding of the M1-M2 relationship, AND in
*      particular no M1-M2 confounder that is itself affected by A
*      (this is the extra assumption that comes with adding a second,
*      causally ordered mediator — $W must include anything that
*      confounds M1 and M2 and is measured pre-treatment)
*   3. Correctly specified mediator and outcome models
*===============================================================================

********* Redefine variables for multiple mediators setting:

**** Pre-treatment baseline confounders 
global W="age sex comorbbin mny asasmri hla pertvt1 ibdbl emmtvt1 comedtvt1 basfitotalt1"

**** Baseline values of mediators
global bM1="crpt1"
global bM2="asdastotalt1_pro"

**** Treatment (binary)
global A="bionew"

**** Mediators (continuous)
global M1="crpt2"
global M2="asdastotalt2_pro"

**** Outcome (continuous)
global Y="basfitotalt2"

order id $W $bM1 $bM2 $A $M1 $M2 $Y 

capture program drop run_gformula_paths
program define run_gformula_paths, rclass

	///////////////// Step 1 — Model the observed data, in causal order

	* --- Outcome model: BASFI ~ bionew + M1 + M2 + W + bM1 + bM2 ---
	regress $Y $A $M1 $M2 $W $bM1 $bM2
	estimates store outcome

	* --- Mediator 2 model: M2 ~ bionew + M1 + W + bM1 + bM2 (CRP causally precedes BASDAI) ---
	regress $M2 $A $M1 $W $bM1 $bM2
	estimates store mediator2

	* --- Mediator 1 model: M1 ~ bionew + W + bM1 + bM2 ---
	regress $M1 $A $W $bM1 $bM2
	estimates store mediator1

	///////////////// Step 2 — Monte Carlo simulation of nested counterfactuals

	*----------------------------------------------------------------------
	* Step 2.1 — Counterfactual M1 under treatment (1) and control (0)
	*----------------------------------------------------------------------
	estimates restore mediator1

	gen M1_1 = _b[_cons] + _b[$A]*1 + ///
	           _b[age]*age + _b[sex]*sex + _b[comorbbin]*comorbbin + ///
	           _b[mny]*mny + _b[asasmri]*asasmri + _b[hla]*hla + ///
	           _b[pertvt1]*pertvt1 + _b[ibdbl]*ibdbl + ///
	           _b[emmtvt1]*emmtvt1 + _b[comedtvt1]*comedtvt1 + ///
	           _b[$bM1]*$bM1 + _b[$bM2]*$bM2 + ///
	           _b[basfitotalt1]*basfitotalt1

	gen M1_0 = _b[_cons] + _b[$A]*0 + ///
	           _b[age]*age + _b[sex]*sex + _b[comorbbin]*comorbbin + ///
	           _b[mny]*mny + _b[asasmri]*asasmri + _b[hla]*hla + ///
	           _b[pertvt1]*pertvt1 + _b[ibdbl]*ibdbl + ///
	           _b[emmtvt1]*emmtvt1 + _b[comedtvt1]*comedtvt1 + ///
	           _b[$bM1]*$bM1 + _b[$bM2]*$bM2 + ///
	           _b[basfitotalt1]*basfitotalt1

	*----------------------------------------------------------------------
	* Step 2.2 — Counterfactual M2 at exactly the 3 combinations needed
	*   M2_A0_M1_0 : bionew=0 in the M2 model, M1 held at its CONTROL value
	*   M2_A0_M1_1 : bionew=0 in the M2 model, M1 held at its TREATED value
	*   M2_A1_M1_1 : bionew=1 in the M2 model, M1 held at its TREATED value
	* (bionew=1, M1=M1_0 is never needed for the Type I decomposition)
	*----------------------------------------------------------------------
	estimates restore mediator2

	gen M2_A0_M1_0 = _b[_cons] + _b[$A]*0 + _b[$M1]*M1_0 + ///
	                 _b[age]*age + _b[sex]*sex + _b[comorbbin]*comorbbin + ///
	                 _b[mny]*mny + _b[asasmri]*asasmri + _b[hla]*hla + ///
	                 _b[pertvt1]*pertvt1 + _b[ibdbl]*ibdbl + ///
	                 _b[emmtvt1]*emmtvt1 + _b[comedtvt1]*comedtvt1 + ///
	                 _b[$bM1]*$bM1 + _b[$bM2]*$bM2 + ///
	                 _b[basfitotalt1]*basfitotalt1

	gen M2_A0_M1_1 = _b[_cons] + _b[$A]*0 + _b[$M1]*M1_1 + ///
	                 _b[age]*age + _b[sex]*sex + _b[comorbbin]*comorbbin + ///
	                 _b[mny]*mny + _b[asasmri]*asasmri + _b[hla]*hla + ///
	                 _b[pertvt1]*pertvt1 + _b[ibdbl]*ibdbl + ///
	                 _b[emmtvt1]*emmtvt1 + _b[comedtvt1]*comedtvt1 + ///
	                 _b[$bM1]*$bM1 + _b[$bM2]*$bM2 + ///
	                 _b[basfitotalt1]*basfitotalt1

	gen M2_A1_M1_1 = _b[_cons] + _b[$A]*1 + _b[$M1]*M1_1 + ///
	                 _b[age]*age + _b[sex]*sex + _b[comorbbin]*comorbbin + ///
	                 _b[mny]*mny + _b[asasmri]*asasmri + _b[hla]*hla + ///
	                 _b[pertvt1]*pertvt1 + _b[ibdbl]*ibdbl + ///
	                 _b[emmtvt1]*emmtvt1 + _b[comedtvt1]*comedtvt1 + ///
	                 _b[$bM1]*$bM1 + _b[$bM2]*$bM2 + ///
	                 _b[basfitotalt1]*basfitotalt1

	*----------------------------------------------------------------------
	* Step 2.3 — Counterfactual BASFI at the 4 nested (a0, M1, M2) combos
	*   g(0,0,0) : everyone untreated, mediators at their untreated values
	*   g(1,0,0) : outcome sees bionew=1, mediators STILL at untreated values
	*              -> isolates the direct effect A->Y
	*   g(1,1,0) : outcome sees bionew=1, M1 now treated, M2 still at its
	*              "A=0,M1=M1_1" value -> the extra move from g(1,0,0) to
	*              g(1,1,0) isolates the path-specific effect via M1
	*              (which may travel on through M2)
	*   g(1,1,1) : outcome sees bionew=1, M1 treated, M2 now fully treated
	*              -> the extra move from g(1,1,0) to g(1,1,1) isolates
	*              the path-specific effect via M2 (direct arrow)
	*----------------------------------------------------------------------
	estimates restore outcome

	gen Y_000 = _b[_cons] + _b[$A]*0 + _b[$M1]*M1_0 + _b[$M2]*M2_A0_M1_0 + ///
	            _b[age]*age + _b[sex]*sex + _b[comorbbin]*comorbbin + ///
	            _b[mny]*mny + _b[asasmri]*asasmri + _b[hla]*hla + ///
	            _b[pertvt1]*pertvt1 + _b[ibdbl]*ibdbl + ///
	            _b[emmtvt1]*emmtvt1 + _b[comedtvt1]*comedtvt1 + ///
	            _b[$bM1]*$bM1 + _b[$bM2]*$bM2 + ///
	            _b[basfitotalt1]*basfitotalt1

	gen Y_100 = _b[_cons] + _b[$A]*1 + _b[$M1]*M1_0 + _b[$M2]*M2_A0_M1_0 + ///
	            _b[age]*age + _b[sex]*sex + _b[comorbbin]*comorbbin + ///
	            _b[mny]*mny + _b[asasmri]*asasmri + _b[hla]*hla + ///
	            _b[pertvt1]*pertvt1 + _b[ibdbl]*ibdbl + ///
	            _b[emmtvt1]*emmtvt1 + _b[comedtvt1]*comedtvt1 + ///
	            _b[$bM1]*$bM1 + _b[$bM2]*$bM2 + ///
	            _b[basfitotalt1]*basfitotalt1

	gen Y_110 = _b[_cons] + _b[$A]*1 + _b[$M1]*M1_1 + _b[$M2]*M2_A0_M1_1 + ///
	            _b[age]*age + _b[sex]*sex + _b[comorbbin]*comorbbin + ///
	            _b[mny]*mny + _b[asasmri]*asasmri + _b[hla]*hla + ///
	            _b[pertvt1]*pertvt1 + _b[ibdbl]*ibdbl + ///
	            _b[emmtvt1]*emmtvt1 + _b[comedtvt1]*comedtvt1 + ///
	            _b[$bM1]*$bM1 + _b[$bM2]*$bM2 + ///
	            _b[basfitotalt1]*basfitotalt1

	gen Y_111 = _b[_cons] + _b[$A]*1 + _b[$M1]*M1_1 + _b[$M2]*M2_A1_M1_1 + ///
	            _b[age]*age + _b[sex]*sex + _b[comorbbin]*comorbbin + ///
	            _b[mny]*mny + _b[asasmri]*asasmri + _b[hla]*hla + ///
	            _b[pertvt1]*pertvt1 + _b[ibdbl]*ibdbl + ///
	            _b[emmtvt1]*emmtvt1 + _b[comedtvt1]*comedtvt1 + ///
	            _b[$bM1]*$bM1 + _b[$bM2]*$bM2 + ///
	            _b[basfitotalt1]*basfitotalt1

	///////////////// Step 3 — Path-specific effects (Type I decomposition)

	qui sum Y_000
	scalar m000 = r(mean)
	qui sum Y_100
	scalar m100 = r(mean)
	qui sum Y_110
	scalar m110 = r(mean)
	qui sum Y_111
	scalar m111 = r(mean)

	return scalar direct = m100 - m000          // bionew -> BASFI
	return scalar pse_m1 = m110 - m100          // bionew -> M1 ~> BASFI
	return scalar pse_m2 = m111 - m110          // bionew -> M2 -> BASFI
	return scalar total  = m111 - m000          // Total effect (should equal direct + pse_m1 + pse_m2)

	///////////////// Clean up before the next bootstrap replication

	capture drop M1_1 M1_0 M2_A0_M1_0 M2_A0_M1_1 M2_A1_M1_1 ///
	             Y_000 Y_100 Y_110 Y_111 _est_outcome _est_mediator1 _est_mediator2

end

* Run the bootstrap across 1000 replications
set seed `seed'
bootstrap Direct=r(direct) PSE_M1=r(pse_m1) PSE_M2=r(pse_m2) Total=r(total), ///
	reps(`sims') nodots: run_gformula_paths

* Display the percentile-based bootstrap 95% Confidence Intervals
noi estat bootstrap, percentile

savebs "S4" "Edge g-formula: main paths without interaction (manual)"
savesheet "S4"



} // close 3


******************************************************************************************************************************************************
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
**** 3.1 - SUPPLEMENTARY BOX S2.1: two causally ordered mediators, full set of path-specific effects (manual)
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
******************************************************************************************************************************************************

if `run' == 3.1 {

*===============================================================================
* METHOD:       Manual extension of the R paths package (type 1 decomposition)
*===============================================================================

*===============================================================================
*
* METHOD:       Parametric edge g-formula for path-specific effects (PSEs)
*               with two causally ordered mediators: M1 -> M2 -> Y
*
* ESTIMAND:     Full set of effects:
*                 - Direct        : A -> Y
*                 - PSE_M1        : A -> M1 ~> Y   (combined, may travel through M2)
*                 - PSE_M1_direct : A -> M1 -> Y   (M1's own edge into Y, skipping M2)
*                 - PSE_M1_via_M2 : A -> M1 -> M2 -> Y (M1's effect that runs through M2)
*                 - PSE_M2        : A -> M2 -> Y   (M2's own edge into Y, direct arrow)
*                 - Total         : A ~> Y
*                 - ATE_M1        : total effect of A on M1 itself (mediator-scale)
*                 - ATE_M2        : total effect of A on M2 itself (mediator-scale,
*                                   combining A's direct edge into M2 AND the part
*                                   that arrives via M1)
*                 - M1_on_M2      : M1's own causal effect on M2 (mediator-to-mediator,
*                                   M2-scale), i.e. ATE_M2 with A's direct edge into
*                                   M2 held at 0 — NOT a path-specific effect on Y
*
* REFERENCE:    Generalizes the single-mediator logic in Box S2 to 2 causally
*               ordered mediators, following the same decomposition logic
*               implemented in the R package "paths" (Zhou & Yamamoto 2020,
*               "Tracing Causal Paths from Experimental and Observational Data")
*               — Type I decomposition — plus the edge g-formula terminology
*               (Shpitser et al.) used to derive the finer PSE_M1 split and the
*               mediator-to-mediator / treatment-to-mediator effects.
*
* IDENTIFYING ASSUMPTIONS (stronger than the single-mediator case):
*   1. No unmeasured confounding of A-Y, A-M1, A-M2, M1-Y, M2-Y (as before)
*   2. No unmeasured confounding of the M1-M2 relationship, AND in
*      particular no M1-M2 confounder that is itself affected by A
*      (this is the extra assumption that comes with adding a second,
*      causally ordered mediator — $W must include anything that
*      confounds M1 and M2 and is measured pre-treatment)
*   3. Correctly specified mediator and outcome models
*
* NOTE: PSE_M1_direct, PSE_M1_via_M2, ATE_M1, ATE_M2, and M1_on_M2 all reuse the
* three fitted models above — no additional model or identifying assumption is
* needed beyond what PSE_M1/PSE_M2/Direct/Total already require.
*===============================================================================


********* Redefine variables for multiple mediators setting:

**** Pre-treatment baseline confounders 
global W="age sex comorbbin mny asasmri hla pertvt1 ibdbl emmtvt1 comedtvt1 basfitotalt1"

**** Baseline values of mediators
global bM1="crpt1"
global bM2="asdastotalt1_pro"

**** Treatment (binary)
global A="bionew"

**** Mediators (continuous)
global M1="crpt2"
global M2="asdastotalt2_pro"

**** Outcome (continuous)
global Y="basfitotalt2"

order id $W $bM1 $bM2 $A $M1 $M2 $Y 

capture program drop run_gformula_paths
program define run_gformula_paths, rclass

	///////////////// Step 1 — Model the observed data, in causal order

	* --- Outcome model: BASFI ~ bionew + M1 + M2 + W + bM1 + bM2 ---
	regress $Y $A $M1 $M2 $W $bM1 $bM2
	estimates store outcome

	* --- Mediator 2 model: M2 ~ bionew + M1 + W + bM1 + bM2 ---
	regress $M2 $A $M1 $W $bM1 $bM2
	estimates store mediator2

	* --- Mediator 1 model: M1 ~ bionew + W + bM1 + bM2 ---
	regress $M1 $A $W $bM1 $bM2
	estimates store mediator1

	///////////////// Step 2 — Monte Carlo simulation of nested counterfactuals

	*----------------------------------------------------------------------
	* Step 2.1 — Counterfactual M1 under treatment (1) and control (0)
	*----------------------------------------------------------------------
	estimates restore mediator1

	gen M1_1 = _b[_cons] + _b[$A]*1 + ///
	           _b[age]*age + _b[sex]*sex + _b[comorbbin]*comorbbin + ///
	           _b[mny]*mny + _b[asasmri]*asasmri + _b[hla]*hla + ///
	           _b[pertvt1]*pertvt1 + _b[ibdbl]*ibdbl + ///
	           _b[emmtvt1]*emmtvt1 + _b[comedtvt1]*comedtvt1 + ///
	           _b[$bM1]*$bM1 + _b[$bM2]*$bM2 + ///
	           _b[basfitotalt1]*basfitotalt1

	gen M1_0 = _b[_cons] + _b[$A]*0 + ///
	           _b[age]*age + _b[sex]*sex + _b[comorbbin]*comorbbin + ///
	           _b[mny]*mny + _b[asasmri]*asasmri + _b[hla]*hla + ///
	           _b[pertvt1]*pertvt1 + _b[ibdbl]*ibdbl + ///
	           _b[emmtvt1]*emmtvt1 + _b[comedtvt1]*comedtvt1 + ///
	           _b[$bM1]*$bM1 + _b[$bM2]*$bM2 + ///
	           _b[basfitotalt1]*basfitotalt1

	*----------------------------------------------------------------------
	* Step 2.2 — Counterfactual M2 at exactly the 3 combinations needed
	*   M2_A0_M1_0 : bionew=0 in the M2 model, M1 held at its CONTROL value
	*   M2_A0_M1_1 : bionew=0 in the M2 model, M1 held at its TREATED value
	*   M2_A1_M1_1 : bionew=1 in the M2 model, M1 held at its TREATED value
	* (bionew=1, M1=M1_0 is never needed for the Type I decomposition)
	*----------------------------------------------------------------------
	estimates restore mediator2

	gen M2_A0_M1_0 = _b[_cons] + _b[$A]*0 + _b[$M1]*M1_0 + ///
	                 _b[age]*age + _b[sex]*sex + _b[comorbbin]*comorbbin + ///
	                 _b[mny]*mny + _b[asasmri]*asasmri + _b[hla]*hla + ///
	                 _b[pertvt1]*pertvt1 + _b[ibdbl]*ibdbl + ///
	                 _b[emmtvt1]*emmtvt1 + _b[comedtvt1]*comedtvt1 + ///
	                 _b[$bM1]*$bM1 + _b[$bM2]*$bM2 + ///
	                 _b[basfitotalt1]*basfitotalt1

	gen M2_A0_M1_1 = _b[_cons] + _b[$A]*0 + _b[$M1]*M1_1 + ///
	                 _b[age]*age + _b[sex]*sex + _b[comorbbin]*comorbbin + ///
	                 _b[mny]*mny + _b[asasmri]*asasmri + _b[hla]*hla + ///
	                 _b[pertvt1]*pertvt1 + _b[ibdbl]*ibdbl + ///
	                 _b[emmtvt1]*emmtvt1 + _b[comedtvt1]*comedtvt1 + ///
	                 _b[$bM1]*$bM1 + _b[$bM2]*$bM2 + ///
	                 _b[basfitotalt1]*basfitotalt1

	gen M2_A1_M1_1 = _b[_cons] + _b[$A]*1 + _b[$M1]*M1_1 + ///
	                 _b[age]*age + _b[sex]*sex + _b[comorbbin]*comorbbin + ///
	                 _b[mny]*mny + _b[asasmri]*asasmri + _b[hla]*hla + ///
	                 _b[pertvt1]*pertvt1 + _b[ibdbl]*ibdbl + ///
	                 _b[emmtvt1]*emmtvt1 + _b[comedtvt1]*comedtvt1 + ///
	                 _b[$bM1]*$bM1 + _b[$bM2]*$bM2 + ///
	                 _b[basfitotalt1]*basfitotalt1

	*----------------------------------------------------------------------
	* Step 2.3 — Counterfactual BASFI at the nested (a0, M1, M2) combos
	*   g(0,0,0)   : everyone untreated, mediators at their untreated values
	*   g(1,0,0)   : outcome sees bionew=1, mediators STILL at untreated values
	*                -> isolates the direct effect A->Y
	*   g_mid      : outcome sees bionew=1, M1 now treated, M2 still at its
	*                "A=0,M1=M1_0" value (i.e. M2 NOT updated for M1's move)
	*                -> the extra move from g(1,0,0) to g_mid isolates M1's
	*                OWN edge into Y, skipping M2 (PSE_M1_direct)
	*   g(1,1,0)   : outcome sees bionew=1, M1 treated, M2 now updated to
	*                reflect M1's move ("A=0,M1=M1_1")
	*                -> the extra move from g_mid to g(1,1,0) isolates the
	*                part of M1's effect that runs through M2 (PSE_M1_via_M2)
	*   g(1,1,1)   : outcome sees bionew=1, M1 treated, M2 now fully treated
	*                -> the extra move from g(1,1,0) to g(1,1,1) isolates
	*                the path-specific effect via M2 (PSE_M2, direct arrow)
	*----------------------------------------------------------------------
	estimates restore outcome

	gen Y_000 = _b[_cons] + _b[$A]*0 + _b[$M1]*M1_0 + _b[$M2]*M2_A0_M1_0 + ///
	            _b[age]*age + _b[sex]*sex + _b[comorbbin]*comorbbin + ///
	            _b[mny]*mny + _b[asasmri]*asasmri + _b[hla]*hla + ///
	            _b[pertvt1]*pertvt1 + _b[ibdbl]*ibdbl + ///
	            _b[emmtvt1]*emmtvt1 + _b[comedtvt1]*comedtvt1 + ///
	            _b[$bM1]*$bM1 + _b[$bM2]*$bM2 + ///
	            _b[basfitotalt1]*basfitotalt1

	gen Y_100 = _b[_cons] + _b[$A]*1 + _b[$M1]*M1_0 + _b[$M2]*M2_A0_M1_0 + ///
	            _b[age]*age + _b[sex]*sex + _b[comorbbin]*comorbbin + ///
	            _b[mny]*mny + _b[asasmri]*asasmri + _b[hla]*hla + ///
	            _b[pertvt1]*pertvt1 + _b[ibdbl]*ibdbl + ///
	            _b[emmtvt1]*emmtvt1 + _b[comedtvt1]*comedtvt1 + ///
	            _b[$bM1]*$bM1 + _b[$bM2]*$bM2 + ///
	            _b[basfitotalt1]*basfitotalt1

	* --- Y_mid: A=1, M1 treated, but M2 held at its "M1 still untreated" value ---
	gen Y_mid = _b[_cons] + _b[$A]*1 + _b[$M1]*M1_1 + _b[$M2]*M2_A0_M1_0 + ///
	            _b[age]*age + _b[sex]*sex + _b[comorbbin]*comorbbin + ///
	            _b[mny]*mny + _b[asasmri]*asasmri + _b[hla]*hla + ///
	            _b[pertvt1]*pertvt1 + _b[ibdbl]*ibdbl + ///
	            _b[emmtvt1]*emmtvt1 + _b[comedtvt1]*comedtvt1 + ///
	            _b[$bM1]*$bM1 + _b[$bM2]*$bM2 + ///
	            _b[basfitotalt1]*basfitotalt1

	gen Y_110 = _b[_cons] + _b[$A]*1 + _b[$M1]*M1_1 + _b[$M2]*M2_A0_M1_1 + ///
	            _b[age]*age + _b[sex]*sex + _b[comorbbin]*comorbbin + ///
	            _b[mny]*mny + _b[asasmri]*asasmri + _b[hla]*hla + ///
	            _b[pertvt1]*pertvt1 + _b[ibdbl]*ibdbl + ///
	            _b[emmtvt1]*emmtvt1 + _b[comedtvt1]*comedtvt1 + ///
	            _b[$bM1]*$bM1 + _b[$bM2]*$bM2 + ///
	            _b[basfitotalt1]*basfitotalt1

	gen Y_111 = _b[_cons] + _b[$A]*1 + _b[$M1]*M1_1 + _b[$M2]*M2_A1_M1_1 + ///
	            _b[age]*age + _b[sex]*sex + _b[comorbbin]*comorbbin + ///
	            _b[mny]*mny + _b[asasmri]*asasmri + _b[hla]*hla + ///
	            _b[pertvt1]*pertvt1 + _b[ibdbl]*ibdbl + ///
	            _b[emmtvt1]*emmtvt1 + _b[comedtvt1]*comedtvt1 + ///
	            _b[$bM1]*$bM1 + _b[$bM2]*$bM2 + ///
	            _b[basfitotalt1]*basfitotalt1

	///////////////// Step 3 — Path-specific effects (Type I decomposition)

	qui sum Y_000
	scalar m000 = r(mean)
	qui sum Y_100
	scalar m100 = r(mean)
	qui sum Y_mid
	scalar mmid = r(mean)
	qui sum Y_110
	scalar m110 = r(mean)
	qui sum Y_111
	scalar m111 = r(mean)

	return scalar direct = m100 - m000          // A -> Y
	return scalar pse_m1 = m110 - m100          // A -> M1 ~> Y (combined)
	return scalar pse_m2 = m111 - m110          // A -> M2 -> Y
	return scalar total  = m111 - m000          // direct + pse_m1 + pse_m2

	* --- Split of pse_m1 into its two constituent paths ---
	return scalar pse_m1_direct = mmid - m100   // A -> M1 -> Y (skips M2)
	return scalar pse_m1_via_m2 = m110 - mmid   // A -> M1 -> M2 -> Y

	///////////////// Step 4 — Effects of treatment on each mediator, and
	///////////////// the mediator-to-mediator effect

	qui sum M1_1
	scalar mean_M1_1 = r(mean)
	qui sum M1_0
	scalar mean_M1_0 = r(mean)
	qui sum M2_A1_M1_1
	scalar mean_M2_A1_M1_1 = r(mean)
	qui sum M2_A0_M1_0
	scalar mean_M2_A0_M1_0 = r(mean)
	qui sum M2_A0_M1_1
	scalar mean_M2_A0_M1_1 = r(mean)

	* Total effect of treatment on M1 itself (mediator-scale, not on Y)
	return scalar ate_m1 = mean_M1_1 - mean_M1_0

	* Total effect of treatment on M2 itself (mediator-scale, not on Y) —
	* combines A's direct edge into M2 AND the part that arrives via M1
	return scalar ate_m2 = mean_M2_A1_M1_1 - mean_M2_A0_M1_0

	* M1's own causal effect on M2 (mediator-to-mediator, M2-scale) —
	* A's direct edge into M2 held at 0 in both terms, so this isolates
	* only the M1 -> M2 edge. NOT a path-specific effect on Y.
	return scalar m1_on_m2 = mean_M2_A0_M1_1 - mean_M2_A0_M1_0

	///////////////// Clean up before the next bootstrap replication

	capture drop M1_1 M1_0 M2_A0_M1_0 M2_A0_M1_1 M2_A1_M1_1 ///
	             Y_000 Y_100 Y_mid Y_110 Y_111 _est_outcome _est_mediator1 _est_mediator2

end

* Run the bootstrap across 1000 replications
set seed `seed'
bootstrap Direct=r(direct) PSE_M1=r(pse_m1) PSE_M2=r(pse_m2) Total=r(total) ///
	PSE_M1_direct=r(pse_m1_direct) PSE_M1_via_M2=r(pse_m1_via_m2) ///
	ATE_M1=r(ate_m1) ATE_M2=r(ate_m2) M1_on_M2=r(m1_on_m2), ///
	reps(`sims') nodots: run_gformula_paths // manuscript: 1000

* Display the percentile-based bootstrap 95% Confidence Intervals
noi estat bootstrap, percentile

savebs "S5" "Edge g-formula: all paths without interaction (manual)"
* the same four main paths also fill Table S4, so the two specifications sit
* side by side there without a table of their own
savebs "S4" "Edge g-formula: main paths without interaction (manual)" "" "Direct PSE_M1 PSE_M2 Total"
savesheet "S5"
savesheet "S4"

*----------------------------------------------------------
* Create a dataset for plotting in R
*----------------------------------------------------------

clear
set obs 4
gen str25 effect=""
gen ATE=.
gen ATE_LL=.
gen ATE_UL=.
matrix b=e(b)
matrix ci=e(ci_percentile)
replace effect="Direct" in 1
replace ATE=b[1,1] in 1
replace ATE_LL=ci[1,1] in 1
replace ATE_UL=ci[2,1] in 1
replace effect="PSE_M1" in 2
replace ATE=b[1,2] in 2
replace ATE_LL=ci[1,2] in 2
replace ATE_UL=ci[2,2] in 2
replace effect="PSE_M2" in 3
replace ATE=b[1,3] in 3
replace ATE_LL=ci[1,3] in 3
replace ATE_UL=ci[2,3] in 3
replace effect="Total" in 4
replace ATE=b[1,4] in 4
replace ATE_LL=ci[1,4] in 4
replace ATE_UL=ci[2,4] in 4
export delimited using "${tables}results_edgegformula.csv", replace


} // close 3.1


******************************************************************************************************************************************************
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
**** 3.2 - SUPPLEMENTARY BOX S2.1: two causally ordered mediators, full set of path-specific effects, with interactions (manual)
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
******************************************************************************************************************************************************

if `run' == 3.2 {

*===============================================================================
* METHOD:       Manual extension of the R paths package (type 1 decomposition)
*               With signficant interactions
*===============================================================================

*** Variables for the two-mediator setting (repeated so the block runs on its own)
**** Pre-treatment baseline confounders 
global W="age sex comorbbin mny asasmri hla pertvt1 ibdbl emmtvt1 comedtvt1 basfitotalt1"

**** Baseline values of mediators
global bM1="crpt1"
global bM2="asdastotalt1_pro"

**** Treatment (binary)
global A="bionew"

**** Mediators (continuous)
global M1="crpt2"
global M2="asdastotalt2_pro"

**** Outcome (continuous)
global Y="basfitotalt2"

order id $W $bM1 $bM2 $A $M1 $M2 $Y 


///////// Test interactions

* The four models are fitted quietly; only the test of each interaction term is shown.

local pcut = 0.20 // significance level for calling an interaction significant

* 1. A x M1 in the outcome model
qui regress $Y $A $M1 $M2 $W $bM1 $bM2 i.$A#c.$M1
qui testparm i.$A#c.$M1
local p1 = r(p)
local v1 "not significant"
if `p1' < `pcut' {
	local v1 "significant"
}

* 2. A x M2 in the outcome model
qui regress $Y $A $M1 $M2 $W $bM1 $bM2 i.$A#c.$M2
qui testparm i.$A#c.$M2
local p2 = r(p)
local v2 "not significant"
if `p2' < `pcut' {
	local v2 "significant"
}

* 3. M1 x M2 in the outcome model
qui regress $Y $A $M1 $M2 $W $bM1 $bM2 c.$M1#c.$M2
qui testparm c.$M1#c.$M2
local p3 = r(p)
local v3 "not significant"
if `p3' < `pcut' {
	local v3 "significant"
}

* 4. A x M1 in the M2 (mediator) model
qui regress $M2 $A $M1 $W $bM1 $bM2 i.$A#c.$M1
qui testparm i.$A#c.$M1
local p4 = r(p)
local v4 "not significant"
if `p4' < `pcut' {
	local v4 "significant"
}

noi di ""
noi di as text "Interaction tests (significant at p < `pcut')"
noi di as text "{hline 70}"
noi di as text "1. $A x $M1 in the outcome model      : p = " as result %6.4f `p1' as text "  `v1'"
noi di as text "2. $A x $M2 in the outcome model      : p = " as result %6.4f `p2' as text "  `v2'"
noi di as text "3. $M1 x $M2 in the outcome model     : p = " as result %6.4f `p3' as text "  `v3'"
noi di as text "4. $A x $M1 in the $M2 mediator model : p = " as result %6.4f `p4' as text "  `v4'"
noi di as text "{hline 70}"
noi di ""

* The two retained interactions, refit together. Both product terms contain
* $M1, so they compete for the same variance and each p-value rises; the joint
* test is what says whether there is an interaction signal at all.
qui regress $Y $A $M1 $M2 $W $bM1 $bM2 c.$A#c.$M1 c.$M1#c.$M2
qui testparm c.$A#c.$M1
local p5 = r(p)
qui testparm c.$M1#c.$M2
local p6 = r(p)
qui testparm c.$A#c.$M1 c.$M1#c.$M2
local p7 = r(p)

noi di as text "Both retained interactions in one model"
noi di as text "{hline 70}"
noi di as text "5. $A x $M1                           : p = " as result %6.4f `p5'
noi di as text "6. $M1 x $M2                          : p = " as result %6.4f `p6'
noi di as text "7. joint test of the two              : p = " as result %6.4f `p7'
noi di as text "{hline 70}"
noi di ""

savediag "S5" "Edge g-formula: all paths with interaction (manual)" "p interaction AxM1 outcome model" `p1'
savediag "S5" "Edge g-formula: all paths with interaction (manual)" "p interaction AxM2 outcome model" `p2'
savediag "S5" "Edge g-formula: all paths with interaction (manual)" "p interaction M1xM2 outcome model" `p3'
savediag "S5" "Edge g-formula: all paths with interaction (manual)" "p interaction AxM1 mediator model" `p4'
savediag "S5" "Edge g-formula: all paths with interaction (manual)" "p interaction AxM1 both in model" `p5'
savediag "S5" "Edge g-formula: all paths with interaction (manual)" "p interaction M1xM2 both in model" `p6'
savediag "S5" "Edge g-formula: all paths with interaction (manual)" "p interaction joint test of the two" `p7'

*===============================================================================
*
* METHOD:       Parametric edge g-formula for path-specific effects (PSEs)
*               with two causally ordered mediators: M1 -> M2 -> Y
*
* ESTIMAND:     Full set of effects discussed:
*                 - Direct        : A -> Y
*                 - PSE_M1        : A -> M1 ~> Y   (combined, may travel through M2)
*                 - PSE_M1_direct : A -> M1 -> Y   (M1's own edge into Y, skipping M2)
*                 - PSE_M1_via_M2 : A -> M1 -> M2 -> Y (M1's effect that runs through M2)
*                 - PSE_M2        : A -> M2 -> Y   (M2's own edge into Y, direct arrow)
*                 - Total         : A ~> Y
*                 - ATE_M1        : total effect of A on M1 itself (mediator-scale)
*                 - ATE_M2        : total effect of A on M2 itself (mediator-scale,
*                                   combining A's direct edge into M2 AND the part
*                                   that arrives via M1)
*                 - M1_on_M2      : M1's own causal effect on M2 (mediator-to-mediator,
*                                   M2-scale), i.e. ATE_M2 with A's direct edge into
*                                   M2 held at 0 — NOT a path-specific effect on Y
*
* REFERENCE:    Generalizes the single-mediator logic in Box S2 to 2 causally
*               ordered mediators, following the same decomposition logic
*               implemented in the R package "paths" (Zhou & Yamamoto 2020,
*               "Tracing Causal Paths from Experimental and Observational Data")
*               — Type I decomposition — plus the edge g-formula terminology
*               (Shpitser et al.) used to derive the finer PSE_M1 split and the
*               mediator-to-mediator / treatment-to-mediator effects.
*
* IDENTIFYING ASSUMPTIONS (stronger than the single-mediator case):
*   1. No unmeasured confounding of A-Y, A-M1, A-M2, M1-Y, M2-Y (as before)
*   2. No unmeasured confounding of the M1-M2 relationship, AND in
*      particular no M1-M2 confounder that is itself affected by A
*      (this is the extra assumption that comes with adding a second,
*      causally ordered mediator — need to make sure $W includes anything that
*      confounds M1 and M2 and is measured pre-treatment)
*   3. Correctly specified mediator and outcome models
*
* NOTE: PSE_M1_direct, PSE_M1_via_M2, ATE_M1, ATE_M2, and M1_on_M2 all reuse the
* three fitted models above — no additional model or identifying assumption is
* needed beyond what PSE_M1/PSE_M2/Direct/Total already require.
*
* THIS VERSION includes two interactions in the outcome model:
*   - bionew x crpt2                (A x M1)
*   - crpt2 x asdastotalt2_pro      (M1 x M2)
* They are the two that pass the one-at-a-time screening at p<0.2; the other
* two candidates (A x M2 in the outcome model, and A x M1 in the M2 model) are
* dropped. Step 0 below runs all of those tests and writes the p-values to the
* diagnostics file.
* Both retained interactions are refit together and addded through every Y_* counterfactual using
* that line's own (A, M1, M2) combination, not the raw observed variables.
*===============================================================================

********* Redefine variables for multiple mediators setting:

**** Pre-treatment baseline confounders 
global W="age sex comorbbin mny asasmri hla pertvt1 ibdbl emmtvt1 comedtvt1 basfitotalt1"

**** Baseline values of mediators
global bM1="crpt1"
global bM2="asdastotalt1_pro"

**** Treatment (binary)
global A="bionew"

**** Mediators (continuous)
global M1="crpt2"
global M2="asdastotalt2_pro"

**** Outcome (continuous)
global Y="basfitotalt2"

order id $W $bM1 $bM2 $A $M1 $M2 $Y 

capture program drop run_gformula_paths
program define run_gformula_paths, rclass

	///////////////// Step 1 — Model the observed data, in causal order

	* --- Outcome model: BASFI ~ bionew + M1 + M2 + W + bM1 + bM2 ---
	* --- PLUS the two interactions retained by Step 0: A x M1 and M1 x M2 ---
	regress $Y $A $M1 $M2 $W $bM1 $bM2 c.$A#c.$M1 c.$M1#c.$M2
	estimates store outcome

	* --- Mediator 2 model: M2 ~ bionew + M1 + W + bM1 + bM2 ---
	regress $M2 $A $M1 $W $bM1 $bM2
	estimates store mediator2

	* --- Mediator 1 model: M1 ~ bionew + W + bM1 + bM2 ---
	regress $M1 $A $W $bM1 $bM2
	estimates store mediator1

	///////////////// Step 2 — Monte Carlo simulation of nested counterfactuals

	*----------------------------------------------------------------------
	* Step 2.1 — Counterfactual M1 under treatment (1) and control (0)
	*----------------------------------------------------------------------
	estimates restore mediator1

	gen M1_1 = _b[_cons] + _b[$A]*1 + ///
	           _b[age]*age + _b[sex]*sex + _b[comorbbin]*comorbbin + ///
	           _b[mny]*mny + _b[asasmri]*asasmri + _b[hla]*hla + ///
	           _b[pertvt1]*pertvt1 + _b[ibdbl]*ibdbl + ///
	           _b[emmtvt1]*emmtvt1 + _b[comedtvt1]*comedtvt1 + ///
	           _b[$bM1]*$bM1 + _b[$bM2]*$bM2 + ///
	           _b[basfitotalt1]*basfitotalt1

	gen M1_0 = _b[_cons] + _b[$A]*0 + ///
	           _b[age]*age + _b[sex]*sex + _b[comorbbin]*comorbbin + ///
	           _b[mny]*mny + _b[asasmri]*asasmri + _b[hla]*hla + ///
	           _b[pertvt1]*pertvt1 + _b[ibdbl]*ibdbl + ///
	           _b[emmtvt1]*emmtvt1 + _b[comedtvt1]*comedtvt1 + ///
	           _b[$bM1]*$bM1 + _b[$bM2]*$bM2 + ///
	           _b[basfitotalt1]*basfitotalt1

	*----------------------------------------------------------------------
	* Step 2.2 — Counterfactual M2 at exactly the 3 combinations needed
	*   M2_A0_M1_0 : bionew=0 in the M2 model, M1 held at its CONTROL value
	*   M2_A0_M1_1 : bionew=0 in the M2 model, M1 held at its TREATED value
	*   M2_A1_M1_1 : bionew=1 in the M2 model, M1 held at its TREATED value
	* (bionew=1, M1=M1_0 is never needed for the Type I decomposition)
	*----------------------------------------------------------------------
	estimates restore mediator2

	gen M2_A0_M1_0 = _b[_cons] + _b[$A]*0 + _b[$M1]*M1_0 + ///
	                 _b[age]*age + _b[sex]*sex + _b[comorbbin]*comorbbin + ///
	                 _b[mny]*mny + _b[asasmri]*asasmri + _b[hla]*hla + ///
	                 _b[pertvt1]*pertvt1 + _b[ibdbl]*ibdbl + ///
	                 _b[emmtvt1]*emmtvt1 + _b[comedtvt1]*comedtvt1 + ///
	                 _b[$bM1]*$bM1 + _b[$bM2]*$bM2 + ///
	                 _b[basfitotalt1]*basfitotalt1

	gen M2_A0_M1_1 = _b[_cons] + _b[$A]*0 + _b[$M1]*M1_1 + ///
	                 _b[age]*age + _b[sex]*sex + _b[comorbbin]*comorbbin + ///
	                 _b[mny]*mny + _b[asasmri]*asasmri + _b[hla]*hla + ///
	                 _b[pertvt1]*pertvt1 + _b[ibdbl]*ibdbl + ///
	                 _b[emmtvt1]*emmtvt1 + _b[comedtvt1]*comedtvt1 + ///
	                 _b[$bM1]*$bM1 + _b[$bM2]*$bM2 + ///
	                 _b[basfitotalt1]*basfitotalt1

	gen M2_A1_M1_1 = _b[_cons] + _b[$A]*1 + _b[$M1]*M1_1 + ///
	                 _b[age]*age + _b[sex]*sex + _b[comorbbin]*comorbbin + ///
	                 _b[mny]*mny + _b[asasmri]*asasmri + _b[hla]*hla + ///
	                 _b[pertvt1]*pertvt1 + _b[ibdbl]*ibdbl + ///
	                 _b[emmtvt1]*emmtvt1 + _b[comedtvt1]*comedtvt1 + ///
	                 _b[$bM1]*$bM1 + _b[$bM2]*$bM2 + ///
	                 _b[basfitotalt1]*basfitotalt1

	*----------------------------------------------------------------------
	* Step 2.3 — Counterfactual BASFI at the nested (a0, M1, M2) combos
	*   g(0,0,0)   : everyone untreated, mediators at their untreated values
	*   g(1,0,0)   : outcome sees bionew=1, mediators STILL at untreated values
	*                -> isolates the direct effect A->Y
	*   g_mid      : outcome sees bionew=1, M1 now treated, M2 still at its
	*                "A=0,M1=M1_0" value (i.e. M2 NOT updated for M1's move)
	*                -> the extra move from g(1,0,0) to g_mid isolates M1's
	*                OWN edge into Y, skipping M2 (PSE_M1_direct)
	*   g(1,1,0)   : outcome sees bionew=1, M1 treated, M2 now updated to
	*                reflect M1's move ("A=0,M1=M1_1")
	*                -> the extra move from g_mid to g(1,1,0) isolates the
	*                part of M1's effect that runs through M2 (PSE_M1_via_M2)
	*   g(1,1,1)   : outcome sees bionew=1, M1 treated, M2 now fully treated
	*                -> the extra move from g(1,1,0) to g(1,1,1) isolates
	*                the path-specific effect via M2 (PSE_M2, direct arrow)
	*----------------------------------------------------------------------
	estimates restore outcome

	gen Y_000 = _b[_cons] + _b[$A]*0 + _b[$M1]*M1_0 + _b[$M2]*M2_A0_M1_0 + ///
	            _b[age]*age + _b[sex]*sex + _b[comorbbin]*comorbbin + ///
	            _b[mny]*mny + _b[asasmri]*asasmri + _b[hla]*hla + ///
	            _b[pertvt1]*pertvt1 + _b[ibdbl]*ibdbl + ///
	            _b[emmtvt1]*emmtvt1 + _b[comedtvt1]*comedtvt1 + ///
	            _b[$bM1]*$bM1 + _b[$bM2]*$bM2 + ///
	            _b[basfitotalt1]*basfitotalt1 + ///
	            _b[c.$A#c.$M1]*(0*M1_0) + ///
	            _b[c.$M1#c.$M2]*(M1_0*M2_A0_M1_0)

	gen Y_100 = _b[_cons] + _b[$A]*1 + _b[$M1]*M1_0 + _b[$M2]*M2_A0_M1_0 + ///
	            _b[age]*age + _b[sex]*sex + _b[comorbbin]*comorbbin + ///
	            _b[mny]*mny + _b[asasmri]*asasmri + _b[hla]*hla + ///
	            _b[pertvt1]*pertvt1 + _b[ibdbl]*ibdbl + ///
	            _b[emmtvt1]*emmtvt1 + _b[comedtvt1]*comedtvt1 + ///
	            _b[$bM1]*$bM1 + _b[$bM2]*$bM2 + ///
	            _b[basfitotalt1]*basfitotalt1 + ///
	            _b[c.$A#c.$M1]*(1*M1_0) + ///
	            _b[c.$M1#c.$M2]*(M1_0*M2_A0_M1_0)

	* --- Y_mid: A=1, M1 treated, but M2 held at its "M1 still untreated" value ---
	gen Y_mid = _b[_cons] + _b[$A]*1 + _b[$M1]*M1_1 + _b[$M2]*M2_A0_M1_0 + ///
	            _b[age]*age + _b[sex]*sex + _b[comorbbin]*comorbbin + ///
	            _b[mny]*mny + _b[asasmri]*asasmri + _b[hla]*hla + ///
	            _b[pertvt1]*pertvt1 + _b[ibdbl]*ibdbl + ///
	            _b[emmtvt1]*emmtvt1 + _b[comedtvt1]*comedtvt1 + ///
	            _b[$bM1]*$bM1 + _b[$bM2]*$bM2 + ///
	            _b[basfitotalt1]*basfitotalt1 + ///
	            _b[c.$A#c.$M1]*(1*M1_1) + ///
	            _b[c.$M1#c.$M2]*(M1_1*M2_A0_M1_0)

	gen Y_110 = _b[_cons] + _b[$A]*1 + _b[$M1]*M1_1 + _b[$M2]*M2_A0_M1_1 + ///
	            _b[age]*age + _b[sex]*sex + _b[comorbbin]*comorbbin + ///
	            _b[mny]*mny + _b[asasmri]*asasmri + _b[hla]*hla + ///
	            _b[pertvt1]*pertvt1 + _b[ibdbl]*ibdbl + ///
	            _b[emmtvt1]*emmtvt1 + _b[comedtvt1]*comedtvt1 + ///
	            _b[$bM1]*$bM1 + _b[$bM2]*$bM2 + ///
	            _b[basfitotalt1]*basfitotalt1 + ///
	            _b[c.$A#c.$M1]*(1*M1_1) + ///
	            _b[c.$M1#c.$M2]*(M1_1*M2_A0_M1_1)

	gen Y_111 = _b[_cons] + _b[$A]*1 + _b[$M1]*M1_1 + _b[$M2]*M2_A1_M1_1 + ///
	            _b[age]*age + _b[sex]*sex + _b[comorbbin]*comorbbin + ///
	            _b[mny]*mny + _b[asasmri]*asasmri + _b[hla]*hla + ///
	            _b[pertvt1]*pertvt1 + _b[ibdbl]*ibdbl + ///
	            _b[emmtvt1]*emmtvt1 + _b[comedtvt1]*comedtvt1 + ///
	            _b[$bM1]*$bM1 + _b[$bM2]*$bM2 + ///
	            _b[basfitotalt1]*basfitotalt1 + ///
	            _b[c.$A#c.$M1]*(1*M1_1) + ///
	            _b[c.$M1#c.$M2]*(M1_1*M2_A1_M1_1)

	///////////////// Step 3 — Path-specific effects (Type I decomposition)

	qui sum Y_000
	scalar m000 = r(mean)
	qui sum Y_100
	scalar m100 = r(mean)
	qui sum Y_mid
	scalar mmid = r(mean)
	qui sum Y_110
	scalar m110 = r(mean)
	qui sum Y_111
	scalar m111 = r(mean)

	return scalar direct = m100 - m000          // A -> Y
	return scalar pse_m1 = m110 - m100          // A -> M1 ~> Y (combined)
	return scalar pse_m2 = m111 - m110          // A -> M2 -> Y
	return scalar total  = m111 - m000          // direct + pse_m1 + pse_m2

	* --- Split of pse_m1 into its two constituent paths ---
	return scalar pse_m1_direct = mmid - m100   // A -> M1 -> Y (skips M2)
	return scalar pse_m1_via_m2 = m110 - mmid   // A -> M1 -> M2 -> Y

	///////////////// Step 4 — Effects of treatment on each mediator, and
	///////////////// the mediator-to-mediator effect

	qui sum M1_1
	scalar mean_M1_1 = r(mean)
	qui sum M1_0
	scalar mean_M1_0 = r(mean)
	qui sum M2_A1_M1_1
	scalar mean_M2_A1_M1_1 = r(mean)
	qui sum M2_A0_M1_0
	scalar mean_M2_A0_M1_0 = r(mean)
	qui sum M2_A0_M1_1
	scalar mean_M2_A0_M1_1 = r(mean)

	* Total effect of treatment on M1 itself (mediator-scale, not on Y)
	return scalar ate_m1 = mean_M1_1 - mean_M1_0

	* Total effect of treatment on M2 itself (mediator-scale, not on Y) —
	* combines A's direct edge into M2 AND the part that arrives via M1
	return scalar ate_m2 = mean_M2_A1_M1_1 - mean_M2_A0_M1_0

	* M1's own causal effect on M2 (mediator-to-mediator, M2-scale) —
	* A's direct edge into M2 held at 0 in both terms, so this isolates
	* only the M1 -> M2 edge. NOT a path-specific effect on Y.
	return scalar m1_on_m2 = mean_M2_A0_M1_1 - mean_M2_A0_M1_0

	///////////////// Clean up before the next bootstrap replication

	capture drop M1_1 M1_0 M2_A0_M1_0 M2_A0_M1_1 M2_A1_M1_1 ///
	             Y_000 Y_100 Y_mid Y_110 Y_111 _est_outcome _est_mediator1 _est_mediator2

end

* Run the bootstrap across 1000 replications
set seed `seed'
bootstrap Direct=r(direct) PSE_M1=r(pse_m1) PSE_M2=r(pse_m2) Total=r(total) ///
	PSE_M1_direct=r(pse_m1_direct) PSE_M1_via_M2=r(pse_m1_via_m2) ///
	ATE_M1=r(ate_m1) ATE_M2=r(ate_m2) M1_on_M2=r(m1_on_m2), ///
	reps(`sims') nodots: run_gformula_paths

* Display the percentile-based bootstrap 95% Confidence Intervals
noi estat bootstrap, percentile

savebs "S5" "Edge g-formula: all paths with interaction (manual)"
savebs "S4" "Edge g-formula: main paths with interaction (manual)" "" "Direct PSE_M1 PSE_M2 Total"
savesheet "S5"
savesheet "S4"



*===============================================================================
* PROGRAM:       Stata medeff
*===============================================================================
* Not supported


*===============================================================================
* PROGRAM:       Stata gformula (mediation syntax) 
*===============================================================================

* Not supported



} // close 3.2


******************************************************************************************************************************************************
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
**** 4 - SUPPLEMENTARY BOX S6: mediation accounting for mediator-outcome confounding (manual)
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
******************************************************************************************************************************************************

if `run' == 4 {

*===============================================================================
* REFERENCE:    Supplementary Box S6 (considering  mediator-outcome confounding)
* METHOD:       Parametric Time-Fixed G-Formula Algorithm for causal mediation
* ESTIMAND:     Causal mediation effects of bDMARDs on BASFI at 6 Months
*===============================================================================

********* Return to original variables for single mediator setting:

**** Pre-treatment baseline confounders 
global W="age sex comorbbin mny asasmri hla pertvt1 ibdbl emmtvt1 comedtvt1 asdastotalt1 basfitotalt1"

**** Treatment (binary)
global A="bionew"

**** Mediator (continuous)
global M="asdastotalt2"

**** Mediator-outcome confounders (MOC) (all binary)
global MOCa ="pertvt2"
global MOCb ="emmtvt2"
global MOCc ="comedtvt2"

**** Outcome (continuous)
global Y="basfitotalt2"

order id $W $A $M $Y 


*===============================================================================
* METHOD:       Manual calculation
*===============================================================================


capture program drop run_gformula
program define run_gformula, rclass // Program for confidence interval (ignore if only interested in point estimate)

///////////////// Step 1 — Model the Observed Data

**** Outcome model with interaction (mandatory: include mediator only if later simulated)
gen inter=$M*$A // This is OLS interaction and not the same as in R mediation (simulation) - both will likely give the same answer 

regress $Y $M $A $MOCa $MOCb $MOCc $W inter
** Save coeficients from the outcome model in a vector
estimates store outcome

**** Outcome model without interaction (mandatory: include mediator only if later simulated)
regress $Y $M $A $MOCa $MOCb $MOCc $W
** Save coeficients from the outcome model in a vector
estimates store outcome

* Note: run the desired outcome model (with or without interaction) last

**** Post-treatment mediator-outcome confounder (MOC) models(mandatory)
logit $MOCa $A $W
estimates store moca
logit $MOCb $A $W
estimates store mocb
logit $MOCc $A $W
estimates store mocc

**** Mediator model (optional: simulation can run without simulating the mediator, mandatory only if mediator included in outcome model)
regress $M $A $MOCa $MOCb $MOCc $W
** Save coeficients from the Mediator model in a vector
estimates store mediator

**** Treatment model (optional: only for diagnostics)
logit $A $W
** Save coeficients from the Mediator model in a vector
estimates store treatment


///////////////// Step 2 — Monte Carlo Simulation (adjust for confounding)

/////////// Step 2.1.  Simulate counterfactual post-treatment mediator-outcome confounder (MOC) at 6 months

****** MOC
foreach x in a b c {
* Restore the MOC model coefficients
estimates restore moc`x'

* Generate MOC`x'1 (Counterfactual MOC at 6 months if everyone were treated, bionew = 1)
gen MOC`x'1 = _b[_cons] + ///
              _b[bionew]*1 + ///
              _b[age]*age + ///
              _b[sex]*sex + ///
              _b[comorbbin]*comorbbin + ///
              _b[mny]*mny + ///
              _b[asasmri]*asasmri + ///
              _b[hla]*hla + ///
              _b[pertvt1]*pertvt1 + ///
              _b[ibdbl]*ibdbl + ///
              _b[emmtvt1]*emmtvt1 + ///
              _b[comedtvt1]*comedtvt1 + ///
              _b[asdastotalt1]*asdastotalt1 + ///
              _b[basfitotalt1]*basfitotalt1 

* Generate MOC`x'0 (Counterfactual MOC at 6 months if everyone were untreated, bionew = 0)
gen MOC`x'0 = _b[_cons] + ///
              _b[bionew]*0 + ///
              _b[age]*age + ///
              _b[sex]*sex + ///
              _b[comorbbin]*comorbbin + ///
              _b[mny]*mny + ///
              _b[asasmri]*asasmri + ///
              _b[hla]*hla + ///
              _b[pertvt1]*pertvt1 + ///
              _b[ibdbl]*ibdbl + ///
              _b[emmtvt1]*emmtvt1 + ///
              _b[comedtvt1]*comedtvt1 + ///
              _b[asdastotalt1]*asdastotalt1 + ///
              _b[basfitotalt1]*basfitotalt1 


replace MOC`x'1= invlogit(MOC`x'1)
*replace MOC`x'1= rbinomial(1, MOC`x'1) // random error added (remove this line to use linear prediction and get stable point estimates)

replace MOC`x'0= invlogit(MOC`x'0)
*replace MOC`x'0= rbinomial(1, MOC`x'0) // random error added (remove this line to use linear prediction and get stable point estimates)


}


/////////// Step 2.2. Simulate counterfactual mediator (ASDAS) at 6 months (differs from BOX S2 because include MOC) 

* Restore the mediator model coefficients
estimates restore mediator

* Generate M1 (Counterfactual ASDAS at 6 months if everyone were treated, bionew = 1 and with all MOC at values by bionew=1)
gen M1 = _b[_cons] + ///
		 _b[pertvt2]*MOCa1 + ///
		 _b[emmtvt2]*MOCb1 + ///
		 _b[comedtvt2]*MOCc1 + ///
         _b[bionew]*1 + ///
         _b[age]*age + ///
         _b[sex]*sex + ///
         _b[comorbbin]*comorbbin + ///
         _b[mny]*mny + ///
         _b[asasmri]*asasmri + ///
         _b[hla]*hla + ///
         _b[pertvt1]*pertvt1 + ///
         _b[ibdbl]*ibdbl + ///
         _b[emmtvt1]*emmtvt1 + ///
         _b[comedtvt1]*comedtvt1 + ///
         _b[asdastotalt1]*asdastotalt1 + ///
         _b[basfitotalt1]*basfitotalt1 

* Generate M0 (Counterfactual ASDAS at 6 months if everyone were untreated, bionew = 0 and with all MOC at values by bionew=0)
gen M0 = _b[_cons] + ///
		 _b[pertvt2]*MOCa0 + ///
		 _b[emmtvt2]*MOCb0 + ///
		 _b[comedtvt2]*MOCc0 + ///
         _b[bionew]*0 + ///
         _b[age]*age + ///
         _b[sex]*sex + ///
         _b[comorbbin]*comorbbin + ///
         _b[mny]*mny + ///
         _b[asasmri]*asasmri + ///
         _b[hla]*hla + ///
         _b[pertvt1]*pertvt1 + ///
         _b[ibdbl]*ibdbl + ///
         _b[emmtvt1]*emmtvt1 + ///
         _b[comedtvt1]*comedtvt1 + ///
         _b[asdastotalt1]*asdastotalt1 + ///
         _b[basfitotalt1]*basfitotalt1 
		 

/////////// Step 2.2. Simulate counterfactual outcome (BASFI) at 6 months (differs from BOX S2 because include MOC)

* Restore the outcome model coefficients
estimates restore outcome

* Generate Y1M1 (Counterfactual BASFI if everyone were treated, with ASDAS and MOC at its treated value M1 and MOC1)
gen Y1M1 = _b[_cons] + ///
           _b[asdastotalt2]*M1 + ///
		   _b[pertvt2]*MOCa1 + ///
		   _b[emmtvt2]*MOCb1 + ///
		   _b[comedtvt2]*MOCc1 + ///
           _b[bionew]*1 + ///
           _b[age]*age + ///
           _b[sex]*sex + ///
           _b[comorbbin]*comorbbin + ///
           _b[mny]*mny + ///
           _b[asasmri]*asasmri + ///
           _b[hla]*hla + ///
           _b[pertvt1]*pertvt1 + ///
           _b[ibdbl]*ibdbl + ///
           _b[emmtvt1]*emmtvt1 + ///
           _b[comedtvt1]*comedtvt1 + ///
           _b[asdastotalt1]*asdastotalt1 + ///
           _b[basfitotalt1]*basfitotalt1 //+ _b[inter]*(M1*1) // Delete last term if no interaction is intended.

* Generate Y1M0 (Counterfactual BASFI if everyone were treated, with ASDAS and MOC at its untreated value M0 and MOC0)
gen Y1M0 = _b[_cons] + ///
           _b[asdastotalt2]*M0 + ///
		   _b[pertvt2]*MOCa0 + ///
		   _b[emmtvt2]*MOCb0 + ///
		   _b[comedtvt2]*MOCc0 + ///
           _b[bionew]*1 + ///
           _b[age]*age + ///
           _b[sex]*sex + ///
           _b[comorbbin]*comorbbin + ///
           _b[mny]*mny + ///
           _b[asasmri]*asasmri + ///
           _b[hla]*hla + ///
           _b[pertvt1]*pertvt1 + ///
           _b[ibdbl]*ibdbl + ///
           _b[emmtvt1]*emmtvt1 + ///
           _b[comedtvt1]*comedtvt1 + ///
           _b[asdastotalt1]*asdastotalt1 + ///
           _b[basfitotalt1]*basfitotalt1 //+ _b[inter]*(M0*1) // Delete last term if no interaction is intended.
		
* Generate Y0M1 (Counterfactual BASFI if everyone were untreated, with ASDAS and MOC at its treated value M1 MOC1)
gen Y0M1 = _b[_cons] + ///
           _b[asdastotalt2]*M1 + ///
		   _b[pertvt2]*MOCa1 + ///
		   _b[emmtvt2]*MOCb1 + ///
		   _b[comedtvt2]*MOCc1 + ///
           _b[bionew]*0 + ///
           _b[age]*age + ///
           _b[sex]*sex + ///
           _b[comorbbin]*comorbbin + ///
           _b[mny]*mny + ///
           _b[asasmri]*asasmri + ///
           _b[hla]*hla + ///
           _b[pertvt1]*pertvt1 + ///
           _b[ibdbl]*ibdbl + ///
           _b[emmtvt1]*emmtvt1 + ///
           _b[comedtvt1]*comedtvt1 + ///
           _b[asdastotalt1]*asdastotalt1 + ///
           _b[basfitotalt1]*basfitotalt1 //+ _b[inter]*(M1*0) // Delete last term if no interaction is intended.
		   
* Generate Y0M0 (Counterfactual BASFI if everyone were untreated, with ASDAS and MOC at its untreated value M0 and MOC0)
gen Y0M0 = _b[_cons] + ///
           _b[asdastotalt2]*M0 + ///
		   _b[pertvt2]*MOCa0 + ///
		   _b[emmtvt2]*MOCb0 + ///
		   _b[comedtvt2]*MOCc0 + ///
           _b[bionew]*0 + ///
           _b[age]*age + ///
           _b[sex]*sex + ///
           _b[comorbbin]*comorbbin + ///
           _b[mny]*mny + ///
           _b[asasmri]*asasmri + ///
           _b[hla]*hla + ///
           _b[pertvt1]*pertvt1 + ///
           _b[ibdbl]*ibdbl + ///
           _b[emmtvt1]*emmtvt1 + ///
           _b[comedtvt1]*comedtvt1 + ///
           _b[asdastotalt1]*asdastotalt1 + ///
           _b[basfitotalt1]*basfitotalt1 //+ _b[inter]*(M0*0) // Delete last term if no interaction is intended.

/////////// Step 3 — Mediation effects (equal to Supplementaty Box S2)

* Calculate mean PO accross all patients and save in scalars for bootstraping

qui sum Y1M1
scalar mean_Y1M1 = r(mean)
 
qui sum Y1M0
scalar mean_Y1M0 = r(mean)

qui sum Y0M1
scalar mean_Y0M1 = r(mean)

qui sum Y0M0
scalar mean_Y0M0 = r(mean)


* Calculate mediation effects accross all patients and save in scalars for bootstraping

gen NIE = Y0M1-Y0M0
qui sum NIE
return scalar mean_NIE = r(mean)
    
gen TIE = Y1M1-Y1M0
qui sum TIE
return scalar mean_TIE = r(mean)

gen AIE = (NIE+TIE)/2
qui sum AIE
return scalar mean_AIE = r(mean)
	
gen NDE = Y1M0-Y0M0
qui sum NDE
return scalar mean_NDE = r(mean)	
 
gen TDE = Y1M1-Y0M1
qui sum TDE
return scalar mean_TDE = r(mean)	

gen ADE = (NDE+TDE)/2
qui sum ADE
return scalar mean_ADE = r(mean)


* Calculate the marginal total treatment effect (ATE) and save in scalar for bootstraping

scalar total_effect = mean_Y1M1 - mean_Y0M0
    
/////////// Return values for the bootstraping
return scalar ate = total_effect
return scalar mean11 = mean_Y1M1
return scalar mean10 = mean_Y1M0
return scalar mean01 = mean_Y0M1
return scalar mean00 = mean_Y0M0


/////////// Optional step. Predict treatment, mediator and outcome under observed data (only for diagnostics)

* Observed treatment
qui sum $A, meanonly 
return scalar mean_$A = r(mean)

* Generate A (predicted treatment under observed data)
estimates restore treatment

gen Ap = _b[_cons] + ///
         _b[age]*age + ///
         _b[sex]*sex + ///
         _b[comorbbin]*comorbbin + ///
         _b[mny]*mny + ///
         _b[asasmri]*asasmri + ///
         _b[hla]*hla + ///
         _b[pertvt1]*pertvt1 + ///
         _b[ibdbl]*ibdbl + ///
         _b[emmtvt1]*emmtvt1 + ///
         _b[comedtvt1]*comedtvt1 + ///
         _b[asdastotalt1]*asdastotalt1 + ///
         _b[basfitotalt1]*basfitotalt1
replace Ap= invlogit(Ap)
*replace Ap= rbinomial(1, Ap) // random error added (remove this line to use linear prediction and get stable point estimates)

qui sum Ap, meanonly 
return scalar mean_Ap = r(mean)

* Observed mediator
qui sum $M, meanonly 
return scalar mean_$M = r(mean)

* Predicted MOC under observed treatment
gen MOCap = Ap*MOCa1 + (1-Ap)*MOCa0
gen MOCbp = Ap*MOCb1 + (1-Ap)*MOCb0
gen MOCcp = Ap*MOCc1 + (1-Ap)*MOCc0

* Restore the mediator model coefficients
estimates restore mediator

* Generate M (Predicted ASDAS at 6 months under observed treatment)
gen Mp = _b[_cons] + ///
         _b[pertvt2]*MOCap + ///
         _b[emmtvt2]*MOCbp + ///
         _b[comedtvt2]*MOCcp + ///
         _b[bionew]*Ap + ///
         _b[age]*age + ///
         _b[sex]*sex + ///
         _b[comorbbin]*comorbbin + ///
         _b[mny]*mny + ///
         _b[asasmri]*asasmri + ///
         _b[hla]*hla + ///
         _b[pertvt1]*pertvt1 + ///
         _b[ibdbl]*ibdbl + ///
         _b[emmtvt1]*emmtvt1 + ///
         _b[comedtvt1]*comedtvt1 + ///
         _b[asdastotalt1]*asdastotalt1 + ///
         _b[basfitotalt1]*basfitotalt1

qui sum Mp, meanonly 
return scalar mean_Mp = r(mean)

* Observed outcome
qui sum $Y, meanonly 
return scalar mean_$Y = r(mean)


* Restore the outcome model coeficients
estimates restore outcome

gen Yp   = _b[_cons] + ///
           _b[asdastotalt2]*Mp + ///
           _b[pertvt2]*MOCap + ///
           _b[emmtvt2]*MOCbp + ///
           _b[comedtvt2]*MOCcp + ///
           _b[bionew]*Ap + ///
           _b[age]*age + ///
           _b[sex]*sex + ///
           _b[comorbbin]*comorbbin + ///
           _b[mny]*mny + ///
           _b[asasmri]*asasmri + ///
           _b[hla]*hla + ///
           _b[pertvt1]*pertvt1 + ///
           _b[ibdbl]*ibdbl + ///
           _b[emmtvt1]*emmtvt1 + ///
           _b[comedtvt1]*comedtvt1 + ///
           _b[asdastotalt1]*asdastotalt1 + ///
           _b[basfitotalt1]*basfitotalt1 //+ _b[inter]*(Mp*Ap) // Delete last term if no interaction is intended.

qui sum Yp, meanonly 
return scalar mean_Yp = r(mean)

/////////// Clean up temporary variables before the next replication

capture drop M1 M0 Y1M1 Y1M0 Y0M1 Y0M0 MOCa1 MOCa0 MOCb1 MOCb0 MOCc1 MOCc0 MOCap MOCbp MOCcp NIE TIE AIE NDE TDE ADE Ap Mp Yp _est_mediator _est_outcome _est_treatment inter

end

* Treatment-mediator interaction, for the Table S10 footnote (diagnostic only: the
* outcome model used is the one stored last inside the program above)
qui gen interdx = $M*$A
qui regress $Y $M $A $MOCa $MOCb $MOCc $W interdx
qui test interdx
local pinterMOC = r(p)
noi di as text "Treatment-mediator interaction ($A x $M): p = " as result %6.4f `pinterMOC'
savediag "S10" "G-formula with MOC (manual)" "p interaction AxM" `pinterMOC'
qui drop interdx

* Run the bootstrap across 1000 replications (second line of bootstrap comman contains only diagnostics and is optional
set seed `seed'
bootstrap ATE=r(ate) Mean_Y1M1=r(mean11) Mean_Y1M0=r(mean10) Mean_Y0M1=r(mean01) Mean_Y0M0=r(mean00) ///
		  NIE=r(mean_NIE) TIE=r(mean_TIE) AIE=r(mean_AIE) NDE=r(mean_NDE) TDE=r(mean_TDE) ADE=r(mean_ADE) ///
          Mean_$A=r(mean_$A) Mean_Ap=r(mean_Ap) Mean_$M=r(mean_$M) Mean_Mp=r(mean_Mp) Mean_$Y=r(mean_$Y) Mean_Yp=r(mean_Yp) ///
		  , reps(`sims') nodots: run_gformula

* Display the percentile-based bootstrap 95% Confidence Intervals
noi estat bootstrap, percentile

savebs "S10" "G-formula with MOC (manual)"
savesheet "S10"



} // close 4


******************************************************************************************************************************************************
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
**** 4.1 - SUPPLEMENTARY BOX S6: mediation accounting for mediator-outcome confounding (Stata gformula, mediation syntax)
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
******************************************************************************************************************************************************

if `run' == 4.1 {

*===============================================================================
* PROGRAM:       Stata medeff
*===============================================================================

* Not supported

*===============================================================================
* PROGRAM:       Stata gformula (mediation syntax) without treatment-mediator interaction (not supported)
*===============================================================================

* Stata gformula simulates the mediator, the MOC and the outcome instead of using their
* predicted means, so its point estimates carry Monte Carlo error on top of the
* sampling error. Two options bring it close to the manual code:
*   minsim:                uses the predicted mean of the outcome instead of drawing it
*   moreMC simulations():  lets the Monte Carlo sample exceed the analysis sample, so
*                          the remaining noise from the mediator and MOC draws shrinks
*                          with the square root of the Monte Carlo sample size

capture matrix drop b se ci_normal ci_percentile ci_bc ci_bca

noi gformula $Y $M $MOCa $MOCb $MOCc $A $W, ///
mediation ex($A) mediator($M) out($Y) post_confs($MOCa $MOCb $MOCc) ///
eq($M: $A $MOCa $MOCb $MOCc $W, $Y: $A $M $MOCa $MOCb $MOCc $W, $MOCa: $A $W, $MOCb: $A $W, $MOCc: $A $W) /// 
com($M:regress, $Y:regress, $MOCa:logit, $MOCb:logit, $MOCc:logit)   ///
obe base_confs($W) ///
minsim moreMC simulations(`mcsims') ///
seed(`seed') samples(`sims') all

savegf "S10" "gformula"
savesheet "S10"




} // close 4.1


******************************************************************************************************************************************************
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
**** 5 - SUPPLEMENTARY BOX S3: marginal structural model with IPTW (manual)
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
******************************************************************************************************************************************************

if `run' == 5 {

*===============================================================================
* REFERENCE:    Supplementary Box S3 
* METHOD:       Marginal structural model (MSM)
* ESTIMAND:     Average Total Effect (ATE) of bDMARDs on BASFI at 6 Months
*===============================================================================

*===============================================================================
* METHOD:       Manual calculation
*===============================================================================

///////////////// Step 1 — Model the Observed Data

**** Treatment model for the denominator
logit $A $W 
estimates store denominator

**** Treatment model for the numerator (intercept only: or marginal probability of treatment)
logit $A 
estimates store numerator

///////////////// Step 2 — Weighting (adjust for confounding)

///////////////// Step 2.1. Calculate the denominator of the IPTW 
estimates restore denominator

** Predict probability of treatment given covariates (propensity score)
predict ps if e(sample), pr

** PS diagnostics
noi sum ps 
local ps_mean=round(r(mean), 0.01) 
local ps_sd=round(r(sd), 0.01) 
local ps_min=round(r(min), 0.01) 
local ps_max=round(r(max), 0.01) 

** Calculate the denominator for treated and untreated
gen denom=ps*$A+(1-ps)*(1-$A) 

** Calculate the unstabilized IPTW (uIPTW) for treated and untreated
qui gen unstabweight= 1/denom 

** uIPTW diagnostics
noi sum unstabweight 
local unstabweightmean=round(r(mean), 0.01) 
local unstabweightsd=round(r(sd), 0.01) 
local unstabweightmin=round(r(min), 0.01) 
local unstabweightmax=round(r(max), 0.01) 


///////////////// Step 2.2. Calculate the numerator of the IPTW 
estimates restore numerator

** Predict marginal probability of treatment
predict num_ps if e(sample), pr 

** Calculate the numerator for treated and untreated (can be calculated non-paramtetric: just tabulate $A)
gen num=num_ps*$A+(1-num_ps)*(1-$A) 

** Calculate the stabilized IPTW (uIPTW) for treated and untreated
gen stabweight= num/denom 

** sIPTW diagnostics
noi sum stabweight 
local stabweightmean=round(r(mean), 0.01) 
local stabweightsd=round(r(sd), 0.01) 
local stabweightmin=round(r(min), 0.01) 
local stabweightmax=round(r(max), 0.01) 

///////////////// Step 3 and — Total Treatment Effect and 95% CI: MSM

noi regress $Y $A [pw=stabweight], robust

///////////////// e-value

* SD of the observed outcome
qui sum $Y
scalar SD_Y = r(sd)

* Compute the e-value
local rr = exp(0.91*abs(_b[$A]/scalar(SD_Y)))
scalar evalue_point = `rr' + sqrt(`rr'*(`rr'-1))
noi di "E-value (point estimate): " %6.3f scalar(evalue_point)

savest "S6" "MSM with sIPTW (manual)" "ATE" `=_b[$A]' `=_se[$A]' `=_b[$A]-1.96*_se[$A]' `=_b[$A]+1.96*_se[$A]' `=scalar(evalue_point)' . $anaN
savediag "S6" "MSM with sIPTW (manual)" "sIPTW mean" `stabweightmean'
savediag "S6" "MSM with sIPTW (manual)" "sIPTW min" `stabweightmin'
savediag "S6" "MSM with sIPTW (manual)" "sIPTW max" `stabweightmax'
savediag "S6" "MSM with sIPTW (manual)" "sIPTW SD" `stabweightsd'


///////////////// PS and IPTW Diagnostics without trimming or truncation

matrix DIAG = J(3, 4, .)

matrix DIAG[1,1] = `ps_mean'
matrix DIAG[1,2] = `ps_sd'
matrix DIAG[1,3] = `ps_min'
matrix DIAG[1,4] = `ps_max'

matrix DIAG[2,1] = `unstabweightmean'
matrix DIAG[2,2] = `unstabweightsd'
matrix DIAG[2,3] = `unstabweightmin'
matrix DIAG[2,4] = `unstabweightmax'

matrix DIAG[3,1] = `stabweightmean'
matrix DIAG[3,2] = `stabweightsd'
matrix DIAG[3,3] = `stabweightmin'
matrix DIAG[3,4] = `stabweightmax'

matrix rownames DIAG = "Propensity score" "uIPTW" "sIPTW"
matrix colnames DIAG = Mean SD Min Max

noi matrix list DIAG, format(%9.3f) title("IPTW Diagnostics")


///////////////// Positivity diagnostics — Area of Common Support (ACS): trimming

* (overlap in PS between treated and unctreated)

///////////////// Step 1 — PS range by treatment group
bysort $A: egen max_ps = max(ps)
bysort $A: egen min_ps = min(ps)

* Store min/max PS separately for treated (A=1) and untreated (A=0)
sum max_ps if $A == 1, meanonly
local max_ps1 = r(mean)
sum max_ps if $A == 0, meanonly
local max_ps0 = r(mean)

sum min_ps if $A == 1, meanonly
local min_ps1 = r(mean)
sum min_ps if $A == 0, meanonly
local min_ps0 = r(mean)

* Define ACS bounds (overlap region)
* Upper bound = min of the two maxima
* Lower bound = max of the two minima
local acsmax = min(`max_ps1', `max_ps0') + 0.0000001
local acsmin = max(`min_ps1', `min_ps0') - 0.0000001

noi di as text "ACS lower bound: `acsmin'"
noi di as text "ACS upper bound: `acsmax'"

* Flag observations within/outside ACS
gen acs = .
replace acs = 1 if ps >= `acsmin' & ps <= `acsmax' & ps !=.
replace acs = 0 if (ps < `acsmin' | ps > `acsmax') & ps !=.

* Count observations outside ACS
sum acs
local n_outside = round(r(N) - (r(N) * r(mean)), 1)
local pct_outside = round((1 - r(mean)) * 100, 0.1)

noi di as text "Observations outside ACS: `n_outside' (`pct_outside'%)"

* Clean up intermediate variables
drop max_ps min_ps

* Total Treatment Effect: MSM within the ACS

noi regress $Y $A [pw=stabweight] if acs==1, robust

savest "S6" "MSM with sIPTW within ACS (manual)" "ATE" `=_b[$A]' `=_se[$A]' `=_b[$A]-1.96*_se[$A]' `=_b[$A]+1.96*_se[$A]' . . $anaN
savediag "S6" "MSM with sIPTW within ACS (manual)" "N outside ACS" `n_outside'
savediag "S6" "MSM with sIPTW within ACS (manual)" "pct outside ACS" `pct_outside'


///////////////// Positivity diagnostics: truncation
foreach x in unstabweight stabweight {
    sum `x', detail
    scalar p1_`x'  = r(p1)
    scalar p99_`x' = r(p99)
}

foreach x in unstabweight stabweight {
    gen `x'_trunc = `x'
    replace `x'_trunc = p1_`x'  if `x' < p1_`x'  & `x' != .
    replace `x'_trunc = p99_`x' if `x' > p99_`x' & `x' != .
    
    * Count truncated observations
    count if `x' < p1_`x'  & `x' != .
    local n_low_`x' = r(N)
    count if `x' > p99_`x' & `x' != .
    local n_high_`x' = r(N)
    local n_trunc_`x' = `n_low_`x'' + `n_high_`x''
    
    * Truncation thresholds
    local p1_`x'  = round(scalar(p1_`x'),  0.001)
    local p99_`x' = round(scalar(p99_`x'), 0.001)

    sum `x'_trunc
    local `x'_truncmean = round(r(mean), 0.01)
    local `x'_truncsd   = round(r(sd),   0.01)
    local `x'_truncmin  = round(r(min),  0.01)
    local `x'_truncmax  = round(r(max),  0.01)
}

noi regress $Y $A [pw=stabweight_trunc], vce(robust)

savest "S6" "MSM with sIPTW truncated P1-P99 (manual)" "ATE" `=_b[$A]' `=_se[$A]' `=_b[$A]-1.96*_se[$A]' `=_b[$A]+1.96*_se[$A]' . . $anaN
savediag "S6" "MSM with sIPTW truncated P1-P99 (manual)" "N truncated sIPTW" `n_trunc_stabweight'
savesheet "S6"

matrix TRUNC = J(2, 6, .)

matrix TRUNC[1,1] = `p1_unstabweight'
matrix TRUNC[1,2] = `p99_unstabweight'
matrix TRUNC[1,3] = `n_low_unstabweight'
matrix TRUNC[1,4] = `n_high_unstabweight'
matrix TRUNC[1,5] = `n_trunc_unstabweight'
matrix TRUNC[1,6] = `unstabweight_truncmean'

matrix TRUNC[2,1] = `p1_stabweight'
matrix TRUNC[2,2] = `p99_stabweight'
matrix TRUNC[2,3] = `n_low_stabweight'
matrix TRUNC[2,4] = `n_high_stabweight'
matrix TRUNC[2,5] = `n_trunc_stabweight'
matrix TRUNC[2,6] = `stabweight_truncmean'

matrix rownames TRUNC = "uIPTW" "sIPTW"
matrix colnames TRUNC = "P1 (lower)" "P99 (upper)" "N low tail" "N high tail" "N truncated" "Mean (post-trunc)"

noi matrix list TRUNC, format(%9.3f) title("Truncation Diagnostics (1st-99th percentile)")


///////////////// Balance before and after weighting (Supplementary Table S8)

* Means of the baseline confounders in the treated and the untreated, in the
* original population and in the pseudo-population created by the sIPTW.
* SMD as in Jackson (Epidemiology 2016): the unweighted and the weighted
* difference in means are divided by the same unweighted pooled SD,
* sqrt((var treated + var untreated)/2), with var = p(1-p) for binary variables.

qui do "${programs}msmbalanceV1.do" // Run supporting program

global balance "basfitotalt1 asdastotalt1 age sex comorbbin mny asasmri hla pertvt1 ibdbl emmtvt1 comedtvt1"

noi msmbalanceV1 $balance, treatvar($A) weight(stabweight) ///
	label("Supplementary Table S8: original population vs pseudo-population (sIPTW)") ///
	save(yes) folder("${path}Tables") filename(Supplementary_Table_S8)

* Group sizes for the column headers: number of patients in the original
* population and sum of the sIPTW in the pseudo-population

qui count if $A == 1
local n_treated = r(N)

qui count if $A == 0
local n_untreated = r(N)

local n_total = `n_treated' + `n_untreated'

qui sum stabweight if $A == 1
local w_treated = r(sum)

qui sum stabweight if $A == 0
local w_untreated = r(sum)

qui sum stabweight
local w_total = r(sum)
local w_mean = r(mean)

noi di as text "Original population:  treated " as result `n_treated' as text ", untreated " as result `n_untreated'
noi di as text "Pseudo-population:    treated " as result %6.1f `w_treated' as text ", untreated " as result %6.1f `w_untreated' ///
	as text ", total " as result %6.1f `w_total' as text " (mean sIPTW " as result %5.3f `w_mean' as text ")"

* Write the group sizes below the balance matrix, in the same sheet
putexcel set "${path}Tables/Supplementary_Table_S8.xlsx", sheet("Supplementary_Table_S8") modify
putexcel A15 = "Group sizes"
putexcel B15 = "Treated"
putexcel C15 = "Untreated"
putexcel D15 = "Total"
putexcel A16 = "N original population"
putexcel B16 = `n_treated'
putexcel C16 = `n_untreated'
putexcel D16 = `n_total'
putexcel A17 = "N pseudo-population (sum of sIPTW)"
putexcel B17 = `w_treated'
putexcel C17 = `w_untreated'
putexcel D17 = `w_total'
putexcel A18 = "Mean sIPTW"
putexcel B18 = `w_mean'
putexcel clear



} // close 5


******************************************************************************************************************************************************
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
**** 5.1 - SUPPLEMENTARY BOX S3: IPTW (Stata teffects ipw)
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
******************************************************************************************************************************************************

if `run' == 5.1 {

*===============================================================================
* PROGRAM:       Stata teffects 
*===============================================================================

noi teffects ipw ($Y) ($A $W, logit), pomeans // potential outcomes 
noi teffects ipw ($Y) ($A $W, logit), ate // total treatment effect

savete "S6" "IPTW (teffects)" "ATE"
savesheet "S6"

* Ony uIPTWs are calculated with teffects (equal to sIPTW here because the MSM is saturated)
* SE are lower with teffects (the program uses influence-function-based SE that accounts for the uncertainty in the propensity score estimation step)




} // close 5.1


******************************************************************************************************************************************************
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
**** 6 - SUPPLEMENTARY BOX S4: targeted maximum likelihood estimation, total effect (manual)
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
******************************************************************************************************************************************************

if `run' == 6 {

*===============================================================================
* REFERENCE:    Supplementary Box S4 
* METHOD:       Targeted Maximum Likelihood Estimation for Total Treatment Effect
*               (ignoring mediator and MOC)
* ESTIMAND:     Average Total Effect (ATE) of bDMARDs on BASFI at 6 Months
*===============================================================================

*===============================================================================
* PROGRAM:       Manual calculation
*===============================================================================

///////////////// Step 1 — Model the Observed Data

**** Outcome model
qui regress $Y $A $W
estimates store outcome

**** Treatment model
qui logit $A $W
estimates store treatment

///////////////// Step 2 — Targeting (adjust for confounding)

///////////////// Step 2.1. Simulate counterfactual outcome (BASFI) at 6 months

* Restore the outcome model coefficients
estimates restore outcome

* Generate Y1 (Counterfactual BASFI at 6 months if everyone were treated, bionew = 1)
gen Y1 = _b[_cons] + ///
         _b[bionew]*1 + ///
         _b[age]*age + ///
         _b[sex]*sex + ///
         _b[comorbbin]*comorbbin + ///
         _b[mny]*mny + ///
         _b[asasmri]*asasmri + ///
         _b[hla]*hla + ///
         _b[pertvt1]*pertvt1 + ///
         _b[ibdbl]*ibdbl + ///
         _b[emmtvt1]*emmtvt1 + ///
         _b[comedtvt1]*comedtvt1 + ///
         _b[asdastotalt1]*asdastotalt1 + ///
         _b[basfitotalt1]*basfitotalt1

* Generate Y0 (Counterfactual BASFI at 6 months if everyone were untreated, bionew = 0)
gen Y0 = _b[_cons] + ///
         _b[bionew]*0 + ///
         _b[age]*age + ///
         _b[sex]*sex + ///
         _b[comorbbin]*comorbbin + ///
         _b[mny]*mny + ///
         _b[asasmri]*asasmri + ///
         _b[hla]*hla + ///
         _b[pertvt1]*pertvt1 + ///
         _b[ibdbl]*ibdbl + ///
         _b[emmtvt1]*emmtvt1 + ///
         _b[comedtvt1]*comedtvt1 + ///
         _b[asdastotalt1]*asdastotalt1 + ///
         _b[basfitotalt1]*basfitotalt1

* Generate Yobs (Predicted BASFI at 6 months under observed treatment)
gen Yobs = _b[_cons] + ///
           _b[bionew]*$A + ///
           _b[age]*age + ///
           _b[sex]*sex + ///
           _b[comorbbin]*comorbbin + ///
           _b[mny]*mny + ///
           _b[asasmri]*asasmri + ///
           _b[hla]*hla + ///
           _b[pertvt1]*pertvt1 + ///
           _b[ibdbl]*ibdbl + ///
           _b[emmtvt1]*emmtvt1 + ///
           _b[comedtvt1]*comedtvt1 + ///
           _b[asdastotalt1]*asdastotalt1 + ///
           _b[basfitotalt1]*basfitotalt1

///////////////// Step 2.2. Predict treatment (bDMARDs) at 6 months

* Restore the treatment model coefficients
estimates restore treatment

* Predict propensity score (PS) and g (PS or 1-PS according to observed treatment)
predict ps_tmle if e(sample), pr
gen g = ps_tmle*$A + (1-ps_tmle)*(1-$A)

///////////////// Step 2.3. Calculate the "clever covariate" (H)

* H = A/PS - (1-A)/(1-PS)
* For the treated (A=1): H = 1/PS
* For the untreated (A=0): H = -1/(1-PS)
gen H = $A/ps_tmle - (1-$A)/(1-ps_tmle)

///////////////// Step 2.4. Targeting step

* Fluctuation model: regress (Y - Yobs) on H with no intercept to fix Yobs coefficient to 1
gen Y_resid = $Y - Yobs
quietly regress Y_resid H, nocons
scalar epsilon = _b[H]

* Compute targeted predictions
gen Y1_target   = Y1   + scalar(epsilon) / ps_tmle
gen Y0_target   = Y0   - scalar(epsilon) / (1-ps_tmle)
gen Yobs_target = Yobs + scalar(epsilon) * H

///////////////// Step 3 — Total Treatment Effect

qui sum Y1_target
scalar mean_Y1 = r(mean)

qui sum Y0_target
scalar mean_Y0 = r(mean)

scalar ATE_tmle = mean_Y1 - mean_Y0

///////////////// Step 4 — 95% Confidence Interval via Efficient Influence Function (EIF)

* EIF for each patient: H*(Yi - Yobs_target) + Y1_target - Y0_target - ATE
gen EIF = H*($Y - Yobs_target) + Y1_target - Y0_target - scalar(ATE_tmle)

qui count
scalar n = r(N)

* Variance of ATE = (1/n²) * sum(EIF²)
gen EIF2 = EIF^2
qui sum EIF2
scalar var_ATE = r(sum) / (scalar(n)^2)
scalar se_ATE  = sqrt(scalar(var_ATE))

scalar lb_ATE = ATE_tmle - 1.96*se_ATE
scalar ub_ATE = ATE_tmle + 1.96*se_ATE

///////////////// Clean up temporary variables

drop Y1 Y0 Yobs g H Y_resid Y1_target Y0_target Yobs_target EIF2

///////////////// Results table

matrix TMLE = J(1, 4, .)
matrix TMLE[1,1] = scalar(mean_Y1)
matrix TMLE[1,2] = scalar(mean_Y0)
matrix TMLE[1,3] = scalar(ATE_tmle)
matrix TMLE[1,4] = scalar(se_ATE)

matrix rownames TMLE = "TMLE (manual)"
matrix colnames TMLE = "E[Y1]" "E[Y0]" "ATE" "SE"

noi matrix list TMLE, format(%9.3f) title("TMLE — Total Treatment Effect of bDMARDs on BASFI at 6 Months")

noi di as text _newline "95% CI (EIF-based): [" %6.3f scalar(lb_ATE) ", " %6.3f scalar(ub_ATE) "]"


///////////////// TMLE diagnostics: % g-truncation and SD of influence curve

* g-truncation: % of patients whose estimated PS hit a bound (commonly [0.025, 0.975] or similar)
* If you did not explicitly truncate ps_tmle, define the bound used as the truncation threshold
local g_trunc_lb = 0.025
local g_trunc_ub = 0.975

gen g_truncated = (ps_tmle < `g_trunc_lb' | ps_tmle > `g_trunc_ub')
qui sum g_truncated, meanonly
local pct_gtrunc = round(100*r(mean), 0.1)

* SD of the influence curve (EIF) — already generated in Step 4 of the manual TMLE code
qui sum EIF
local sd_IC = round(r(sd), 0.001)

noi di as text _newline "TMLE diagnostics"
noi di as text "% g-truncated: `pct_gtrunc'%"
noi di as text "SD of influence curve: `sd_IC'"

savediag "S7" "TMLE (manual)" "pct g-truncated" `pct_gtrunc'
savediag "S7" "TMLE (manual)" "SD influence curve" `sd_IC'
savediag "S7" "TMLE (manual)" "Mean Y1" `=scalar(mean_Y1)'
savediag "S7" "TMLE (manual)" "Mean Y0" `=scalar(mean_Y0)'

drop ps_tmle EIF g_truncated


///////////////// e-value

* SD of the observed outcome (added: it was being taken from the Box S1 block)
qui sum $Y
scalar SD_Y = r(sd)
local rr = exp(0.91*abs(ATE_tmle/SD_Y))
scalar evalue_point = `rr' + sqrt(`rr'*(`rr'-1))
noi di "E-value (point estimate): " %6.3f scalar(evalue_point)

savest "S7" "TMLE (manual)" "ATE" `=scalar(ATE_tmle)' `=scalar(se_ATE)' `=scalar(lb_ATE)' `=scalar(ub_ATE)' `=scalar(evalue_point)' . $anaN
savesheet "S7"




} // close 6


******************************************************************************************************************************************************
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
**** 6.1 - SUPPLEMENTARY BOX S4: augmented IPW, total effect (Stata teffects aipw)
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
******************************************************************************************************************************************************

if `run' == 6.1 {

*===============================================================================
* PROGRAM:       Stata teffects
*===============================================================================

noi teffects aipw ($Y $W) ($A $W, logit), ate // Augmented IPW (AIPW) estimator: similar (double robust) but not the same as TMLE

savete "S7" "AIPW (teffects)" "ATE"
savesheet "S7"



} // close 6.1


******************************************************************************************************************************************************
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
**** 7 - SUPPLEMENTARY TABLE S9: mediation within levels of baseline characteristics (manual)
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
******************************************************************************************************************************************************

if `run' == 7 {

*===============================================================================
* REFERENCE:    Supplementary Table S9
* METHOD:       Parametric Time-Fixed G-Formula Algorithm for causal mediation
* ESTIMAND:     Causal mediation effects of bDMARDs on BASFI at 6 Months
*               within levels of baseline characteristics
*===============================================================================


*===============================================================================
* METHOD:       Manual calculation
*===============================================================================


capture program drop run_gformula
program define run_gformula, rclass // Program for confidence interval (ignore if only interested in point estimate)

///////////////// Step 1 — Model the Observed Data

**** Outcome model with interaction (mandatory: include mediator only if later simulated)
gen inter=$M*$A // This is OLS interaction and not the same as in R mediation (simulation) - both will likely give the same answer 

regress $Y $M $A $W inter
** Save coeficients from the outcome model in a vector
estimates store outcome

**** Outcome model without interaction (mandatory: include mediator only if later simulated)
regress $Y $M $A $W
** Save coeficients from the outcome model in a vector
estimates store outcome

* Note: run the desired outcome model (with or without interaction) last

**** Mediator model (optional: simulation can run without simulating the mediator, mandatory only if mediator included in outcome model)
regress $M $A $W
** Save coeficients from the Mediator model in a vector
estimates store mediator

**** Treatment model (optional: only for diagnostics)
logit $A $W
** Save coeficients from the Mediator model in a vector
estimates store treatment

///////////////// Step 2 — Monte Carlo Simulation (adjust for confounding)

/////////// Step 2.1. Simulate counterfactual mediator (ASDAS) at 6 months (Optional, only if included in the outcome model)

* Restore the mediator model coefficients
estimates restore mediator

* Generate M1 (Counterfactual ASDAS at 6 months if everyone were treated, bionew = 1)
gen M1 = _b[_cons] + ///
         _b[bionew]*1 + ///
         _b[age]*age + ///
         _b[sex]*sex + ///
         _b[comorbbin]*comorbbin + ///
         _b[mny]*mny + ///
         _b[asasmri]*asasmri + ///
         _b[hla]*hla + ///
         _b[pertvt1]*pertvt1 + ///
         _b[ibdbl]*ibdbl + ///
         _b[emmtvt1]*emmtvt1 + ///
         _b[comedtvt1]*comedtvt1 + ///
         _b[asdastotalt1]*asdastotalt1 + ///
         _b[basfitotalt1]*basfitotalt1 

* Generate M0 (Counterfactual ASDAS at 6 months if everyone were untreated, bionew = 0)
gen M0 = _b[_cons] + ///
         _b[bionew]*0 + ///
         _b[age]*age + ///
         _b[sex]*sex + ///
         _b[comorbbin]*comorbbin + ///
         _b[mny]*mny + ///
         _b[asasmri]*asasmri + ///
         _b[hla]*hla + ///
         _b[pertvt1]*pertvt1 + ///
         _b[ibdbl]*ibdbl + ///
         _b[emmtvt1]*emmtvt1 + ///
         _b[comedtvt1]*comedtvt1 + ///
         _b[asdastotalt1]*asdastotalt1 + ///
         _b[basfitotalt1]*basfitotalt1 
		 
/////////// Step 2.2. Simulate counterfactual outcome (BASFI) at 6 months (differs from BOX S1 only from this point onward)

* Restore the outcome model coefficients
estimates restore outcome

* Generate Y1M1 (Counterfactual BASFI if everyone were treated, with ASDAS at its treated value M1)
gen Y1M1 = _b[_cons] + ///
           _b[asdastotalt2]*M1 + ///
           _b[bionew]*1 + ///
           _b[age]*age + ///
           _b[sex]*sex + ///
           _b[comorbbin]*comorbbin + ///
           _b[mny]*mny + ///
           _b[asasmri]*asasmri + ///
           _b[hla]*hla + ///
           _b[pertvt1]*pertvt1 + ///
           _b[ibdbl]*ibdbl + ///
           _b[emmtvt1]*emmtvt1 + ///
           _b[comedtvt1]*comedtvt1 + ///
           _b[asdastotalt1]*asdastotalt1 + ///
           _b[basfitotalt1]*basfitotalt1 //+ _b[inter]*(M1*1) // Delete last term if no interaction is intended.

* Generate Y1M0 (Counterfactual BASFI if everyone were treated, with ASDAS at its untreated value M0)
gen Y1M0 = _b[_cons] + ///
           _b[asdastotalt2]*M0 + ///
           _b[bionew]*1 + ///
           _b[age]*age + ///
           _b[sex]*sex + ///
           _b[comorbbin]*comorbbin + ///
           _b[mny]*mny + ///
           _b[asasmri]*asasmri + ///
           _b[hla]*hla + ///
           _b[pertvt1]*pertvt1 + ///
           _b[ibdbl]*ibdbl + ///
           _b[emmtvt1]*emmtvt1 + ///
           _b[comedtvt1]*comedtvt1 + ///
           _b[asdastotalt1]*asdastotalt1 + ///
           _b[basfitotalt1]*basfitotalt1 //+ _b[inter]*(M0*1) // Delete last term if no interaction is intended.
		
* Generate Y0M1 (Counterfactual BASFI if everyone were untreated, with ASDAS at its treated value M1)
gen Y0M1 = _b[_cons] + ///
           _b[asdastotalt2]*M1 + ///
           _b[bionew]*0 + ///
           _b[age]*age + ///
           _b[sex]*sex + ///
           _b[comorbbin]*comorbbin + ///
           _b[mny]*mny + ///
           _b[asasmri]*asasmri + ///
           _b[hla]*hla + ///
           _b[pertvt1]*pertvt1 + ///
           _b[ibdbl]*ibdbl + ///
           _b[emmtvt1]*emmtvt1 + ///
           _b[comedtvt1]*comedtvt1 + ///
           _b[asdastotalt1]*asdastotalt1 + ///
           _b[basfitotalt1]*basfitotalt1 //+ _b[inter]*(M1*0) // Delete last term if no interaction is intended.
		   
* Generate Y0M0 (Counterfactual BASFI if everyone were untreated, with ASDAS at its untreated value M0)
gen Y0M0 = _b[_cons] + ///
           _b[asdastotalt2]*M0 + ///
           _b[bionew]*0 + ///
           _b[age]*age + ///
           _b[sex]*sex + ///
           _b[comorbbin]*comorbbin + ///
           _b[mny]*mny + ///
           _b[asasmri]*asasmri + ///
           _b[hla]*hla + ///
           _b[pertvt1]*pertvt1 + ///
           _b[ibdbl]*ibdbl + ///
           _b[emmtvt1]*emmtvt1 + ///
           _b[comedtvt1]*comedtvt1 + ///
           _b[asdastotalt1]*asdastotalt1 + ///
           _b[basfitotalt1]*basfitotalt1 //+ _b[inter]*(M0*0) // Delete last term if no interaction is intended.

/////////// Step 3 — Mediation effects

* Calculate mean PO accross all patients and save in scalars for bootstraping

qui sum Y1M1
scalar mean_Y1M1 = r(mean)
 
qui sum Y1M0
scalar mean_Y1M0 = r(mean)

qui sum Y0M1
scalar mean_Y0M1 = r(mean)

qui sum Y0M0
scalar mean_Y0M0 = r(mean)


* Calculate mediation effects accross all patients and save in scalars for bootstraping

gen NIE = Y0M1-Y0M0
qui sum NIE
return scalar mean_NIE = r(mean)
    
gen TIE = Y1M1-Y1M0
qui sum TIE
return scalar mean_TIE = r(mean)

gen AIE = (NIE+TIE)/2
qui sum AIE
return scalar mean_AIE = r(mean)
	
gen NDE = Y1M0-Y0M0
qui sum NDE
return scalar mean_NDE = r(mean)	
 
gen TDE = Y1M1-Y0M1
qui sum TDE
return scalar mean_TDE = r(mean)	

gen ADE = (NDE+TDE)/2
qui sum ADE
return scalar mean_ADE = r(mean)


* Calculate the marginal total treatment effect (ATE) and save in scalar for bootstraping

scalar total_effect = mean_Y1M1 - mean_Y0M0
    
/////////// Return values for the bootstraping
return scalar ate = total_effect
return scalar mean11 = mean_Y1M1
return scalar mean10 = mean_Y1M0
return scalar mean01 = mean_Y0M1
return scalar mean00 = mean_Y0M0


/////////// Optional step. Predict treatment, mediator and outcome under observed data (only for diagnostics)

* Observed treatment
qui sum $A, meanonly 
return scalar mean_$A = r(mean)

* Generate A (predicted treatment under observed data)
estimates restore treatment

gen Ap = _b[_cons] + ///
         _b[age]*age + ///
         _b[sex]*sex + ///
         _b[comorbbin]*comorbbin + ///
         _b[mny]*mny + ///
         _b[asasmri]*asasmri + ///
         _b[hla]*hla + ///
         _b[pertvt1]*pertvt1 + ///
         _b[ibdbl]*ibdbl + ///
         _b[emmtvt1]*emmtvt1 + ///
         _b[comedtvt1]*comedtvt1 + ///
         _b[asdastotalt1]*asdastotalt1 + ///
         _b[basfitotalt1]*basfitotalt1
replace Ap= invlogit(Ap)
replace Ap= rbinomial(1, Ap) // random error added

qui sum Ap, meanonly 
return scalar mean_Ap = r(mean)

* Observed mediator
qui sum $M, meanonly 
return scalar mean_$M = r(mean)

* Restore the mediator model coefficients
estimates restore mediator

* Generate M (Predicted ASDAS at 6 months under observed treatment)
gen Mp = _b[_cons] + ///
         _b[bionew]*Ap + ///
         _b[age]*age + ///
         _b[sex]*sex + ///
         _b[comorbbin]*comorbbin + ///
         _b[mny]*mny + ///
         _b[asasmri]*asasmri + ///
         _b[hla]*hla + ///
         _b[pertvt1]*pertvt1 + ///
         _b[ibdbl]*ibdbl + ///
         _b[emmtvt1]*emmtvt1 + ///
         _b[comedtvt1]*comedtvt1 + ///
         _b[asdastotalt1]*asdastotalt1 + ///
         _b[basfitotalt1]*basfitotalt1

qui sum Mp, meanonly 
return scalar mean_Mp = r(mean)

* Observed outcome
qui sum $Y, meanonly 
return scalar mean_$Y = r(mean)


* Restore the outcome model coeficients
estimates restore outcome

gen Yp   = _b[_cons] + ///
           _b[asdastotalt2]*Mp + ///
           _b[bionew]*Ap + ///
           _b[age]*age + ///
           _b[sex]*sex + ///
           _b[comorbbin]*comorbbin + ///
           _b[mny]*mny + ///
           _b[asasmri]*asasmri + ///
           _b[hla]*hla + ///
           _b[pertvt1]*pertvt1 + ///
           _b[ibdbl]*ibdbl + ///
           _b[emmtvt1]*emmtvt1 + ///
           _b[comedtvt1]*comedtvt1 + ///
           _b[asdastotalt1]*asdastotalt1 + ///
           _b[basfitotalt1]*basfitotalt1 //+ _b[inter]*(Mp*Ap) // Delete last term if no interaction is intended.

qui sum Yp, meanonly 
return scalar mean_Yp = r(mean)



/////////// Clean up temporary variables before the next replication

capture drop M1 M0 Y1M1 Y1M0 Y0M1 Y0M0 NIE TIE AIE NDE TDE ADE Ap Mp Yp _est_mediator _est_outcome _est_treatment inter

end

*====================================================================
* Subgroup analyses with bootstrap percentile CIs
*====================================================================

tempname posth

postfile `posth' ///
    str20 subgroup ///
    byte level ///
    long N ///
    long REPS ///
    double ATE ATE_LCL ATE_UCL ///
    double ADE ADE_LCL ADE_UCL ///
    double AIE AIE_LCL AIE_UCL ///
    double P_INTER ///
    using "${tables}Supplementary_Table_S9.dta", replace

foreach subgroup in crpelevatedt1 asasmri mny { 

    * Interaction between the baseline characteristic and treatment on the total effect
    * (outcome model without the mediator)
    local extra ""
    if strpos(" $W ", " `subgroup' ") == 0 local extra "`subgroup'"
    qui regress $Y $A $W `extra' i.$A#i.`subgroup'
    qui testparm i.$A#i.`subgroup'
    local pinter = r(p)
    noi di as text "Interaction `subgroup' x $A: p = " as result %6.4f `pinter'

    foreach level in 0 1 {

        preserve

        keep if `subgroup' == `level'

        qui count
        local nsub = r(N)

        noi di as text "Running subgroup: `subgroup' = `level' (n = `nsub')"

        set seed `seed'

        noi bootstrap ///
            ATE=r(ate) ///
            AIE=r(mean_AIE) ///
            ADE=r(mean_ADE), ///
            reps(`sims') nodots : run_gformula // manuscript: 1000

        * Extract point estimates, percentile CIs and completed replications
        local nreps = e(N_reps)
        matrix B  = e(b)
        matrix CI = e(ci_percentile)

        local ate     = B[1,1]
        local aie     = B[1,2]
        local ade     = B[1,3]

        local ate_lcl = CI[1,1]
        local ate_ucl = CI[2,1]

        local aie_lcl = CI[1,2]
        local aie_ucl = CI[2,2]

        local ade_lcl = CI[1,3]
        local ade_ucl = CI[2,3]

        * Store results
        post `posth' ///
            ("`subgroup'") ///
            (`level') ///
            (`nsub') ///
            (`nreps') ///
            (`ate') (`ate_lcl') (`ate_ucl') ///
            (`ade') (`ade_lcl') (`ade_ucl') ///
            (`aie') (`aie_lcl') (`aie_ucl') ///
            (`pinter')

        restore
    }
}

postclose `posth'

* keep the analysis data: the results table below replaces it in memory
tempfile anadata
save `anadata'

use "${tables}Supplementary_Table_S9.dta", clear

noi list, clean noobs

* Excel shows the decimals of the Stata display format; the p-value is written as text (3 decimals)
format ATE ATE_LCL ATE_UCL ADE ADE_LCL ADE_UCL AIE AIE_LCL AIE_UCL %12.2f
format level N REPS %12.0f
gen str8 P_TXT = strtrim(string(P_INTER, "%9.3f"))
drop P_INTER
rename P_TXT P_INTER

export excel using "${tables}Supplementary_Table_S9.xlsx", ///
    firstrow(variables) replace

savesheet "S9"

use `anadata', clear // back to the analysis data for Table S9.1


*====================================================================
* Supplementary Table S9.1: two causally ordered mediators (CRP, then
* ASDAS-PRO) within levels of baseline characteristics. Same models as
* run 3 (main paths, without interaction), fitted within each stratum.
*====================================================================

**** Pre-treatment baseline confounders (baseline ASDAS replaced by the baseline mediators)
global W="age sex comorbbin mny asasmri hla pertvt1 ibdbl emmtvt1 comedtvt1 basfitotalt1"

**** Baseline values of mediators
global bM1="crpt1"
global bM2="asdastotalt1_pro"

**** Mediators (continuous)
global M1="crpt2"
global M2="asdastotalt2_pro"

capture program drop run_gformula_paths_sub
program define run_gformula_paths_sub, rclass

	///////////////// Step 1 — Model the observed data, in causal order

	* --- Outcome model: BASFI ~ bionew + M1 + M2 + W + bM1 + bM2 ---
	regress $Y $A $M1 $M2 $W $bM1 $bM2
	estimates store outcome

	* --- Mediator 2 model: M2 ~ bionew + M1 + W + bM1 + bM2 (CRP causally precedes BASDAI) ---
	regress $M2 $A $M1 $W $bM1 $bM2
	estimates store mediator2

	* --- Mediator 1 model: M1 ~ bionew + W + bM1 + bM2 ---
	regress $M1 $A $W $bM1 $bM2
	estimates store mediator1

	///////////////// Step 2 — Monte Carlo simulation of nested counterfactuals

	*----------------------------------------------------------------------
	* Step 2.1 — Counterfactual M1 under treatment (1) and control (0)
	*----------------------------------------------------------------------
	estimates restore mediator1

	gen M1_1 = _b[_cons] + _b[$A]*1 + ///
	           _b[age]*age + _b[sex]*sex + _b[comorbbin]*comorbbin + ///
	           _b[mny]*mny + _b[asasmri]*asasmri + _b[hla]*hla + ///
	           _b[pertvt1]*pertvt1 + _b[ibdbl]*ibdbl + ///
	           _b[emmtvt1]*emmtvt1 + _b[comedtvt1]*comedtvt1 + ///
	           _b[$bM1]*$bM1 + _b[$bM2]*$bM2 + ///
	           _b[basfitotalt1]*basfitotalt1

	gen M1_0 = _b[_cons] + _b[$A]*0 + ///
	           _b[age]*age + _b[sex]*sex + _b[comorbbin]*comorbbin + ///
	           _b[mny]*mny + _b[asasmri]*asasmri + _b[hla]*hla + ///
	           _b[pertvt1]*pertvt1 + _b[ibdbl]*ibdbl + ///
	           _b[emmtvt1]*emmtvt1 + _b[comedtvt1]*comedtvt1 + ///
	           _b[$bM1]*$bM1 + _b[$bM2]*$bM2 + ///
	           _b[basfitotalt1]*basfitotalt1

	*----------------------------------------------------------------------
	* Step 2.2 — Counterfactual M2 at exactly the 3 combinations needed
	*   M2_A0_M1_0 : bionew=0 in the M2 model, M1 held at its CONTROL value
	*   M2_A0_M1_1 : bionew=0 in the M2 model, M1 held at its TREATED value
	*   M2_A1_M1_1 : bionew=1 in the M2 model, M1 held at its TREATED value
	* (bionew=1, M1=M1_0 is never needed for the Type I decomposition)
	*----------------------------------------------------------------------
	estimates restore mediator2

	gen M2_A0_M1_0 = _b[_cons] + _b[$A]*0 + _b[$M1]*M1_0 + ///
	                 _b[age]*age + _b[sex]*sex + _b[comorbbin]*comorbbin + ///
	                 _b[mny]*mny + _b[asasmri]*asasmri + _b[hla]*hla + ///
	                 _b[pertvt1]*pertvt1 + _b[ibdbl]*ibdbl + ///
	                 _b[emmtvt1]*emmtvt1 + _b[comedtvt1]*comedtvt1 + ///
	                 _b[$bM1]*$bM1 + _b[$bM2]*$bM2 + ///
	                 _b[basfitotalt1]*basfitotalt1

	gen M2_A0_M1_1 = _b[_cons] + _b[$A]*0 + _b[$M1]*M1_1 + ///
	                 _b[age]*age + _b[sex]*sex + _b[comorbbin]*comorbbin + ///
	                 _b[mny]*mny + _b[asasmri]*asasmri + _b[hla]*hla + ///
	                 _b[pertvt1]*pertvt1 + _b[ibdbl]*ibdbl + ///
	                 _b[emmtvt1]*emmtvt1 + _b[comedtvt1]*comedtvt1 + ///
	                 _b[$bM1]*$bM1 + _b[$bM2]*$bM2 + ///
	                 _b[basfitotalt1]*basfitotalt1

	gen M2_A1_M1_1 = _b[_cons] + _b[$A]*1 + _b[$M1]*M1_1 + ///
	                 _b[age]*age + _b[sex]*sex + _b[comorbbin]*comorbbin + ///
	                 _b[mny]*mny + _b[asasmri]*asasmri + _b[hla]*hla + ///
	                 _b[pertvt1]*pertvt1 + _b[ibdbl]*ibdbl + ///
	                 _b[emmtvt1]*emmtvt1 + _b[comedtvt1]*comedtvt1 + ///
	                 _b[$bM1]*$bM1 + _b[$bM2]*$bM2 + ///
	                 _b[basfitotalt1]*basfitotalt1

	*----------------------------------------------------------------------
	* Step 2.3 — Counterfactual BASFI at the 4 nested (a0, M1, M2) combos
	*   g(0,0,0) : everyone untreated, mediators at their untreated values
	*   g(1,0,0) : outcome sees bionew=1, mediators STILL at untreated values
	*              -> isolates the direct effect A->Y
	*   g(1,1,0) : outcome sees bionew=1, M1 now treated, M2 still at its
	*              "A=0,M1=M1_1" value -> the extra move from g(1,0,0) to
	*              g(1,1,0) isolates the path-specific effect via M1
	*              (which may travel on through M2)
	*   g(1,1,1) : outcome sees bionew=1, M1 treated, M2 now fully treated
	*              -> the extra move from g(1,1,0) to g(1,1,1) isolates
	*              the path-specific effect via M2 (direct arrow)
	*----------------------------------------------------------------------
	estimates restore outcome

	gen Y_000 = _b[_cons] + _b[$A]*0 + _b[$M1]*M1_0 + _b[$M2]*M2_A0_M1_0 + ///
	            _b[age]*age + _b[sex]*sex + _b[comorbbin]*comorbbin + ///
	            _b[mny]*mny + _b[asasmri]*asasmri + _b[hla]*hla + ///
	            _b[pertvt1]*pertvt1 + _b[ibdbl]*ibdbl + ///
	            _b[emmtvt1]*emmtvt1 + _b[comedtvt1]*comedtvt1 + ///
	            _b[$bM1]*$bM1 + _b[$bM2]*$bM2 + ///
	            _b[basfitotalt1]*basfitotalt1

	gen Y_100 = _b[_cons] + _b[$A]*1 + _b[$M1]*M1_0 + _b[$M2]*M2_A0_M1_0 + ///
	            _b[age]*age + _b[sex]*sex + _b[comorbbin]*comorbbin + ///
	            _b[mny]*mny + _b[asasmri]*asasmri + _b[hla]*hla + ///
	            _b[pertvt1]*pertvt1 + _b[ibdbl]*ibdbl + ///
	            _b[emmtvt1]*emmtvt1 + _b[comedtvt1]*comedtvt1 + ///
	            _b[$bM1]*$bM1 + _b[$bM2]*$bM2 + ///
	            _b[basfitotalt1]*basfitotalt1

	gen Y_110 = _b[_cons] + _b[$A]*1 + _b[$M1]*M1_1 + _b[$M2]*M2_A0_M1_1 + ///
	            _b[age]*age + _b[sex]*sex + _b[comorbbin]*comorbbin + ///
	            _b[mny]*mny + _b[asasmri]*asasmri + _b[hla]*hla + ///
	            _b[pertvt1]*pertvt1 + _b[ibdbl]*ibdbl + ///
	            _b[emmtvt1]*emmtvt1 + _b[comedtvt1]*comedtvt1 + ///
	            _b[$bM1]*$bM1 + _b[$bM2]*$bM2 + ///
	            _b[basfitotalt1]*basfitotalt1

	gen Y_111 = _b[_cons] + _b[$A]*1 + _b[$M1]*M1_1 + _b[$M2]*M2_A1_M1_1 + ///
	            _b[age]*age + _b[sex]*sex + _b[comorbbin]*comorbbin + ///
	            _b[mny]*mny + _b[asasmri]*asasmri + _b[hla]*hla + ///
	            _b[pertvt1]*pertvt1 + _b[ibdbl]*ibdbl + ///
	            _b[emmtvt1]*emmtvt1 + _b[comedtvt1]*comedtvt1 + ///
	            _b[$bM1]*$bM1 + _b[$bM2]*$bM2 + ///
	            _b[basfitotalt1]*basfitotalt1

	///////////////// Step 3 — Path-specific effects (Type I decomposition)

	qui sum Y_000
	scalar m000 = r(mean)
	qui sum Y_100
	scalar m100 = r(mean)
	qui sum Y_110
	scalar m110 = r(mean)
	qui sum Y_111
	scalar m111 = r(mean)

	return scalar direct = m100 - m000          // bionew -> BASFI
	return scalar pse_m1 = m110 - m100          // bionew -> M1 ~> BASFI
	return scalar pse_m2 = m111 - m110          // bionew -> M2 -> BASFI
	return scalar total  = m111 - m000          // Total effect (should equal direct + pse_m1 + pse_m2)

	* Effects of treatment on the two mediators themselves (mediator scale)
	qui sum M1_1
	scalar mean_M1_1 = r(mean)
	qui sum M1_0
	scalar mean_M1_0 = r(mean)
	qui sum M2_A1_M1_1
	scalar mean_M2_A1_M1_1 = r(mean)
	qui sum M2_A0_M1_0
	scalar mean_M2_A0_M1_0 = r(mean)

	return scalar ate_m1 = mean_M1_1 - mean_M1_0               // bionew -> CRP
	return scalar ate_m2 = mean_M2_A1_M1_1 - mean_M2_A0_M1_0   // bionew -> ASDAS-PRO

	///////////////// Clean up before the next bootstrap replication

	capture drop M1_1 M1_0 M2_A0_M1_0 M2_A0_M1_1 M2_A1_M1_1 ///
	             Y_000 Y_100 Y_110 Y_111 _est_outcome _est_mediator1 _est_mediator2

end

tempname posth2

postfile `posth2' ///
    str20 subgroup ///
    byte level ///
    long N ///
    long REPS ///
    double TOTAL TOTAL_LCL TOTAL_UCL ///
    double DIRECT DIRECT_LCL DIRECT_UCL ///
    double PSE_CRP PSE_CRP_LCL PSE_CRP_UCL ///
    double PSE_PRO PSE_PRO_LCL PSE_PRO_UCL ///
    double ATE_CRP ATE_CRP_LCL ATE_CRP_UCL ///
    double ATE_PRO ATE_PRO_LCL ATE_PRO_UCL ///
    using "${tables}Supplementary_Table_S9_1.dta", replace

foreach subgroup in crpelevatedt1 asasmri mny { 

    foreach level in 0 1 {

        preserve

        keep if `subgroup' == `level'

        qui count
        local nsub = r(N)

        noi di as text "Running subgroup, two mediators: `subgroup' = `level' (n = `nsub')"

        set seed `seed'

        noi bootstrap ///
            Total=r(total) ///
            Direct=r(direct) ///
            PSE_M1=r(pse_m1) ///
            PSE_M2=r(pse_m2) ///
            ATE_M1=r(ate_m1) ///
            ATE_M2=r(ate_m2), ///
            reps(`sims') nodots : run_gformula_paths_sub

        local nreps = e(N_reps)
        matrix B  = e(b)
        matrix CI = e(ci_percentile)

        * Store results (columns of B in the order of the bootstrap call)
        post `posth2' ///
            ("`subgroup'") ///
            (`level') ///
            (`nsub') ///
            (`nreps') ///
            (B[1,1]) (CI[1,1]) (CI[2,1]) ///
            (B[1,2]) (CI[1,2]) (CI[2,2]) ///
            (B[1,3]) (CI[1,3]) (CI[2,3]) ///
            (B[1,4]) (CI[1,4]) (CI[2,4]) ///
            (B[1,5]) (CI[1,5]) (CI[2,5]) ///
            (B[1,6]) (CI[1,6]) (CI[2,6])

        restore
    }
}

postclose `posth2'

use "${tables}Supplementary_Table_S9_1.dta", clear

noi list, clean noobs

* Excel shows the decimals of the Stata display format
format TOTAL TOTAL_LCL TOTAL_UCL DIRECT DIRECT_LCL DIRECT_UCL PSE_CRP PSE_CRP_LCL PSE_CRP_UCL PSE_PRO PSE_PRO_LCL PSE_PRO_UCL ATE_CRP ATE_CRP_LCL ATE_CRP_UCL ATE_PRO ATE_PRO_LCL ATE_PRO_UCL %12.2f
format level N REPS %12.0f

export excel using "${tables}Supplementary_Table_S9_1.xlsx", ///
    firstrow(variables) replace

savesheet "S9_1"


} // close 7




******************************************************************************************************************************************************
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
**** 8 - DESCRIPTIVE: patients with complete vs missing 6-month data, and positivity of being observed (Supplementary Table S11)
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
******************************************************************************************************************************************************

if `run' == 8 {

*===============================================================================
* PURPOSE:      Describe the patients lost to missing 6-month data and check that
*               every patient had a non-zero probability of being observed.
*
* POPULATION:   Eligible patients with complete baseline covariates ($W): n = 481.
*               Of these, 419 have ASDAS and BASFI at 6 months (the analysis
*               cohort) and 62 do not.
*
* NOTE:         No estimate is produced here. With linear mediator and outcome
*               models without interactions, the mediation effects do not depend
*               on the covariate distribution, so standardising to the 481 returns
*               the same point estimates as the 419. The 419 estimates therefore
*               apply to all 481 if missingness depends only on treatment and
*               baseline characteristics. This run describes how plausible that is.
*===============================================================================

assert `cohort' == 800 // defined for the two-visit cohort only

*-------------------------------------------------------------------------------
* Eligible patients with complete baseline covariates, one row per patient
*-------------------------------------------------------------------------------

import delimited "${datasets}originalfulllong800.csv", clear varnames(1) case(preserve)
keep if t == 6 // one row per patient

local nW : word count $W
egen byte nbase = rownonmiss($W)
gen byte basecomplete = 0
replace basecomplete = 1 if nbase == `nW'
keep if basecomplete == 1

qui count
assert r(N) == 481 // stops here if the source file or the selection changes
local neligible = r(N)
global anaN = r(N)

*-------------------------------------------------------------------------------
* Complete 6-month data (ASDAS and BASFI): the analysis cohort
*-------------------------------------------------------------------------------

gen byte complete6 = 0
replace complete6 = 1 if !missing(asdastotalt2) & !missing(basfitotalt2)

qui count if complete6 == 1
assert r(N) == 419 // the analysis cohort
local ncomplete = r(N)

qui count if complete6 == 0
local nmissing = r(N)

* The two-mediator analysis also needs CRP at 6 months: no completer lacks it
qui count if complete6 == 1 & missing(crpt2)
assert r(N) == 0

label define complete6lbl 1 "Complete 6M" 0 "Missing 6M"
label values complete6 complete6lbl

noi di ""
noi di as text "Eligible with complete baseline covariates: " as result `neligible'
noi di as text "  complete ASDAS and BASFI at 6 months:     " as result `ncomplete'
noi di as text "  missing ASDAS or BASFI at 6 months:       " as result `nmissing'
noi di ""
noi tab bionew complete6, row

*-------------------------------------------------------------------------------
* Baseline characteristics, complete vs missing 6-month data (my program)
*-------------------------------------------------------------------------------

qui do "${programs}smdtableV2.do" // Run supporting program

global descriptive="bionew age sex symptomsduration comorbbin mny asasmri hla crpelevatedt1 arthtvt1 enthtvt1 psotvt1 ibdbl aautvt1 asdastotalt1 i.asdascatt1 basfitotalt1 nsaidnewt1 csdmardnewt1 gcnewt1"
global groupvars "basecomplete==1"
noi smdtableV2, selector($groupvars) visit(6) analysis(table1) variable($descriptive) treatvar(complete6) roundvar(0) roundsmd(2) difference(0.1) save(yes) ///
    folder($tables) filename(Supplementary_Table_S11)

*-------------------------------------------------------------------------------
* Positivity of being observed: P(complete 6-month data | treatment, baseline)
*-------------------------------------------------------------------------------
* The same model a censoring weight would use. Every patient needs a predicted
* probability clearly above zero, in the treated and in the untreated.

noi logit complete6 bionew $W
* Baseline IBD predicts complete data perfectly (all IBD patients are complete):
* logit drops ibdbl and those patients from the model. The option rules gives them
* their predicted probability of 1 instead of a missing value.
predict double pcomplete, pr rules

qui count if complete6 == 0 & bionew == 1
local nmiss1 = r(N)
qui count if complete6 == 0 & bionew == 0
local nmiss0 = r(N)

qui sum pcomplete if bionew == 1
local pmin1 = r(min)
local pmax1 = r(max)
qui sum pcomplete if bionew == 0
local pmin0 = r(min)
local pmax0 = r(max)

noi di ""
noi di as text "Predicted probability of complete 6-month data"
noi di as text "  treated:   min " as result %5.3f `pmin1' as text "  max " as result %5.3f `pmax1'
noi di as text "  untreated: min " as result %5.3f `pmin0' as text "  max " as result %5.3f `pmax0'
noi di ""

savediag "T2" "Missing 6-month data (descriptive)" "N eligible with complete baseline" `neligible'
savediag "T2" "Missing 6-month data (descriptive)" "N complete 6-month data" `ncomplete'
savediag "T2" "Missing 6-month data (descriptive)" "N missing 6-month data" `nmissing'
savediag "T2" "Missing 6-month data (descriptive)" "N missing among treated" `nmiss1'
savediag "T2" "Missing 6-month data (descriptive)" "N missing among untreated" `nmiss0'
savediag "T2" "Missing 6-month data (descriptive)" "Min P(complete) treated" `pmin1'
savediag "T2" "Missing 6-month data (descriptive)" "Min P(complete) untreated" `pmin0'
savesheet "T2"


} // close 8


******************************************************************************************************************************************************
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
**** 8.1 - SENSITIVITY (Table 2): treatment defined as at least 1 month of bDMARD exposure (manual, single and two mediators)
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
******************************************************************************************************************************************************

if `run' == 8.1 {

*===============================================================================
* REFERENCE:    Boxes S2 and S2.1, main specification, same cohort (n = 419)
* CHANGE:       Treatment = bionew1m (at least 30 days of bDMARD exposure between
*               baseline and 6 months) instead of bionew (any bDMARD recorded).
*
* SPECIFICATION: As in the main analysis: no treatment-mediator interaction
*               (single mediator) and no interaction (two mediators). The
*               specification is fixed, not re-tested, so only the treatment
*               definition changes.
*===============================================================================

assert `cohort' == 800 // defined for the two-visit cohort only

**** Treatment (binary)
global A="bionew1m"

qui count if missing($A)
assert r(N) == 0
qui count if $A == 1
assert r(N) == 112 // stops here if the source file or the definition changes
local ntreat = r(N)
qui count if $A == 0
local nuntreat = r(N)

noi di ""
noi di as text "Treatment $A: treated " as result `ntreat' as text ", untreated " as result `nuntreat'
noi di ""


*===============================================================================
* A. Single mediator (ASDAS): Box S2, manual
*===============================================================================

**** Pre-treatment baseline confounders
global W="age sex comorbbin mny asasmri hla pertvt1 ibdbl emmtvt1 comedtvt1 asdastotalt1 basfitotalt1"

**** Mediator (continuous)
global M="asdastotalt2"

**** Outcome (continuous)
global Y="basfitotalt2"

capture program drop run_gformula_t2
program define run_gformula_t2, rclass

///////////////// Step 1 — Model the Observed Data

**** Outcome model
regress $Y $M $A $W
estimates store outcome

**** Mediator model
regress $M $A $W
estimates store mediator

///////////////// Step 2 — Monte Carlo Simulation

/////////// Step 2.1. Counterfactual mediator (ASDAS) at 6 months

estimates restore mediator

* M1: ASDAS at 6 months if everyone were treated
gen M1 = _b[_cons] + ///
         _b[$A]*1 + ///
         _b[age]*age + ///
         _b[sex]*sex + ///
         _b[comorbbin]*comorbbin + ///
         _b[mny]*mny + ///
         _b[asasmri]*asasmri + ///
         _b[hla]*hla + ///
         _b[pertvt1]*pertvt1 + ///
         _b[ibdbl]*ibdbl + ///
         _b[emmtvt1]*emmtvt1 + ///
         _b[comedtvt1]*comedtvt1 + ///
         _b[asdastotalt1]*asdastotalt1 + ///
         _b[basfitotalt1]*basfitotalt1

* M0: ASDAS at 6 months if everyone were untreated
gen M0 = _b[_cons] + ///
         _b[$A]*0 + ///
         _b[age]*age + ///
         _b[sex]*sex + ///
         _b[comorbbin]*comorbbin + ///
         _b[mny]*mny + ///
         _b[asasmri]*asasmri + ///
         _b[hla]*hla + ///
         _b[pertvt1]*pertvt1 + ///
         _b[ibdbl]*ibdbl + ///
         _b[emmtvt1]*emmtvt1 + ///
         _b[comedtvt1]*comedtvt1 + ///
         _b[asdastotalt1]*asdastotalt1 + ///
         _b[basfitotalt1]*basfitotalt1

/////////// Step 2.2. Counterfactual outcome (BASFI) at 6 months

estimates restore outcome

* Y1M1: treated, ASDAS at its treated value
gen Y1M1 = _b[_cons] + ///
           _b[$M]*M1 + ///
           _b[$A]*1 + ///
           _b[age]*age + ///
           _b[sex]*sex + ///
           _b[comorbbin]*comorbbin + ///
           _b[mny]*mny + ///
           _b[asasmri]*asasmri + ///
           _b[hla]*hla + ///
           _b[pertvt1]*pertvt1 + ///
           _b[ibdbl]*ibdbl + ///
           _b[emmtvt1]*emmtvt1 + ///
           _b[comedtvt1]*comedtvt1 + ///
           _b[asdastotalt1]*asdastotalt1 + ///
           _b[basfitotalt1]*basfitotalt1

* Y1M0: treated, ASDAS at its untreated value
gen Y1M0 = _b[_cons] + ///
           _b[$M]*M0 + ///
           _b[$A]*1 + ///
           _b[age]*age + ///
           _b[sex]*sex + ///
           _b[comorbbin]*comorbbin + ///
           _b[mny]*mny + ///
           _b[asasmri]*asasmri + ///
           _b[hla]*hla + ///
           _b[pertvt1]*pertvt1 + ///
           _b[ibdbl]*ibdbl + ///
           _b[emmtvt1]*emmtvt1 + ///
           _b[comedtvt1]*comedtvt1 + ///
           _b[asdastotalt1]*asdastotalt1 + ///
           _b[basfitotalt1]*basfitotalt1

* Y0M1: untreated, ASDAS at its treated value
gen Y0M1 = _b[_cons] + ///
           _b[$M]*M1 + ///
           _b[$A]*0 + ///
           _b[age]*age + ///
           _b[sex]*sex + ///
           _b[comorbbin]*comorbbin + ///
           _b[mny]*mny + ///
           _b[asasmri]*asasmri + ///
           _b[hla]*hla + ///
           _b[pertvt1]*pertvt1 + ///
           _b[ibdbl]*ibdbl + ///
           _b[emmtvt1]*emmtvt1 + ///
           _b[comedtvt1]*comedtvt1 + ///
           _b[asdastotalt1]*asdastotalt1 + ///
           _b[basfitotalt1]*basfitotalt1

* Y0M0: untreated, ASDAS at its untreated value
gen Y0M0 = _b[_cons] + ///
           _b[$M]*M0 + ///
           _b[$A]*0 + ///
           _b[age]*age + ///
           _b[sex]*sex + ///
           _b[comorbbin]*comorbbin + ///
           _b[mny]*mny + ///
           _b[asasmri]*asasmri + ///
           _b[hla]*hla + ///
           _b[pertvt1]*pertvt1 + ///
           _b[ibdbl]*ibdbl + ///
           _b[emmtvt1]*emmtvt1 + ///
           _b[comedtvt1]*comedtvt1 + ///
           _b[asdastotalt1]*asdastotalt1 + ///
           _b[basfitotalt1]*basfitotalt1

/////////// Step 3 — Mediation effects

qui sum Y1M1
scalar mean_Y1M1 = r(mean)
qui sum Y1M0
scalar mean_Y1M0 = r(mean)
qui sum Y0M1
scalar mean_Y0M1 = r(mean)
qui sum Y0M0
scalar mean_Y0M0 = r(mean)

gen NIE = Y0M1-Y0M0
gen TIE = Y1M1-Y1M0
gen AIE = (NIE+TIE)/2
qui sum AIE
return scalar mean_AIE = r(mean)

gen NDE = Y1M0-Y0M0
gen TDE = Y1M1-Y0M1
gen ADE = (NDE+TDE)/2
qui sum ADE
return scalar mean_ADE = r(mean)

return scalar ate    = mean_Y1M1 - mean_Y0M0
return scalar mean11 = mean_Y1M1
return scalar mean00 = mean_Y0M0

/////////// Clean up temporary variables before the next replication

capture drop M1 M0 Y1M1 Y1M0 Y0M1 Y0M0 NIE TIE AIE NDE TDE ADE _est_mediator _est_outcome

end

set seed `seed'
bootstrap ATE=r(ate) AIE=r(mean_AIE) ADE=r(mean_ADE) Mean_Y1M1=r(mean11) Mean_Y0M0=r(mean00) ///
          , reps(`sims') nodots: run_gformula_t2

noi estat bootstrap, percentile

savebs "T2" "Single mediator - bDMARD at least 1 month (manual)"
savediag "T2" "Single mediator - bDMARD at least 1 month (manual)" "N treated" `ntreat'
savediag "T2" "Single mediator - bDMARD at least 1 month (manual)" "N untreated" `nuntreat'


*===============================================================================
* B. Two causally ordered mediators (CRP, then ASDAS-PRO): Box S2.1, manual,
*    main paths, without interaction
*===============================================================================

**** Pre-treatment baseline confounders (baseline ASDAS replaced by the baseline mediators)
global W="age sex comorbbin mny asasmri hla pertvt1 ibdbl emmtvt1 comedtvt1 basfitotalt1"

**** Baseline values of mediators
global bM1="crpt1"
global bM2="asdastotalt1_pro"

**** Mediators (continuous)
global M1="crpt2"
global M2="asdastotalt2_pro"

**** Outcome (continuous)
global Y="basfitotalt2"

capture program drop run_gformula_paths_t2
program define run_gformula_paths_t2, rclass

	///////////////// Step 1 — Model the observed data, in causal order

	regress $Y $A $M1 $M2 $W $bM1 $bM2
	estimates store outcome

	regress $M2 $A $M1 $W $bM1 $bM2
	estimates store mediator2

	regress $M1 $A $W $bM1 $bM2
	estimates store mediator1

	///////////////// Step 2 — Monte Carlo simulation of nested counterfactuals

	* Step 2.1 — Counterfactual M1 under treatment (1) and control (0)
	estimates restore mediator1

	gen M1_1 = _b[_cons] + _b[$A]*1 + ///
	           _b[age]*age + _b[sex]*sex + _b[comorbbin]*comorbbin + ///
	           _b[mny]*mny + _b[asasmri]*asasmri + _b[hla]*hla + ///
	           _b[pertvt1]*pertvt1 + _b[ibdbl]*ibdbl + ///
	           _b[emmtvt1]*emmtvt1 + _b[comedtvt1]*comedtvt1 + ///
	           _b[$bM1]*$bM1 + _b[$bM2]*$bM2 + ///
	           _b[basfitotalt1]*basfitotalt1

	gen M1_0 = _b[_cons] + _b[$A]*0 + ///
	           _b[age]*age + _b[sex]*sex + _b[comorbbin]*comorbbin + ///
	           _b[mny]*mny + _b[asasmri]*asasmri + _b[hla]*hla + ///
	           _b[pertvt1]*pertvt1 + _b[ibdbl]*ibdbl + ///
	           _b[emmtvt1]*emmtvt1 + _b[comedtvt1]*comedtvt1 + ///
	           _b[$bM1]*$bM1 + _b[$bM2]*$bM2 + ///
	           _b[basfitotalt1]*basfitotalt1

	* Step 2.2 — Counterfactual M2 at the 3 combinations needed
	estimates restore mediator2

	gen M2_A0_M1_0 = _b[_cons] + _b[$A]*0 + _b[$M1]*M1_0 + ///
	                 _b[age]*age + _b[sex]*sex + _b[comorbbin]*comorbbin + ///
	                 _b[mny]*mny + _b[asasmri]*asasmri + _b[hla]*hla + ///
	                 _b[pertvt1]*pertvt1 + _b[ibdbl]*ibdbl + ///
	                 _b[emmtvt1]*emmtvt1 + _b[comedtvt1]*comedtvt1 + ///
	                 _b[$bM1]*$bM1 + _b[$bM2]*$bM2 + ///
	                 _b[basfitotalt1]*basfitotalt1

	gen M2_A0_M1_1 = _b[_cons] + _b[$A]*0 + _b[$M1]*M1_1 + ///
	                 _b[age]*age + _b[sex]*sex + _b[comorbbin]*comorbbin + ///
	                 _b[mny]*mny + _b[asasmri]*asasmri + _b[hla]*hla + ///
	                 _b[pertvt1]*pertvt1 + _b[ibdbl]*ibdbl + ///
	                 _b[emmtvt1]*emmtvt1 + _b[comedtvt1]*comedtvt1 + ///
	                 _b[$bM1]*$bM1 + _b[$bM2]*$bM2 + ///
	                 _b[basfitotalt1]*basfitotalt1

	gen M2_A1_M1_1 = _b[_cons] + _b[$A]*1 + _b[$M1]*M1_1 + ///
	                 _b[age]*age + _b[sex]*sex + _b[comorbbin]*comorbbin + ///
	                 _b[mny]*mny + _b[asasmri]*asasmri + _b[hla]*hla + ///
	                 _b[pertvt1]*pertvt1 + _b[ibdbl]*ibdbl + ///
	                 _b[emmtvt1]*emmtvt1 + _b[comedtvt1]*comedtvt1 + ///
	                 _b[$bM1]*$bM1 + _b[$bM2]*$bM2 + ///
	                 _b[basfitotalt1]*basfitotalt1

	* Step 2.3 — Counterfactual BASFI at the 4 nested combinations
	estimates restore outcome

	gen Y_000 = _b[_cons] + _b[$A]*0 + _b[$M1]*M1_0 + _b[$M2]*M2_A0_M1_0 + ///
	            _b[age]*age + _b[sex]*sex + _b[comorbbin]*comorbbin + ///
	            _b[mny]*mny + _b[asasmri]*asasmri + _b[hla]*hla + ///
	            _b[pertvt1]*pertvt1 + _b[ibdbl]*ibdbl + ///
	            _b[emmtvt1]*emmtvt1 + _b[comedtvt1]*comedtvt1 + ///
	            _b[$bM1]*$bM1 + _b[$bM2]*$bM2 + ///
	            _b[basfitotalt1]*basfitotalt1

	gen Y_100 = _b[_cons] + _b[$A]*1 + _b[$M1]*M1_0 + _b[$M2]*M2_A0_M1_0 + ///
	            _b[age]*age + _b[sex]*sex + _b[comorbbin]*comorbbin + ///
	            _b[mny]*mny + _b[asasmri]*asasmri + _b[hla]*hla + ///
	            _b[pertvt1]*pertvt1 + _b[ibdbl]*ibdbl + ///
	            _b[emmtvt1]*emmtvt1 + _b[comedtvt1]*comedtvt1 + ///
	            _b[$bM1]*$bM1 + _b[$bM2]*$bM2 + ///
	            _b[basfitotalt1]*basfitotalt1

	gen Y_110 = _b[_cons] + _b[$A]*1 + _b[$M1]*M1_1 + _b[$M2]*M2_A0_M1_1 + ///
	            _b[age]*age + _b[sex]*sex + _b[comorbbin]*comorbbin + ///
	            _b[mny]*mny + _b[asasmri]*asasmri + _b[hla]*hla + ///
	            _b[pertvt1]*pertvt1 + _b[ibdbl]*ibdbl + ///
	            _b[emmtvt1]*emmtvt1 + _b[comedtvt1]*comedtvt1 + ///
	            _b[$bM1]*$bM1 + _b[$bM2]*$bM2 + ///
	            _b[basfitotalt1]*basfitotalt1

	gen Y_111 = _b[_cons] + _b[$A]*1 + _b[$M1]*M1_1 + _b[$M2]*M2_A1_M1_1 + ///
	            _b[age]*age + _b[sex]*sex + _b[comorbbin]*comorbbin + ///
	            _b[mny]*mny + _b[asasmri]*asasmri + _b[hla]*hla + ///
	            _b[pertvt1]*pertvt1 + _b[ibdbl]*ibdbl + ///
	            _b[emmtvt1]*emmtvt1 + _b[comedtvt1]*comedtvt1 + ///
	            _b[$bM1]*$bM1 + _b[$bM2]*$bM2 + ///
	            _b[basfitotalt1]*basfitotalt1

	///////////////// Step 3 — Path-specific effects (Type I decomposition)

	qui sum Y_000
	scalar m000 = r(mean)
	qui sum Y_100
	scalar m100 = r(mean)
	qui sum Y_110
	scalar m110 = r(mean)
	qui sum Y_111
	scalar m111 = r(mean)

	return scalar direct = m100 - m000          // bDMARD -> BASFI
	return scalar pse_m1 = m110 - m100          // bDMARD -> CRP ~> BASFI
	return scalar pse_m2 = m111 - m110          // bDMARD -> ASDAS-PRO -> BASFI
	return scalar total  = m111 - m000          // total effect

	///////////////// Clean up before the next bootstrap replication

	capture drop M1_1 M1_0 M2_A0_M1_0 M2_A0_M1_1 M2_A1_M1_1 ///
	             Y_000 Y_100 Y_110 Y_111 _est_outcome _est_mediator1 _est_mediator2

end

set seed `seed'
bootstrap Direct=r(direct) PSE_M1=r(pse_m1) PSE_M2=r(pse_m2) Total=r(total), ///
	reps(`sims') nodots: run_gformula_paths_t2

noi estat bootstrap, percentile

savebs "T2" "Two mediators - bDMARD at least 1 month (manual)"
savediag "T2" "Two mediators - bDMARD at least 1 month (manual)" "N treated" `ntreat'
savediag "T2" "Two mediators - bDMARD at least 1 month (manual)" "N untreated" `nuntreat'

savesheet "T2"


} // close 8.1


******************************************************************************************************************************************************
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
**** 8.2 - SENSITIVITY (Table 2): treatment defined as at least 3 months of bDMARD exposure (manual, single and two mediators)
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
******************************************************************************************************************************************************

if `run' == 8.2 {

*===============================================================================
* REFERENCE:    Boxes S2 and S2.1, main specification, same cohort (n = 419)
* CHANGE:       Treatment = bionew3m (at least 90 days of bDMARD exposure between
*               baseline and 6 months) instead of bionew (any bDMARD recorded).
*               Same definition as in the target trial emulation.
* SPECIFICATION: As in the main analysis: no treatment-mediator interaction
*               (single mediator) and no interaction (two mediators). The
*               specification is fixed, not re-tested, so only the treatment
*               definition changes.
*===============================================================================

assert `cohort' == 800 // defined for the two-visit cohort only

**** Treatment (binary)
global A="bionew3m"

qui count if missing($A)
assert r(N) == 0
qui count if $A == 1
assert r(N) == 81 // stops here if the source file or the definition changes
local ntreat = r(N)
qui count if $A == 0
local nuntreat = r(N)

noi di ""
noi di as text "Treatment $A: treated " as result `ntreat' as text ", untreated " as result `nuntreat'
noi di ""


*===============================================================================
* A. Single mediator (ASDAS): Box S2, manual
*===============================================================================

**** Pre-treatment baseline confounders
global W="age sex comorbbin mny asasmri hla pertvt1 ibdbl emmtvt1 comedtvt1 asdastotalt1 basfitotalt1"

**** Mediator (continuous)
global M="asdastotalt2"

**** Outcome (continuous)
global Y="basfitotalt2"

capture program drop run_gformula_t2
program define run_gformula_t2, rclass

///////////////// Step 1 — Model the Observed Data

**** Outcome model
regress $Y $M $A $W
estimates store outcome

**** Mediator model
regress $M $A $W
estimates store mediator

///////////////// Step 2 — Monte Carlo Simulation

/////////// Step 2.1. Counterfactual mediator (ASDAS) at 6 months

estimates restore mediator

* M1: ASDAS at 6 months if everyone were treated
gen M1 = _b[_cons] + ///
         _b[$A]*1 + ///
         _b[age]*age + ///
         _b[sex]*sex + ///
         _b[comorbbin]*comorbbin + ///
         _b[mny]*mny + ///
         _b[asasmri]*asasmri + ///
         _b[hla]*hla + ///
         _b[pertvt1]*pertvt1 + ///
         _b[ibdbl]*ibdbl + ///
         _b[emmtvt1]*emmtvt1 + ///
         _b[comedtvt1]*comedtvt1 + ///
         _b[asdastotalt1]*asdastotalt1 + ///
         _b[basfitotalt1]*basfitotalt1

* M0: ASDAS at 6 months if everyone were untreated
gen M0 = _b[_cons] + ///
         _b[$A]*0 + ///
         _b[age]*age + ///
         _b[sex]*sex + ///
         _b[comorbbin]*comorbbin + ///
         _b[mny]*mny + ///
         _b[asasmri]*asasmri + ///
         _b[hla]*hla + ///
         _b[pertvt1]*pertvt1 + ///
         _b[ibdbl]*ibdbl + ///
         _b[emmtvt1]*emmtvt1 + ///
         _b[comedtvt1]*comedtvt1 + ///
         _b[asdastotalt1]*asdastotalt1 + ///
         _b[basfitotalt1]*basfitotalt1

/////////// Step 2.2. Counterfactual outcome (BASFI) at 6 months

estimates restore outcome

* Y1M1: treated, ASDAS at its treated value
gen Y1M1 = _b[_cons] + ///
           _b[$M]*M1 + ///
           _b[$A]*1 + ///
           _b[age]*age + ///
           _b[sex]*sex + ///
           _b[comorbbin]*comorbbin + ///
           _b[mny]*mny + ///
           _b[asasmri]*asasmri + ///
           _b[hla]*hla + ///
           _b[pertvt1]*pertvt1 + ///
           _b[ibdbl]*ibdbl + ///
           _b[emmtvt1]*emmtvt1 + ///
           _b[comedtvt1]*comedtvt1 + ///
           _b[asdastotalt1]*asdastotalt1 + ///
           _b[basfitotalt1]*basfitotalt1

* Y1M0: treated, ASDAS at its untreated value
gen Y1M0 = _b[_cons] + ///
           _b[$M]*M0 + ///
           _b[$A]*1 + ///
           _b[age]*age + ///
           _b[sex]*sex + ///
           _b[comorbbin]*comorbbin + ///
           _b[mny]*mny + ///
           _b[asasmri]*asasmri + ///
           _b[hla]*hla + ///
           _b[pertvt1]*pertvt1 + ///
           _b[ibdbl]*ibdbl + ///
           _b[emmtvt1]*emmtvt1 + ///
           _b[comedtvt1]*comedtvt1 + ///
           _b[asdastotalt1]*asdastotalt1 + ///
           _b[basfitotalt1]*basfitotalt1

* Y0M1: untreated, ASDAS at its treated value
gen Y0M1 = _b[_cons] + ///
           _b[$M]*M1 + ///
           _b[$A]*0 + ///
           _b[age]*age + ///
           _b[sex]*sex + ///
           _b[comorbbin]*comorbbin + ///
           _b[mny]*mny + ///
           _b[asasmri]*asasmri + ///
           _b[hla]*hla + ///
           _b[pertvt1]*pertvt1 + ///
           _b[ibdbl]*ibdbl + ///
           _b[emmtvt1]*emmtvt1 + ///
           _b[comedtvt1]*comedtvt1 + ///
           _b[asdastotalt1]*asdastotalt1 + ///
           _b[basfitotalt1]*basfitotalt1

* Y0M0: untreated, ASDAS at its untreated value
gen Y0M0 = _b[_cons] + ///
           _b[$M]*M0 + ///
           _b[$A]*0 + ///
           _b[age]*age + ///
           _b[sex]*sex + ///
           _b[comorbbin]*comorbbin + ///
           _b[mny]*mny + ///
           _b[asasmri]*asasmri + ///
           _b[hla]*hla + ///
           _b[pertvt1]*pertvt1 + ///
           _b[ibdbl]*ibdbl + ///
           _b[emmtvt1]*emmtvt1 + ///
           _b[comedtvt1]*comedtvt1 + ///
           _b[asdastotalt1]*asdastotalt1 + ///
           _b[basfitotalt1]*basfitotalt1

/////////// Step 3 — Mediation effects

qui sum Y1M1
scalar mean_Y1M1 = r(mean)
qui sum Y1M0
scalar mean_Y1M0 = r(mean)
qui sum Y0M1
scalar mean_Y0M1 = r(mean)
qui sum Y0M0
scalar mean_Y0M0 = r(mean)

gen NIE = Y0M1-Y0M0
gen TIE = Y1M1-Y1M0
gen AIE = (NIE+TIE)/2
qui sum AIE
return scalar mean_AIE = r(mean)

gen NDE = Y1M0-Y0M0
gen TDE = Y1M1-Y0M1
gen ADE = (NDE+TDE)/2
qui sum ADE
return scalar mean_ADE = r(mean)

return scalar ate    = mean_Y1M1 - mean_Y0M0
return scalar mean11 = mean_Y1M1
return scalar mean00 = mean_Y0M0

/////////// Clean up temporary variables before the next replication

capture drop M1 M0 Y1M1 Y1M0 Y0M1 Y0M0 NIE TIE AIE NDE TDE ADE _est_mediator _est_outcome

end

set seed `seed'
bootstrap ATE=r(ate) AIE=r(mean_AIE) ADE=r(mean_ADE) Mean_Y1M1=r(mean11) Mean_Y0M0=r(mean00) ///
          , reps(`sims') nodots: run_gformula_t2

noi estat bootstrap, percentile

savebs "T2" "Single mediator - bDMARD at least 3 months (manual)"
savediag "T2" "Single mediator - bDMARD at least 3 months (manual)" "N treated" `ntreat'
savediag "T2" "Single mediator - bDMARD at least 3 months (manual)" "N untreated" `nuntreat'


*===============================================================================
* B. Two causally ordered mediators (CRP, then ASDAS-PRO): Box S2.1, manual,
*    main paths, without interaction
*===============================================================================

**** Pre-treatment baseline confounders (baseline ASDAS replaced by the baseline mediators)
global W="age sex comorbbin mny asasmri hla pertvt1 ibdbl emmtvt1 comedtvt1 basfitotalt1"

**** Baseline values of mediators
global bM1="crpt1"
global bM2="asdastotalt1_pro"

**** Mediators (continuous)
global M1="crpt2"
global M2="asdastotalt2_pro"

**** Outcome (continuous)
global Y="basfitotalt2"

capture program drop run_gformula_paths_t2
program define run_gformula_paths_t2, rclass

	///////////////// Step 1 — Model the observed data, in causal order

	regress $Y $A $M1 $M2 $W $bM1 $bM2
	estimates store outcome

	regress $M2 $A $M1 $W $bM1 $bM2
	estimates store mediator2

	regress $M1 $A $W $bM1 $bM2
	estimates store mediator1

	///////////////// Step 2 — Monte Carlo simulation of nested counterfactuals

	* Step 2.1 — Counterfactual M1 under treatment (1) and control (0)
	estimates restore mediator1

	gen M1_1 = _b[_cons] + _b[$A]*1 + ///
	           _b[age]*age + _b[sex]*sex + _b[comorbbin]*comorbbin + ///
	           _b[mny]*mny + _b[asasmri]*asasmri + _b[hla]*hla + ///
	           _b[pertvt1]*pertvt1 + _b[ibdbl]*ibdbl + ///
	           _b[emmtvt1]*emmtvt1 + _b[comedtvt1]*comedtvt1 + ///
	           _b[$bM1]*$bM1 + _b[$bM2]*$bM2 + ///
	           _b[basfitotalt1]*basfitotalt1

	gen M1_0 = _b[_cons] + _b[$A]*0 + ///
	           _b[age]*age + _b[sex]*sex + _b[comorbbin]*comorbbin + ///
	           _b[mny]*mny + _b[asasmri]*asasmri + _b[hla]*hla + ///
	           _b[pertvt1]*pertvt1 + _b[ibdbl]*ibdbl + ///
	           _b[emmtvt1]*emmtvt1 + _b[comedtvt1]*comedtvt1 + ///
	           _b[$bM1]*$bM1 + _b[$bM2]*$bM2 + ///
	           _b[basfitotalt1]*basfitotalt1

	* Step 2.2 — Counterfactual M2 at the 3 combinations needed
	estimates restore mediator2

	gen M2_A0_M1_0 = _b[_cons] + _b[$A]*0 + _b[$M1]*M1_0 + ///
	                 _b[age]*age + _b[sex]*sex + _b[comorbbin]*comorbbin + ///
	                 _b[mny]*mny + _b[asasmri]*asasmri + _b[hla]*hla + ///
	                 _b[pertvt1]*pertvt1 + _b[ibdbl]*ibdbl + ///
	                 _b[emmtvt1]*emmtvt1 + _b[comedtvt1]*comedtvt1 + ///
	                 _b[$bM1]*$bM1 + _b[$bM2]*$bM2 + ///
	                 _b[basfitotalt1]*basfitotalt1

	gen M2_A0_M1_1 = _b[_cons] + _b[$A]*0 + _b[$M1]*M1_1 + ///
	                 _b[age]*age + _b[sex]*sex + _b[comorbbin]*comorbbin + ///
	                 _b[mny]*mny + _b[asasmri]*asasmri + _b[hla]*hla + ///
	                 _b[pertvt1]*pertvt1 + _b[ibdbl]*ibdbl + ///
	                 _b[emmtvt1]*emmtvt1 + _b[comedtvt1]*comedtvt1 + ///
	                 _b[$bM1]*$bM1 + _b[$bM2]*$bM2 + ///
	                 _b[basfitotalt1]*basfitotalt1

	gen M2_A1_M1_1 = _b[_cons] + _b[$A]*1 + _b[$M1]*M1_1 + ///
	                 _b[age]*age + _b[sex]*sex + _b[comorbbin]*comorbbin + ///
	                 _b[mny]*mny + _b[asasmri]*asasmri + _b[hla]*hla + ///
	                 _b[pertvt1]*pertvt1 + _b[ibdbl]*ibdbl + ///
	                 _b[emmtvt1]*emmtvt1 + _b[comedtvt1]*comedtvt1 + ///
	                 _b[$bM1]*$bM1 + _b[$bM2]*$bM2 + ///
	                 _b[basfitotalt1]*basfitotalt1

	* Step 2.3 — Counterfactual BASFI at the 4 nested combinations
	estimates restore outcome

	gen Y_000 = _b[_cons] + _b[$A]*0 + _b[$M1]*M1_0 + _b[$M2]*M2_A0_M1_0 + ///
	            _b[age]*age + _b[sex]*sex + _b[comorbbin]*comorbbin + ///
	            _b[mny]*mny + _b[asasmri]*asasmri + _b[hla]*hla + ///
	            _b[pertvt1]*pertvt1 + _b[ibdbl]*ibdbl + ///
	            _b[emmtvt1]*emmtvt1 + _b[comedtvt1]*comedtvt1 + ///
	            _b[$bM1]*$bM1 + _b[$bM2]*$bM2 + ///
	            _b[basfitotalt1]*basfitotalt1

	gen Y_100 = _b[_cons] + _b[$A]*1 + _b[$M1]*M1_0 + _b[$M2]*M2_A0_M1_0 + ///
	            _b[age]*age + _b[sex]*sex + _b[comorbbin]*comorbbin + ///
	            _b[mny]*mny + _b[asasmri]*asasmri + _b[hla]*hla + ///
	            _b[pertvt1]*pertvt1 + _b[ibdbl]*ibdbl + ///
	            _b[emmtvt1]*emmtvt1 + _b[comedtvt1]*comedtvt1 + ///
	            _b[$bM1]*$bM1 + _b[$bM2]*$bM2 + ///
	            _b[basfitotalt1]*basfitotalt1

	gen Y_110 = _b[_cons] + _b[$A]*1 + _b[$M1]*M1_1 + _b[$M2]*M2_A0_M1_1 + ///
	            _b[age]*age + _b[sex]*sex + _b[comorbbin]*comorbbin + ///
	            _b[mny]*mny + _b[asasmri]*asasmri + _b[hla]*hla + ///
	            _b[pertvt1]*pertvt1 + _b[ibdbl]*ibdbl + ///
	            _b[emmtvt1]*emmtvt1 + _b[comedtvt1]*comedtvt1 + ///
	            _b[$bM1]*$bM1 + _b[$bM2]*$bM2 + ///
	            _b[basfitotalt1]*basfitotalt1

	gen Y_111 = _b[_cons] + _b[$A]*1 + _b[$M1]*M1_1 + _b[$M2]*M2_A1_M1_1 + ///
	            _b[age]*age + _b[sex]*sex + _b[comorbbin]*comorbbin + ///
	            _b[mny]*mny + _b[asasmri]*asasmri + _b[hla]*hla + ///
	            _b[pertvt1]*pertvt1 + _b[ibdbl]*ibdbl + ///
	            _b[emmtvt1]*emmtvt1 + _b[comedtvt1]*comedtvt1 + ///
	            _b[$bM1]*$bM1 + _b[$bM2]*$bM2 + ///
	            _b[basfitotalt1]*basfitotalt1

	///////////////// Step 3 — Path-specific effects (Type I decomposition)

	qui sum Y_000
	scalar m000 = r(mean)
	qui sum Y_100
	scalar m100 = r(mean)
	qui sum Y_110
	scalar m110 = r(mean)
	qui sum Y_111
	scalar m111 = r(mean)

	return scalar direct = m100 - m000          // bDMARD -> BASFI
	return scalar pse_m1 = m110 - m100          // bDMARD -> CRP ~> BASFI
	return scalar pse_m2 = m111 - m110          // bDMARD -> ASDAS-PRO -> BASFI
	return scalar total  = m111 - m000          // total effect

	///////////////// Clean up before the next bootstrap replication

	capture drop M1_1 M1_0 M2_A0_M1_0 M2_A0_M1_1 M2_A1_M1_1 ///
	             Y_000 Y_100 Y_110 Y_111 _est_outcome _est_mediator1 _est_mediator2

end

set seed `seed'
bootstrap Direct=r(direct) PSE_M1=r(pse_m1) PSE_M2=r(pse_m2) Total=r(total), ///
	reps(`sims') nodots: run_gformula_paths_t2

noi estat bootstrap, percentile

savebs "T2" "Two mediators - bDMARD at least 3 months (manual)"
savediag "T2" "Two mediators - bDMARD at least 3 months (manual)" "N treated" `ntreat'
savediag "T2" "Two mediators - bDMARD at least 3 months (manual)" "N untreated" `nuntreat'

savesheet "T2"


} // close 8.2



******************************************************************************************************************************************************
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
**** 8.3 - SENSITIVITY (Figure 3B): total effect with MSM, stabilised IPTW x stabilised censoring weights (IPCW) (manual), n = 481
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
******************************************************************************************************************************************************

if `run' == 8.3 {

*===============================================================================
* REFERENCE:    Supplementary Box S7 (Box S3 with censoring weights added)
* METHOD:       Marginal structural model (MSM) weighted by the product of the
*               stabilised IPTW and the stabilised inverse probability of
*               censoring weight (IPCW)
* ESTIMAND:     Average Total Effect (ATE) of bDMARDs on BASFI at 6 months in
*               the eligible patients with complete baseline data
*
* POPULATION:   481 eligible patients with complete baseline covariates ($W).
*               Censored = ASDAS or BASFI missing at 6 months (62 patients), so
*               the 419 observed are exactly the patients of the main analysis:
*               the censoring weights reweight them to represent all 481.
*
* ASSUMPTION:   Missing ASDAS or BASFI at 6 months depends only on treatment and the
*               baseline confounders (the variables in the censoring model).
*
* CENSORING MODEL: $A and $W without ibdbl. All patients with baseline IBD have
*               ASDAS and BASFI at 6 months, so ibdbl predicts being observed perfectly and
*               logit would drop it (and the patients) from the model. It is left out on
*               purpose so every patient keeps a predicted probability.
*===============================================================================

assert `cohort' == 800 // defined for the two-visit cohort only

*-------------------------------------------------------------------------------
* Eligible patients with complete baseline covariates, one row per patient
*-------------------------------------------------------------------------------

import delimited "${datasets}originalfulllong800.csv", clear varnames(1) case(preserve)
keep if t == 6 // one row per patient

**** Pre-treatment baseline confounders, treatment and outcome (as in the main analysis)
global W="age sex comorbbin mny asasmri hla pertvt1 ibdbl emmtvt1 comedtvt1 asdastotalt1 basfitotalt1"
global A="bionew"
global Y="basfitotalt2"

**** Covariates of the censoring model ($W without ibdbl, see above)
global Wc="age sex comorbbin mny asasmri hla pertvt1 emmtvt1 comedtvt1 asdastotalt1 basfitotalt1"

local nW : word count $W
egen byte nbase = rownonmiss($W)
gen byte basecomplete = 0
replace basecomplete = 1 if nbase == `nW'
keep if basecomplete == 1

qui count
assert r(N) == 481 // stops here if the source file or the selection changes
global anaN = r(N)

*-------------------------------------------------------------------------------
* Censoring indicator: 1 = ASDAS and BASFI observed at 6 months (the main analysis)
*-------------------------------------------------------------------------------

gen byte uncensored = 0
replace uncensored = 1 if !missing($Y) & !missing(asdastotalt2)

qui count if uncensored == 0
assert r(N) == 62
local ncens = r(N)
qui count if uncensored == 0 & $A == 1
local ncens1 = r(N)
qui count if uncensored == 0 & $A == 0
local ncens0 = r(N)

* ibdbl can only be left out of the censoring model while no IBD patient is censored
qui count if uncensored == 0 & ibdbl == 1
assert r(N) == 0

noi di ""
noi di as text "Censored (ASDAS or BASFI missing at 6 months): " as result `ncens' as text " (treated " as result `ncens1' as text ", untreated " as result `ncens0' as text ")"
noi di ""

///////////////// Step 1 — Model the Observed Data

**** Treatment model for the denominator (all 481)
logit $A $W
estimates store denominator

**** Treatment model for the numerator
logit $A
estimates store numerator

**** Censoring model for the denominator: P(uncensored | treatment, baseline)
logit uncensored $A $Wc
estimates store censoring

**** Censoring model for the numerator: P(uncensored | treatment)
logit uncensored $A
estimates store censoring_num

///////////////// Step 2 — Weighting

///////////////// Step 2.1. Stabilised IPTW (as in Box S3)

estimates restore denominator
predict ps, pr
gen denom = ps*$A + (1-ps)*(1-$A)

estimates restore numerator
predict num_ps, pr
gen num = num_ps*$A + (1-num_ps)*(1-$A)

gen stabweight = num/denom

///////////////// Step 2.2. Stabilised IPCW

estimates restore censoring

* pc1: probability of being observed if treated
gen pc1 = invlogit(_b[_cons] + ///
                   _b[$A]*1 + ///
                   _b[age]*age + ///
                   _b[sex]*sex + ///
                   _b[comorbbin]*comorbbin + ///
                   _b[mny]*mny + ///
                   _b[asasmri]*asasmri + ///
                   _b[hla]*hla + ///
                   _b[pertvt1]*pertvt1 + ///
                   _b[emmtvt1]*emmtvt1 + ///
                   _b[comedtvt1]*comedtvt1 + ///
                   _b[asdastotalt1]*asdastotalt1 + ///
                   _b[basfitotalt1]*basfitotalt1)

* pc0: probability of being observed if untreated
gen pc0 = invlogit(_b[_cons] + ///
                   _b[$A]*0 + ///
                   _b[age]*age + ///
                   _b[sex]*sex + ///
                   _b[comorbbin]*comorbbin + ///
                   _b[mny]*mny + ///
                   _b[asasmri]*asasmri + ///
                   _b[hla]*hla + ///
                   _b[pertvt1]*pertvt1 + ///
                   _b[emmtvt1]*emmtvt1 + ///
                   _b[comedtvt1]*comedtvt1 + ///
                   _b[asdastotalt1]*asdastotalt1 + ///
                   _b[basfitotalt1]*basfitotalt1)

* pc: probability of being observed under the treatment actually received
gen pc = pc1*$A + pc0*(1-$A)

estimates restore censoring_num
predict num_c, pr

gen stabweight_c = num_c/pc

///////////////// Step 2.3. Combined weight (used only for the uncensored)

gen weight_tc = stabweight*stabweight_c // not 'weight': the dataset already has a variable weight (body weight)

** Diagnostics
qui sum pc if $A == 1
local pcmin1 = r(min)
qui sum pc if $A == 0
local pcmin0 = r(min)
qui sum stabweight_c if uncensored == 1
local swcmean = r(mean)
local swcmin  = r(min)
local swcmax  = r(max)
qui sum weight_tc if uncensored == 1
local wmean = r(mean)
local wmin  = r(min)
local wmax  = r(max)

noi di as text "Min P(observed): treated " as result %5.3f `pcmin1' as text ", untreated " as result %5.3f `pcmin0'
noi di as text "sIPCW: mean " as result %5.3f `swcmean' as text ", min " as result %5.3f `swcmin' as text ", max " as result %5.3f `swcmax'
noi di as text "sIPTW x sIPCW: mean " as result %5.3f `wmean' as text ", min " as result %5.3f `wmin' as text ", max " as result %5.3f `wmax'

///////////////// Step 3 — Total Treatment Effect and 95% CI: MSM in the uncensored

noi regress $Y $A [pw=weight_tc] if uncensored == 1, robust

savest "T2" "MSM with sIPTW and sIPCW (manual)" "ATE" `=_b[$A]' `=_se[$A]' `=_b[$A]-1.96*_se[$A]' `=_b[$A]+1.96*_se[$A]' . . $anaN
savediag "T2" "MSM with sIPTW and sIPCW (manual)" "N censored" `ncens'
savediag "T2" "MSM with sIPTW and sIPCW (manual)" "N censored treated" `ncens1'
savediag "T2" "MSM with sIPTW and sIPCW (manual)" "N censored untreated" `ncens0'
savediag "T2" "MSM with sIPTW and sIPCW (manual)" "Min P(observed) treated" `pcmin1'
savediag "T2" "MSM with sIPTW and sIPCW (manual)" "Min P(observed) untreated" `pcmin0'
savediag "T2" "MSM with sIPTW and sIPCW (manual)" "sIPCW mean" `swcmean'
savediag "T2" "MSM with sIPTW and sIPCW (manual)" "sIPCW max" `swcmax'
savediag "T2" "MSM with sIPTW and sIPCW (manual)" "sIPTW x sIPCW max" `wmax'
savesheet "T2"


} // close 8.3


******************************************************************************************************************************************************
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
**** 8.4 - SENSITIVITY (Figure 3B): total effect with TMLE including censoring (manual), n = 481
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
******************************************************************************************************************************************************

if `run' == 8.4 {

*===============================================================================
* REFERENCE:    Supplementary Box S7 (Box S4 with censoring added)
* METHOD:       TMLE in which the treatment mechanism g is the product of the
*               probability of the treatment received and the probability of
*               being observed (uncensored) given treatment and baseline
* ESTIMAND:     Average Total Effect (ATE) of bDMARDs on BASFI at 6 months in
*               the eligible patients with complete baseline data
*
* POPULATION, CENSORING AND CENSORING MODEL: as in run 8.3
*
* CHANGES FROM BOX S4:
*   - the outcome model is fitted in the uncensored
*   - the clever covariate is H = C x [A/(PS x pc1) - (1-A)/((1-PS) x pc0)],
*     with C = 1 if uncensored: it is zero for the censored
*   - the fluctuation is fitted in the uncensored; the counterfactual
*     predictions and the effect are averaged over all 481
*===============================================================================

assert `cohort' == 800 // defined for the two-visit cohort only

*-------------------------------------------------------------------------------
* Eligible patients with complete baseline covariates, one row per patient
*-------------------------------------------------------------------------------

import delimited "${datasets}originalfulllong800.csv", clear varnames(1) case(preserve)
keep if t == 6 // one row per patient

**** Pre-treatment baseline confounders, treatment and outcome (as in the main analysis)
global W="age sex comorbbin mny asasmri hla pertvt1 ibdbl emmtvt1 comedtvt1 asdastotalt1 basfitotalt1"
global A="bionew"
global Y="basfitotalt2"

**** Covariates of the censoring model ($W without ibdbl, see run 8.3)
global Wc="age sex comorbbin mny asasmri hla pertvt1 emmtvt1 comedtvt1 asdastotalt1 basfitotalt1"

local nW : word count $W
egen byte nbase = rownonmiss($W)
gen byte basecomplete = 0
replace basecomplete = 1 if nbase == `nW'
keep if basecomplete == 1

qui count
assert r(N) == 481 // stops here if the source file or the selection changes
global anaN = r(N)

gen byte uncensored = 0
replace uncensored = 1 if !missing($Y) & !missing(asdastotalt2)

qui count if uncensored == 0
assert r(N) == 62
local ncens = r(N)

qui count if uncensored == 0 & ibdbl == 1
assert r(N) == 0

///////////////// Step 1 — Model the Observed Data

**** Outcome model, in the uncensored
qui regress $Y $A $W if uncensored == 1
estimates store outcome

**** Treatment model (all 481)
qui logit $A $W
estimates store treatment

**** Censoring model (all 481)
qui logit uncensored $A $Wc
estimates store censoring

///////////////// Step 2 — Targeting

///////////////// Step 2.1. Counterfactual outcome (BASFI) at 6 months, all 481

estimates restore outcome

* Y1: BASFI at 6 months if everyone were treated
gen Y1 = _b[_cons] + ///
         _b[$A]*1 + ///
         _b[age]*age + ///
         _b[sex]*sex + ///
         _b[comorbbin]*comorbbin + ///
         _b[mny]*mny + ///
         _b[asasmri]*asasmri + ///
         _b[hla]*hla + ///
         _b[pertvt1]*pertvt1 + ///
         _b[ibdbl]*ibdbl + ///
         _b[emmtvt1]*emmtvt1 + ///
         _b[comedtvt1]*comedtvt1 + ///
         _b[asdastotalt1]*asdastotalt1 + ///
         _b[basfitotalt1]*basfitotalt1

* Y0: BASFI at 6 months if everyone were untreated
gen Y0 = _b[_cons] + ///
         _b[$A]*0 + ///
         _b[age]*age + ///
         _b[sex]*sex + ///
         _b[comorbbin]*comorbbin + ///
         _b[mny]*mny + ///
         _b[asasmri]*asasmri + ///
         _b[hla]*hla + ///
         _b[pertvt1]*pertvt1 + ///
         _b[ibdbl]*ibdbl + ///
         _b[emmtvt1]*emmtvt1 + ///
         _b[comedtvt1]*comedtvt1 + ///
         _b[asdastotalt1]*asdastotalt1 + ///
         _b[basfitotalt1]*basfitotalt1

* Yobs: BASFI at 6 months under the treatment actually received
gen Yobs = _b[_cons] + ///
           _b[$A]*$A + ///
           _b[age]*age + ///
           _b[sex]*sex + ///
           _b[comorbbin]*comorbbin + ///
           _b[mny]*mny + ///
           _b[asasmri]*asasmri + ///
           _b[hla]*hla + ///
           _b[pertvt1]*pertvt1 + ///
           _b[ibdbl]*ibdbl + ///
           _b[emmtvt1]*emmtvt1 + ///
           _b[comedtvt1]*comedtvt1 + ///
           _b[asdastotalt1]*asdastotalt1 + ///
           _b[basfitotalt1]*basfitotalt1

///////////////// Step 2.2. Propensity score

estimates restore treatment
predict ps_tmle, pr

///////////////// Step 2.3. Probability of being observed if treated (pc1) and if untreated (pc0)

estimates restore censoring

gen pc1 = invlogit(_b[_cons] + ///
                   _b[$A]*1 + ///
                   _b[age]*age + ///
                   _b[sex]*sex + ///
                   _b[comorbbin]*comorbbin + ///
                   _b[mny]*mny + ///
                   _b[asasmri]*asasmri + ///
                   _b[hla]*hla + ///
                   _b[pertvt1]*pertvt1 + ///
                   _b[emmtvt1]*emmtvt1 + ///
                   _b[comedtvt1]*comedtvt1 + ///
                   _b[asdastotalt1]*asdastotalt1 + ///
                   _b[basfitotalt1]*basfitotalt1)

gen pc0 = invlogit(_b[_cons] + ///
                   _b[$A]*0 + ///
                   _b[age]*age + ///
                   _b[sex]*sex + ///
                   _b[comorbbin]*comorbbin + ///
                   _b[mny]*mny + ///
                   _b[asasmri]*asasmri + ///
                   _b[hla]*hla + ///
                   _b[pertvt1]*pertvt1 + ///
                   _b[emmtvt1]*emmtvt1 + ///
                   _b[comedtvt1]*comedtvt1 + ///
                   _b[asdastotalt1]*asdastotalt1 + ///
                   _b[basfitotalt1]*basfitotalt1)

///////////////// Step 2.4. Clever covariate (zero for the censored)

* Uncensored treated:   H =  1/(PS x pc1)
* Uncensored untreated: H = -1/((1-PS) x pc0)
* Censored:             H =  0
gen H = 0
replace H = $A/(ps_tmle*pc1) - (1-$A)/((1-ps_tmle)*pc0) if uncensored == 1

///////////////// Step 2.5. Targeting step, in the uncensored

gen Y_resid = $Y - Yobs // used in the uncensored only
quietly regress Y_resid H if uncensored == 1, nocons
scalar epsilon = _b[H]

gen Y1_target   = Y1   + scalar(epsilon)/(ps_tmle*pc1)
gen Y0_target   = Y0   - scalar(epsilon)/((1-ps_tmle)*pc0)
gen Yobs_target = Yobs + scalar(epsilon)*H

///////////////// Step 3 — Total Treatment Effect, averaged over all 481

qui sum Y1_target
scalar mean_Y1 = r(mean)

qui sum Y0_target
scalar mean_Y0 = r(mean)

scalar ATE_tmle = mean_Y1 - mean_Y0

///////////////// Step 4 — 95% CI via the efficient influence function (EIF)

* EIF = H x (Y - Yobs_target) + Y1_target - Y0_target - ATE
* The first term is zero for the censored (H = 0, Y not observed)
gen EIF = Y1_target - Y0_target - scalar(ATE_tmle)
replace EIF = EIF + H*($Y - Yobs_target) if uncensored == 1

qui count
scalar n = r(N)

gen EIF2 = EIF^2
qui sum EIF2
scalar var_ATE = r(sum)/(scalar(n)^2)
scalar se_ATE  = sqrt(scalar(var_ATE))

scalar lb_ATE = ATE_tmle - 1.96*se_ATE
scalar ub_ATE = ATE_tmle + 1.96*se_ATE

noi di ""
noi di as text "TMLE with censoring (n = " as result $anaN as text ", censored " as result `ncens' as text ")"
noi di as text "  E[Y1] = " as result %6.3f scalar(mean_Y1) as text "   E[Y0] = " as result %6.3f scalar(mean_Y0)
noi di as text "  ATE   = " as result %6.3f scalar(ATE_tmle) as text "   SE " as result %6.3f scalar(se_ATE) ///
       as text "   95% CI (" as result %6.3f scalar(lb_ATE) as text "; " as result %6.3f scalar(ub_ATE) as text ")"

///////////////// Diagnostics

* g = PS x pc for the treated, (1-PS) x pc for the untreated: smallest value
gen gc = ps_tmle*pc1*$A + (1-ps_tmle)*pc0*(1-$A)
qui sum gc
local gcmin = r(min)

qui sum EIF
local sd_IC = r(sd)

noi di as text "  smallest g x P(observed): " as result %6.4f `gcmin' as text "   SD influence curve: " as result %6.3f `sd_IC'

savest "T2" "TMLE with censoring (manual)" "ATE" `=scalar(ATE_tmle)' `=scalar(se_ATE)' `=scalar(lb_ATE)' `=scalar(ub_ATE)' . . $anaN
savediag "T2" "TMLE with censoring (manual)" "N censored" `ncens'
savediag "T2" "TMLE with censoring (manual)" "Mean Y1" `=scalar(mean_Y1)'
savediag "T2" "TMLE with censoring (manual)" "Mean Y0" `=scalar(mean_Y0)'
savediag "T2" "TMLE with censoring (manual)" "Min g x P(observed)" `gcmin'
savediag "T2" "TMLE with censoring (manual)" "SD influence curve" `sd_IC'
savesheet "T2"


} // close 8.4



} // close overall qui
