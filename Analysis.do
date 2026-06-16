******************************************************************************************************************************************************************
**************************************************************** CAUSAL Mediation Analysis in METEOR *************************************************************
******************************************************************************************************************************************************************

global path your path
use "${path}causalmediation.dta", clear


*===============================================================================
* Dataset: wide format
* Visits:  t1 (baseline), t2 (6 months)
*===============================================================================


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
set seed 1234
bootstrap ATE=r(ate) Mean_Y1M1=r(mean11) Mean_Y0M0=r(mean00) ///
          Mean_$A=r(mean_$A) Mean_Ap=r(mean_Ap) Mean_$M=r(mean_$M) Mean_Mp=r(mean_Mp) Mean_$Y=r(mean_$Y) Mean_Yp=r(mean_Yp) ///
		  , reps(1000) nodots: run_gformula

* Display the percentile-based bootstrap 95% Confidence Intervals
estat bootstrap, percentile


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
di "E-value (point estimate): " %6.3f scalar(evalue_point)


*===============================================================================
* PROGRAM:       Stata medeff
*===============================================================================

medeff (regress $M $A $W) (regress $Y $A $M $W), mediate($M) treat($A) vce(bootstrap, reps(1000))

*===============================================================================
* PROGRAM:       Stata gformula (mediation syntax)
*===============================================================================

gformula $Y $M $A $W, ///
mediation ex($A) mediator($M) out($Y) ///
eq($M: $A $W, $Y: $A $M $W) /// 
com($M:regress, $Y:regress)   ///
obe base_confs($W) ///
seed(1) samples(1000) all



*===============================================================================
* REFERENCE:    Supplementary Box S2 (ignoring  mediator-outcome confounding)
* METHOD:       Parametric Time-Fixed G-Formula Algorithm for causal mediation
* ESTIMAND:     Causal mediation effects of bDMARDs on BASFI at 6 Months
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

* Run the bootstrap across 1000 replications (second line of bootstrap comman contains only diagnostics and is optional
set seed 1234
bootstrap ATE=r(ate) Mean_Y1M1=r(mean11) Mean_Y1M0=r(mean10) Mean_Y0M1=r(mean01) Mean_Y0M0=r(mean00) ///
		  NIE=r(mean_NIE) TIE=r(mean_TIE) AIE=r(mean_AIE) NDE=r(mean_NDE) TDE=r(mean_TDE) ADE=r(mean_ADE) ///
          Mean_$A=r(mean_$A) Mean_Ap=r(mean_Ap) Mean_$M=r(mean_$M) Mean_Mp=r(mean_Mp) Mean_$Y=r(mean_$Y) Mean_Yp=r(mean_Yp) ///
		  , reps(100) nodots: run_gformula

* Display the percentile-based bootstrap 95% Confidence Intervals
estat bootstrap, percentile

*===============================================================================
* PROGRAM:       Stata medeff
*===============================================================================
gen inter=$M*$A 
medeff (regress $M $A $W) (regress $Y $A $M $W inter), mediate($M) treat($A) vce(bootstrap, reps(1000))  interact(inter)


*===============================================================================
* PROGRAM:       Stata gformula (mediation syntax) 
*===============================================================================

* Not supported



*===============================================================================
* REFERENCE:    Supplementary Box S3 (considering  mediator-outcome confounding)
* METHOD:       Parametric Time-Fixed G-Formula Algorithm for causal mediation
* ESTIMAND:     Causal mediation effects of bDMARDs on BASFI at 6 Months
*===============================================================================

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

capture drop M1 M0 Y1M1 Y1M0 Y0M1 Y0M0 MOCa1 MOCa0 MOCb1 MOCb0 MOCc1 MOCc0 NIE TIE AIE NDE TDE ADE Ap Mp Yp _est_mediator _est_outcome _est_treatment inter

end

* Run the bootstrap across 1000 replications (second line of bootstrap comman contains only diagnostics and is optional
set seed 123
bootstrap ATE=r(ate) Mean_Y1M1=r(mean11) Mean_Y1M0=r(mean10) Mean_Y0M1=r(mean01) Mean_Y0M0=r(mean00) ///
		  NIE=r(mean_NIE) TIE=r(mean_TIE) AIE=r(mean_AIE) NDE=r(mean_NDE) TDE=r(mean_TDE) ADE=r(mean_ADE) ///
          Mean_$A=r(mean_$A) Mean_Ap=r(mean_Ap) Mean_$M=r(mean_$M) Mean_Mp=r(mean_Mp) Mean_$Y=r(mean_$Y) Mean_Yp=r(mean_Yp) ///
		  , reps(1000) nodots: run_gformula // higher number of reps needed because of random error added by MOC

* Display the percentile-based bootstrap 95% Confidence Intervals
estat bootstrap, percentile


*===============================================================================
* PROGRAM:       Stata medeff
*===============================================================================

* Not supported

*===============================================================================
* PROGRAM:       Stata gformula (mediation syntax) without treatment-mediator interaction (not supported)
*===============================================================================


gformula $Y $M $MOCa $MOCb $MOCc $A $W, ///
mediation ex($A) mediator($M) out($Y) post_confs($MOCa $MOCb $MOCc) ///
eq($M: $A $MOCa $MOCb $MOCc $W, $Y: $A $M $MOCa $MOCb $MOCc $W, $MOCa: $A $W, $MOCb: $A $W, $MOCc: $A $W) /// 
com($M:regress, $Y:regress, $MOCa:logit, $MOCb:logit, $MOCc:logit)   ///
obe base_confs($W) ///
seed(12345) samples(2000) all




*===============================================================================
* REFERENCE:    Supplementary Box S4 
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
sum ps 
local ps_mean=round(r(mean), 0.01) 
local ps_sd=round(r(sd), 0.01) 
local ps_min=round(r(min), 0.01) 
local ps_max=round(r(max), 0.01) 

** Calculate the denominator for treated and untreated
gen denom=ps*$A+(1-ps)*(1-$A) 

** Calculate the unstabilized IPTW (uIPTW) for treated and untreated
qui gen unstabweight= 1/denom 

** uIPTW diagnostics
sum unstabweight 
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
sum stabweight 
local stabweightmean=round(r(mean), 0.01) 
local stabweightsd=round(r(sd), 0.01) 
local stabweightmin=round(r(min), 0.01) 
local stabweightmax=round(r(max), 0.01) 

///////////////// Step 3 and — Total Treatment Effect and 95% CI: MSM

regress $Y $A [pw=stabweight], robust


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

matrix list DIAG, format(%9.3f) title("IPTW Diagnostics")


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

di as text "ACS lower bound: `acsmin'"
di as text "ACS upper bound: `acsmax'"

* Flag observations within/outside ACS
gen acs = .
replace acs = 1 if ps >= `acsmin' & ps <= `acsmax' & ps !=.
replace acs = 0 if (ps < `acsmin' | ps > `acsmax') & ps !=.

* Count observations outside ACS
sum acs
local n_outside = round(r(N) - (r(N) * r(mean)), 1)
local pct_outside = round((1 - r(mean)) * 100, 0.1)

di as text "Observations outside ACS: `n_outside' (`pct_outside'%)"

* Clean up intermediate variables
drop max_ps min_ps

* Total Treatment Effect: MSM within the ACS

regress $Y $A [pw=stabweight] if acs==1, robust


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

regress $Y $A [pw=stabweight_trunc], vce(robust)

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

matrix list TRUNC, format(%9.3f) title("Truncation Diagnostics (1st-99th percentile)")



*===============================================================================
* PROGRAM:       Stata teffects 
*===============================================================================

teffects ipw ($Y) ($A $W, logit), pomeans // potential outcomes 
teffects ipw ($Y) ($A $W, logit), ate // total treatment effect

* Ony uIPTWs are calculated with teffects (equal to sIPTW here because the MSM is saturated)
* SE are lower with teffects (the program uses influence-function-based SE that accounts for the uncertainty in the propensity score estimation step)



*===============================================================================
* REFERENCE:    Supplementary Box S5 
* METHOD:       Targeted Maximum Likelihood Estimation for Total Treatment Effect
*               (ignoring mediator and MOC)
* ESTIMAND:     Average Total Effect (ATE) of bDMARDs on BASFI at 6 Months
*===============================================================================

*===============================================================================
* PROGRAM:       Manual calculation
*===============================================================================

///////////////// Step 1 — Model the Observed Data

**** Outcome model
quietly regress $Y $A $W
estimates store outcome

**** Treatment model
quietly logit $A $W
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

matrix list TMLE, format(%9.3f) title("TMLE — Total Treatment Effect of bDMARDs on BASFI at 6 Months")

di as text _newline "95% CI (EIF-based): [" %6.3f scalar(lb_ATE) ", " %6.3f scalar(ub_ATE) "]"


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

di as text _newline "TMLE diagnostics"
di as text "% g-truncated: `pct_gtrunc'%"
di as text "SD of influence curve: `sd_IC'"

drop ps_tmle EIF g_truncated


///////////////// e-value
local rr = exp(0.91*abs(ATE_tmle/SD_Y))
scalar evalue_point = `rr' + sqrt(`rr'*(`rr'-1))
di "E-value (point estimate): " %6.3f scalar(evalue_point)



*===============================================================================
* PROGRAM:       Stata teffects
*===============================================================================

teffects aipw ($Y $W) ($A $W, logit), ate // Augmented IPW (AIPW) estimator: similar (double robust) but not the same as TMLE


