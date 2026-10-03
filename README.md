# Causal-mediation-in-METEOR

This repository accompanies the manuscript:

**"Biological DMARDs improve function by suppressing disease activity in patients with axial spondyloarthritis: a causal mediation analysis in the METEOR registry."**

The code provided here is intended **for illustrative purposes only**. Its purpose is to help readers understand the intuition and computational steps behind the causal inference methods described in the manuscript and supplementary material.

The scripts reproduce the illustrative algorithms presented in:

* Supplementary Box S1: Parametric g-formula for the total treatment effect
* Supplementary Box S2: Causal mediation analysis with one mediator (ASDAS)
* Supplementary Box S2.1: Causal mediation analysis with two causally ordered mediators (CRP, then ASDAS-PRO; edge g-formula)
* Supplementary Box S3: Marginal structural models (MSM) with inverse probability of treatment weights
* Supplementary Box S4: Targeted maximum likelihood estimation (TMLE)
* Supplementary Box S5: E-value
* Supplementary Box S6: Causal mediation analysis accounting for post-treatment mediator–outcome confounding
* Supplementary Box S7: Censoring weights in MSM and TMLE (missing data at 6 months)

## Files

* `Analysis.do`: Stata script (version 15.1). Each analysis is selected with `local run` at the top of the file.
* `Analysis.R`: R script (version 4.5.1), with the same analyses and the code for Figures 2 and 3.
* `Stata Programs/`: supporting Stata programs called by `Analysis.do`
  * `saveestV1.do`: collects the estimates and diagnostics for the supplementary tables
  * `msmbalanceV1.do`: covariate balance before and after weighting (Supplementary Table S8)
  * `smdtableV2.do`: descriptive tables (Table 1, Supplementary Tables S1 and S11)

Before running, set the working folder at the top of each script (`C:\mypath\` in Stata, `C:/mypath/` in R). The scripts expect the subfolders `Datasets`, `Tables`, `Figures` and `Stata Programs` inside that folder.

The METEOR data are not included in this repository. They are available upon reasonable request, as described in the data availability statement of the manuscript.

## Software

Each analysis is implemented manually, following the algorithms in the Supplementary Boxes, and with existing programs where available.

These scripts are **not validated statistical software** and should **not be used as implementations for applied research**.

For actual analyses, readers are referred to available software:

* **Parametric g-formula and mediation analysis**

  * Stata: `gformula` and `medeff` (user-written programs)
  * R: `mediation`

* **Path-specific effects with two causally ordered mediators**

  * R: `paths`

* **Marginal structural models**

  * Stata: `teffects ipw`
  * R: `ipw` and `survey`

* **Targeted maximum likelihood estimation**

  * R: `ltmle`
  * Stata: `teffects aipw` (augmented inverse probability weighting, a doubly robust estimator closely related to TMLE)

* **Tables and figures in R**

  * `writexl` and `ggplot2`

The methodological references, software references, and links to the original implementations are provided in the supplementary material (section "Main references for each method and program").

If you use this repository, please cite the accompanying manuscript and the original methodological references.
