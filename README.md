# Causal-mediation-in-METEOR

This repository accompanies the manuscript:

**"Biological DMARDs improve function by suppressing inflammation in patients with axial spondyloarthritis: a causal mediation analysis in the METEOR registry."**

The code provided here is intended **for educational and illustrative purposes only**. Its purpose is to help readers understand the intuition and computational steps behind the causal inference methods described in the manuscript and supplementary material.

The scripts reproduce the illustrative algorithms presented in:

* Supplementary Box S1: Parametric g-formula for the total treatment effect
* Supplementary Box S2: Mediation analysis using the mediation formula
* Supplementary Box S3: Mediation analysis accounting for post-treatment mediator–outcome confounding
* Supplementary Box S4: Marginal structural models (MSM)
* Supplementary Box S5: Targeted maximum likelihood estimation (TMLE)

These scripts are **not validated statistical software** and should **not be used as reference implementations for applied research or production analyses**.

For actual analyses, readers should use established software implementations:

* **Parametric g-formula and mediation analysis**

  * Stata: `gformula` and `medeff`
  * R: `mediation`

* **Marginal structural models**

  * Stata: `teffects ipw`

* **Targeted maximum likelihood estimation**

  * R: `ltmle`
  * Stata: `teffects aipw` (augmented inverse probability weighting, a doubly robust estimator closely related to TMLE)

The methodological references, software references, and links to the original implementations are provided in the supplementary material (Main References section).

If you use this repository, please cite the accompanying manuscript and the original methodological references.
