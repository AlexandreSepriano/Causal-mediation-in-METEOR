# Causal-mediation-in-METEOR

## Overview

This repository accompanies the manuscript:

**"Biological DMARDs improve function by suppressing inflammation in patients with axial spondyloarthritis: a causal mediation analysis in the METEOR registry."**

The purpose of this repository is to provide  annotated examples of the causal inference methods used in the study, including:

* Parametric g-formula
* Causal mediation analysis
* Marginal structural models (MSMs)
* Targeted maximum likelihood estimation (TMLE)
* Diagnostic procedures and sensitivity analyses

The code is intended to help readers understand the intuition behind these methods and the steps involved in their implementation.

---

## Important Disclaimer

The code in this repository was developed primarily as an educational resource to accompany the manuscript and its supplementary material.

It is intended to:

* Illustrate the logic of the algorithms.
* Demonstrate how causal estimands are identified and estimated.
* Facilitate understanding of the methodological concepts described in the paper.
* Allow interested readers to reproduce the analyses presented in the manuscript.

### Not a validated software implementation

The scripts provided here **should not be considered validated statistical software** and **should not be used as a reference implementation for applied research, clinical decision-making, regulatory submissions, or production analyses**.

While every effort was made to ensure that the code reproduces the analyses reported in the manuscript, the repository was not developed, tested, or maintained according to the standards of dedicated methodological software packages.

Users are responsible for verifying all results independently.

---

## Recommended Software for Applied Analyses

For actual research applications, users should rely on established and validated software implementations maintained by the methodological community.

### Parametric g-formula and mediation analysis

**Stata**

* `gformula` (Daniel et al., 2011)
* `medeff` (Hicks & Tingley, 2011)

**R**

* `mediation` package (Tingley et al., 2014)

### Marginal Structural Models (IPTW)

**Stata**

* `teffects ipw`

### Targeted Maximum Likelihood Estimation (TMLE)

**R**

* `ltmle` package (Lendle et al., 2017)

**Stata**

* `teffects aipw` provides an augmented inverse probability weighting (AIPW) estimator, which is closely related to TMLE but is not a TMLE implementation.

---

## Relation to the Manuscript

The algorithms implemented in this repository correspond to the methodological descriptions provided in the supplementary material of the manuscript.

In particular, the repository contains illustrative implementations of:

* Parametric g-formula for total treatment effects
* Mediation analysis using the mediation formula
* Parametric g-formula mediation with and without post-treatment mediator–outcome confounding
* Marginal structural models using inverse probability weighting
* TMLE-style estimation procedures
* Diagnostic and sensitivity analyses

These implementations were developed to enhance transparency and facilitate learning rather than to replace established software packages.

---

## Reproducibility

The code is provided to facilitate reproducibility of the analyses reported in the manuscript.

Researchers wishing to reproduce the published results should carefully consult:

1. The manuscript.
2. The supplementary material.
3. The original methodological references listed below.

Small numerical differences between implementations may arise because different software packages use different estimation procedures, simulation approaches, bootstrap settings, or variance estimators.

---

## Key References

### Parametric g-formula

Robins JM. A new approach to causal inference in mortality studies with a sustained exposure period—application to control of the healthy worker survivor effect. *Mathematical Modelling*. 1986;7:1393–1512.

### Causal mediation

Robins JM, Greenland S. Identifiability and exchangeability for direct and indirect effects. *Epidemiology*. 1992;3:143–155.

Pearl J. The causal mediation formula—a guide to the assessment of pathways and mechanisms. *Prevention Science*. 2012;13:426–436.

### Marginal Structural Models

Robins JM, Hernán MA, Brumback B. Marginal structural models and causal inference in epidemiology. *Epidemiology*. 2000;11:550–560.

### TMLE

van der Laan MJ, Rubin D. Targeted maximum likelihood learning. *International Journal of Biostatistics*. 2006;2(1):Article 11.

### Practical Reference

Hernán MA, Robins JM. *Causal Inference: What If*. Chapman & Hall/CRC; 2020.

---

## Citation

If you use material from this repository, please cite the accompanying manuscript and the original methodological references listed above.

---

## Contact

For questions, corrections, or suggestions, please open an issue in this repository.
