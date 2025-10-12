Bayesian Latent Severity LCA
================

A **Bayesian latent class model with continuous severity** for
diagnostic test evaluation.

This project extends classical Hui–Walter and random-effects LCA
frameworks by introducing a continuous *latent severity* parameter
$S_i$.  
The approach enables efficient estimation of diagnostic accuracies and
prevalence using **fewer than four assays**, addressing common
identifiability limitations in small-panel diagnostic evaluations.

------------------------------------------------------------------------

## Overview

Traditional latent class models assume conditional independence between
binary tests given a latent disease status.  
The **Latent Severity LCA (S-LCA)** relaxes this by allowing
within-class variability in disease expression:

$$
T_{ij}^* = \beta_j S_i + \gamma_j + \varepsilon_{ij}, 
\qquad T_{ij} = \mathbb{I}(T_{ij}^* > 0)
$$

Here, $S_i$ captures the strength or intensity of infection for diseased
subjects ($D_i = 1$), producing more realistic inference on test
sensitivity, specificity, and prevalence.

------------------------------------------------------------------------

### 1. `bayesian_severity_LCA()`

A unified function implementing two prior choices for the latent
severity parameter $S_i$:

- **Gamma Severity Model**  
  $S_i \mid D_i = 1 \sim \mathrm{Gamma}(\alpha_S, \beta_S)$, 
  $\beta_S \sim \mathrm{Gamma}(a_\beta, b_\beta)$

- **Normal Moment (NM⁺) Severity Model**  
  $p(S_i \mid D_i = 1) \propto S_i^2 \exp[-(S_i - \mu_0)^2 / (2\tau^2)]$, 
  $S_i > 0$,  $\tau^2 \sim \mathrm{Inv\text{-}Gamma}(a_\tau, b_\tau)$

Both models employ MCMC with truncated-normal data augmentation for
latent test variables, and include posterior updates for: - Class
indicator $D_i$ - Test parameters $(\beta_j, \gamma_j, \sigma_j^2)$ -
Prevalence $\rho$ - Severity priors (Gamma or NM+)

------------------------------------------------------------------------

### Candidate Priors for Severity

The following figure compares the Gamma and NM+ priors for $S_i$ when
each has variance = 1:

![](README_files/unnamed-chunk-1-1.png)<!-- -->
