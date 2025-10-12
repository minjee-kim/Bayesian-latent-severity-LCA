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

![](README_files/unnamed-chunk-1-1.png)<!-- -->
