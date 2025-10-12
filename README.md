---
title: "Bayesian Latent Severity LCA"
output: html_document
date: "2025-10-12"
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

A **Bayesian latent class model with continuous severity** for diagnostic test evaluation.

This project extends classical Hui–Walter and random-effects LCA frameworks by introducing a continuous *latent severity* parameter \(S_i\).  
The approach enables efficient estimation of diagnostic accuracies and prevalence using **fewer than four assays**, addressing common identifiability limitations in small-panel diagnostic evaluations.
