
# Joseph et al. (1995) Strongyloides Data
# ---------------------------------------------------
# This dataset is reconstructed from published frequency counts in:
# Joseph, L., Gyorkos, T. W., & Coupal, L. (1995).
# Bayesian estimation of disease prevalence and the parameters of diagnostic tests
# in the absence of a gold standard.
# American Journal of Epidemiology, 141(3), 263–272.
#
# The study estimated the prevalence of Strongyloides infection among
# Cambodian refugees in Canada using two imperfect diagnostic tests.
# The data are available only in aggregated form in the original article.
# This script programmatically reconstructs an expanded binary dataset
# for research and reproducibility purposes only. Do not redistribute without citation.
set.seed(123)

patterns <- c("00", "10", "01", "11")
freq     <- c(38, 2, 87, 35)

expanded <- rep(patterns, freq)
k <- nchar(patterns[1])
data <- do.call(rbind, strsplit(expanded, ""))
data <- as.data.frame(apply(data, 2, as.numeric))
table(data)


##################################################################################
######### Reproducing the results on Joseph, Dendukuri 2001 page 159 #############
##################################################################################
prior_CI = list(
  prev = c(1,1),  # Beta(1,1) on prevalence
  tests = list(
    list(sens=c(4.44, 13.31), spec=c(71.25, 3.75)),  # Test 1
    list(sens=c(21.96, 5.49), spec=c(4.1, 1.76))   # Test 2
  )
)
fit_CI <- bayes_2LCR(data = data, model="CI",
                        iterations = 200000, burnin=50000, thin=1,
                        prior_input=prior_CI)
# saveRDS(fit_CI, "Strongyloides_CI.RDS")

prior_rand <- list(
  prev = c(1,1),   # Beta(a,b) for prevalence
  tests = list(
    # Test 1 priors (means & SDs)
    list(a1 = list(mean = -0.811, sd = 0.380),   
         a0 = list(mean =  -2.171, sd = 0.261),  ## in the paper spec was pnorm(a0/ sqrt(1+b0^2))
         b1 = list(mean =  0.668, sd = 0.5),  
         b0 = list(mean =  0.861, sd = 0.5)), 
    # Test 2 priors
    list(a1 = list(mean =  1.012, sd = 0.268),
         a0 = list(mean =  -0.692, sd = 0.560),
         b1 = list(mean =  0.668, sd = 0.5),    
         b0 = list(mean =  0.861, sd = 0.5))
  )
)

fit_rand <- bayes_2LCR(data = data, model="2LCR1", 
                       common_slopes = FALSE,
                       iterations = 200000, burnin=50000, thin=1,
                       prior_input=prior_rand)

# saveRDS(fit_rand, "Strongyloides_random.RDS")



##################################################################################
######### Running our model on the same data #####################################
##################################################################################
ranges <- list(
  list(sens=c(0.07, 0.47), spec=c(0.89, 0.99)),   # test 1
  list(sens=c(0.63, 0.92), spec=c(0.31, 0.96))   # test 2
)

pr_CI <- build_priors_from_ranges(ranges, severity="CI")
fitBLS_CI <- Bayesian_LCA_severity(
  data       = data,
  iterations = 200000,
  burnin     = 100000,
  thin       = 100,
  severity   = "CI",   
  mu_beta    = pr_CI$mu_beta,
  sd_beta    = pr_CI$sd_beta,
  m_gamma    = pr_CI$m_gamma,
  sd_gamma   = pr_CI$sd_gamma,
  rho_beta   = c(1,1)
)
saveRDS(fitBLS_CI, "Strongyloides_BLS_CI.RDS")

pr_gamma <- build_priors_from_ranges(ranges, severity="gamma", aS=3, bS=sqrt(3))
fitBLS_Gamma <- Bayesian_LCA_severity(
  data       = data,
  iterations = 200000,
  burnin     = 100000,
  thin       = 100,
  severity   = "gamma",   
  mu_beta    = pr_gamma$mu_beta,
  sd_beta    = pr_gamma$sd_beta,
  m_gamma    = pr_gamma$m_gamma,
  sd_gamma   = pr_gamma$sd_gamma,
  rho_beta   = c(1,1)
)
saveRDS(fitBLS_Gamma, "Strongyloides_BLS_Gamma3.RDS")





