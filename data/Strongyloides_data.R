
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
prior_input = list(
  prev = c(1,1),  # Beta(1,1) on prevalence
  tests = list(
    list(sens=c(4.44, 13.31), spec=c(71.25, 3.75)),  # Test 1
    list(sens=c(21.96, 5.49), spec=c(4.1, 1.76))   # Test 2
  )
)
fit_fixed <- bayes_2LCR(data = data, model="fixed",
                        iterations=20000, burnin=3000, thin=2,
                        prior_input=prior_input)
# saveRDS(fit_fixed, "Strongyloides_fixed.RDS")

prior_input <- list(
  prev = c(1,1),   # Beta(a,b) for prevalence
  tests = list(
    # Test 1 priors (means & SDs)
    list(a1 = list(mean = -0.811, sd = 0.380),   # class 0
         a0 = list(mean =  2.171, sd = 0.261),   # class 1
         b1 = list(mean =  0.668, sd = 0.5),   # class 0 slope
         b0 = list(mean =  0.861, sd = 0.5)),   # class 1 slope
    # Test 2 priors
    list(a1 = list(mean =  1.012, sd = 0.268),
         a0 = list(mean =  0.692, sd = 0.560),
         b1 = list(mean =  0.668, sd = 0.5),           # class 0 slope
         b0 = list(mean =  0.861, sd = 0.5))
  ),
  # for 2LCR1 you can set this to TRUE or specify in the model argument 
  common_slopes = FALSE  # TRUE: b0,b1 shared across tests
)

fit_rand <- bayes_2LCR(data = data, model="2LCR1", 
                       iterations = 200000, burnin=30000, thin=2,
                       common_slopes = TRUE,
                       prior_input=prior_input)
# saveRDS(fit_rand, "Strongyloides_random.RDS")

