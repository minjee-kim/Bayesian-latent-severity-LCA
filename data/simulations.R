


set.seed(123)
source("Bayesian_severity_LCA.R")
source("bayes_2LCR.R")

##################################################################
########## Conditional Independence ##############################
##################################################################

N <- 4000
J <- 4
rho_true <- 0.35
D <- rbinom(N, 1, rho_true)
gamma_true <- c(qnorm(0.95), qnorm(0.80), qnorm(0.60), qnorm(0.99))
beta_true  <- c(1.4, 2.4, 3.0, 0.7)
S <- ifelse(D==1, 1, 0)
Mu <- outer(D*S, beta_true) - matrix(gamma_true, N, J, byrow=TRUE)
Tij <- matrix( rbinom(N*J, 1, pnorm(Mu)), nrow=N, ncol=J)

# Per-test prior ranges for Se/Sp (interpreted as 95% intervals)

prior_2LCR <- list(
  prev = c(1,1),
  ranges = list(
    list(sens=c(0.20, 0.99), spec=c(0.30, 0.999)),
    list(sens=c(0.20, 0.99), spec=c(0.30, 0.999)),
    list(sens=c(0.20, 0.99), spec=c(0.30, 0.999)),
    list(sens=c(0.20, 0.99), spec=c(0.30, 0.999))
  ),
  range_ci = 0.95
)
CI_bayes2LCR_CI <- bayes_2LCR(data = Tij, model="CI",
                     iterations = 200000 + 200000, burnin=200000, thin=1,
                     prior_input = prior_2LCR )
# saveRDS(CI_bayes2LCR_CI, "CIsim_bayes2LCR_CI.RDS")

CI_bayes2LCR_RE <- bayes_2LCR(
  data = Tij, model = "random",
  iterations = 200000 + 200000,  # total iters
  burnin = 200000,
  thin = 20,
  prior_input = prior_2LCR,
  common_slopes = FALSE   # set TRUE if you intend shared slopes, else leave FALSE
)
# saveRDS(CI_bayes2LCR_RE, "CIsim_bayes2LCR_RE.RDS")

## Our models
ranges <- list(
  list(sens=c(0.20, 0.99), spec=c(0.30, 0.999)), # test 1
  list(sens=c(0.20, 0.99), spec=c(0.30, 0.999)), # test 2
  list(sens=c(0.20, 0.99), spec=c(0.30, 0.999)), # test 3
  list(sens=c(0.20, 0.99), spec=c(0.30, 0.999))  # test 4
)
CI_fitBLS_CI <- Bayesian_LCA_severity(
  data       = Tij,
  iterations = 200000,
  burnin     = 200000,
  thin       = 2,
  severity_prior = list(type="ci"),
  ranges     = ranges,
  rho_prior  = list(beta = c(1,1))
)
saveRDS(CI_fitBLS_CI, "CIsim_fitBLS_CI.RDS")

CI_fitBLS_Gamma <- Bayesian_LCA_severity(
  data       = Tij,
  iterations = 200000,
  burnin     = 200000,
  thin       = 2,
  severity_prior = list(type="gamma", aS=3, bS=sqrt(3)),
  ranges     = ranges,
  rho_prior  = list(beta = c(1,1))
)
saveRDS(CI_fitBLS_Gamma, "CIsim_fitBLS_Gamma.RDS")

CI_fitBLS_NM <- Bayesian_LCA_severity(
  data       = Tij,
  iterations = 200000,
  burnin     = 200000,
  thin       = 2,
  severity_prior = list(type="normal moment", mu0=0, tau=1.48495),
  ranges     = ranges,
  rho_prior  = list(beta = c(1,1))
)
saveRDS(CI_fitBLS_NM, "CIsim_fitBLS_NM.RDS")

##################################################################
########## Dependence ############################################
##################################################################

N <- 4000
J <- 4
rho_true <- 0.35
D <- rbinom(N, 1, rho_true)
gamma_true <- c(qnorm(0.95), qnorm(0.80), qnorm(0.60), qnorm(0.99))
beta_true  <- c(1.4, 2.4, 3.0, 0.7)
S <- ifelse(D==1, rgamma(N, shape=4.5, rate=sqrt(4.5)), 0)
Mu <- outer(D*S, beta_true) - matrix(gamma_true, N, J, byrow=TRUE)
Tij <- matrix( rbinom(N*J, 1, pnorm(Mu)), nrow=N, ncol=J)



RE_bayes2LCR_CI <- bayes_2LCR(data = Tij, model="CI",
                              iterations = 200000 + 200000, burnin=200000, thin=1,
                              prior_input = prior_2LCR)

RE_bayes2LCR_RE <- bayes_2LCR(data = Tij, model="random",
                              iterations = 200000 + 200000, burnin=200000, thin=1,
                              prior_input = prior_2LCR)

# Per-test prior ranges for Se/Sp (interpreted as 95% intervals)
ranges <- list(
  list(sens=c(0.20, 0.99), spec=c(0.30, 0.999)), # test 1
  list(sens=c(0.20, 0.99), spec=c(0.30, 0.999)), # test 2
  list(sens=c(0.20, 0.99), spec=c(0.30, 0.999)), # test 3
  list(sens=c(0.20, 0.99), spec=c(0.30, 0.999))  # test 4
)

# Run Models
fitBLS_CI <- Bayesian_LCA_severity(
  data       = Tij,
  iterations = 200000,
  burnin     = 200000,
  thin       = 2,
  severity_prior = list(type="ci"),
  ranges     = ranges,
  rho_prior  = list(beta = c(1,1))
)

fitBLS_Gamma <- Bayesian_LCA_severity(
  data       = Tij,
  iterations = 200000,
  burnin     = 200000,
  thin       = 2,
  severity_prior = list(type="gamma", aS=3, bS=sqrt(3)),
  ranges     = ranges,
  rho_prior  = list(beta = c(1,1))
)

fitBLS_NM <- Bayesian_LCA_severity(
  data       = Tij,
  iterations = 200000,
  burnin     = 200000,
  thin       = 2,
  severity_prior = list(type="normal moment", mu0=0, tau=1.48495),
  ranges     = ranges,
  rho_prior  = list(beta = c(1,1))
)


