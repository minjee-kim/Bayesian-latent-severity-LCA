


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

################## Bayesian RE Models########################################
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
                     iterations = 200000, burnin=20000, thin=100,
                     prior_input = prior_2LCR )
saveRDS(CI_bayes2LCR_CI, "CIsim_bayes2LCR_CI.RDS")

CI_bayes2LCR_RE <- bayes_2LCR(
  data = Tij, model = "random",
  iterations = 200000, 
  burnin = 20000,
  thin = 100,
  prior_input = prior_2LCR,
  common_slopes = FALSE   # set TRUE if you intend shared slopes, else leave FALSE
)
saveRDS(CI_bayes2LCR_RE, "CIsim_bayes2LCR_RE.RDS")

################## Our Models########################################
ranges <- list(
  list(sens=c(0.20, 0.99), spec=c(0.30, 0.999)), # test 1
  list(sens=c(0.20, 0.99), spec=c(0.30, 0.999)), # test 2
  list(sens=c(0.20, 0.99), spec=c(0.30, 0.999)), # test 3
  list(sens=c(0.20, 0.99), spec=c(0.30, 0.999))  # test 4
)

pr_CI <- build_priors_from_ranges(ranges, severity="CI")
CI_fitBLS_CI <- Bayesian_LCA_severity(
  data       = Tij,
  iterations = 200000,
  burnin     = 20000,
  thin       = 100,
  severity   = "CI",   
  mu_beta    = pr_CI$mu_beta,
  sd_beta    = pr_CI$sd_beta,
  m_gamma    = pr_CI$m_gamma,
  sd_gamma   = pr_CI$sd_gamma,
  rho_beta   = c(1,1)
)
saveRDS(CI_fitBLS_CI, "CIsim_fitBLS_CI.RDS")

pr_gamma <- build_priors_from_ranges(ranges, severity="gamma", aS=3, bS=sqrt(3))
CI_fitBLS_Gamma <- Bayesian_LCA_severity(
  data       = Tij,
  iterations = 200000,
  burnin     = 20000,
  thin       = 100,
  severity   = "gamma",   
  mu_beta    = pr_gamma$mu_beta,
  sd_beta    = pr_gamma$sd_beta,
  m_gamma    = pr_gamma$m_gamma,
  sd_gamma   = pr_gamma$sd_gamma,
  rho_beta   = c(1,1)
)
saveRDS(CI_fitBLS_Gamma, "CIsim_fitBLS_Gamma.RDS")

pr_nm <- build_priors_from_ranges(ranges, severity="nm+", mu0 = 0, tau = 1.48495)
CI_fitBLS_NM <- Bayesian_LCA_severity(
  data       = Tij,
  iterations = 200000,
  burnin     = 20000,
  thin       = 100,
  severity   = "nm+",   
  mu_beta    = pr_nm$mu_beta,
  sd_beta    = pr_nm$sd_beta,
  m_gamma    = pr_nm$m_gamma,
  sd_gamma   = pr_nm$sd_gamma,
  rho_beta   = c(1,1)
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

################## Bayesian RE Models#######################################
RE_bayes2LCR_CI <- bayes_2LCR(data = Tij, model="CI",
                              iterations = 200000, burnin=20000, thin=100,
                              prior_input = prior_2LCR )

RE_bayes2LCR_RE <- bayes_2LCR(data = Tij, model="random",
                              iterations = 200000, burnin=20000, thin=100,
                              prior_input = prior_2LCR)

################## Our Models########################################

fitBLS_CI <- Bayesian_LCA_severity(
  data       = Tij,
  iterations = 200000,
  burnin     = 20000,
  thin       = 100,
  severity   = "CI",   
  mu_beta    = pr_CI$mu_beta,
  sd_beta    = pr_CI$sd_beta,
  m_gamma    = pr_CI$m_gamma,
  sd_gamma   = pr_CI$sd_gamma,
  rho_beta   = c(1,1)
)
saveRDS(fitBLS_CI, "REsim_fitBLS_CI.RDS")

fitBLS_Gamma <- Bayesian_LCA_severity(
  data       = Tij,
  iterations = 200000,
  burnin     = 20000,
  thin       = 100,
  severity   = "gamma",   
  mu_beta    = pr_gamma$mu_beta,
  sd_beta    = pr_gamma$sd_beta,
  m_gamma    = pr_gamma$m_gamma,
  sd_gamma   = pr_gamma$sd_gamma,
  rho_beta   = c(1,1)
)
saveRDS(fitBLS_Gamma, "REsim_fitBLS_Gamma.RDS")

fitBLS_NM <- Bayesian_LCA_severity(
  data       = Tij,
  iterations = 200000,
  burnin     = 20000,
  thin       = 100,
  severity   = "nm+",   
  mu_beta    = pr_nm$mu_beta,
  sd_beta    = pr_nm$sd_beta,
  m_gamma    = pr_nm$m_gamma,
  sd_gamma   = pr_nm$sd_gamma,
  rho_beta   = c(1,1)
)
saveRDS(fitBLS_NM, "REsim_fitBLS_NM.RDS")

