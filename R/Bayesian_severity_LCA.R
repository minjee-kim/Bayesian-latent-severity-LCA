

source("Gamma_LCA_severity.R")
source("CI_LCA_probit.R")
source("build_priors_from_ranges.R")

Bayesian_LCA_severity <- function(
    data, iterations, burnin, thin = 1,
    severity = c("CI", "gamma"),
    # per-test priors (required): vectors of length J
    mu_beta,  sd_beta,    # beta_j ~ N+(mu_beta[j], sd_beta[j]^2) (truncated at 0)
    m_gamma,  sd_gamma,   # gamma_j ~ N(m_gamma[j], sd_gamma[j]^2)
    # prevalence prior
    rho_beta = c(1, 1),   # rho ~ Beta(a,b)
    # severity hyper (only if not "CI")
    aS = 3, bS = sqrt(3) # for gamma S|D=1
){
  library(truncnorm)
  severity <- match.arg(tolower(severity), c("ci","gamma"))
  data = as.matrix(data)
  if (!all(data %in% c(0,1))) stop("data must be 0/1.")
  
  # Input sanity
  J <- ncol(data)
  req_len <- function(x, nm) if (length(x) != J) stop(sprintf("%s must have length %d", nm, J))
  req_len(mu_beta, "mu_beta"); req_len(sd_beta, "sd_beta")
  req_len(m_gamma, "m_gamma"); req_len(sd_gamma, "sd_gamma")
  
  a_rho <- rho_beta[1]; b_rho <- rho_beta[2]
  if (any(!is.finite(c(a_rho, b_rho))) || any(c(a_rho, b_rho) <= 0))
    stop("rho_beta must be positive.")
  
  fit <- switch(severity,
                "ci" = CI_LCA_probit(data       = data,
                                     iterations = iterations,
                                     burnin     = burnin,
                                     thin       = thin,
                                     mu_beta     = mu_beta,
                                     sd_beta    = sd_beta,
                                     m_gamma    = m_gamma,
                                     sd_gamma   = sd_gamma,
                                     a_rho      = a_rho,
                                     b_rho      = b_rho),
                "gamma" = Gamma_LCA_severity(data       = data,
                                             iterations = iterations,
                                             burnin     = burnin,
                                             thin       = thin,
                                             aS         = aS,
                                             bS         = bS,
                                             mu_beta    = mu_beta,
                                             sd_beta    = sd_beta,
                                             m_gamma    = m_gamma,
                                             sd_gamma   = sd_gamma,
                                             a_rho      = a_rho,
                                             b_rho      = b_rho)
  )
  
  priors <- list(
    severity = tolower(severity),
    aS = aS, bS = bS,
    rho_ab = c(a = rho_beta[1], b = rho_beta[2]),
    per_test = Map(function(mb, sdb, mg, sdg)
      list(m_beta = mb, sd_beta = sdb, m_gamma = mg, sd_gamma = sdg),
      mu_beta, sd_beta, m_gamma, sd_gamma)
  )
  attach(priors)
  fit$priors <- priors
  fit
}
