

source("Gamma_LCA_severity.R")
source("CI_LCA_probit.R")
source("build_priors_from_ranges.R")

Bayesian_LCA_severity <- function(
    data, iterations, burnin, thin = 1,
    severity = c("CI", "gamma"),
    # per-test priors (required): vectors of length J (or scalars to recycle)
    mu_beta,  sd_beta,    # beta_j ~ N+(mu_beta[j], sd_beta[j]^2) (truncated at 0)
    m_gamma,  sd_gamma,   # gamma_j ~ N(m_gamma[j], sd_gamma[j]^2)
    # prevalence prior
    rho_beta = c(1, 1),   # rho ~ Beta(a,b)
    # severity hyper (only if not "CI")
    aS = 3, bS = sqrt(3)  # for gamma S|D=1
){
  severity <- match.arg(tolower(severity), c("ci","gamma"))
  
  # data checks
  data <- as.matrix(data)
  if (!all(data %in% c(0,1))) stop("data must be 0/1.")
  if (!is.numeric(iterations) || !is.numeric(burnin) || !is.numeric(thin)) stop("iterations/burnin/thin must be numeric.")
  if (iterations <= burnin) stop("iterations must be > burnin.")
  if (thin < 1) stop("thin must be >= 1.")
  if (iterations != as.integer(iterations) || burnin != as.integer(burnin) || thin != as.integer(thin))
    warning("iterations/burnin/thin will be coerced to integers by subroutines.")
  
  J <- ncol(data)
  
  # accept scalars or length-J vectors; recycle if needed
  as_lenJ <- function(x, nm){
    if (length(x) == 1) rep(x, J)
    else if (length(x) == J) x
    else stop(sprintf("%s must be length 1 or %d", nm, J))
  }
  mu_beta  <- as_lenJ(mu_beta,  "mu_beta")
  sd_beta  <- as_lenJ(sd_beta,  "sd_beta")
  m_gamma  <- as_lenJ(m_gamma,  "m_gamma")
  sd_gamma <- as_lenJ(sd_gamma, "sd_gamma")
  
  # prevalence prior
  a_rho <- rho_beta[1]; b_rho <- rho_beta[2]
  if (any(!is.finite(c(a_rho, b_rho))) || any(c(a_rho, b_rho) <= 0))
    stop("rho_beta must be positive and finite (Beta shape params).")
  
  fit <- switch(
    severity,
    "ci" = CI_LCA_probit(
      data       = data,
      iterations = iterations,
      burnin     = burnin,
      thin       = thin,
      mu_beta    = mu_beta,
      sd_beta    = sd_beta,
      m_gamma    = m_gamma,
      sd_gamma   = sd_gamma,
      a_rho      = a_rho,
      b_rho      = b_rho
    ),
    "gamma" = Gamma_LCA_severity(
      data       = data,
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
      b_rho      = b_rho
    )
  )
  
  # attach priors to fit object
  priors <- list(
    severity = severity,
    aS = aS, bS = bS,
    rho_ab = c(a = a_rho, b = b_rho),
    per_test = Map(function(mb, sdb, mg, sdg)
      list(m_beta = mb, sd_beta = sdb, m_gamma = mg, sd_gamma = sdg),
      mu_beta, sd_beta, m_gamma, sd_gamma)
  )
  fit$priors <- priors
  fit
}




