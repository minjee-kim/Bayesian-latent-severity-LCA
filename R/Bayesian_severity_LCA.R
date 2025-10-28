

source("Gamma_LCA_severity.R")
source("CI_LCA_probit.R")

Bayesian_LCA_severity <- function(
    data, iterations, burnin, thin = 1,
    severity = c("CI", "gamma"),
    # CI
    mu_beta = NULL, sd_beta = NULL,   # beta_j ~ TN+(mu_beta[j], sd_beta[j]^2)
    mu_gamma = NULL, sd_gamma = NULL,  # gamma_j ~ N(mu_gamma[j], sd_gamma[j]^2)
    # gamma
    aB = NULL, bB = NULL, # beta_j ~ Gamma(aB[j], bB[j])  (rate)
    aS = 3, bS = sqrt(3), # S|D=1 ~ Gamma(aS,bS)
    # prevalence prior
    rho_beta = c(1, 1)
){
  severity <- match.arg(tolower(severity), c("ci","gamma"))
  data <- as.matrix(data)
  if (!all(data %in% c(0,1))) stop("data must be 0/1.")
  if (!is.numeric(iterations) || !is.numeric(burnin) || !is.numeric(thin)) stop("iterations/burnin/thin must be numeric.")
  if (iterations <= burnin) stop("iterations must be > burnin.")
  if (thin < 1) stop("thin must be >= 1.")
  
  J <- ncol(data)
  as_lenJ <- function(x, nm){
    if (length(x) == 1) rep(x, J)
    else if (length(x) == J) x
    else stop(sprintf("%s must be length 1 or %d", nm, J))
  }
  
  # prevalence prior
  a_rho <- rho_beta[1]; b_rho <- rho_beta[2]
  if (any(!is.finite(c(a_rho, b_rho))) || any(c(a_rho, b_rho) <= 0))
    stop("rho_beta must be positive and finite (Beta shape params).")
  
  if (severity == "ci") {
    # CI
    if (is.null(mu_beta) || is.null(sd_beta) || is.null(mu_gamma) || is.null(sd_gamma))
      stop("For severity='CI', provide mu_beta, sd_beta, mu_gamma, sd_gamma.")
    mu_beta <- as_lenJ(mu_beta, "mu_beta")
    sd_beta <- as_lenJ(sd_beta, "sd_beta")
    mu_gamma <- as_lenJ(mu_gamma, "mu_gamma")
    sd_gamma <- as_lenJ(sd_gamma, "sd_gamma")
    if (any(!is.finite(mu_beta)) || any(!is.finite(sd_beta)) ||
        any(!is.finite(mu_gamma)) || any(!is.finite(sd_gamma)))
      stop("All CI priors must be finite.")
    if (any(sd_beta <= 0) || any(sd_gamma <= 0))
      stop("sd_beta and sd_gamma must be > 0.")
    
    fit <- CI_LCA_probit(
      data       = data,
      iterations = iterations,
      burnin     = burnin,
      thin       = thin,
      mu_beta    = mu_beta,
      sd_beta    = sd_beta,
      mu_gamma    = mu_gamma,
      sd_gamma   = sd_gamma,
      a_rho      = a_rho,
      b_rho      = b_rho
    )
    
    priors <- list(
      severity = "ci",
      rho_ab = c(a = a_rho, b = b_rho),
      per_test = Map(function(mb, sdb, mg, sdg)
        list(beta_prior = "TN+",
             m_beta = mb, sd_beta = sdb,
             mu_gamma = mg, sd_gamma = sdg),
        mu_beta, sd_beta, mu_gamma, sd_gamma)
    )
    
  } else { 
    # Gamma-severity
    if (is.null(aB) || is.null(bB) || is.null(mu_gamma) || is.null(sd_gamma))
      stop("For severity='gamma', provide aB, bB, mu_gamma, sd_gamma.")
    aB <- as_lenJ(aB, "aB")
    bB <- as_lenJ(bB, "bB")
    mu_gamma <- as_lenJ(mu_gamma, "mu_gamma")
    sd_gamma <- as_lenJ(sd_gamma, "sd_gamma")
    if (any(!is.finite(aB)) || any(!is.finite(bB)) ||
        any(!is.finite(mu_gamma)) || any(!is.finite(sd_gamma)))
      stop("All gamma-severity priors must be finite.")
    if (any(aB <= 0) || any(bB <= 0) || any(sd_gamma <= 0))
      stop("Require aB>0, bB>0, sd_gamma>0.")
    if (!is.finite(aS) || !is.finite(bS) || aS <= 1 || bS <= 0)
      stop("S prior requires aS>1 and bS>0 (Gamma rate parameterization).")
    
    fit <- Gamma_LCA_severity(
      data        = data,
      iterations  = iterations,
      burnin      = burnin,
      thin        = thin,
      aS          = aS,
      bS          = bS,
      aB          = aB,
      bB          = bB,
      mu_gamma    = mu_gamma,
      sd_gamma    = sd_gamma,
      a_rho       = a_rho,
      b_rho       = b_rho
    )
    
    priors <- list(
      severity = "gamma",
      aS = aS, bS = bS,
      rho_ab = c(a = a_rho, b = b_rho),
      per_test = Map(function(sh, rt, mg, sdg)
        list(beta_prior = "Gamma(rate)",
             aB = sh, bB = rt,
             mu_gamma = mg, sd_gamma = sdg),
        aB, bB, mu_gamma, sd_gamma)
    )
  }
  
  fit$priors <- priors
  fit
}
