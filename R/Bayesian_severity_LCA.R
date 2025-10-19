

source("Gamma_LCA_severity.R")
source("NM_LCA_severity.R")
source("CI_LCA_probit.R")

Bayesian_LCA_severity <- function(
    data, iterations, burnin, thin = 1,
    severity = c("CI", "gamma", "nm+"),
    # per-test priors (required): vectors of length J
    mu_beta,  sd_beta,    # β_j ~ N+(mu_beta[j], sd_beta[j]^2) (truncated at 0)
    m_gamma,  sd_gamma,   # γ_j ~ N(m_gamma[j], sd_gamma[j]^2)
    # prevalence prior
    rho_beta = c(1, 1),   # ρ ~ Beta(a,b)
    # severity hyper (only if not "CI")
    aS = NULL, bS = NULL, # for gamma S|D=1
    mu0 = NULL, tau = NULL # for NM+ S|D=1
){
  library(truncnorm)
  severity <- match.arg(tolower(severity), c("CI","gamma","nm+"))
  if (!all(data %in% c(0,1))) stop("data must be 0/1.")
  
  # Input sanity
  J <- ncol(data)
  req_len <- function(x, nm) if (length(x) != J) stop(sprintf("%s must have length %d", nm, J))
  req_len(mu_beta, "mu_beta"); req_len(sd_beta, "sd_beta")
  req_len(m_gamma, "m_gamma"); req_len(sd_gamma, "sd_gamma")
  
  a_rho <- rho_beta[1]; b_rho <- rho_beta[2]
  if (any(!is.finite(c(a_rho, b_rho))) || any(c(a_rho, b_rho) <= 0))
    stop("rho_beta must be positive.")
  
  # build per-test priors from simple Se/Sp ranges 
  build_priors_from_ranges <- function(ranges, ci_level = 0.95, clamp_sd = 1.5) {
    # ranges: list of length J; each [[j]] has $sens=c(L,U), $spec=c(L,U)
    J <- length(ranges)
    z <- qnorm((1 + ci_level)/2)
    
    clamp <- function(x, lo, hi) pmin(pmax(x, lo), hi)
    
    m_gamma <- sd_gamma <- mu_beta <- sd_beta <- numeric(J)
    
    for (j in seq_len(J)) {
      rg <- ranges[[j]]
      
      # gamma from specificity: Sp_j = Phi(gamma_j)
      Lsp <- clamp(min(rg$spec), 1e-4, 1 - 1e-4)
      Usp <- clamp(max(rg$spec), 1e-4, 1 - 1e-4)
      m_g  <- (qnorm(Lsp) + qnorm(Usp)) / 2
      sd_g <- (qnorm(Usp) - qnorm(Lsp)) / (2 * z)
      sd_g <- clamp(sd_g, 0.02, clamp_sd)
      
      # beta (approx) from sensitivity probit: Se ≈ Phi(beta - gamma)
      Lse <- clamp(min(rg$sens), 1e-4, 1 - 1e-4)
      Use <- clamp(max(rg$sens), 1e-4, 1 - 1e-4)
      m_eta  <- (qnorm(Lse) + qnorm(Use)) / 2
      sd_eta <- (qnorm(Use) - qnorm(Lse)) / (2 * z)
      sd_eta <- clamp(sd_eta, 0.02, clamp_sd)
      
      m_b  <- m_eta + m_g
      sd_b <- sqrt(sd_eta^2 + sd_g^2) # rough delta method
      sd_b <- clamp(sd_b, 0.05, clamp_sd)
      
      m_gamma[j] <- m_g
      sd_gamma[j] <- sd_g
      mu_beta[j] <- max(m_b, 0)   # truncate mean at 0
      sd_beta[j] <- sd_b
    }
    
    list(mu_beta = mu_beta, sd_beta = sd_beta,
         m_gamma = m_gamma, sd_gamma = sd_gamma)
  }
  
  
  if (severity == "CI") {
    return(CI_LCA_probit(
      data       = data,
      iterations = iterations,
      burnin     = burnin,
      thin       = thin,
      m_beta     = mu_beta,
      sd_beta    = sd_beta,
      m_gamma    = m_gamma,
      sd_gamma   = sd_gamma,
      a_rho      = a_rho,
      b_rho      = b_rho
    ))
  }
  
  if (severity == "gamma") {
    return(Gamma_LCA_severity(
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
    ))
  }
  
  # NM+ severity
  if (severity == "nm+") {
    return(NM_LCA_severity(
    data       = data,
    iterations = iterations,
    burnin     = burnin,
    thin       = thin,
    mu0        = mu0,
    tau        = tau,
    mu_beta    = mu_beta,
    sd_beta    = sd_beta,
    m_gamma    = m_gamma,
    sd_gamma   = sd_gamma,
    a_rho      = a_rho,
    b_rho      = b_rho
    ))
    }
}
