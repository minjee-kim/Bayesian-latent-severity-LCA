

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