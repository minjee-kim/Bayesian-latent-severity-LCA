


build_priors_from_ranges <- function(
    ranges,
    severity = list(type = "gamma")
){
  aS = 3
  bS = sqrt(3)
  ci_level = 0.95
  clamp_sd = 1.5
  Sm  <- S_moments(severity_prior)
  muS <- Sm$mean
  vS  <- Sm$var
  z   <- qnorm((1 + ci_level) / 2)
  
  # Solve β from target p = Φ( (β μS − mγ) / sqrt(1 + β^2 vS) )
  solve_beta_from_Se <- function(p, m_gamma, muS, vS){
    p <- min(max(p, 1e-9), 1 - 1e-9)
    t <- qnorm(p)
    A <- muS^2 - (t^2) * vS
    B <- -2 * muS * m_gamma
    C <- m_gamma^2 - t^2
    if (abs(A) < 1e-10) {
      beta <- -C / B
      return(max(beta, 1e-9))
    } else {
      disc <- B*B - 4*A*C
      disc <- max(disc, 0) 
      b1 <- (-B + sqrt(disc)) / (2*A)
      b2 <- (-B - sqrt(disc)) / (2*A)
      cand <- c(b1, b2)
      cand <- cand[is.finite(cand) & cand > 0]
      if (length(cand) == 0) {
        return(1e-6)
      } else {
        best <- cand[1]
        best_err <- abs(p - pnorm((best*muS - m_gamma)/sqrt(1 + best*best*vS)))
        for (bb in cand[-1]) {
          err <- abs(p - pnorm((bb*muS - m_gamma)/sqrt(1 + bb*bb*vS)))
          if (err < best_err) { best <- bb; best_err <- err }
        }
        return(best)
      }
    }
  }
  
  map_spec_to_gamma <- function(Lsp, Usp, z, clamp_sd){
    Lsp <- min(max(Lsp, 1e-6), 1 - 1e-6)
    Usp <- min(max(Usp, 1e-6), 1 - 1e-6)
    zL  <- qnorm(Lsp); zU <- qnorm(Usp)
    m_gamma  <- 0.5 * (zL + zU)
    sd_gamma <- (zU - zL) / (2 * z)
    sd_gamma <- min(max(sd_gamma, 0.05), clamp_sd)
    c(m_gamma = m_gamma, sd_gamma = sd_gamma)
  }
  
  out <- lapply(ranges, function(rg){
    Lsp <- min(rg$spec); Usp <- max(rg$spec)
    gg  <- map_spec_to_gamma(Lsp, Usp, z, clamp_sd)
    m_gamma <- gg["m_gamma"]; sd_gamma <- gg["sd_gamma"]
    
    Lse <- min(rg$sens); Use <- max(rg$sens)
    Lse <- min(max(Lse, 1e-6), 1 - 1e-6)
    Use <- min(max(Use, 1e-6), 1 - 1e-6)
    mse <- 0.5 * (Lse + Use)
    
    m_beta <- solve_beta_from_Se(mse, m_gamma, muS, vS)
    
    bL <- solve_beta_from_Se(Lse, m_gamma, muS, vS)
    bU <- solve_beta_from_Se(Use, m_gamma, muS, vS)
    sd_beta <- (bU - bL) / (2 * z)
    sd_beta <- max(sd_beta, 0.05)
    
    if (tolower(severity_prior$type) == "ci") {
      list(
        m_gamma = as.numeric(m_gamma),
        sd_gamma = as.numeric(sd_gamma),
        m_beta  = as.numeric(m_beta),
        sd_beta = as.numeric(sd_beta),
        beta_prior = "TN+"
      )
    } else {
      aB <- (m_beta / sd_beta)^2
      bB <-  m_beta / (sd_beta^2)
      aB <- max(aB, 1e-6)
      bB <- max(bB, 1e-6)
      list(
        m_gamma = as.numeric(m_gamma),
        sd_gamma = as.numeric(sd_gamma),
        aB = as.numeric(aB),
        bB = as.numeric(bB),
        beta_prior = "Gamma(rate)"
      )
    }
  })
  
  out
}



