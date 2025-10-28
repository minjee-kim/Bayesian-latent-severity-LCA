
build_priors_from_ranges <- function(
    ranges,
    severity_type = c("CI","gamma"),
    aS = 3, bS = sqrt(3),     
    ci_level = 0.95,
    clamp_sd = 1.5
){
  severity_type <- match.arg(tolower(severity_type), c("ci","gamma"))
  J <- length(ranges); stopifnot(J > 0)
  
  # moments of S
  if (severity_type == "ci") {
    muS <- 1
    vS  <- 0
  } else {
    if (!is.finite(aS) || !is.finite(bS) || aS <= 0 || bS <= 0)
      stop("For gamma severity: aS>0 and bS>0 (rate parameterization).")
    muS <- aS / bS
    vS  <- aS / (bS^2)
  }
  
  z <- stats::qnorm((1 + ci_level)/2)
  clamp01 <- function(p) pmin(pmax(p, 1e-6), 1 - 1e-6)
  
  # Sp range - Normal prior for gamma
  map_spec_to_gamma <- function(Lsp, Usp){
    Lsp <- clamp01(min(Lsp, Usp)); Usp <- clamp01(max(Lsp, Usp))
    zL  <- stats::qnorm(Lsp); zU <- stats::qnorm(Usp)
    m_gamma  <- 0.5 * (zL + zU)
    sd_gamma <- (zU - zL) / (2*z)
    sd_gamma <- min(max(sd_gamma, 0.05), clamp_sd)
    c(m_gamma = m_gamma, sd_gamma = sd_gamma)
  }
  
  # Solve for beta 
  solve_beta_from_Se <- function(p, m_gamma, muS, vS){
    p <- clamp01(p)
    t <- stats::qnorm(p)
    A <- muS^2 - (t^2) * vS
    B <- -2 * muS * m_gamma
    C <- m_gamma^2 - t^2
    if (abs(A) < 1e-12) {
      beta <- -C / B
      return(max(beta, 1e-9))
    } else {
      disc <- B*B - 4*A*C
      disc <- max(disc, 0)
      b1 <- (-B + sqrt(disc)) / (2*A)
      b2 <- (-B - sqrt(disc)) / (2*A)
      cand <- c(b1, b2)
      cand <- cand[is.finite(cand) & cand > 0]
      if (!length(cand)) return(1e-6)
      errs <- sapply(cand, function(bb)
        abs(p - stats::pnorm((bb*muS - m_gamma)/sqrt(1 + bb*bb*vS))))
      cand[which.min(errs)]
    }
  }
  
  # accumulate per-test
  mu_gamma <- sd_gamma <- numeric(J)
  if (severity_type == "ci") {
    mu_beta <- sd_beta <- numeric(J)
  } else {
    aB <- bB <- numeric(J)
  }
  
  for (j in seq_len(J)){
    rj <- ranges[[j]]
    gg <- map_spec_to_gamma(rj$spec[1], rj$spec[2])
    m_g  <- gg["m_gamma"]; s_g <- gg["sd_gamma"]
    mu_gamma[j] <- m_g; sd_gamma[j] <- s_g
    
    Lse <- clamp01(min(rj$sens)); Use <- clamp01(max(rj$sens))
    mse <- 0.5*(Lse + Use)
    
    m_beta <- solve_beta_from_Se(mse, m_g, muS, vS)
    bL     <- solve_beta_from_Se(Lse, m_g, muS, vS)
    bU     <- solve_beta_from_Se(Use, m_g, muS, vS)
    sd_b   <- max((bU - bL) / (2*z), 0.05)
    
    if (severity_type == "ci") {
      mu_beta[j] <- m_beta
      sd_beta[j] <- sd_b
    } else {
      # match mean/sd to Gamma(rate): mean = a/b, var = a/b^2
      aB[j] <- (m_beta / sd_b)^2
      bB[j] <-  m_beta / (sd_b^2)
      aB[j] <- max(aB[j], 1e-6); bB[j] <- max(bB[j], 1e-6)
    }
  }
  
  if (severity_type == "ci") {
    return(list(
      mu_beta  = mu_beta,
      sd_beta  = sd_beta,
      mu_gamma = mu_gamma,
      sd_gamma = sd_gamma
    ))
  } else {
    return(list(
      aB       = aB,
      bB       = bB,
      mu_gamma = mu_gamma,
      sd_gamma = sd_gamma,
      aS = aS, bS = bS
    ))
  }
}
