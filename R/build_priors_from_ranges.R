

# Build per-test priors (β_j ~ N^+(mu_beta[j], sd_beta[j]^2), γ_j ~ N(m_gamma[j], sd_gamma[j]^2))
# from user-specified Se/Sp ranges, accounting for severity shape.
#
# ranges: list of length J; each [[j]] is a list with $sens=c(L,U), $spec=c(L,U)
# severity: "CI" (constant S=1), "gamma" (S~Gamma(aS, rate=bS)), or "nm+" (S~NM+(mu0,tau))
# Returns a list with vectors mu_beta, sd_beta, m_gamma, sd_gamma (length J).
build_priors_from_ranges <- function(
    ranges,
    severity = c("CI","gamma","nm+"),
    aS = 3, bS = sqrt(3),
    mu0 = 0, tau = 1.48495
){
  severity <- match.arg(tolower(severity), c("ci","gamma","nm+"))
  J <- length(ranges); if (J == 0L) stop("ranges must be a non-empty list.")
  ci_level        <- 0.95
  clamp_sd_beta   <- 1.5
  clamp_sd_gamma  <- 1.5
  nm_mc           <- 40000L
  z               <- qnorm((1 + ci_level) / 2)
  
  clamp01  <- function(p) pmin(pmax(p, 1e-4), 1 - 1e-4)
  clamp    <- function(x, lo, hi) pmin(pmax(x, lo), hi)
  probit_mid_sd <- function(L, U, z) {
    L <- clamp01(L); U <- clamp01(U)
    m  <- (qnorm(L) + qnorm(U)) / 2
    sd <- (qnorm(U) - qnorm(L)) / (2 * z)
    list(m = m, sd = sd)
  }
  
  # NM+ sampler 
  sample_NM <- function(n, mu0, tau, grid_size = 1e4) {
    x_vals <- seq(0, mu0 + 6 * tau, length.out = grid_size)
    dens   <- 2 * (x_vals^2) / (sqrt(2*pi) * tau^3) * exp(-((x_vals - mu0)^2) / (2*tau^2))
    dx     <- x_vals[2] - x_vals[1]
    dens   <- dens / sum(dens * dx)
    cdf    <- cumsum(dens * dx)
    u      <- runif(n)
    approx(cdf, x_vals, xout = u, rule = 2)$y
  }
  
  # --- moments of S|D=1 for the chosen severity prior ---
  sev <- switch(
    severity,
    "ci"    = list(mu = 1,           var = 0),
    "gamma" = { if (is.null(aS) || is.null(bS)) stop("Need aS,bS for severity='gamma'.")
      list(mu = aS / bS,   var = aS / (bS^2)) },
    "nm+"   = { if (is.null(mu0) || is.null(tau)) stop("Need mu0,tau for severity='nm+'.")
      S <- sample_NM(nm_mc, mu0, tau); list(mu = mean(S), var = var(S)) }
  )
  muS <- sev$mu
  sig2 <- sev$var
  
  # Solve mean(beta) from: m_eta = (beta*muS - m_g) / sqrt(1 + beta^2 * sig2)
  solve_beta_mean <- function(m_eta, m_g, muS, sig2) {
    if (sig2 <= 0) return(max(0, m_eta + m_g))  # CI case
    A <- (muS^2 - (m_eta^2) * sig2)
    B <- -2 * muS * m_g
    C <- (m_g^2 - m_eta^2)
    if (abs(A) < 1e-12) {
      f <- function(b) ((b*muS - m_g) / sqrt(1 + b^2*sig2)) - m_eta
      root <- try(uniroot(f, c(0, 10))$root, silent = TRUE)
      return(if (inherits(root, "try-error")) max(0, m_eta + m_g) else max(0, root))
    }
    disc <- B*B - 4*A*C
    if (disc < 0) return(max(0, m_eta + m_g))
    r1 <- (-B + sqrt(disc)) / (2*A)
    r2 <- (-B - sqrt(disc)) / (2*A)
    cand <- c(r1, r2); cand <- cand[is.finite(cand) & cand >= 0]
    if (length(cand)) min(cand) else max(0, m_eta + m_g)
  }
  
  m_gamma <- sd_gamma <- mu_beta <- sd_beta <- numeric(J)
  
  
  for (j in seq_len(J)) {
    rg <- ranges[[j]]
    if (is.null(rg$spec) || is.null(rg$sens))
      stop(sprintf("ranges[[%d]] must have $sens and $spec.", j))
    
    # gamma prior from specificity: Sp = Phi(gamma)
    g <- probit_mid_sd(min(rg$spec), max(rg$spec), z)
    m_g  <- g$m
    sd_g <- clamp(g$sd, 0.02, clamp_sd_gamma)
    
    # target probit for sensitivity
    e <- probit_mid_sd(min(rg$sens), max(rg$sens), z)
    m_eta  <- e$m
    sd_eta <- clamp(e$sd, 0.02, clamp_sd_beta)
    
    # beta mean & sd (delta method)
    m_b <- solve_beta_mean(m_eta, m_g, muS, sig2)
    D   <- sqrt(1 + m_b^2 * sig2)
    dAdb <- ( muS*(1 + m_b^2*sig2) - (m_b*muS - m_g)*m_b*sig2 ) / (D^3)
    varA <- sd_eta^2 + (sd_g^2) / (D^2)
    sd_b <- clamp(sqrt(varA / max(dAdb^2, 1e-12)), 0.05, clamp_sd_beta)
    
    m_gamma[j] <- m_g
    sd_gamma[j] <- sd_g
    mu_beta[j] <- m_b
    sd_beta[j] <- sd_b
  }
  
  list(
    mu_beta  = mu_beta,
    sd_beta  = sd_beta,
    m_gamma  = m_gamma,
    sd_gamma = sd_gamma
  )
}
