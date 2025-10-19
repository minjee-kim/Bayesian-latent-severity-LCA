

# Build per-test priors from Se/Sp ranges, with severity-aware Se calibration.
# Returns a list of equal-length vectors: mu_beta, sd_beta, m_gamma, sd_gamma.
build_priors_from_ranges <- function(
    ranges,
    severity = c("CI","gamma","nm+"),
    aS = 3, bS = sqrt(3),       # for severity="gamma"
    mu0 = 0, tau = 1.48495,     # for severity="nm+"
    ci_level = 0.95
){
  nsim_cal = 40000 # MC size for calibration (gamma/NM+)
  clamp_sd_beta = 1.5
  clamp_sd_gamma = 1.5
  
  severity <- match.arg(tolower(severity), c("ci","gamma","nm+"))
  J <- length(ranges); if (J == 0L) stop("ranges must be a non-empty list.")
  z <- qnorm((1 + ci_level) / 2)
  
  clamp01 <- function(p) pmin(pmax(p, 1e-4), 1 - 1e-4)
  clamp   <- function(x, lo, hi) pmin(pmax(x, lo), hi)
  
  probit_mid_sd <- function(L, U, z) {
    L <- clamp01(L); U <- clamp01(U)
    m  <- (qnorm(L) + qnorm(U)) / 2
    sd <- (qnorm(U) - qnorm(L)) / (2 * z)
    list(m = m, sd = sd)
  }
  
  # S samplers for severity models
  .sample_S <- function(n, severity, aS, bS, mu0, tau){
    if (severity == "ci")    return(rep(1, n))
    if (severity == "gamma") return(rgamma(n, shape = aS, rate = bS))
    # NM+ inverse-CDF sampler
    x  <- seq(0, mu0 + 6*tau, length.out = 10000L)
    dx <- x[2] - x[1]
    dens <- 2*(x^2)/(sqrt(2*pi)*tau^3) * exp(- (x-mu0)^2/(2*tau^2))
    dens <- dens / sum(dens*dx)
    cdf  <- cumsum(dens*dx)
    approx(cdf, x, xout = runif(n), rule = 2)$y
  }
  
  # Calibrate (m_beta, sd_beta) so Se quantiles match [Lse, Use]
  calibrate_beta_to_Se <- function(Lse, Use, m_gamma, sd_gamma,
                                   severity, aS, bS, mu0, tau,
                                   nsim = nsim_cal){
    stopifnot(0 < Lse, Lse < Use, Use < 1)
    S <- .sample_S(nsim, severity, aS, bS, mu0, tau)
    target <- c(L = Lse, U = Use)
    
    obj <- function(par){                # par = log(m), log(sd)
      m  <- exp(par[1]); sd <- exp(par[2])
      if (!is.finite(m) || !is.finite(sd) || m < 1e-4 || m > 15 || sd < 1e-3 || sd > 5) return(1e6)
      beta  <- truncnorm::rtruncnorm(nsim, a = 0, b = Inf, mean = m, sd = sd)
      gamma <- rnorm(nsim, m_gamma, sd_gamma)
      se    <- pnorm(beta * S - gamma)
      q     <- quantile(se, c(0.025, 0.975), names = FALSE)
      sum((q - target)^2)
    }
    o <- optim(c(log(1), log(0.6)), obj, method = "Nelder-Mead",
               control = list(maxit = 300, reltol = 1e-8))
    list(m_beta = exp(o$par[1]), sd_beta = exp(o$par[2]), ok = (o$convergence == 0))
  }
  
  m_gamma <- sd_gamma <- mu_beta <- sd_beta <- numeric(J)
  
  for (j in seq_len(J)) {
    rg <- ranges[[j]]
    if (is.null(rg$spec) || is.null(rg$sens))
      stop(sprintf("ranges[[%d]] must have $sens and $spec.", j))
    
    # From specificity range -> gamma prior (independent of severity)
    g     <- probit_mid_sd(min(rg$spec), max(rg$spec), z)
    m_g   <- g$m
    sd_g  <- clamp(g$sd, 0.02, clamp_sd_gamma)
  
    # From sensitivity range -> beta prior (depends on severity)
    if (severity == "ci") {
      # CI: Se = Phi(beta - gamma). Use eta = beta - gamma.
      e     <- probit_mid_sd(min(rg$sens), max(rg$sens), z)
      m_eta <- e$m
      sd_eta<- clamp(e$sd, 0.02, clamp_sd_beta)
      m_b   <- max(0, m_eta + m_g)
      sd_b  <- clamp(sqrt(sd_eta^2 + sd_g^2), 0.05, clamp_sd_beta)  # simple delta
    } else {
      cf  <- calibrate_beta_to_Se(min(rg$sens), max(rg$sens),
                                  m_gamma = m_g, sd_gamma = sd_g,
                                  severity = severity, aS = aS, bS = bS, mu0 = mu0, tau = tau)
      m_b  <- cf$m_beta
      sd_b <- cf$sd_beta
    }
    
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
