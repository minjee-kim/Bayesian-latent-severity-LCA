

build_priors_from_ranges <- function(
    ranges,
    severity = c("CI","gamma","nm+"),
    aS = 3, bS = sqrt(3),
    mu0 = 0, tau = 1.48495,
    nsim = 100000,
    verbose = TRUE
){
  set.seed(123)
  severity <- match.arg(tolower(severity), c("ci","gamma","nm+"))
  J <- length(ranges)
  
  clamp01 <- function(p) pmin(pmax(p, 1e-6), 1 - 1e-6)
  qtrunc <- function(u, mean, sd) truncnorm::qtruncnorm(u, a = 0, b = Inf, mean = mean, sd = sd)
  
  map_spec_to_gamma <- function(Lsp, Usp){
    z   <- qnorm(0.975)
    Lsp <- clamp01(Lsp); Usp <- clamp01(Usp)
    m_g  <- 0.5 * (qnorm(Lsp) + qnorm(Usp))
    sd_g <- (qnorm(Usp) - qnorm(Lsp)) / (2*z)
    c(m_g = m_g, sd_g = sd_g)
  }
  
  sample_S_once <- function(sev, n, aS, bS, mu0, tau){
    if (sev == "ci")    return(rep(1, n))
    if (sev == "gamma") return(rgamma(n, shape = aS, rate = bS))
    x  <- seq(0, mu0 + 6*tau, length.out = 10000)
    dx <- x[2] - x[1]
    dens <- 2*(x^2)/(sqrt(2*pi)*tau^3) * exp(- (x - mu0)^2/(2*tau^2))
    dens <- dens / sum(dens*dx)
    cdf  <- cumsum(dens*dx)
    approx(cdf, x, xout = runif(n), rule = 2)$y
  }
  
  # Common random numbers
  U_beta  <- runif(nsim)
  Z_gamma <- rnorm(nsim)
  S_vec <- sample_S_once(severity, nsim, aS, bS, mu0, tau)
  
  make_calibrator <- function(S_vec, U_beta, Z_gamma){
    function(Lse, Use, m_g, sd_g){
      zL <- qnorm(clamp01(Lse))
      zU <- qnorm(clamp01(Use))
      obj <- function(par){
        m <- exp(par[1]); s <- exp(par[2])
        if (!is.finite(m) || !is.finite(s) || s <= 1e-8) return(1e12)
        beta  <- qtrunc(U_beta, mean = m, sd = s)
        gamma <- m_g + sd_g * Z_gamma
        eta   <- beta * S_vec - gamma
        qs    <- as.numeric(quantile(eta, c(0.025, 0.5, 0.975), names = FALSE))
        # match tails strongly; nudge median toward the midpoint
        (qs[1] - zL)^2 + (qs[3] - zU)^2 + 0.1*(qs[2] - 0.5*(zL+zU))^2
      }
      inits <- rbind(c(log(1.0), log(0.6)),
                     c(log(0.4), log(0.4)),
                     c(log(2.0), log(0.8)),
                     c(log(3.0), log(1.2)),
                     c(log(0.5), log(1.5)))
      best <- list(val = Inf, par = inits[1,])
      for (k in seq_len(nrow(inits))){
        o <- optim(inits[k,], obj, method = "Nelder-Mead",
                   control = list(maxit = 800, reltol = 1e-9))
        if (is.finite(o$value) && o$value < best$val) best <- list(val = o$value, par = o$par)
      }
      c(m_beta = exp(best$par[1]), sd_beta = exp(best$par[2]))
    }
  }
  calib <- make_calibrator(S_vec, U_beta, Z_gamma)
  
  mu_beta <- sd_beta <- m_gamma <- sd_gamma <- numeric(J)
  for (j in seq_len(J)){
    rj  <- ranges[[j]]
    Lsp <- min(rj$spec); Usp <- max(rj$spec)
    Lse <- min(rj$sens); Use <- max(rj$sens)
    gg <- map_spec_to_gamma(Lsp, Usp)
    bb <- calib(Lse, Use, m_g = gg["m_g"], sd_g = gg["sd_g"])
    m_gamma[j] <- gg["m_g"];  sd_gamma[j] <- gg["sd_g"]
    mu_beta[j] <- bb["m_beta"]; sd_beta[j] <- bb["sd_beta"]
  }
  
  if (verbose) {
    # Reuse the SAME S_vec used in calibration
    rows <- lapply(seq_len(J), function(j) {
      beta  <- truncnorm::rtruncnorm(nsim, a = 0, b = Inf, mean = mu_beta[j], sd = sd_beta[j])
      gamma <- stats::rnorm(nsim, mean = m_gamma[j], sd = sd_gamma[j])
      Se <- stats::pnorm(beta * S_vec - gamma)
      Sp <- stats::pnorm(gamma)
      qSe <- stats::quantile(Se, c(0.025, 0.50, 0.975), names = FALSE)
      qSp <- stats::quantile(Sp, c(0.025, 0.50, 0.975), names = FALSE)
      data.frame(
        Test=j,
        beta_mu=round(mu_beta[j],3), beta_sd=round(sd_beta[j],3),
        gamma_mu=round(m_gamma[j],3), gamma_sd=round(sd_gamma[j],3),
        Se_q025=round(qSe[1],3), Se_q50=round(qSe[2],3), Se_q975=round(qSe[3],3),
        Sp_q025=round(qSp[1],3), Sp_q50=round(qSp[2],3), Sp_q975=round(qSp[3],3)
      )
    })
    cat("\n================ PRIOR CHECK (", toupper(severity), ") =================\n", sep = "")
    print(do.call(rbind, rows))
  }
  
  list(mu_beta = mu_beta, sd_beta = sd_beta,
       m_gamma = m_gamma, sd_gamma = sd_gamma)
}



