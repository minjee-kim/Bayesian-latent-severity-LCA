
summarize_bayesian_lca_priors <- function(fit, ci = 0.95, nsim = 20000) {
  `%||%` <- function(x, y) if (is.null(x)) y else x
  pr <- attr(fit, "priors")
  if (is.null(pr)) stop("No prior info found on object (attr(fit, 'priors') is NULL).")
  
  sp  <- pr$severity_prior
  rp  <- pr$rho_prior 
  pls <- pr$priors   
  
  stype <- tolower(sp$type %||% "gamma")
  if (stype %in% c("nmplus","nm_plus","normal moment")) stype <- "nm+"
  
  sample_S <- switch(
    stype,
    "gamma" = {
      aS <- sp$aS; bS <- sp$bS
      function(n) rgamma(n, shape = aS, rate = bS)
    },
    "nm+" = {
      mu0 <- sp$mu0 %||% 0; tau <- sp$tau %||% 1.48495
      function(n, grid_size = 10000) {
        # inverse-CDF sampler for NM+ prior
        x  <- seq(0, mu0 + 6*tau, length.out = grid_size)
        dx <- x[2] - x[1]
        dens <- 2 * (x^2) / (sqrt(2*pi) * tau^3) * exp(- (x - mu0)^2 / (2*tau^2))
        dens <- dens / sum(dens * dx)
        cdf  <- cumsum(dens * dx)
        u    <- runif(n)
        approx(cdf, x, xout = u, rule = 2)$y
      }
    },
    "ci" = { function(n) rep(NA_real_, n) }  # not used for CI model
  )
  
  # per-test param names differ by model type
  is_CI <- (stype == "ci")
  
  # Monte Carlo to get implied Se/Sp priors
  z <- qnorm((1 + ci)/2)
  draw_beta <- function(m, s, n) {
    # beta is truncated at 0 for severity models
    rnorm(n, mean = m, sd = s)
  }
  
  test_tbl <- do.call(rbind, lapply(seq_along(pls), function(j) {
    pj <- pls[[j]]
    if (is_CI) {
      # CI model: Se prior via beta ~ N(m_beta, sd_beta), Sp via gamma ~ N(m_gamma, sd_gamma)
      m_beta <- pj$m_beta; sd_beta <- pj$sd_beta
      m_g    <- pj$m_gamma; sd_g   <- pj$sd_gamma
      
      beta_draw <- rnorm(nsim, m_beta, sd_beta)
      gamma_draw <- rnorm(nsim, m_g, sd_g)
      
      se_draw  <- pnorm(beta_draw)
      sp_draw  <- pnorm(gamma_draw)
      
      data.frame(
        test = j,
        prior_type = "CI(probit)",
        par1 = sprintf("beta ~ N(%.3f, %.3f^2)", m_beta, sd_beta),
        par2 = sprintf("gamma ~ N(%.3f, %.3f^2)", m_g, sd_g),
        Se_mean = mean(se_draw), Se_L = quantile(se_draw, (1-ci)/2), Se_U = quantile(se_draw, 1-(1-ci)/2),
        Sp_mean = mean(sp_draw), Sp_L = quantile(sp_draw, (1-ci)/2), Sp_U = quantile(sp_draw, 1-(1-ci)/2)
      )
    } else {
      # Severity models: beta ~ N^(m_beta, sd_beta^2), gamma ~ N(m_gamma, sd_gamma^2)
      m_b <- pj$m_beta; sd_b <- pj$sd_beta
      m_g <- pj$m_gamma; sd_g <- pj$sd_gamma
      
      beta_draw  <- draw_beta(m_b, sd_b, nsim)
      gamma_draw <- rnorm(nsim, m_g, sd_g)
      S_draw     <- sample_S(nsim)
      
      se_draw <- pnorm(beta_draw * S_draw - gamma_draw)
      sp_draw <- pnorm(gamma_draw)
      
      data.frame(
        test = j,
        prior_type = if (stype == "gamma") "Gamma-severity" else "NM+-severity",
        par1 = sprintf("beta ~ TN(%.3f, %.3f^2; [0,∞))", m_b, sd_b),
        par2 = sprintf("gamma ~ N(%.3f, %.3f^2)", m_g, sd_g),
        Se_mean = mean(se_draw), Se_L = quantile(se_draw, (1-ci)/2), Se_U = quantile(se_draw, 1-(1-ci)/2),
        Sp_mean = mean(sp_draw), Sp_L = quantile(sp_draw, (1-ci)/2), Sp_U = quantile(sp_draw, 1-(1-ci)/2)
      )
    }
  }))
  
  a <- unname(rp["a"])
  b <- unname(rp["b"])
  prev_mean <- a / (a + b)
  prev_L <- qbeta((1-ci)/2, a, b)
  prev_U <- qbeta(1 - (1-ci)/2, a, b)
  
  list(
    severity_type = stype,
    severity_prior = sp,
    prevalence_prior = list(
      beta = c(a = a, b = b),
      mean = prev_mean,
      ci = c(L = prev_L, U = prev_U)
    ),
    tests = test_tbl
  )
}




