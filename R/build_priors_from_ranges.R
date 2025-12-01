

build_priors_from_ranges <- function(
    ranges,
    severity_type = c("CI","gamma"),
    aS = 3, bS = sqrt(3),
    ci_level = 0.95,
    nsim_cal = 100000,
    seed = 123,
    spec_level = 0.95,
    n_implied = 100000
){
  severity_type <- match.arg(tolower(severity_type), c("ci","gamma"))
  stopifnot(is.list(ranges), length(ranges) > 0)
  J <- length(ranges)
  
  clamp01 <- function(p) pmin(pmax(p, 1e-6), 1 - 1e-6)
  q3 <- function(x) as.numeric(stats::quantile(x, c(0.025, 0.5, 0.975), names = FALSE))
  
  # Convert an equal-tail Beta CI [L,U] at given level -> Beta(a,b)
  beta_from_interval <- function(L, U, level = 0.95) {
    L <- clamp01(min(L, U))
    U <- clamp01(max(L, U))
    alpha <- (1 - level)/2
    probs <- c(alpha, 1 - alpha)
    
    obj <- function(par) {
      a <- exp(par[1]); b <- exp(par[2])
      q <- stats::qbeta(probs, a, b)
      sum((q - c(L, U))^2)
    }
    
    o <- stats::optim(
      par    = c(0, 0),
      fn     = obj,
      method = "BFGS"
    )
    a <- exp(o$par[1]); b <- exp(o$par[2])
    list(a = a, b = b, opt = o)
  }
  
  # Map specificity interval -> probit(Sp) Normal prior via Beta + MC
  map_spec_to_gamma <- function(Lsp, Usp, level = spec_level, nspec = nsim_cal){
    bf <- beta_from_interval(Lsp, Usp, level)
    Sp_draw <- stats::rbeta(nspec, bf$a, bf$b)
    gam <- stats::qnorm(clamp01(Sp_draw))
    mu <- mean(gam)
    sd <- sd(gam)
    sd <- min(max(sd, 0.001), 2)
    c(mu_gamma = mu, sd_gamma = sd, a_sp = bf$a, b_sp = bf$b)
  }
  
  alpha <- (1 - ci_level)/2
  quantile_probs <- c(alpha, 0.5, 1 - alpha)
  taus <- quantile_probs   # target CDF levels at the elicited Se quantiles
  
  sd_beta_bounds <- c(0.001, 2)
  mu_gamma <- sd_gamma <- numeric(J)
  
  # store underlying Beta priors for Se/Sp
  se_beta_ab <- matrix(NA_real_, nrow = J, ncol = 2,
                       dimnames = list(NULL, c("a","b")))
  sp_beta_ab <- matrix(NA_real_, nrow = J, ncol = 2,
                       dimnames = list(NULL, c("a","b")))
  
  if (severity_type == "ci"){
    set.seed(seed)
    U_B <- stats::runif(nsim_cal)
    Zg  <- stats::rnorm(nsim_cal)
    
    mu_beta <- sd_beta <- numeric(J)
    opt_log <- vector("list", J)
    
    for(j in seq_len(J)){
      r <- ranges[[j]]
      
      # Specificity -> probit Normal via Beta + MC
      gs <- map_spec_to_gamma(r$spec[1], r$spec[2],
                              level = spec_level, nspec = nsim_cal)
      mu_g <- unname(gs["mu_gamma"])
      sd_g <- unname(gs["sd_gamma"])
      mu_gamma[j] <- mu_g
      sd_gamma[j] <- sd_g
      sp_beta_ab[j, ] <- c(gs["a_sp"], gs["b_sp"])
      
      # Sensitivity Beta prior from interval
      se_b <- beta_from_interval(r$sens[1], r$sens[2], level = ci_level)
      se_beta_ab[j, ] <- c(se_b$a, se_b$b)
      se_targets <- stats::qbeta(quantile_probs, se_b$a, se_b$b)
      
      # crude initial guess ignoring variability, just to start optim
      mb  <- stats::qnorm(se_targets[2]) + mu_g
      sdb <- (stats::qnorm(se_targets[3]) - stats::qnorm(se_targets[1])) /
        (2 * stats::qnorm(1 - alpha))
      sdb <- max(min(sdb, sd_beta_bounds[2]), sd_beta_bounds[1])
      
      obj <- function(par){
        mu <- par[1]
        sd <- exp(par[2])
        beta  <- truncnorm::qtruncnorm(U_B, a = 0, b = Inf, mean = mu, sd = sd)
        gamma <- mu_g + sd_g * Zg
        Se <- stats::pnorm(beta - gamma)
        
        F_model <- vapply(se_targets, function(qt) mean(Se <= qt), numeric(1))
        sum((F_model - taus)^2)
      }
      
      par_init <- c(mb, log(sdb))
      o <- stats::optim(
        par     = par_init,
        fn      = obj,
        method  = "L-BFGS-B",
        lower   = c(-10, log(sd_beta_bounds[1])),
        upper   = c( 10, log(sd_beta_bounds[2])),
        control = list(maxit = 600, factr = 1e6)
      )
      
      mu_beta[j] <- o$par[1]
      sd_beta[j] <- exp(o$par[2])
      opt_log[[j]] <- list(convergence = o$convergence, value = o$value, message = o$message)
    }
    
    out <- list(
      mu_beta    = mu_beta,
      sd_beta    = sd_beta,
      mu_gamma   = mu_gamma,
      sd_gamma   = sd_gamma,
      severity   = "ci",
      se_beta_ab = se_beta_ab,
      sp_beta_ab = sp_beta_ab,
      optim      = opt_log
    )
    
  } else {
    if (!is.finite(aS) || !is.finite(bS) || aS <= 0 || bS <= 0)
      stop("For gamma severity: aS>0 and bS>0 (rate parameterization).")
    
    set.seed(seed)
    Zg  <- stats::rnorm(nsim_cal)
    U_S <- stats::runif(nsim_cal)
    U_B <- stats::runif(nsim_cal)
    
    mu_beta <- sd_beta <- numeric(J)
    opt_log <- vector("list", J)
    
    # numeric solver for median Se under prior -> beta "scale" b
    approx_beta <- function(p, mu_g, sd_g, aS, bS, nsim = 4000) {
      p <- clamp01(p)
      S <- stats::rgamma(nsim, shape = aS, rate = bS)
      gamma <- stats::rnorm(nsim, mean = mu_g, sd = sd_g)
      f_b <- function(b) {
        mean(stats::pnorm(b * S - gamma)) - p
      }
      grid <- exp(seq(log(1e-3), log(50), length.out = 40))
      vals <- sapply(grid, f_b)
      if (all(!is.finite(vals))) return(1)
      sgn <- sign(vals)
      if (!(any(sgn > 0) && any(sgn < 0))) {
        return(grid[which.min(abs(vals))])
      }
      idx <- which(sgn[-length(sgn)] * sgn[-1L] < 0)[1]
      uniroot(f_b, interval = c(grid[idx], grid[idx + 1]))$root
    }
    
    mu_bounds <- c(0, 20)
    sd_bounds <- c(0.005, 5)
    
    for (j in seq_len(J)){
      r  <- ranges[[j]]
      
      # Specificity -> probit Normal via Beta + MC
      gs <- map_spec_to_gamma(r$spec[1], r$spec[2],
                              level = spec_level, nspec = nsim_cal)
      mu_g <- unname(gs["mu_gamma"])
      sd_g <- unname(gs["sd_gamma"])
      mu_gamma[j] <- mu_g
      sd_gamma[j] <- sd_g
      sp_beta_ab[j, ] <- c(gs["a_sp"], gs["b_sp"])
      
      # Sensitivity Beta prior from interval
      se_b <- beta_from_interval(r$sens[1], r$sens[2], level = ci_level)
      se_beta_ab[j, ] <- c(se_b$a, se_b$b)
      se_targets <- stats::qbeta(quantile_probs, se_b$a, se_b$b)
      
      # initial guess for beta location via numeric prior-predictive match of median
      b_mid <- approx_beta(se_targets[2], mu_g, sd_g, aS, bS)
      mu0   <- max(min(b_mid, mu_bounds[2]), mu_bounds[1])
      sd0   <- max(min(mu0 / 2, sd_bounds[2]), sd_bounds[1])
      
      obj <- function(par){
        mu <- par[1]
        sd <- exp(par[2])
        beta  <- truncnorm::qtruncnorm(U_B, a = 0, b = Inf, mean = mu, sd = sd)
        S     <- stats::qgamma(U_S, shape = aS, rate = bS)
        gamma <- mu_g + sd_g * Zg
        Se <- stats::pnorm(beta * S - gamma)
        
        F_model <- vapply(se_targets, function(qt) mean(Se <= qt), numeric(1))
        sum((F_model - taus)^2)
      }
      
      par_init <- c(mu0, log(sd0))
      o <- stats::optim(
        par     = par_init,
        fn      = obj,
        method  = "L-BFGS-B",
        lower   = c(mu_bounds[1], log(sd_bounds[1])),
        upper   = c(mu_bounds[2], log(sd_bounds[2])),
        control = list(maxit = 800, factr = 1e6)
      )
      
      mu_beta[j] <- o$par[1]
      sd_beta[j] <- exp(o$par[2])
      opt_log[[j]] <- list(convergence = o$convergence, value = o$value, message = o$message)
    }
    
    out <- list(
      mu_beta    = mu_beta,
      sd_beta    = sd_beta,
      mu_gamma   = mu_gamma,
      sd_gamma   = sd_gamma,
      aS         = aS,
      bS         = bS,
      severity   = "gamma",
      se_beta_ab = se_beta_ab,
      sp_beta_ab = sp_beta_ab,
      optim      = opt_log
    )
  }
  
  # Implied Se/Sp summary under the calibrated priors
  implied_rows <- vector("list", J)
  for (j in seq_len(J)){
    gam <- stats::rnorm(n_implied, mean = mu_gamma[j], sd = sd_gamma[j])
    if (severity_type == "ci"){
      bet <- truncnorm::rtruncnorm(n_implied, a = 0, b = Inf,
                                   mean = out$mu_beta[j], sd = out$sd_beta[j])
      S   <- 1
    } else {
      bet <- truncnorm::rtruncnorm(n_implied, a = 0, b = Inf,
                                   mean = out$mu_beta[j], sd = out$sd_beta[j])
      S   <- stats::rgamma(n_implied, shape = out$aS, rate = out$bS)
    }
    Se <- stats::pnorm(bet * S - gam)
    Sp <- stats::pnorm(gam)
    qsSe <- q3(Se)
    qsSp <- q3(Sp)
    implied_rows[[j]] <- data.frame(
      test   = j,
      Se_q025 = qsSe[1], Se_q50 = qsSe[2], Se_q975 = qsSe[3],
      Sp_q025 = qsSp[1], Sp_q50 = qsSp[2], Sp_q975 = qsSp[3],
      Se_L = ranges[[j]]$sens[1], Se_U = ranges[[j]]$sens[2],
      Sp_L = ranges[[j]]$spec[1], Sp_U = ranges[[j]]$spec[2]
    )
  }
  out$implied <- do.call(rbind, implied_rows)
  
  cat("\n================ Implied Se/Sp Prior Quantiles ================\n")
  print(out$implied, row.names = FALSE)
  cat("===============================================================\n")
  
  out
}


