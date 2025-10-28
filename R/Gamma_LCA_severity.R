


Gamma_LCA_severity <- function(data, iterations, burnin, thin = 1,
                               aS = 3, bS = sqrt(3),
                               aB = 2, bB = 1,  
                               mu_gamma = 2, sd_gamma = 1.5,
                               a_rho = 1, b_rho = 1) {    
  library(truncnorm)
  
  Tij <- as.matrix(data)
  if (!all(Tij %in% c(0,1))) stop("data must be 0/1.")
  N <- nrow(Tij)
  
  J <- ncol(Tij)
  
  as_lenJ <- function(x, nm){
    if (length(x) == 1) rep(x, J)
    else if (length(x) == J) x
    else stop(sprintf("%s must be length 1 or %d", nm, J))
  }
  aB <- as_lenJ(aB, "aB")
  bB <- as_lenJ(bB, "bB")
  mu_gamma  <- as_lenJ(mu_gamma,  "mu_gamma")
  sd_gamma  <- as_lenJ(sd_gamma,  "sd_gamma")
  
  if (aS <= 1 || bS <= 0) stop("Gamma(aS,bS) needs aS>1, bS>0.")
  if (any(aB <= 0) || any(bB <= 0)) stop("Gamma(aB,bB) needs aB>0, bB>0.")
  if (any(sd_gamma <= 0)) stop("sd_gamma must be > 0.")
  if (iterations <= burnin || thin < 1) stop("iterations>burnin and thin>=1.")
  
  iterations_tot <- iterations + burnin
  n_keep <- as.integer((iterations_tot - burnin) / thin)
  
  # initialize
  rho <- rbeta(1, a_rho, b_rho)
  Di  <- rbinom(N, 1, rho)
  Si  <- numeric(N)
  if (any(Di == 1)) Si[Di == 1] <- rgamma(sum(Di == 1), shape = aS, rate = bS)
  
  beta  <- rgamma(J, shape = aB, rate = bB)
  gamma <- rnorm(J, mean = mu_gamma, sd = sd_gamma)
  
  # storage
  rho_samples  <- rep(NA, n_keep)
  D_samples    <- matrix(NA, n_keep, N)
  S_samples    <- matrix(NA, n_keep, N)
  beta_samples <- matrix(NA, n_keep, J)
  gam_samples  <- matrix(NA, n_keep, J)
  sens_samples <- matrix(NA, n_keep, J)
  spec_samples <- matrix(NA, n_keep, J)
  
  # latent V sampler
  sample_Vij <- function(Tij, Di, Si, beta, gamma) {
    mu <- outer(Di * Si, beta) - matrix(gamma, N, J, byrow = TRUE)
    a  <- ifelse(Tij == 1, 0, -Inf)
    b  <- ifelse(Tij == 1, Inf, 0)
    v  <- rtruncnorm(length(mu), a = as.vector(a), b = as.vector(b),
                     mean = as.vector(mu), sd = 1)
    matrix(v, nrow = N, ncol = J)
  }
  
  # gamma_j | V, (D,S), beta  (conjugate Normal)
  update_gamma <- function(V, Di, Si, beta, mu_gamma, sd_gamma) {
    x <- Di * Si
    out <- numeric(length(beta))
    for (j in seq_along(beta)) {
      prec0_j <- 1 / (sd_gamma[j]^2)
      s_num   <- prec0_j * mu_gamma[j] + sum(beta[j] * x - V[, j])
      prec    <- prec0_j + length(x)
      out[j]  <- rnorm(1, mean = s_num / prec, sd = sqrt(1 / prec))
    }
    out
  }
  
  # log-posterior in theta = log(beta_j) up to additive constant
  logpost_logbeta_j <- function(theta, Vj, x, gamma_j, aB, bB) {
    bj <- exp(theta)
    # likelihood: Vj + gamma_j - x*bj ~ N(0,1)
    res <- Vj + gamma_j - x * bj
    ll  <- sum(dnorm(res, 0, 1, log = TRUE))
    # prior: bj ~ Gamma(aB,bB) with jacobian (from beta to theta)
    lp  <- (aB - 1) * log(bj) - bB * bj + theta  # +theta is |d beta / d theta| = exp(theta)
    ll + lp
  }
  
  # RW-MH on log(beta_j)
  update_beta_logmh <- function(V, Di, Si, beta, gamma, aB, bB, prop_sd_logbeta = 0.25) {
    x <- Di * Si
    for (j in 1:length(beta)) {
      theta_curr <- log(beta[j])
      # log posterior at theta = log β_j
      logpost <- function(theta){
        bj <- exp(theta)
        res <- V[, j] + gamma[j] - x * bj
        ll  <- sum(dnorm(res, 0, 1, log = TRUE))
        # prior β_j ~ Gamma(aB[j], bB[j])  + Jacobian term (theta)
        lp  <- (aB[j] - 1) * log(bj) - bB[j] * bj + theta
        ll + lp
      }
      lp_curr  <- logpost(theta_curr)
      theta_pr <- rnorm(1, mean = theta_curr, sd = prop_sd_logbeta)
      lp_prop  <- logpost(theta_pr)
      if (is.finite(lp_prop - lp_curr) && log(runif(1)) < (lp_prop - lp_curr)) {
        beta[j] <- exp(theta_pr)
      }
    }
    beta
  }
  
  
  # log posterior for log S_i
  logpost_logS_i <- function(logS, Vij_row, beta, gamma, aS, bS) {
    S <- exp(logS)
    if (!is.finite(S) || S <= 0) return(-Inf)
    A <- sum(beta^2)
    B <- sum(beta * (Vij_row + gamma))
    # prior Gamma(aS,bS): (aS-1)log S - bS*S ; plus likelihood -0.5*A S^2 + B S
    (aS - 1) * logS - bS * S - 0.5 * A * S^2 + B * S
  }
  
  slice_sample_logS <- function(logS0, Vij_row, beta, gamma, aS, bS, w = 0.5, m = 50) {
    y0 <- logpost_logS_i(logS0, Vij_row, beta, gamma, aS, bS) - rexp(1)
    L <- logS0 - runif(1, 0, w); R <- L + w
    JL <- floor(runif(1, 0, m)); KR <- (m - 1) - JL
    while (JL > 0 && logpost_logS_i(L, Vij_row, beta, gamma, aS, bS) > y0) { L <- L - w; JL <- JL - 1 }
    while (KR > 0 && logpost_logS_i(R, Vij_row, beta, gamma, aS, bS) > y0) { R <- R + w; KR <- KR - 1 }
    repeat {
      logS <- runif(1, L, R)
      if (logpost_logS_i(logS, Vij_row, beta, gamma, aS, bS) >= y0) return(logS)
      if (logS < logS0) L <- logS else R <- logS
    }
  }
  
  # joint D,S flip (unchanged)
  DS_flip <- function(i, Tij_row, Di, Si, beta, gamma, aS, bS, a_rho, b_rho) {
    Sd <- sum(Di) - Di[i]
    N  <- length(Di)
    lpodds <- log((a_rho + Sd) / (b_rho + (N - 1) - Sd))
    
    eps <- 1e-12
    p0 <- pnorm(-gamma); p0 <- pmin(pmax(p0, eps), 1 - eps)
    ll0 <- sum(Tij_row * log(p0) + (1 - Tij_row) * log1p(-p0))
    
    if (Di[i] == 0) {
      Sprop <- rgamma(1, shape = aS, rate = bS)
      p1 <- pnorm(beta * Sprop - gamma); p1 <- pmin(pmax(p1, eps), 1 - eps)
      ll1 <- sum(Tij_row * log(p1) + (1 - Tij_row) * log1p(-p1))
      if (log(runif(1)) < lpodds + (ll1 - ll0)) return(list(D = 1L, S = Sprop))
      else return(list(D = 0L, S = 0))
    } else {
      Scurr <- max(Si[i], .Machine$double.eps)
      p1 <- pnorm(beta * Scurr - gamma); p1 <- pmin(pmax(p1, eps), 1 - eps)
      ll1 <- sum(Tij_row * log(p1) + (1 - Tij_row) * log1p(-p1))
      if (log(runif(1)) < (-lpodds) + (ll0 - ll1)) return(list(D = 0L, S = 0))
      else return(list(D = 1L, S = Scurr))
    }
  }
  
  keep_idx <- 0L
  for (iter in 1:iterations_tot) {
    # D,S
    for (i in 1:N) {
      upd <- DS_flip(i, Tij[i, ], Di, Si, beta, gamma, aS, bS, a_rho, b_rho)
      Di[i] <- upd$D
      Si[i] <- upd$S
    }
    
    V <- sample_Vij(Tij, Di, Si, beta, gamma)
    
    idx1 <- which(Di == 1L)
    if (length(idx1)) {
      for (ii in idx1) {
        logS0 <- log(max(Si[ii], .Machine$double.eps))
        Si[ii] <- exp(slice_sample_logS(logS0, V[ii, ], beta, gamma, aS, bS))
      }
    }
    Si[Di == 0L] <- 0
    
    beta <- update_beta_logmh(V, Di, Si, beta, gamma, aB, bB)
    gamma <- update_gamma(V, Di, Si, beta, mu_gamma, sd_gamma)
    rho   <- rbeta(1, a_rho + sum(Di), b_rho + (N - sum(Di)))
    
    if (iter > burnin && ((iter - burnin) %% thin == 0)) {
      keep_idx <- keep_idx + 1L
      rho_samples[keep_idx]    <- rho
      D_samples[keep_idx, ]    <- Di
      S_samples[keep_idx, ]    <- Si
      beta_samples[keep_idx, ] <- beta
      gam_samples[keep_idx, ]  <- gamma
      if (length(idx1)) {
        sens_samples[keep_idx, ] <- sapply(1:J, function(j) mean(pnorm(beta[j] * Si[idx1] - gamma[j])))
      } else {
        sens_samples[keep_idx, ] <- rep(NA, J)
      }
      spec_samples[keep_idx, ] <- pnorm(gamma)
    }
  }
  
  list(
    rho_Samples         = rho_samples,
    D_Samples           = D_samples,
    S_Samples           = S_samples,
    beta_Samples        = beta_samples,
    gamma_Samples       = gam_samples,
    sensitivity_Samples = sens_samples,
    specificity_Samples = spec_samples
  )
}

