


Gamma_LCA_severity <- function(data, iterations, burnin, thin = 1,
                               aS = 3, bS = sqrt(3),
                               mu_beta = 1, sd_beta = 5,
                               m_gamma = 2, sd_gamma = 1.5,
                               a_rho = 1, b_rho = 1) {
  library(truncnorm)
  
  Tij   <- as.matrix(data)
  N     <- nrow(Tij)
  J <- ncol(Tij)

  as_lenJ <- function(x, name){
    if (length(x) == 1) rep(x, J)
    else if (length(x) == J) x
    else stop(sprintf("%s must be length 1 or %d", name, J))
  }
  mu_beta  <- as_lenJ(mu_beta,  "mu_beta")
  sd_beta  <- as_lenJ(sd_beta,  "sd_beta")
  m_gamma  <- as_lenJ(m_gamma,  "m_gamma")
  sd_gamma <- as_lenJ(sd_gamma, "sd_gamma")
  
  iterations_tot = iterations + burnin
  n_keep <- (iterations_tot - burnin) / thin
  n_keep <- as.integer(n_keep)
  
  rho <- rbeta(1, a_rho, b_rho)
  Di <- rbinom(N, 1, rho)
  Si  <- numeric(N)
  if (any(Di == 1)) Si[Di == 1] <- rgamma(sum(Di==1), shape = aS, rate = bS)
  
  beta  <- rtruncnorm(J, a=0, b=Inf, mean = mu_beta, sd = sd_beta)
  gamma <- rnorm(J, m_gamma, sd_gamma)

  rho_samples  <- rep(NA, n_keep)
  D_samples    <- matrix(NA, n_keep, N)
  S_samples    <- matrix(NA, n_keep, N)
  beta_samples <- matrix(NA, n_keep, J)
  gam_samples  <- matrix(NA, n_keep, J)
  sens_samples <- matrix(NA, n_keep, J)
  spec_samples <- matrix(NA, n_keep, J)
  
  sample_Vij <- function(Tij, Di, Si, beta, gamma) {
    mu <- outer(Di*Si, beta) - matrix(gamma, N, J, byrow=TRUE)
    a  <- ifelse(Tij==1, 0, -Inf)
    b  <- ifelse(Tij==1, Inf, 0)
    v  <- rtruncnorm(n = length(mu),
                     a = as.vector(a), b = as.vector(b),
                     mean = as.vector(mu), sd = 1)
    matrix(v, nrow = N, ncol = J, byrow = FALSE)
  }
  
  update_beta <- function(V, Di, Si, gamma, mu_beta, sd_beta) {
    x  <- Di*Si
    x2 <- sum(x^2)
    out <- numeric(J)
    for (j in 1:J) {
      vtil  <- V[, j] + gamma[j]
      s1    <- sum(x * vtil)
      prec0 <- 1/(sd_beta[j]^2)
      prec  <- prec0 + x2
      meanj <- (prec0*mu_beta[j] + s1)/prec
      out[j] <- rtruncnorm(1, a=0, b=Inf, mean=meanj, sd=1/sqrt(prec))
    }
    out
  }
  
  update_gamma <- function(V, Di, Si, beta, m_gamma, sd_gamma) {
    x <- Di*Si
    out <- numeric(J)
    for (j in 1:J) {
      prec0 <- 1/(sd_gamma[j]^2)
      s_num <- prec0*m_gamma[j] + sum(beta[j]*x - V[, j]) 
      prec  <- prec0 + N
      out[j] <- rnorm(1, s_num/prec, sqrt(1/prec))
    }
    out
  }
  
  logpost_logS_i <- function(logS, Vij_row, beta, gamma, aS, bS) {
    S <- exp(logS); if (!is.finite(S) || S <= 0) return(-Inf)
    A <- sum(beta^2)
    B <- sum(beta * (Vij_row + gamma))
    aS*logS - bS*S - 0.5*A*S^2 + B*S
  }
  slice_sample_logS <- function(logS0, Vij_row, beta, gamma, aS, bS, w=1.0, m=50) {
    y0 <- logpost_logS_i(logS0, Vij_row, beta, gamma, aS, bS) - rexp(1)
    L <- logS0 - runif(1, 0, w); R <- L + w
    JL <- floor(runif(1, 0, m)); KR <- (m-1) - JL
    while (JL>0 && logpost_logS_i(L, Vij_row, beta, gamma, aS, bS) > y0) { L <- L - w; JL <- JL - 1 }
    while (KR>0 && logpost_logS_i(R, Vij_row, beta, gamma, aS, bS) > y0) { R <- R + w; KR <- KR - 1 }
    repeat {
      logS <- runif(1, L, R)
      if (logpost_logS_i(logS, Vij_row, beta, gamma, aS, bS) >= y0) return(logS)
      if (logS < logS0) L <- logS else R <- logS
    }
  }
  
  DS_flip <- function(i, Tij_row, Di, Si, beta, gamma, aS, bS, a_rho, b_rho) {
    Sd <- sum(Di) - Di[i]
    N  <- length(Di)
    lpodds <- log( (a_rho + Sd) / (b_rho + (N - 1) - Sd) )
    
    eps <- 1e-12
    p0 <- pnorm(-gamma)
    p0 <- pmin(pmax(p0, eps), 1 - eps)
    ll0 <- sum(Tij_row * log(p0) + (1 - Tij_row) * log1p(-p0))
    
    if (Di[i] == 0) {
      Sprop <- rgamma(1, shape = aS, rate = bS)
      eta1  <- beta * Sprop - gamma
      p1    <- pnorm(eta1)
      p1 <- pmin(pmax(p1, eps), 1 - eps)
      ll1   <- sum(Tij_row * log(p1) + (1 - Tij_row) * log1p(-p1))
      log_acc <- lpodds + (ll1 - ll0)
      if (is.finite(log_acc) && log(runif(1)) < log_acc) {
        return(list(D = 1L, S = Sprop))
      } else return(list(D = 0L, S = 0))
    } else {
      Scurr <- max(Si[i], .Machine$double.eps)
      eta1  <- beta * Scurr - gamma
      p1    <- pnorm(eta1)
      p1 <- pmin(pmax(p1, eps), 1 - eps)
      ll1   <- sum(Tij_row * log(p1) + (1 - Tij_row) * log1p(-p1))
      log_acc <- (-lpodds) + (ll0 - ll1)
      if (is.finite(log_acc) && log(runif(1)) < log_acc) {
        return(list(D = 0L, S = 0))
      } else return(list(D = 1L, S = Scurr))
    }
  }
  
  keep_idx <- 0L
  for (iter in 1:iterations_tot) {
    
    for (i in 1:N) {
      upd <- DS_flip(i, Tij[i,], Di, Si, beta, gamma, aS, bS, a_rho, b_rho)
      Di[i] <- upd$D
      Si[i] <- upd$S
    }
    
    V <- sample_Vij(Tij, Di, Si, beta, gamma)
    
    for (i in which(Di==1L)) {
      logS0 <- log(Si[i])
      logS1 <- slice_sample_logS(logS0, V[i,], beta, gamma, aS, bS)
      Si[i] <- exp(logS1)
    }
    Si[Di==0L] <- 0
    
    beta  <- update_beta(V, Di, Si, gamma, mu_beta, sd_beta)
    gamma <- update_gamma(V, Di, Si, beta, m_gamma, sd_gamma)
    
    rho <- rbeta(1, a_rho + sum(Di), b_rho + (N - sum(Di)))
    
    # Save
    if (iter > burnin && ((iter - burnin) %% thin == 0)) {
      keep_idx <- keep_idx + 1
      rho_samples[keep_idx]   <- rho
      D_samples[keep_idx, ]   <- Di
      S_samples[keep_idx, ]   <- Si
      beta_samples[keep_idx,] <- beta
      gam_samples[keep_idx, ] <- gamma

      idx <- which(Di == 1)
      sens <- if (length(idx)) {
        sapply(1:J, function(j) mean(pnorm(beta[j]*Si[idx] - gamma[j])))
      } else rep(NA, J)
      spec <- pnorm(gamma)
      sens_samples[keep_idx, ] <- sens
      spec_samples[keep_idx, ] <- spec
    }
  }
  
  list(
    rho_Samples = rho_samples,
    D_Samples   = D_samples,
    S_Samples   = S_samples,
    beta_Samples= beta_samples,
    gamma_Samples= gam_samples,
    sensitivity_Samples = sens_samples,
    specificity_Samples = spec_samples
  )
}
