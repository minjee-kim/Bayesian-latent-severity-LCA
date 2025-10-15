


NM_LCA_severity <- function(
    data, iterations, burnin, thin = 1,
    mu_beta = 0,  sd_beta = 5,       
    m_gamma = 0, sd_gamma = 1.5,    
    mu0 = 0, tau = 1.48495,
    a_rho = 1, b_rho = 1
){
  library(truncnorm)
  
  Tij <- as.matrix(data)
  N  <- nrow(Tij)
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

  n_keep <- (iterations - burnin) / thin
  n_keep <- as.integer(n_keep)
  
  sample_NM <- function(n, mu0, tau, grid_size = 1e4) {
    x_vals <- seq(0, mu0 + 6 * tau, length.out = grid_size)
    dens_vals <- 2 * (x_vals^2) / (sqrt(2 * pi) * tau^3) * 
      exp(-((x_vals - mu0)^2) / (2 * tau^2))
    dx <- x_vals[2] - x_vals[1]
    dens_vals <- dens_vals / sum(dens_vals * dx)
    
    cdf_vals <- cumsum(dens_vals * dx)
    
    # Sample from uniform and interpolate inverse CDF
    u <- runif(n)
    samples <- approx(cdf_vals, x_vals, xout = u, rule = 2)$y
    
    return(samples)
  }
  
  rho <- rbeta(1, a_rho, b_rho)
  Di <- rbinom(N, 1, rho)
  Si <- numeric(N)
  if (any(Di==1)) {
    Si[Di==1] <- sample_NM(sum(Di==1), mu0, tau)
  }
  
  beta  <- rnorm(J, mu_beta, sd_beta)
  gamma <- rnorm(J, m_gamma, sd_gamma)
  
  rho_samples   <- numeric(n_keep)
  D_samples     <- matrix(NA, n_keep, N)
  S_samples     <- matrix(NA, n_keep, N)
  beta_samples  <- matrix(NA, n_keep, J)
  gamma_samples <- matrix(NA, n_keep, J)
  sens_samples <- matrix(NA, n_keep, J)
  spec_samples <- matrix(NA, n_keep, J)
  
  sample_Vij <- function(Tij, Di, Si, beta, gamma) {
    mu <- outer(Di * Si, beta) - matrix(gamma, N, J, byrow = TRUE)
    a  <- ifelse(Tij == 1, 0, -Inf)
    b  <- ifelse(Tij == 1, Inf, 0)
    v  <- rtruncnorm(length(mu), a = as.vector(a), b = as.vector(b),
                                mean = as.vector(mu), sd = 1)
    matrix(v, nrow = N, ncol = J, byrow = FALSE)
  }
  
  # Beta | V, D, S, gamma 
  update_beta <- function(V, D, S, gamma, mu_beta, sd_beta) {
    x  <- D * S
    x2 <- sum(x^2)
    out <- numeric(J)
    for (j in 1:J) {
      vtil   <- V[, j] + gamma[j]
      s1     <- sum(x * vtil)
      prec0  <- 1 / (sd_beta[j]^2)      
      prec   <- prec0 + x2
      meanj  <- (prec0 * mu_beta[j] + s1) / prec
      sdj    <- 1 / sqrt(prec)
      out[j] <- rtruncnorm(1, a = 0, b = Inf, mean = meanj, sd = sdj)
    }
    out
  }
  
  # Gamma | V, D, S, beta 
  update_gamma <- function(V, D, S, beta, m_gamma, sd_gamma) {
    x <- D * S
    out <- numeric(J)
    for (j in 1:J) {
      prec0  <- 1 / (sd_gamma[j]^2)         
      s_num  <- prec0 * m_gamma[j] + sum(beta[j] * x - V[, j])  
      prec   <- prec0 + N
      mu     <- s_num / prec
      sd     <- sqrt(1 / prec)
      out[j] <- rnorm(1, mean = mu, sd = sd)
    }
    out
  }
  
  logpost_logS_i <- function(logS, Vij_row, beta, gamma, mu0, tau) {
    s <- exp(logS); if (!is.finite(s) || s <= 0) return(-Inf)
    A <- sum(beta^2)
    B <- sum(beta*(Vij_row + gamma))
    3*logS - (s - mu0)^2/(2*tau^2) - 0.5*A*s^2 + B*s
  }
  
  slice_logS <- function(logS0, Vij_row, beta, gamma, mu0, tau, w=1, m=50) {
    y0 <- logpost_logS_i(logS0, Vij_row, beta, gamma, mu0, tau) - rexp(1)
    L <- logS0 - runif(1,0,w); R <- L + w
    JL <- sample.int(m,1)-1; KR <- (m-1)-JL
    while (JL>0 && logpost_logS_i(L, Vij_row, beta, gamma, mu0, tau) > y0) { L <- L - w; JL <- JL - 1 }
    while (KR>0 && logpost_logS_i(R, Vij_row, beta, gamma, mu0, tau) > y0) { R <- R + w; KR <- KR - 1 }
    repeat {
      ls <- runif(1, L, R)
      if (logpost_logS_i(ls, Vij_row, beta, gamma, mu0, tau) >= y0) return(ls)
      if (ls < logS0) L <- ls else R <- ls
    }
  }

  DS_flip_NM <- function(i, Tij_row, Di, Si, beta, gamma, mu0, tau, a_rho, b_rho) {
    Sd <- sum(Di) - Di[i]
    Nn <- length(Di)
    lpodds <- log( (a_rho + Sd) / (b_rho + (Nn - 1) - Sd) )
    eps <- 1e-12
    p0 <- pnorm(-gamma)
    p0 <- pmin(pmax(p0, eps), 1 - eps)
    ll0 <- sum(Tij_row*log(p0) + (1-Tij_row)*log1p(-p0))
    
    if (Di[i]==0) {
      Sprop <- sample_NM(1, mu0, tau)
      p1 <- pnorm(beta*Sprop - gamma)
      p1 <- pmin(pmax(p1, eps), 1 - eps)
      ll1 <- sum(Tij_row*log(p1) + (1-Tij_row)*log1p(-p1))
      log_acc <- lpodds + (ll1 - ll0)
      if (is.finite(log_acc) && log(runif(1)) < log_acc) list(D=1, S=Sprop) else list(D=0, S=0)
    } else {
      Scurr <- max(Si[i], .Machine$double.eps)
      p1 <- pnorm(beta*Scurr - gamma)
      p1 <- pmin(pmax(p1, eps), 1 - eps)
      ll1 <- sum(Tij_row*log(p1) + (1-Tij_row)*log1p(-p1))
      log_acc <- (-lpodds) + (ll0 - ll1)
      if (is.finite(log_acc) && log(runif(1)) < log_acc) list(D=0, S=0) else list(D=1, S=Scurr)
    }
  }
  
  keep <- 0
  for (iter in 1:iterations) {
    
    for (i in 1:N) {
      upd <- DS_flip_NM(i, Tij[i,], Di, Si, beta, gamma, mu0, tau, a_rho, b_rho)
      Di[i] <- upd$D
      Si[i] <- upd$S
    }
    
    V <- sample_Vij(Tij, Di, Si, beta, gamma)
    
    idx1 <- which(Di==1)
    if (length(idx1)){
      for (i in idx1) {
        logS0 <- log(max(Si[i], .Machine$double.eps)) 
        logS1 <- slice_logS(logS0, V[i,], beta, gamma, mu0, tau)
        Si[i] <- exp(logS1)
      }
    }
    Si[Di==0] <- 0
    
    beta  <- update_beta(V, Di, Si, gamma, mu_beta, sd_beta)
    gamma <- update_gamma(V, Di, Si, beta, m_gamma, sd_gamma)
  
    rho <- rbeta(1, a_rho + sum(Di), b_rho + (N - sum(Di)))
    
    # save
    if (iter > burnin && ((iter - burnin) %% thin == 0)) {
      keep <- keep + 1
      rho_samples[keep]    <- rho
      D_samples[keep, ]    <- Di
      S_samples[keep, ]    <- Si
      beta_samples[keep, ] <- beta
      gamma_samples[keep, ]<- gamma
      
      idx <- which(Di == 1)
      sens <- if (length(idx)) {
        sapply(1:J, function(j) mean(pnorm(beta[j]*Si[idx] - gamma[j])))
      } else rep(NA, J)
      spec <- pnorm(gamma)
      
      sens_samples[keep, ] <- sens
      spec_samples[keep, ] <- spec
    }
  }
  
  list(
    rho_Samples = rho_samples,
    D_Samples   = D_samples,
    S_Samples   = S_samples,
    beta_Samples= beta_samples,
    gamma_Samples= gamma_samples,
    sensitivity_Samples = sens_samples,
    specificity_Samples = spec_samples
  )
}
