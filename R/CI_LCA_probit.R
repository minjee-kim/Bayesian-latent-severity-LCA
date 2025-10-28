

CI_LCA_probit <- function(data, iterations, burnin, thin=1,
                          mu_beta, sd_beta,  # beta_j: probit Se_j
                          mu_gamma, sd_gamma, # gamma_j: probit Sp_j
                          a_rho=1, b_rho=1){
  library(truncnorm)
  Tij = as.matrix(data)
  N = nrow(Tij)
  T_col = ncol(Tij)
  rho = runif(1) 
  Di = rbinom(N, 1, rho)
  gamma_j = rnorm(T_col, mu_gamma, sd_gamma)
  beta_j = rtruncnorm(T_col, a = 0, b = Inf, mean = mu_beta, sd = sd_beta)
  
  ## initializing samples
  iterations_tot = iterations + burnin
  rho_samples <- rep(NA, (iterations_tot - burnin)/thin)
  D_samples <- matrix(NA, (iterations_tot - burnin)/thin, N)
  gamma_samples <- matrix(NA, (iterations_tot - burnin)/thin, T_col)
  beta_samples <- matrix(NA, (iterations_tot - burnin)/thin, T_col)
  sensitivity_samples <- matrix(NA, (iterations_tot - burnin)/thin, T_col)
  specificity_samples <- matrix(NA, (iterations_tot - burnin)/thin, T_col)
  
  sample_Vij <- function(Tij, Di, beta_j, gamma_j){
    N = length(Di)
    T_col = length(beta_j)
    Vij <- matrix(NA, nrow = N, ncol = T_col)
    Vij_mean <- outer(Di, beta_j, "*")
    
    for(i in 1:N){
      for(j in 1:T_col){
        tij = Tij[i, j]
        v_mean = Vij_mean[i, j] 
        gamma = gamma_j[j] 
        
        if(tij == 1){
          ## Draw Vij from the right of gamma_j if Tij = 1
          Vij[i, j] = rtruncnorm(n = 1, a = gamma, b = Inf, mean = v_mean, sd = 1)
        } else {
          ## Draw Vij from the left of gamma_j if Tij = 0
          Vij[i, j] = rtruncnorm(n = 1, a = -Inf, b = gamma, mean = v_mean, sd = 1)
        }
      }
    }
    
    return(Vij)
  }
  
  log_likelihood <- function(Tij, Di, beta_j, gamma_j){
    N <- nrow(Tij)  
    T_col <- ncol(Tij)  
    
    likelihood_mat <- matrix(NA, nrow = N, ncol = T_col)
    Vij_mean <- outer(Di, beta_j, "*")
    
    for(i in 1:N){
      for(j in 1:T_col){
        if(Tij[i,j] == 1){
          likelihood_mat[i,j] = log(pmax(1 - pnorm((gamma_j[j] - Vij_mean[i, j])), 1e-10))
        }else{
          likelihood_mat[i,j] =  log(pmax(pnorm((gamma_j[j] - Vij_mean[i, j])), 1e-10))
        }
      }
    }
    
    return(likelihood_mat)
  }
  
  
  gamma_update <- function(Vij, Di, mu_gamma, sd_gamma, beta_j) {
    N <- length(Di)
    J <- ncol(Vij)
    out <- numeric(J)
    for (j in 1:J) {
      Zj <- Vij[, j] - gamma_j[j]
      y  <- Zj - beta_j[j] * Di   
      s2_post <- 1 / (N + 1/(sd_gamma[j]^2))
      m_post  <- s2_post * (-sum(y) + mu_gamma[j]/(sd_gamma[j]^2))
      out[j]  <- rnorm(1, m_post, sqrt(s2_post))
    }
    out
  }
  
  
  
  D_update <- function(rho, Tij, beta_j, gamma_j) {
    N <- nrow(Tij)
    D_new <- rep(0, N)
    
    for (i in 1:N) {
      log_prior_1 <- log(rho)
      log_prior_0 <- log(1 - rho)
      
      loglik_1 <- sum(log_likelihood(matrix(Tij[i, ], nrow = 1),
                                     1, beta_j, gamma_j))
      
      loglik_0 <- sum(log_likelihood(matrix(Tij[i, ], nrow = 1),
                                     0, beta_j, gamma_j))
      
      log_posterior_1 <- log_prior_1 + loglik_1
      log_posterior_0 <- log_prior_0 + loglik_0
      
      logit <- (log_prior_1 + loglik_1) - (log_prior_0 + loglik_0)
      p1 <- 1 / (1 + exp(-logit))
      
      D_new[i] <- rbinom(1, 1, p1)
    }
    
    return(D_new)
  }
  
  beta_mh <- function(Di, beta_j, Vij, proposal_sd, mu_beta, sd_beta) {
    T_col <- length(beta_j)
    beta_new <- beta_j
    idx <- which(Di == 1L)
    for (j in 1:T_col) {
      if (length(idx) == 0) next
      beta_curr <- beta_j[j]
      resid_curr <- Vij[idx, j] - beta_curr  # mu = beta for D=1
      ll_curr <- sum(dnorm(resid_curr, 0, 1, log = TRUE))
      
      m_bj  <- if (length(mu_beta)  > 1) mu_beta[j]  else mu_beta
      sd_bj <- if (length(sd_beta) > 1) sd_beta[j] else sd_beta
      
      log_prior_curr <- log(dtruncnorm(beta_curr, a = 0, b = Inf, mean = m_bj, sd = sd_bj))
      
      beta_prop <- rtruncnorm(1, a = 0, b = Inf, mean = beta_curr, sd = proposal_sd)
      resid_prop <- Vij[idx, j] - beta_prop
      ll_prop <- sum(dnorm(resid_prop, 0, 1, log = TRUE))
      log_prior_prop <- log(dtruncnorm(beta_prop, a = 0, b = Inf, mean = m_bj, sd = sd_bj))
      
      log_q_cp <- log(dtruncnorm(beta_prop, a = 0, b = Inf, mean = beta_curr, sd = proposal_sd))
      log_q_pc <- log(dtruncnorm(beta_curr, a = 0, b = Inf, mean = beta_prop, sd = proposal_sd))
      
      log_alpha <- (ll_prop + log_prior_prop + log_q_pc) -
        (ll_curr + log_prior_curr + log_q_cp)
      
      if (is.finite(log_alpha) && log(runif(1)) < log_alpha) beta_new[j] <- beta_prop
    }
    beta_new
  }
  
  ### MCMC 
  for(iter in 1:iterations_tot){
    
    ## Di update 
    Di = D_update(rho, Tij, beta_j, gamma_j)
    
    ## sample Vij 
    Vij <- sample_Vij(Tij, Di, beta_j, gamma_j)
    
    ## gamma_j update 
    gamma_j <- gamma_update(Vij, Di, mu_gamma, sd_gamma, beta_j)
    
    
    ## beta_j update
    beta_j <- beta_mh(Di, beta_j, Vij, proposal_sd = 0.05, mu_beta = mu_beta, sd_beta = sd_beta)    
    
    ## rho update
    rho <- rbeta(1, a_rho + sum(Di), b_rho + (N - sum(Di)))
    
    if(iter > burnin && ((iter - burnin) %% thin == 0)){
      index <- (iter - burnin) / thin
      rho_samples[index] <- rho
      D_samples[index, ] <- Di
      gamma_samples[index, ] <- gamma_j
      beta_samples[index, ] <- beta_j
      sens <- pnorm(beta_j - gamma_j) 
      spec <- pnorm(gamma_j)
      
      sensitivity_samples[index, ] <- sens
      specificity_samples[index, ] <- spec
    }
  }
  
  return(list(
    rho_Samples = rho_samples, 
    D_Samples = D_samples, 
    gamma_Samples = gamma_samples, 
    beta_Samples = beta_samples, 
    sensitivity_Samples = sensitivity_samples,
    specificity_Samples = specificity_samples
  )) 
}
