

CI_LCA_probit <- function(data, iterations, burnin, thin=1,
                          mu_beta, sd_beta,  # beta_j: probit Se_j
                          m_gamma, sd_gamma, # gamma_j: probit Sp_j
                          a_rho=1, b_rho=1){
  
  Tij = as.matrix(data)
  N = nrow(Tij)
  T_col = ncol(Tij)
  rho = runif(1) 
  Di = rbinom(N, 1, rho)
  gamma_j = rnorm(T_col, m_gamma, sd_gamma)
  beta_j = rtruncnorm(T_col, a = 0, b = Inf, mean = 0, sd = 1) # beta_j ~ norm(0, 2^2)
  
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
  
  
  gamma_update <- function(Tij, Di, beta_j, gamma_j, Vij, proposal_sd = 0.1,
                           m_gamma, sd_gamma) {
    N <- nrow(Tij); T_col <- ncol(Tij)
    for (j in 1:T_col) {
      gamma_curr <- gamma_j[j]
      
      # bounds from truncated-normal augmentation
      Vj <- Vij[, j]
      V0 <- Vj[Tij[, j] == 0]
      V1 <- Vj[Tij[, j] == 1]
      if (length(V0) == 0L || length(V1) == 0L) next
      
      a <- max(V0) - 0.05; b <- min(V1) + 0.05
      if (a >= b) next
      
      gamma_prop <- rtruncnorm(1, a = a, b = b, mean = gamma_curr, sd = proposal_sd)
      
      mu_ij  <- beta_j[j] * Di
      z_curr <- gamma_curr - mu_ij
      z_prop <- gamma_prop - mu_ij
      
      loglik_curr <- sum(ifelse(Tij[, j] == 1, log1p(-pnorm(z_curr)),
                                log(pmax(pnorm(z_curr), 1e-12))))
      loglik_prop <- sum(ifelse(Tij[, j] == 1, log1p(-pnorm(z_prop)),
                                log(pmax(pnorm(z_prop), 1e-12))))
      
      # per-test hyperparameters (support scalar or vector inputs)
      m_gj  <- if (length(m_gamma)  > 1) m_gamma[j]  else m_gamma
      sd_gj <- if (length(sd_gamma) > 1) sd_gamma[j] else sd_gamma
      
      log_prior_curr <- dnorm(gamma_curr, mean = m_gj,  sd = sd_gj, log = TRUE)
      log_prior_prop <- dnorm(gamma_prop, mean = m_gj,  sd = sd_gj, log = TRUE)
      
      log_q_curr_to_prop <- log(dtruncnorm(gamma_prop, a = a, b = b, mean = gamma_curr, sd = proposal_sd) + 1e-12)
      log_q_prop_to_curr <- log(dtruncnorm(gamma_curr, a = a, b = b, mean = gamma_prop, sd = proposal_sd) + 1e-12)
      
      log_alpha <- (loglik_prop + log_prior_prop + log_q_prop_to_curr) -
        (loglik_curr + log_prior_curr + log_q_curr_to_prop)
      
      # both comparisons are scalars now
      if (is.finite(log_alpha) && (log(runif(1)) < log_alpha) && (gamma_prop <= 5)) {
        gamma_j[j] <- gamma_prop
      }
    }
    gamma_j
  }
  
  
  ### MCMC 
  for(iter in 1:iterations_tot){
    
    ## Di update 
    Di = D_update(rho, Tij, beta_j, gamma_j)
    
    ## sample Vij 
    Vij <- sample_Vij(Tij, Di, beta_j, gamma_j)
    
    ## gamma_j update 
    gamma_j <- gamma_update(Tij, Di, beta_j, gamma_j, Vij, proposal_sd = 0.01,
                            m_gamma = m_gamma, sd_gamma = sd_gamma)
    
    ## beta_j update
    beta_j <- beta_mh(Di, beta_j, Vij, proposal_sd = 0.01, mu_beta = mu_beta, sd_beta = sd_beta)    
    
    ## rho update
    rho = rbeta(1, 1 + sum(Di), 1 + (N - sum(Di)))
    
    if(iter > burnin && ((iter - burnin) %% thin == 0)){
      index <- (iter - burnin) / thin
      rho_samples[index] <- rho
      D_samples[index, ] <- Di
      gamma_samples[index, ] <- gamma_j
      beta_samples[index, ] <- beta_j
      
      sens = sapply(1:T_col, function(j) mean(1 - pnorm(gamma_j[j], mean = beta_j[j] * Di[Di == 1])) )
      spec = pnorm(gamma_j)
      
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
