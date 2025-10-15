

CI_LCA_probit <- function(data, iterations, burnin, thin=1,
                          m_beta, sd_beta,  # beta_j: probit Se_j
                          m_gamma, sd_gamma, # gamma_j: probit Sp_j
                          a_rho=1, b_rho=1){
  Tij <- as.matrix(data); N <- nrow(Tij); J <- ncol(Tij)
  as_lenJ <- function(x,nm){ if(length(x)==1) rep(x,J) else if(length(x)==J) x else stop(sprintf("%s len",nm)) }
  m_beta <- as_lenJ(m_beta,"m_beta"); sd_beta <- as_lenJ(sd_beta,"sd_beta")
  m_gamma<- as_lenJ(m_gamma,"m_gamma"); sd_gamma<- as_lenJ(sd_gamma,"sd_gamma")
  
  # init
  rho <- rbeta(1, a_rho, b_rho)
  D   <- rbinom(N,1,rho)
  beta <- rnorm(J, m_beta, sd_beta)
  gamma<- rnorm(J, m_gamma, sd_gamma)
  
  n_keep <- (iterations - burnin) %/% thin
  rho_s  <- numeric(n_keep)
  D_s    <- matrix(NA,n_keep,N)
  beta_s <- matrix(NA,n_keep,J)
  gamma_s<- matrix(NA,n_keep,J)
  Se_s   <- matrix(NA,n_keep,J)
  Sp_s   <- matrix(NA,n_keep,J)
  
  log_odds_collapsed <- function(D,i){
    Sdi <- sum(D) - D[i]
    log( (a_rho + Sdi) / (b_rho + (N-1) - Sdi) )
  }
  
  for (it in 1:iterations){
    Z <- matrix(NA,N,J)
    mu1 <- matrix(rep(beta, each=N), N, J)
    mu0 <- matrix(rep(-gamma, each=N), N, J)
    Mu  <- ifelse(D==1, mu1, mu0)
    A   <- ifelse(Tij==1, 0, -Inf)
    B   <- ifelse(Tij==1, Inf, 0)
    Z[] <- truncnorm::rtruncnorm(N*J, a=as.vector(A), b=as.vector(B),
                                 mean=as.vector(Mu), sd=1)
    
    for (j in 1:J){
      idx1 <- which(D==1)
      n1   <- length(idx1)
      prec0 <- 1/(sd_beta[j]^2)
      s_num <- prec0*m_beta[j] + sum(Z[idx1,j])
      prec  <- prec0 + n1
      beta[j] <- rnorm(1, s_num/prec, sqrt(1/prec))
      
      idx0 <- which(D==0)
      n0   <- length(idx0)
      prec0g <- 1/(sd_gamma[j]^2)
      s_numg <- prec0g*m_gamma[j] - sum(Z[idx0,j])
      precg  <- prec0g + n0
      gamma[j] <- rnorm(1, s_numg/precg, sqrt(1/precg))
    }
    
    for (i in 1:N){
      lp <- log_odds_collapsed(D,i)
      # ll difference for D=1 vs 0 for subject i
      mu1_i <- beta
      mu0_i <- -gamma
      # compute likelihood contribution with current observed Tij via latent Z density
      # Using prob of T given class under probit: P(T=1|class)=Phi(mu), T=0: 1-Phi(mu)
      p1 <- pnorm(mu1_i); p0 <- pnorm(mu0_i)
      ll1 <- sum( ifelse(Tij[i,]==1, log(p1), log1p(-p1)) )
      ll0 <- sum( ifelse(Tij[i,]==1, log(p0), log1p(-p0)) )
      log_acc <- lp + (ll1 - ll0)
      if (is.finite(log_acc) && log(runif(1)) < log_acc) D[i] <- 1 else D[i] <- 0
    }
    
    rho <- rbeta(1, a_rho + sum(D), b_rho + (N - sum(D)))
    
    # save
    if (it > burnin && ((it-burnin) %% thin == 0)){
      k <- (it-burnin)%/%thin
      rho_s[k]   <- rho
      D_s[k,]    <- D
      beta_s[k,] <- beta
      gamma_s[k,]<- gamma
      Se_s[k,]   <- pnorm(beta)
      Sp_s[k,]   <- pnorm(gamma)
    }
  }
  list(rho_Samples=rho_s, D_Samples=D_s, beta_Samples=beta_s, gamma_Samples=gamma_s,
       sensitivity_Samples=Se_s, specificity_Samples=Sp_s)
}
