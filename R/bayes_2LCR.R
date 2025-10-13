

library(truncnorm)

bayes_2LCR <- function(data, model = c("2LC","2LCR1","2LCR"), 
                       iterations = 6000, burnin = 3000, thin = 1,
                       a_rho = 0.5, b_rho = 0.5, # rho ~ Beta(a_rho, b_rho)
                       mu_a0 = 1, s2_a0 = 5^2, # a_{j,0} ~ N(mu_a0, s2_a0)
                       mu_a1 = 1, s2_a1 = 5^2, # a_{j,1} ~ N(mu_a1, s2_a1)
                       mu_b0 = 1, s2_b0 = 5^2, # b_{j,0} ~ N(mu_b0, s2_b0)
                       mu_b1 = 1, s2_b1 = 5^2 # b_{j,1} ~ N(mu_b1, s2_b1)
                       ){
  
  Y <- as.matrix(data)
  N <- nrow(Y); J <- ncol(Y)
  model <- match.arg(model)
  
  # identifiability
  if (model == "2LC"   && J <= 3) stop("2LC requires p > 3 tests.")
  if (model == "2LCR1" && J <= 4) stop("2LCR1 requires p > 4 tests.")
  
  vJ <- function(x) if (length(x) == 1) rep(x, J) else x
  mu_a0 <- vJ(mu_a0)
  s2_a0 <- vJ(s2_a0)
  mu_a1 <- vJ(mu_a1)
  s2_a1 <- vJ(s2_a1)
  mu_b0 <- vJ(mu_b0)
  s2_b0 <- vJ(s2_b0)
  mu_b1 <- vJ(mu_b1)
  s2_b1 <- vJ(s2_b1)
  
  # Storage
  n_keep <- floor((iterations - burnin)/thin)
  RHO <- numeric(n_keep)
  A0 <- matrix(NA, n_keep, J)  # class 0 intercepts
  A1 <- matrix(NA, n_keep, J)  # class 1 intercepts
  B0 <- matrix(NA, n_keep, J)  # class 0 slopes (for 2LCR)
  B1 <- matrix(NA, n_keep, J)  # class 1 slopes (for 2LCR)
  SENS <- matrix(NA, n_keep, J)
  SPEC <- matrix(NA, n_keep, J)
  D_keep <- matrix(NA, n_keep, N)

  # initialize
  rho <- 0.1
  D <- rbinom(N, 1, rho)
  a0 <- rnorm(J, mean = mu_a0, sd = sqrt(s2_a0))
  a1 <- rnorm(J, mean = mu_a1, sd = sqrt(s2_a1))

  if (model == "2LC") { # CI model: no random effects
    b0 <- rep(0, J)  # slopes are zero
    b1 <- rep(0, J)
    I  <- rep(0, N)  # no subject effect
  } else if (model == "2LCR1") { ## Qu & Kutner 1996 (and Dendukuri & Joseph 2001)
    # one slope per class
    b0 <- rnorm(1, mu_b0[1], sqrt(s2_b0[1]))
    b1 <- rnorm(1, mu_b1[1], sqrt(s2_b1[1]))
    I  <- rnorm(N, 0, 1)  # subject effects ~ N(0,1)
  } else if (model == "2LCR"){
    # 2LCR: per-test slopes 
    b0 <- rnorm(J, mu_b0, sqrt(s2_b0))
    b1 <- rnorm(J, mu_b1, sqrt(s2_b1))
    I  <- rnorm(N, 0, 1)
  }
  
  keep <- 0
  cap <- 8 ## to use for capping 
  
  for (it in 1:iterations){
  
    # 1) latent Z
    if (model == "2LC") {
      eta0 <- matrix(a0, N, J, byrow = TRUE)
      eta1 <- matrix(a1, N, J, byrow = TRUE)
    } else if (model == "2LCR1") {
      eta0 <- matrix(a0, N, J, byrow = TRUE) + matrix(b0 * I, N, J, byrow = FALSE)
      eta1 <- matrix(a1, N, J, byrow = TRUE) + matrix(b1 * I, N, J, byrow = FALSE)
    } else if (model == "2LCR"){
      eta0 <- matrix(a0, N, J, byrow = TRUE) + tcrossprod(I, b0)
      eta1 <- matrix(a1, N, J, byrow = TRUE) + tcrossprod(I, b1)
    }
    
    eta0 <- pmin(pmax(eta0, -cap), cap)
    eta1 <- pmin(pmax(eta1, -cap), cap)
    
    MU <- eta0
    if (any(D == 1)) MU[D == 1, ] <- eta1[D == 1, ]
    
    Z <- matrix(0, N, J)
    m1 <- (Y == 1)
    m1[is.na(m1)] <- FALSE
    m0 <- (Y == 0)
    m0[is.na(m0)] <- FALSE
    
    s1 <- sum(m1)
    if (s1 > 0) Z[m1] <- rtruncnorm(s1, a = 0,    b = Inf, mean = MU[m1], sd = 1)
    s0 <- sum(m0)
    if (s0 > 0) Z[m0] <- rtruncnorm(s0, a = -Inf, b = 0,   mean = MU[m0], sd = 1)

    # 2) update a0_j, a1_j conjugate normals
    idx0_sub <- which(D == 0)
    idx1_sub <- which(D == 1)
    n0_sub <- length(idx0_sub)
    n1_sub <- length(idx1_sub)
    
    # for class 0
    if (n0_sub > 0){
      Z0 <- Z[idx0_sub, , drop = FALSE]
      I0 <- I[idx0_sub]
      
      S0 <- switch(model,
                   "2LC"   = Z0,
                   "2LCR1" = Z0 - matrix(I0 * b0, n0_sub, J, byrow = FALSE),
                   "2LCR"  = Z0 - (matrix(I0, n0_sub, J, byrow = FALSE) * matrix(b0, n0_sub, J, byrow = TRUE))
      )
      post_var0  <- 1 / (n0_sub + 1 / s2_a0)
      post_mean0 <- post_var0 * (colSums(S0) + mu_a0 / s2_a0)
      a0 <- rnorm(J, post_mean0, sqrt(post_var0))
    } else {
      a0 <- rnorm(J, mu_a0, sqrt(s2_a0))
    }
    
    # for class 1
    if (n1_sub > 0){
      Z1 <- Z[idx1_sub, , drop = FALSE]; I1 <- I[idx1_sub]
      S1 <- switch(model,
                   "2LC"   = Z1,
                   "2LCR1" = Z1 - matrix(I1 * b1, n1_sub, J, byrow = FALSE),
                   "2LCR"  = Z1 - (matrix(I1, n1_sub, J, byrow = FALSE) * matrix(b1, n1_sub, J, byrow = TRUE))
      )
      post_var1  <- 1 / (n1_sub + 1 / s2_a1)
      post_mean1 <- post_var1 * (colSums(S1) + mu_a1 / s2_a1)
      a1 <- rnorm(J, post_mean1, sqrt(post_var1))
    } else {
      a1 <- rnorm(J, mu_a1, sqrt(s2_a1))
    }
    
    # 3) Update random effects I (2LCR only)
    if (model != "2LC") {
      if (n0_sub > 0) {
        diff0 <- Z[idx0_sub, , drop = FALSE] - matrix(a0, n0_sub, J, byrow = TRUE)
        b0_vec <- if (model == "2LCR1") rep(b0, J) else b0
        R0 <- rowSums(diff0 * matrix(b0_vec, n0_sub, J, byrow = TRUE))
        varI0  <- 1 / (sum(b0_vec^2) + 1)  # prior Var(I)=1
        I[idx0_sub] <- rnorm(n0_sub, varI0 * R0, sqrt(varI0))
      }
      if (n1_sub > 0) {
        diff1 <- Z[idx1_sub, , drop = FALSE] - matrix(a1, n1_sub, J, byrow = TRUE)
        b1_vec <- if (model == "2LCR1") rep(b1, J) else b1
        R1 <- rowSums(diff1 * matrix(b1_vec, n1_sub, J, byrow = TRUE))
        varI1  <- 1 / (sum(b1_vec^2) + 1)
        I[idx1_sub] <- rnorm(n1_sub, varI1 * R1, sqrt(varI1))
      }
    }
    
    # 4) Update slopes b0, b1 (2LCR only)
    if (model == "2LCR1") {
      # equal slopes (scalars)
      if (n0_sub > 0) {
        I0 <- I[idx0_sub]
        diff0 <- Z[idx0_sub, , drop = FALSE] - matrix(a0, n0_sub, J, byrow = TRUE)
        sum_term <- sum(I0 * rowSums(diff0))
        var_b0  <- 1 / (J * sum(I0^2) + 1 / s2_b0[1])
        b0 <- rnorm(1, var_b0 * (sum_term + mu_b0[1]/s2_b0[1]), sqrt(var_b0))
      } else b0 <- rnorm(1, mu_b0[1], sqrt(s2_b0[1]))
      
      if (n1_sub > 0) {
        I1 <- I[idx1_sub]
        diff1 <- Z[idx1_sub, , drop = FALSE] - matrix(a1, n1_sub, J, byrow = TRUE)
        sum_term <- sum(I1 * rowSums(diff1))
        var_b1  <- 1 / (J * sum(I1^2) + 1 / s2_b1[1])
        b1 <- rnorm(1, var_b1 * (sum_term + mu_b1[1]/s2_b1[1]), sqrt(var_b1))
      } else b1 <- rnorm(1, mu_b1[1], sqrt(s2_b1[1]))
      
    } else if (model == "2LCR") {
      # per-test slopes (vectorized)
      if (n0_sub > 0) {
        I0 <- I[idx0_sub]; I0I0 <- sum(I0^2)
        diff0 <- Z[idx0_sub, , drop = FALSE] - matrix(a0, n0_sub, J, byrow = TRUE)
        IR0 <- matrix(I0, n0_sub, J, byrow = FALSE) * diff0
        post_var_b0  <- 1 / (I0I0 + 1 / s2_b0)
        post_mean_b0 <- post_var_b0 * (colSums(IR0) + mu_b0 / s2_b0)
        b0 <- rnorm(J, post_mean_b0, sqrt(post_var_b0))
      } else b0 <- rnorm(J, mu_b0, sqrt(s2_b0))
      
      if (n1_sub > 0) {
        I1 <- I[idx1_sub]; I1I1 <- sum(I1^2)
        diff1 <- Z[idx1_sub, , drop = FALSE] - matrix(a1, n1_sub, J, byrow = TRUE)
        IR1 <- matrix(I1, n1_sub, J, byrow = FALSE) * diff1
        post_var_b1  <- 1 / (I1I1 + 1 / s2_b1)
        post_mean_b1 <- post_var_b1 * (colSums(IR1) + mu_b1 / s2_b1)
        b1 <- rnorm(J, post_mean_b1, sqrt(post_var_b1))
      } else b1 <- rnorm(J, mu_b1, sqrt(s2_b1))
    }
    
    if (model == "2LC") {
      eta0 <- matrix(a0, N, J, byrow = TRUE)
      eta1 <- matrix(a1, N, J, byrow = TRUE)
    } else if (model == "2LCR1") {
      eta0 <- matrix(a0, N, J, byrow = TRUE) + matrix(b0 * I, N, J, byrow = FALSE)
      eta1 <- matrix(a1, N, J, byrow = TRUE) + matrix(b1 * I, N, J, byrow = FALSE)
    } else {
      eta0 <- matrix(a0, N, J, byrow = TRUE) + tcrossprod(I, b0)
      eta1 <- matrix(a1, N, J, byrow = TRUE) + tcrossprod(I, b1)
    }
    eta0 <- pmin(pmax(eta0, -cap), cap)
    eta1 <- pmin(pmax(eta1, -cap), cap)
    
    ll1 <- ifelse(Y == 1, pnorm(eta1, log.p = TRUE),
                  pnorm(eta1, lower.tail = FALSE, log.p = TRUE))
    ll0 <- ifelse(Y == 1, pnorm(eta0, log.p = TRUE),
                  pnorm(eta0, lower.tail = FALSE, log.p = TRUE))
    ll1[is.na(Y)] <- 0
    ll0[is.na(Y)] <- 0
    
    loglik1 <- rowSums(ll1)
    loglik0 <- rowSums(ll0)
    logit_p1 <- (loglik1 - loglik0) + (log(rho) - log1p(-rho))
    D <- rbinom(N, 1, plogis(logit_p1))
    
    
    # prevalence
    rho <- rbeta(1, a_rho + sum(D), b_rho + (N - sum(D)))    
    
    
    # save
    if (it > burnin && ((it - burnin) %% thin == 0)) {
      keep <- keep + 1
      RHO[keep]  <- rho
      A0[keep, ] <- a0
      A1[keep, ] <- a1
      if (model == "2LCR1") {
        B0[keep, ] <- rep(b0, J)
        B1[keep, ] <- rep(b1, J)
      } else {
        B0[keep, ] <- b0
        B1[keep, ] <- b1
      }
      SENS[keep, ] <- pnorm( a1 / sqrt(1 + B1[keep, ]^2) )
      SPEC[keep, ] <- pnorm( -a0 / sqrt(1 + B0[keep, ]^2) )
      D_keep[keep, ] <- D
    }
  }
  
  out <- list(rho  = RHO, 
              a0 = A0, 
              a1 = A1, 
              b0 = B0, 
              b1 = B1,
              sens = SENS, 
              spec = SPEC, 
              D = D_keep)
  
  return(out)
}
