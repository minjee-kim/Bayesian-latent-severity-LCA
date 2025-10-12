

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
  N <- nrow(Y)
  J <- ncol(Y)
  model <- match.arg(model)
  
  # identifiability checks
  if (model == "2LC" && J <= 3)  stop("2LC (CI) requires p > 3 tests for identifiability.")
  if (model == "2LCR1" && J <= 4) stop("2LCR1 (equal slopes) requires p > 4 tests for identifiability.")
  
  
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
      # b0 and b1 are scalars (length 1)
      eta0 <- matrix(a0, N, J, byrow = TRUE) + matrix(b0 * I, N, J, byrow = FALSE)
      eta1 <- matrix(a1, N, J, byrow = TRUE) + matrix(b1 * I, N, J, byrow = FALSE)
    } else if (model == "2LCR"){
      # per-test slopes: b0 and b1 are vectors length J
      eta0 <- matrix(a0, N, J, byrow = TRUE) + (I %*% t(b0))
      eta1 <- matrix(a1, N, J, byrow = TRUE) + (I %*% t(b1))
    }
    eta0 <- pmin(pmax(eta0, -cap), cap)
    eta1 <- pmin(pmax(eta1, -cap), cap)
    
    idx1 <- (Y == 1)
    idx0 <- (Y == 0)

    MU <- eta0
    if (any(D == 1L)) MU[D == 1L, ] <- eta1[D == 1L, ]
    
    Z <- matrix(0, N, J)
    n1 <- sum(idx1)
    n0 <- sum(idx0)

    if (n1 > 0) Z[idx1] <- rtruncnorm(n1, a = 0,    b = Inf, mean = MU[idx1], sd = 1)
    if (n0 > 0) Z[idx0] <- rtruncnorm(n0, a = -Inf, b = 0,   mean = MU[idx0], sd = 1)
    
    # 2) update a0_j, a1_j conjugate normals
    idx0_sub <- which(D == 0)
    idx1_sub <- which(D == 1)
    n0_sub <- length(idx0_sub)
    n1_sub <- length(idx1_sub)
    
    # for class 0
    if (n0_sub > 0){
      Z0 <- Z[idx0_sub, , drop = FALSE]
      I0 <- I[idx0_sub]
      
      if (model == "2LC") {
        S0 <- Z0  # no slope
      } else if (model == "2LCR1") {
        S0 <- Z0 - matrix(I0 * b0, n0_sub, J, byrow = FALSE)
      } else if (model == "2LCR") {
        S0 <- Z0 - (matrix(I0, n0_sub, J, byrow = FALSE) * matrix(b0, n0_sub, J, byrow = TRUE))
      }
      post_var0  <- 1 / (n0_sub + 1 / s2_a0)
      post_mean0 <- post_var0 * (colSums(S0) + mu_a0 / s2_a0)
      a0 <- rnorm(J, post_mean0, sqrt(post_var0))
    } else {
      a0 <- rnorm(J, mu_a0, sqrt(s2_a0))
    }
    
    # for class 1
    if (n1_sub > 0){
      Z1 <- Z[idx1_sub, , drop = FALSE]
      I1 <- I[idx1_sub]
      
      if (model == "2LC") {
        S1 <- Z1
      } else if (model == "2LCR1") {
        S1 <- Z1 - matrix(I1 * b1, n1_sub, J, byrow = FALSE)
      } else {
        S1 <- Z1 - (matrix(I1, n1_sub, J, byrow = FALSE) * matrix(b1, n1_sub, J, byrow = TRUE))
      }
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
        varI0  <- 1 / (sum(b0_vec^2) + 1)
        meanI0 <- varI0 * R0
        I[idx0_sub] <- rnorm(n0_sub, meanI0, sqrt(varI0))
      }
      if (n1_sub > 0) {
        diff1 <- Z[idx1_sub, , drop = FALSE] - matrix(a1, n1_sub, J, byrow = TRUE)
        b1_vec <- if (model == "2LCR1") rep(b1, J) else b1
        R1 <- rowSums(diff1 * matrix(b1_vec, n1_sub, J, byrow = TRUE))
        varI1  <- 1 / (sum(b1_vec^2) + 1)
        meanI1 <- varI1 * R1
        I[idx1_sub] <- rnorm(n1_sub, meanI1, sqrt(varI1))
      }
    }
    
    # 4) Update slopes b0, b1 (2LCR only)
    if (model == "2LCR1") {
      if (n0_sub > 0) {
        I0 <- I[idx0_sub]
        diff0 <- Z[idx0_sub, , drop = FALSE] - matrix(a0, n0_sub, J, byrow = TRUE)
        sum_term <- sum(I0 * rowSums(diff0))
        var_b0  <- 1 / (J * sum(I0^2) + 1 / s2_b0[1])
        mean_b0 <- var_b0 * (sum_term + mu_b0[1] / s2_b0[1])
        b0 <- rnorm(1, mean_b0, sqrt(var_b0))
      } else {
        b0 <- rnorm(1, mu_b0[1], sqrt(s2_b0[1]))
      }
      if (n1_sub > 0) {
        I1 <- I[idx1_sub]
        diff1 <- Z[idx1_sub, , drop = FALSE] - matrix(a1, n1_sub, J, byrow = TRUE)
        sum_term <- sum(I1 * rowSums(diff1))
        var_b1  <- 1 / (J * sum(I1^2) + 1 / s2_b1[1])
        mean_b1 <- var_b1 * (sum_term + mu_b1[1] / s2_b1[1])
        b1 <- rnorm(1, mean_b1, sqrt(var_b1))
      } else {
        b1 <- rnorm(1, mu_b1[1], sqrt(s2_b1[1]))
      }
      
    } else if (model == "2LCR") {
      # Per-test slopes: b0 and b1 are length J
      # Class 0 slopes b0_j
      if (n0_sub > 0) {
        I0 <- I[idx0_sub]
        diff0 <- Z[idx0_sub, , drop = FALSE] - matrix(a0, n0_sub, J, byrow = TRUE)
        sum_term <- sum(I0 * rowSums(diff0))
        var_b0  <- 1 / (J * sum(I0^2) + 1 / s2_b0[1])
        mean_b0 <- var_b0 * (sum_term + mu_b0[1] / s2_b0[1])
        b0 <- rnorm(1, mean_b0, sqrt(var_b0))
      } else {
        b0 <- rnorm(1, mu_b0[1], sqrt(s2_b0[1]))
      }
      if (n1_sub > 0) {
        I1 <- I[idx1_sub]
        diff1 <- Z[idx1_sub, , drop = FALSE] - matrix(a1, n1_sub, J, byrow = TRUE)
        sum_term <- sum(I1 * rowSums(diff1))
        var_b1  <- 1 / (J * sum(I1^2) + 1 / s2_b1[1])
        mean_b1 <- var_b1 * (sum_term + mu_b1[1] / s2_b1[1])
        b1 <- rnorm(1, mean_b1, sqrt(var_b1))
      } else {
        b1 <- rnorm(1, mu_b1[1], sqrt(s2_b1[1]))
      }
    }
    
    # Update D
    if (model == "2LC") {
      eta1 <- matrix(a1, N, J, byrow = TRUE)
      eta0 <- matrix(a0, N, J, byrow = TRUE)
    } else if (model == "2LCR1") {
      eta1 <- matrix(a1, N, J, byrow = TRUE) + matrix(b1 * I, N, J, byrow = FALSE)
      eta0 <- matrix(a0, N, J, byrow = TRUE) + matrix(b0 * I, N, J, byrow = FALSE)
    } else if (model == "2LCR"){
      eta0 <- matrix(a0, N, J, byrow = TRUE) + (I %*% t(b0))
      eta1 <- matrix(a1, N, J, byrow = TRUE) + (I %*% t(b1))
    }
    
    eta1 <- pmin(pmax(eta1, -cap), cap)
    eta0 <- pmin(pmax(eta0, -cap), cap)
    
    ll1_cells <- ifelse(Y == 1, pnorm(eta1, log.p = TRUE),
                        pnorm(eta1, lower.tail = FALSE, log.p = TRUE))
    ll0_cells <- ifelse(Y == 1, pnorm(eta0, log.p = TRUE),
                        pnorm(eta0, lower.tail = FALSE, log.p = TRUE))
    ll1_cells[is.na(Y)] <- 0
    ll0_cells[is.na(Y)] <- 0
    loglik1 <- rowSums(ll1_cells)
    loglik0 <- rowSums(ll0_cells)
    
    # log posterior odds for D=1
    logit_p1 <- (loglik1 - loglik0) + (log(rho) - log1p(-rho))
    p1 <- plogis(logit_p1)
    D  <- rbinom(N, 1, p1)
    
    # prevalence
    rho <- rbeta(1, a_rho + sum(D), b_rho + (N - sum(D)))
    
    # save
    if (it > burnin && ((it - burnin) %% thin == 0)){
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
      SENS[keep, ] <- pnorm( a1 / sqrt(1 + b1^2) )
      SPEC[keep, ] <- pnorm( -a0 / sqrt(1 + b0^2) )
      D_keep[keep, ] <- D
    }
  }
  
  out <-   list(
    rho  = RHO,
    a0   = A0,
    a1   = A1,
    b0   = B0,
    b1   = B1,
    sens = SENS,
    spec = SPEC,
    D    = D_keep)
  
  return(out)
}
