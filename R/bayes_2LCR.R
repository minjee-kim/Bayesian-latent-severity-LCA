

# LCM based on Dendukuri & Joseph 2001 (probit link)
# -----------------------------------------------------------
# Models:
#   - model = "CI"   : Conditional-independence (CI / Hui–Walter). Parameters are:
#       prevalence  ρ ~ Beta(aρ,bρ)
#       per test j: Sens_j ~ Beta(aSe_j, bSe_j),  Spec_j ~ Beta(aSp_j, bSp_j)
#     Priors are provided via prior_input$prev = c(aρ,bρ) and
#     prior_input$tests[[j]]$sens = c(aSe_j,bSe_j), prior_input$tests[[j]]$spec = c(aSp_j,bSp_j).
#     DEFAULT if omitted: Beta(1,1) for all (ρ, Se_j, Sp_j).
#
#   - model = "random"  : Probit random-effects with per-test slopes:
#       P(Y_ij=1 | D_i=0, I_i) = Φ(a0_j + b0_j I_i),  I_i ~ N(0,1)
#       P(Y_ij=1 | D_i=1, I_i) = Φ(a1_j + b1_j I_i)
#     Priors are Normal on a0_j, a1_j, b0_j, b1_j with user-specified means/sds
#     via prior_input$tests[[j]]$a0/a1/b0/b1 = list(mean=..., sd=...).
#     (Equivalently you can pass scalars or length-J vectors:
#        prior_input$mu_a0, s2_a0, mu_a1, s2_a1, mu_b0, s2_b0, mu_b1, s2_b1.)
#     DEFAULT if omitted: a0_j,a1_j,b0_j,b1_j ~ Normal(mean=0, sd=1).
#
#   - model = "2LCR1"   : Same as random but **common slopes** across tests
#       b0_j ≡ b0,  b1_j ≡ b1
#     For priors, use the first element (or supply tests[[1]]); defaults mean=0, sd=1.
#
# Identifiability guidance (rules of thumb):
#   - CI model needs p > 3 tests; 2LCR1 needs p > 4 tests.
#
# Outputs:
#   - rho: MCMC draws of prevalence
#   - sens, spec: per-draw implied Se/Sp (for random/2LCR1 computed via probit mapping)
#   - a0,a1,b0,b1: probit parameters (NA in CI model)
#   - D: sampled latent disease indicators per iteration
#
# Notes:
#   - This version does NOT use matrix masking; it assumes complete data. If NAs are present,
#     their contribution is safely zeroed in the likelihood, but no per-cell masking is used
#     in parameter updates (use complete data for best results).
#
library(truncnorm)

bayes_2LCR <- function(data, model = c("CI","random","2LCR1"),
                       prior_input = NULL, common_slopes = FALSE, 
                       iterations = 5000, burnin = 2000, thin = 1){
  
  
  # Helper: fill defaults for priors
  .fill_default_tests_CI <- function(J, tests = NULL) {
    out <- vector("list", J)
    for (j in 1:J) {
      tj <- if (!is.null(tests) && length(tests) >= j) tests[[j]] else NULL
      se <- if (!is.null(tj) && !is.null(tj$sens) && length(tj$sens) == 2) tj$sens else c(1,1)
      sp <- if (!is.null(tj) && !is.null(tj$spec) && length(tj$spec) == 2) tj$spec else c(1,1)
      out[[j]] <- list(sens = se, spec = sp)
    }
    out
  }
  
  .fill_default_tests_random <- function(J) {
    one <- list(
      a0 = list(mean = 0, sd = 1),
      a1 = list(mean = 0, sd = 1),
      b0 = list(mean = 0, sd = 1),
      b1 = list(mean = 0, sd = 1)
    )
    rep(list(one), J)
  }
  
  # Parse per-test normals into mu/s2 vectors
  parse_dj_priors <- function(prior_input, J, model, common_slopes) {
    # prevalence prior
    if (is.null(prior_input$prev)) prior_input$prev <- c(1,1)
    if (!(is.numeric(prior_input$prev) && length(prior_input$prev)==2 &&
          all(is.finite(prior_input$prev)) && all(prior_input$prev>0))) {
      stop("prior_input$prev must be c(a,b) with positive, finite entries.")
    }
    out <- list(a_rho = prior_input$prev[1], b_rho = prior_input$prev[2])
    
    if (model == "CI") return(out)
    
    # defaults
    tests_list <- prior_input$tests
    if (is.null(tests_list)) tests_list <- .fill_default_tests_random(J)
    if (!is.list(tests_list) || length(tests_list)!=J)
      stop("prior_input$tests must be a list of length J for random/2LCR1.")
    
    # helper
    coerce_to_list <- function(x) if (is.data.frame(x)) as.list(x) else x
    mv <- function(L, lab, j) {
      L <- coerce_to_list(L)
      if (is.null(L$mean) || is.null(L$sd))
        stop(sprintf("tests[[%d]]$%s must have mean and sd.", j, lab))
      m <- as.numeric(L$mean); s <- as.numeric(L$sd)
      if (!is.finite(m) || !is.finite(s) || s<0)
        stop(sprintf("Bad %s prior at test %d.", lab, j))
      c(m, s^2)
    }
    
    # pull b priors first (needed if mapping se/sp -> a0/a1)
    mu_b0 <- mu_b1 <- numeric(J)
    s2_b0 <- s2_b1 <- numeric(J)
    
    for (j in 1:J) {
      pj <- coerce_to_list(tests_list[[j]])
      if (!is.null(pj$b) && is.null(pj$b0) && is.null(pj$b1)) { pj$b0 <- pj$b; pj$b1 <- pj$b }
      for (req in c("b0","b1")) if (is.null(pj[[req]]))
        stop(sprintf("tests[[%d]] missing '%s'", j, req))
      bb0 <- mv(pj$b0,"b0",j); bb1 <- mv(pj$b1,"b1",j)
      mu_b0[j] <- bb0[1]; s2_b0[j] <- bb0[2]
      mu_b1[j] <- bb1[1]; s2_b1[j] <- bb1[2]
    }
    
    mu_a0 <- mu_a1 <- numeric(J)
    s2_a0 <- s2_a1 <- numeric(J)
    
    for (j in 1:J) {
      pj <- coerce_to_list(tests_list[[j]])
      for (req in c("a0","a1")) if (is.null(pj[[req]]))
        stop(sprintf("tests[[%d]] missing '%s'", j, req))
      aa0 <- mv(pj$a0,"a0",j); aa1 <- mv(pj$a1,"a1",j)
      mu_a0[j] <- aa0[1]; s2_a0[j] <- aa0[2]
      mu_a1[j] <- aa1[1]; s2_a1[j] <- aa1[2]
    }
    
    # allow vector overrides if provided (mu_*/s2_*)
    vJ <- function(x, key) {
      if (length(x)==1) rep(x, J) else if (length(x)==J) x
      else stop(sprintf("%s must have length 1 or J.", key))
    }
    apply_vec <- function(name_mu, name_s2, cur_mu, cur_s2) {
      if (!is.null(prior_input[[name_mu]])) {
        cur_mu <- vJ(as.numeric(prior_input[[name_mu]]), name_mu)
        if (!all(is.finite(cur_mu))) stop(name_mu)
      }
      if (!is.null(prior_input[[name_s2]])) {
        cur_s2 <- vJ(as.numeric(prior_input[[name_s2]]), name_s2)
        if (!all(is.finite(cur_s2)) || any(cur_s2<=0)) stop(name_s2)
      }
      list(mu=cur_mu, s2=cur_s2)
    }
    tmp <- apply_vec("mu_a0","s2_a0", mu_a0, s2_a0); mu_a0 <- tmp$mu; s2_a0 <- tmp$s2
    tmp <- apply_vec("mu_a1","s2_a1", mu_a1, s2_a1); mu_a1 <- tmp$mu; s2_a1 <- tmp$s2
    tmp <- apply_vec("mu_b0","s2_b0", mu_b0, s2_b0); mu_b0 <- tmp$mu; s2_b0 <- tmp$s2
    tmp <- apply_vec("mu_b1","s2_b1", mu_b1, s2_b1); mu_b1 <- tmp$mu; s2_b1 <- tmp$s2
    
    list(a_rho=out$a_rho, b_rho=out$b_rho,
         mu_a0=mu_a0, s2_a0=s2_a0,
         mu_a1=mu_a1, s2_a1=s2_a1,
         mu_b0=mu_b0, s2_b0=s2_b0,
         mu_b1=mu_b1, s2_b1=s2_b1)
  }
  
  
  
  Y <- as.matrix(data)
  N <- nrow(Y); J <- ncol(Y)
  model <- match.arg(model)
  
  # Defaults container if user passes NULL
  if (is.null(prior_input)) prior_input <- list()
  
  vJ <- function(x) if (length(x) == 1) rep(x, J) else x
  
  # Storage
  n_keep <- floor((iterations - burnin)/thin)
  RHO <- numeric(n_keep)
  A0 <- matrix(NA, n_keep, J)
  A1 <- matrix(NA, n_keep, J)
  B0 <- matrix(NA, n_keep, J)
  B1 <- matrix(NA, n_keep, J)
  SENS <- matrix(NA, n_keep, J)
  SPEC <- matrix(NA, n_keep, J)
  D_keep <- matrix(NA, n_keep, N)
  
  keep <- 0
  
  # =========================
  # (A) Condiitonal Independence MODEL
  # =========================
  if (model == "CI") {
    if (is.null(prior_input$prev)) prior_input$prev <- c(1,1)
    a_rho <- prior_input$prev[1]; b_rho <- prior_input$prev[2]
    
    if (is.null(prior_input$tests)) prior_input$tests <- .fill_default_tests_CI(J)
    if (length(prior_input$tests) != J)
      stop("prior_input$tests must be length J (or omit to use defaults).")
    
    rho <- rbeta(1, a_rho, b_rho)
    sens <- numeric(J); spec <- numeric(J)
    for (j in 1:J) {
      pj <- prior_input$tests[[j]]
      if (length(pj$sens)!=2 || length(pj$spec)!=2)
        stop(sprintf("tests[[%d]] must have sens/spec as c(a,b).", j))
      sens[j] <- rbeta(1, pj$sens[1], pj$sens[2])
      spec[j] <- rbeta(1, pj$spec[1], pj$spec[2])
    }
    D <- rbinom(N, 1, rho)
    
    for (it in 1:iterations) {
      # D | params
      log_pY1 <- rep(0, N)
      log_pY0 <- rep(0, N)
      for (j in 1:J) {
        yj <- Y[,j]
        prob1 <- rep(sens[j], N)        # P(Y=1|D=1)
        prob0 <- rep(1 - spec[j], N)    # P(Y=1|D=0)
        lj1 <- dbinom(yj, 1, prob1, log=TRUE)
        lj0 <- dbinom(yj, 1, prob0, log=TRUE)
        # ignore accidental NA
        lj1[is.na(lj1)] <- 0; lj0[is.na(lj0)] <- 0
        log_pY1 <- log_pY1 + lj1
        log_pY0 <- log_pY0 + lj0
      }
      logit_pi <- (log_pY1 - log_pY0) + (log(rho) - log1p(-rho))
      D <- rbinom(N, 1, plogis(logit_pi))
      
      # prevalence
      rho <- rbeta(1, prior_input$prev[1] + sum(D),
                   prior_input$prev[2] + (N - sum(D)))
      
      # Se/Sp per test
      for (j in 1:J) {
        pj <- prior_input$tests[[j]]
        yj <- Y[,j]; Dj <- D
        tp <- sum(yj==1 & Dj==1, na.rm=TRUE)
        fn <- sum(yj==0 & Dj==1, na.rm=TRUE)
        tn <- sum(yj==0 & Dj==0, na.rm=TRUE)
        fp <- sum(yj==1 & Dj==0, na.rm=TRUE)
        sens[j] <- rbeta(1, pj$sens[1] + tp, pj$sens[2] + fn)
        spec[j] <- rbeta(1, pj$spec[1] + tn, pj$spec[2] + fp)
      }
      
      if (it > burnin && ((it - burnin) %% thin == 0)) {
        keep <- keep + 1
        RHO[keep] <- rho
        SENS[keep,] <- sens
        SPEC[keep,] <- spec
        D_keep[keep,] <- D
      }
    }
    
    return(list(
      rho = RHO,
      sens = SENS, spec = SPEC,
      a0 = A0, a1 = A1, b0 = B0, b1 = B1,  # NA for CI
      D = D_keep, model = model
    ))
  }
  
  # =========================
  # (B) RANDOM / 2LCR1
  # =========================
  pri <- parse_dj_priors(prior_input, J, model, common_slopes)
  a_rho <- pri$a_rho; b_rho <- pri$b_rho
  mu_a0 <- pri$mu_a0; s2_a0 <- pri$s2_a0
  mu_a1 <- pri$mu_a1; s2_a1 <- pri$s2_a1
  mu_b0 <- pri$mu_b0; s2_b0 <- pri$s2_b0
  mu_b1 <- pri$mu_b1; s2_b1 <- pri$s2_b1
  
  cap <- 8
  rho <- rbeta(1, a_rho, b_rho)
  D   <- rbinom(N, 1, rho)
  a0  <- rnorm(J, mu_a0, sqrt(s2_a0))
  a1  <- rnorm(J, mu_a1, sqrt(s2_a1))
  
  if (model == "2LCR1") {
    if (common_slopes) {
      # single slope shared across tests AND classes
      b  <- rnorm(1, mean = pri$mu_b0[1], sd = sqrt(pri$s2_b0[1]))  # pick b0 prior (or average of b0/b1)
    } else {
      # one slope per class, shared across tests (your current 2LCR1)
      b0 <- rnorm(1, pri$mu_b0[1], sqrt(pri$s2_b0[1]))
      b1 <- rnorm(1, pri$mu_b1[1], sqrt(pri$s2_b1[1]))
    }
  } else {
    if (common_slopes) {
      # one slope per test, shared across classes
      b  <- rnorm(J, pri$mu_b0, sqrt(pri$s2_b0))  # use b0 prior (or average)
    } else {
      # DJ-general: per-test, per-class
      b0 <- rnorm(J, pri$mu_b0, sqrt(pri$s2_b0))
      b1 <- rnorm(J, pri$mu_b1, sqrt(pri$s2_b1))
    }
  }
  I  <- rnorm(N, 0, 1)
  
  for (it in 1:iterations) {
    # 1) latent Z
    if (model == "2LCR1" && common_slopes) {
      eta0 <- matrix(a0, N, J, byrow=TRUE) + matrix(b * I, N, J)
      eta1 <- matrix(a1, N, J, byrow=TRUE) + matrix(b * I, N, J)
    } else if (model == "2LCR1" && !common_slopes) {
      eta0 <- matrix(a0, N, J, byrow=TRUE) + matrix(b0 * I, N, J)
      eta1 <- matrix(a1, N, J, byrow=TRUE) + matrix(b1 * I, N, J)
    } else if (model != "2LCR1" && common_slopes) {
      eta0 <- matrix(a0, N, J, byrow=TRUE) + tcrossprod(I, b)
      eta1 <- matrix(a1, N, J, byrow=TRUE) + tcrossprod(I, b)
    } else { # DJ-general
      eta0 <- matrix(a0, N, J, byrow=TRUE) + tcrossprod(I, b0)
      eta1 <- matrix(a1, N, J, byrow=TRUE) + tcrossprod(I, b1)
    }
    eta0 <- pmin(pmax(eta0, -cap), cap)
    eta1 <- pmin(pmax(eta1, -cap), cap)
    
    MU <- eta0
    if (any(D==1)) MU[D==1,] <- eta1[D==1,]
    
    Z <- matrix(0, N, J)
    m1 <- (Y == 1); m1[is.na(m1)] <- FALSE
    m0 <- (Y == 0); m0[is.na(m0)] <- FALSE
    s1 <- sum(m1); if (s1 > 0) Z[m1] <- rtruncnorm(s1, a=0,    b=Inf, mean=MU[m1], sd=1)
    s0 <- sum(m0); if (s0 > 0) Z[m0] <- rtruncnorm(s0, a=-Inf, b=0,   mean=MU[m0], sd=1)
    
    # 2) update a0, a1  (complete-data formulas; per-test n = n0/n1)
    idx0 <- which(D==0); n0 <- length(idx0)
    idx1 <- which(D==1); n1 <- length(idx1)
    
    if (n0 > 0) {
      Z0 <- Z[idx0,, drop=FALSE]; I0 <- I[idx0]
      S0 <- if (model=="2LCR1" && common_slopes) {
        Z0 - matrix(I0 * b,  n0, J, byrow=FALSE)
      } else if (model=="2LCR1") {
        Z0 - matrix(I0 * b0, n0, J, byrow=FALSE)
      } else if (common_slopes) {
        Z0 - (matrix(I0, n0, J, byrow=FALSE) * matrix(b,  n0, J, byrow=TRUE))
      } else {
        Z0 - (matrix(I0, n0, J, byrow=FALSE) * matrix(b0, n0, J, byrow=TRUE))
      }
      post_var0  <- 1 / (n0 + 1 / s2_a0)
      post_mean0 <- post_var0 * (colSums(S0) + mu_a0 / s2_a0)
      a0 <- rnorm(J, post_mean0, sqrt(post_var0))
    } else {
      a0 <- rnorm(J, mu_a0, sqrt(s2_a0))
    }
    
    if (n1 > 0) {
      Z1 <- Z[idx1,, drop=FALSE]; I1 <- I[idx1]
      S1 <- if (model=="2LCR1" && common_slopes) {
        Z1 - matrix(I1 * b,  n1, J, byrow=FALSE)
      } else if (model=="2LCR1") {
        Z1 - matrix(I1 * b1, n1, J, byrow=FALSE)
      } else if (common_slopes) {
        Z1 - (matrix(I1, n1, J, byrow=FALSE) * matrix(b,  n1, J, byrow=TRUE))
      } else {
        Z1 - (matrix(I1, n1, J, byrow=FALSE) * matrix(b1, n1, J, byrow=TRUE))
      }
      post_var1  <- 1 / (n1 + 1 / s2_a1)
      post_mean1 <- post_var1 * (colSums(S1) + mu_a1 / s2_a1)
    
      a1 <- rnorm(J, mean = post_mean1, sd = sqrt(post_var1))
    } else {
      a1 <- rnorm(J, mean = mu_a1, sd = sqrt(s2_a1))
    }
    
    # 3) update I (complete-data formulas)
    if (n0 > 0) {
      diff0 <- Z[idx0,,drop=FALSE] - matrix(a0, n0, J, byrow=TRUE)
      b0v <- if (model=="2LCR1") {
        if (common_slopes) rep(b, J) else rep(b0, J)
      } else {
        if (common_slopes) b else b0
      }
      R0 <- rowSums(diff0 * matrix(b0v, n0, J, byrow=TRUE))
      varI0 <- 1 / (sum(b0v^2) + 1)
      I[idx0] <- rnorm(n0, varI0 * R0, sqrt(varI0))
    }
    if (n1 > 0) {
      diff1 <- Z[idx1,,drop=FALSE] - matrix(a1, n1, J, byrow=TRUE)
      b1v <- if (model=="2LCR1") {
        if (common_slopes) rep(b, J) else rep(b1, J)
      } else {
        if (common_slopes) b else b1
      }
      R1 <- rowSums(diff1 * matrix(b1v, n1, J, byrow=TRUE))
      varI1 <- 1 / (sum(b1v^2) + 1)
      I[idx1] <- rnorm(n1, varI1 * R1, sqrt(varI1))
    }
    
    # 4) update slopes
    if (model=="2LCR1" && common_slopes) {
      diff_all <- Z - ifelse(matrix(D, N, J, byrow=FALSE)==1,
                             matrix(a1, N, J, byrow=TRUE),
                             matrix(a0, N, J, byrow=TRUE))
      sumI2 <- sum(I^2)
      sum_term <- sum(I * rowSums(diff_all))
      var_b <- 1 / (J * sumI2 + 1 / pri$s2_b0[1])   # use b0 prior (or average b0/b1)
      b <- rnorm(1, var_b * (sum_term + pri$mu_b0[1]/pri$s2_b0[1]), sqrt(var_b))
    }
    if (model == "2LCR1" && common_slopes == FALSE) {
      if (n0 > 0) {
        I0 <- I[idx0]
        diff0 <- Z[idx0,, drop=FALSE] - matrix(a0, n0, J, byrow=TRUE)
        sum_term0 <- sum(I0 * rowSums(diff0))
        var_b0 <- 1 / (J * sum(I0^2) + 1 / s2_b0[1])
        b0 <- rnorm(1, var_b0 * (sum_term0 + mu_b0[1]/s2_b0[1]), sqrt(var_b0))
      } else b0 <- rnorm(1, mu_b0[1], sqrt(s2_b0[1]))
      
      if (n1 > 0) {
        I1 <- I[idx1]
        diff1 <- Z[idx1,, drop=FALSE] - matrix(a1, n1, J, byrow=TRUE)
        sum_term1 <- sum(I1 * rowSums(diff1))
        var_b1 <- 1 / (J * sum(I1^2) + 1 / s2_b1[1])
        b1 <- rnorm(1, var_b1 * (sum_term1 + mu_b1[1]/s2_b1[1]), sqrt(var_b1))
      } else b1 <- rnorm(1, mu_b1[1], sqrt(s2_b1[1]))
    } else if(model == "random" && common_slopes == FALSE) {
      if (n0 > 0) {
        I0 <- I[idx0]; I0I0 <- sum(I0^2)
        diff0 <- Z[idx0,, drop=FALSE] - matrix(a0, n0, J, byrow=TRUE)
        IR0 <- matrix(I0, n0, J, byrow=FALSE) * diff0
        post_var_b0  <- 1 / (I0I0 + 1 / s2_b0)
        post_mean_b0 <- post_var_b0 * (colSums(IR0) + mu_b0 / s2_b0)
        b0 <- rnorm(J, post_mean_b0, sqrt(post_var_b0))
      } else b0 <- rnorm(J, mu_b0, sqrt(s2_b0))
      
      if (n1 > 0) {
        I1 <- I[idx1]; I1I1 <- sum(I1^2)
        diff1 <- Z[idx1,, drop=FALSE] - matrix(a1, n1, J, byrow=TRUE)
        IR1 <- matrix(I1, n1, J, byrow=FALSE) * diff1
        post_var_b1  <- 1 / (I1I1 + 1 / s2_b1)
        post_mean_b1 <- post_var_b1 * (colSums(IR1) + mu_b1 / s2_b1)
        b1 <- rnorm(J, post_mean_b1, sqrt(post_var_b1))
      } else b1 <- rnorm(J, mu_b1, sqrt(s2_b1))
    }else if(model == "random" && common_slopes){
      II <- sum(I^2)
      # pooled regression per test
      diff_all <- Z - ifelse(matrix(D, N, J, byrow=FALSE)==1,
                             matrix(a1, N, J, byrow=TRUE),
                             matrix(a0, N, J, byrow=TRUE))
      IR <- colSums(matrix(I, N, J, byrow=FALSE) * diff_all)
      post_var_b  <- 1 / (II + 1 / pri$s2_b0)
      post_mean_b <- post_var_b * (IR + pri$mu_b0 / pri$s2_b0)
      b <- rnorm(J, post_mean_b, sqrt(post_var_b))
    }
    
    # 5) update D
    ll1 <- ifelse(Y == 1, pnorm(eta1, log.p=TRUE),
                  pnorm(eta1, lower.tail=FALSE, log.p=TRUE))
    ll0 <- ifelse(Y == 1, pnorm(eta0, log.p=TRUE),
                  pnorm(eta0, lower.tail=FALSE, log.p=TRUE))
    ll1[is.na(ll1)] <- 0
    ll0[is.na(ll0)] <- 0
    
    loglik1 <- rowSums(ll1)
    loglik0 <- rowSums(ll0)
    logit_p1 <- (loglik1 - loglik0) + (log(rho) - log1p(-rho))
    D <- rbinom(N, 1, plogis(logit_p1))
    
    # 6) prevalence
    rho <- rbeta(1, a_rho + sum(D), b_rho + (N - sum(D)))
    
    # save
    if (it > burnin && ((it - burnin) %% thin == 0)) {
      keep <- keep + 1
      RHO[keep] <- rho
      A0[keep,] <- a0
      A1[keep,] <- a1
      if (model=="2LCR1" && common_slopes) {
        B0[keep,] <- rep(b, J)
        B1[keep,] <- rep(b, J)
      } else if (model=="2LCR1" && !common_slopes) {
        B0[keep,] <- rep(b0, J)
        B1[keep,] <- rep(b1, J)
      } else if (model!="2LCR1" && common_slopes) {
        B0[keep,] <- b
        B1[keep,] <- b
      } else { # DJ-general
        B0[keep,] <- b0
        B1[keep,] <- b1
      }
      SENS[keep,] <- pnorm(A1[keep,] / sqrt(1 + B1[keep,]^2))
      SPEC[keep,] <- pnorm(-A0[keep,] / sqrt(1 + B0[keep,]^2))
      D_keep[keep,] <- D
    }
  }
  
  list(
    rho = RHO,
    a0  = A0, a1 = A1,
    b0  = B0, b1 = B1,
    sens = SENS, spec = SPEC,
    D = D_keep,
    model = model
  )
}
