
# LCM based on Dendukuri & Joseph 2001 (probit link)
# -----------------------------------------------------------
# Models:
#   - model = "CI"   : Conditional-independence (CI / Hui–Walter).
#       prevalence  ρ ~ Beta(aρ,bρ)
#       per test j: Sens_j ~ Beta(aSe_j, bSe_j),  Spec_j ~ Beta(aSp_j, bSp_j)
#     Priors: prior_input$prev = c(aρ,bρ)
#             prior_input$tests[[j]]$sens = c(aSe_j,bSe_j)
#             prior_input$tests[[j]]$spec = c(aSp_j,bSp_j)
#     Or provide prior_input$ranges[[j]]$sens/spec = c(L,U) and we’ll convert
#     those ranges to Beta(a,b) via a CI match (default 95%).
#
#   - model = "random"  : Probit random-effects with per-test slopes:
#       P(Y_ij=1 | D_i=0, I_i) = Φ(a0_j + b0_j I_i),  I_i ~ N(0,1)
#       P(Y_ij=1 | D_i=1, I_i) = Φ(a1_j + b1_j I_i)
#     Priors are Normal on a0_j, a1_j, b0_j, b1_j via prior_input$tests[[j]]:
#       a0/a1/b0/b1 = list(mean=..., sd=...).
#     Or provide Se/Sp ranges in prior_input$ranges (95% by default); we turn
#     those into Normal priors on a0/a1, given the slope priors b0/b1 using
#     the identity Se = Φ(a1/sqrt(1+b1^2)), Sp = Φ(-a0/sqrt(1+b0^2)).
#
#   - model = "2LCR1"   : Same as random but common slopes across tests
#       b0_j ≡ b0,  b1_j ≡ b1  (we take slope priors from tests[[1]])
#
# Outputs:
#   - rho: MCMC draws of prevalence
#   - sens, spec: implied Se/Sp per draw (for RE/2LCR1 computed via probit mapping)
#   - a0,a1,b0,b1: probit parameters (NA in CI model)
#   - D: sampled latent disease indicators per iteration
#
# Notes:
#   - Assumes complete data. NAs are tolerated but not masked in updates.
#
library(truncnorm)

bayes_2LCR <- function(data, model = c("CI","random","2LCR1"),
                       prior_input = NULL, common_slopes = FALSE, 
                       iterations = 5000, burnin = 2000, thin = 1){
  
  `%||%` <- function(x, y) if (is.null(x)) y else x
  .clamp01 <- function(x) pmin(pmax(x, 1e-6), 1-1e-6)
  
  # Beta(a,b) whose central CI matches a range [L,U]
  .beta_from_range <- function(range, ci = 0.95) {
    stopifnot(is.numeric(range), length(range)==2)
    L <- .clamp01(min(range)); U <- .clamp01(max(range))
    if (!(L < U)) stop("Range lower must be < upper.")
    pL <- (1 - ci)/2; pU <- 1 - pL
    o <- optim(c(log(2), log(2)),
               function(par){
                 a <- exp(par[1]); b <- exp(par[2])
                 (qbeta(pL, a, b) - L)^2 + (qbeta(pU, a, b) - U)^2
               },
               method="Nelder-Mead",
               control=list(reltol=1e-10, maxit=2000))
    c(a = exp(o$par[1]), b = exp(o$par[2]))
  }
  
  # For CI: build per-test Beta priors from explicit Beta and/or Se/Sp ranges
  .build_CI_tests_priors <- function(J,
                                     tests_beta = NULL,   # list[[j]]$sens/spec = c(a,b)
                                     ranges     = NULL,   # list[[j]]$sens/spec = c(L,U)
                                     range_ci   = 0.95,
                                     ranges_override = TRUE,
                                     default_beta = c(1,1)) {
    out <- vector("list", J)
    for (j in 1:J) out[[j]] <- list(sens = default_beta, spec = default_beta)
    
    # overlay explicit Beta if given
    if (!is.null(tests_beta)) {
      if (!is.list(tests_beta) || length(tests_beta) != J)
        stop(sprintf("'tests' must be a list of length %d.", J))
      for (j in 1:J) {
        tj <- tests_beta[[j]]
        if (!is.null(tj$sens)) {
          stopifnot(is.numeric(tj$sens), length(tj$sens)==2)
          out[[j]]$sens <- as.numeric(tj$sens)
        }
        if (!is.null(tj$spec)) {
          stopifnot(is.numeric(tj$spec), length(tj$spec)==2)
          out[[j]]$spec <- as.numeric(tj$spec)
        }
      }
    }
    
    # overlay ranges (converted to Beta) if provided
    if (!is.null(ranges)) {
      if (!is.list(ranges)) stop("'ranges' must be a list.")
      if (length(ranges) != J) {
        if (length(ranges) == 1L) ranges <- rep(ranges, J)
        else stop(sprintf("'ranges' must have length 1 or %d.", J))
      }
      for (j in 1:J) {
        rj <- ranges[[j]]
        if (!is.null(rj$sens)) {
          ab_se <- .beta_from_range(rj$sens, ci = range_ci)
          if (ranges_override || all(out[[j]]$sens == default_beta))
            out[[j]]$sens <- c(ab_se["a"], ab_se["b"])
        }
        if (!is.null(rj$spec)) {
          ab_sp <- .beta_from_range(rj$spec, ci = range_ci)
          if (ranges_override || all(out[[j]]$spec == default_beta))
            out[[j]]$spec <- c(ab_sp["a"], ab_sp["b"])
        }
      }
    }
    out
  }
  
  # Turn Se/Sp ranges into Normal priors on probit-scale etas
  .eta_norm_from_ranges <- function(sens_range, spec_range, ci = 0.95) {
    z <- qnorm((1+ci)/2)
    Lse <- .clamp01(min(sens_range)); Use <- .clamp01(max(sens_range))
    m_se  <- qnorm((Lse + Use)/2)
    sd_se <- (qnorm(Use) - qnorm(Lse)) / (2*z); sd_se <- max(sd_se, 0.02) # stability floor
    Lsp <- .clamp01(min(spec_range)); Usp <- .clamp01(max(spec_range))
    m_sp  <- qnorm((Lsp + Usp)/2)
    sd_sp <- (qnorm(Usp) - qnorm(Lsp)) / (2*z); sd_sp <- max(sd_sp, 0.02)
    list(m_se=m_se, sd_se=sd_se, m_sp=m_sp, sd_sp=sd_sp)
  }
  
  # For RE/2LCR1: build Normal priors on a0/a1 from Se/Sp ranges, given b priors
  .build_RE_priors_from_ranges <- function(J, tests_norm = NULL, ranges = NULL,
                                           range_ci = 0.95, common_slopes = FALSE,
                                           default_norm = list(mean=0, sd=1),
                                           ndraw = 200000L) {
    # containers initialized to defaults (used if no overrides)
    mu_a0 <- rep(default_norm$mean, J); s2_a0 <- rep(default_norm$sd^2, J)
    mu_a1 <- rep(default_norm$mean, J); s2_a1 <- rep(default_norm$sd^2, J)
    mu_b0 <- rep(default_norm$mean, J); s2_b0 <- rep(default_norm$sd^2, J)
    mu_b1 <- rep(default_norm$mean, J); s2_b1 <- rep(default_norm$sd^2, J)
    
    # overlay explicit normals if provided
    if (!is.null(tests_norm)) {
      stopifnot(is.list(tests_norm), length(tests_norm)==J)
      getm <- function(x, key) { if (is.null(x[[key]])) default_norm$mean else as.numeric(x[[key]]$mean %||% default_norm$mean) }
      gets <- function(x, key) { if (is.null(x[[key]])) default_norm$sd   else as.numeric(x[[key]]$sd   %||% default_norm$sd) }
      for (j in 1:J) {
        tj <- tests_norm[[j]]
        mu_a0[j] <- getm(tj,"a0"); s2_a0[j] <- gets(tj,"a0")^2
        mu_a1[j] <- getm(tj,"a1"); s2_a1[j] <- gets(tj,"a1")^2
        mu_b0[j] <- getm(tj,"b0"); s2_b0[j] <- gets(tj,"b0")^2
        mu_b1[j] <- getm(tj,"b1"); s2_b1[j] <- gets(tj,"b1")^2
      }
    }
    
    # If ranges supplied: compute a0/a1 to match them (marginalizing over b)
    if (!is.null(ranges)) {
      # common slopes: draw slope once, shared across tests
      if (common_slopes) {
        mu_b0_c <- mu_b0[1]; sd_b0_c <- sqrt(s2_b0[1])
        mu_b1_c <- mu_b1[1]; sd_b1_c <- sqrt(s2_b1[1])
        b0_draw <- rnorm(ndraw, mu_b0_c, sd_b0_c); s0_draw <- sqrt(1 + b0_draw^2)
        b1_draw <- rnorm(ndraw, mu_b1_c, sd_b1_c); s1_draw <- sqrt(1 + b1_draw^2)
      }
      if (length(ranges) != J) {
        if (length(ranges) == 1L) ranges <- rep(ranges, J)
        else stop(sprintf("'ranges' must have length 1 or %d.", J))
      }
      for (j in 1:J) {
        if (!is.null(ranges[[j]])) {
          er <- .eta_norm_from_ranges(ranges[[j]]$sens, ranges[[j]]$spec, ci = range_ci)
          if (!common_slopes) {
            b0_draw <- rnorm(ndraw, mu_b0[j], sqrt(s2_b0[j])); s0_draw <- sqrt(1+b0_draw^2)
            b1_draw <- rnorm(ndraw, mu_b1[j], sqrt(s2_b1[j])); s1_draw <- sqrt(1+b1_draw^2)
          }
          eta_se <- rnorm(ndraw, er$m_se, er$sd_se)
          eta_sp <- rnorm(ndraw, er$m_sp, er$sd_sp)
          a1_draw <- eta_se * s1_draw
          a0_draw <- - eta_sp * s0_draw
          mu_a1[j] <- mean(a1_draw); s2_a1[j] <- var(a1_draw)
          mu_a0[j] <- mean(a0_draw); s2_a0[j] <- var(a0_draw)
        }
      }
    }
    
    list(mu_a0=mu_a0, s2_a0=s2_a0,
         mu_a1=mu_a1, s2_a1=s2_a1,
         mu_b0=mu_b0, s2_b0=s2_b0,
         mu_b1=mu_b1, s2_b1=s2_b1)
  }
  
  # Fill defaults for explicit CI priors if needed
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
  
  # Defaults for RE normals if user provides none
  .fill_default_tests_random <- function(J) {
    one <- list(
      a0 = list(mean = 0, sd = 1),
      a1 = list(mean = 0, sd = 1),
      b0 = list(mean = 0, sd = 1),
      b1 = list(mean = 0, sd = 1)
    )
    rep(list(one), J)
  }
  
  # Parse priors for RE/2LCR1, with optional ranges->(a0,a1) mapping
  parse_dj_priors <- function(prior_input, J, model, common_slopes) {
    # prevalence prior
    if (is.null(prior_input$prev)) prior_input$prev <- c(1,1)
    if (!(is.numeric(prior_input$prev) && length(prior_input$prev)==2 &&
          all(is.finite(prior_input$prev)) && all(prior_input$prev>0))) {
      stop("prior_input$prev must be c(a,b) with positive, finite entries.")
    }
    out <- list(a_rho = prior_input$prev[1], b_rho = prior_input$prev[2])
    if (model == "CI") return(out)
    
    if (!is.null(prior_input$ranges)) {
      re_from_ranges <- .build_RE_priors_from_ranges(
        J,
        tests_norm = prior_input$tests %||% NULL,     # use any user normals for b0/b1
        ranges     = prior_input$ranges,
        range_ci   = prior_input$range_ci %||% 0.95,
        common_slopes = common_slopes,
        default_norm = list(mean = 0, sd = 1),
        ndraw = 200000L
      )
      mu_a0 <- re_from_ranges$mu_a0; s2_a0 <- re_from_ranges$s2_a0
      mu_a1 <- re_from_ranges$mu_a1; s2_a1 <- re_from_ranges$s2_a1
      # keep mu_b0/s2_b0, mu_b1/s2_b1 from explicit tests_norm or defaults
    }
    
    
    # defaults for normals
    tests_list <- prior_input$tests
    if (is.null(tests_list)) tests_list <- .fill_default_tests_random(J)
    if (!is.list(tests_list) || length(tests_list)!=J)
      stop("prior_input$tests must be a list of length J for random/2LCR1.")
    
    # pull b priors (means/vars)
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
    
    # a priors (will be overridden by ranges if provided)
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
    
    # allow scalar/vector overrides on mu_/s2_ if present in prior_input
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
    
    # If ranges are provided: override a0/a1 using mapping via b priors
    if (!is.null(prior_input$ranges)) {
      re <- .build_RE_priors_from_ranges(
        J = J,
        tests_norm   = tests_list,                    # read b0/b1 means/sds if set
        ranges       = prior_input$ranges,
        range_ci     = prior_input$range_ci %||% 0.95,
        common_slopes= common_slopes,
        default_norm = list(mean=0, sd=1),
        ndraw        = prior_input$range_ndraw %||% 200000L
      )
      mu_a0 <- re$mu_a0; s2_a0 <- re$s2_a0
      mu_a1 <- re$mu_a1; s2_a1 <- re$s2_a1
      # keep mu_b*/s2_b* as above (from tests_list or overrides)
    }
    
    list(a_rho=out$a_rho, b_rho=out$b_rho,
         mu_a0=mu_a0, s2_a0=s2_a0,
         mu_a1=mu_a1, s2_a1=s2_a1,
         mu_b0=mu_b0, s2_b0=s2_b0,
         mu_b1=mu_b1, s2_b1=s2_b1)
  }
  
  Y <- as.matrix(data)
  N <- nrow(Y); J <- ncol(Y)
  model <- match.arg(model)
  
  if (is.null(prior_input)) prior_input <- list()
  
  # Storage
  n_keep <- floor((iterations)/thin)
  RHO <- numeric(n_keep)
  A0 <- matrix(NA, n_keep, J)
  A1 <- matrix(NA, n_keep, J)
  B0 <- matrix(NA, n_keep, J)
  B1 <- matrix(NA, n_keep, J)
  SENS <- matrix(NA, n_keep, J)
  SPEC <- matrix(NA, n_keep, J)
  D_keep <- matrix(NA, n_keep, N)
  keep <- 0
  
  # ========== CI model ==========
  if (model == "CI") {
    # prevalence prior
    if (is.null(prior_input$prev)) prior_input$prev <- c(1,1)
    a_rho <- prior_input$prev[1]; b_rho <- prior_input$prev[2]
    
    # build tests Beta priors from explicit Beta and/or ranges
    if (!is.null(prior_input$ranges)) {
      prior_input$tests <- .build_CI_tests_priors(
        J,
        tests_beta = prior_input$tests %||% NULL,
        ranges     = prior_input$ranges,
        range_ci   = prior_input$range_ci %||% 0.95,
        ranges_override = TRUE
      )
    }
    
    tests_beta <- prior_input$tests
    tests_rng  <- prior_input$ranges
    range_ci   <- prior_input$range_ci %||% 0.95
    if (is.null(tests_beta) && is.null(tests_rng)) {
      tests_beta <- .fill_default_tests_CI(J)
    } else if (!is.null(tests_rng)) {
      tests_beta <- .build_CI_tests_priors(
        J, tests_beta = tests_beta, ranges = tests_rng, range_ci = range_ci,
        ranges_override = TRUE, default_beta = c(1,1)
      )
    } else {
      if (length(tests_beta) != J) stop("prior_input$tests must be length J.")
    }
    
    rho <- rbeta(1, a_rho, b_rho)
    sens <- numeric(J); spec <- numeric(J)
    for (j in 1:J) {
      pj <- tests_beta[[j]]
      sens[j] <- rbeta(1, pj$sens[1], pj$sens[2])
      spec[j] <- rbeta(1, pj$spec[1], pj$spec[2])
    }
    D <- rbinom(N, 1, rho)
    
    for (it in 1:(iterations + burnin)) {
      # D | params
      log_pY1 <- rep(0, N)
      log_pY0 <- rep(0, N)
      for (j in 1:J) {
        yj <- Y[,j]
        prob1 <- rep(sens[j], N)        # P(Y=1|D=1)
        prob0 <- rep(1 - spec[j], N)    # P(Y=1|D=0)
        lj1 <- dbinom(yj, 1, prob1, log=TRUE)
        lj0 <- dbinom(yj, 1, prob0, log=TRUE)
        lj1[is.na(lj1)] <- 0; lj0[is.na(lj0)] <- 0
        log_pY1 <- log_pY1 + lj1
        log_pY0 <- log_pY0 + lj0
      }
      logit_pi <- (log_pY1 - log_pY0) + (log(rho) - log1p(-rho))
      D <- rbinom(N, 1, plogis(logit_pi))
      
      # prevalence
      rho <- rbeta(1, a_rho + sum(D), b_rho + (N - sum(D)))
      
      # Se/Sp per test (conjugate updates)
      for (j in 1:J) {
        pj <- tests_beta[[j]]
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
        RHO[keep]  <- rho
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
  
  # ========== RANDOM / 2LCR1 ==========
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
      b  <- rnorm(1, mean = mu_b0[1], sd = sqrt(s2_b0[1]))  # you could average b0/b1 if desired
    } else {
      # one slope per class, shared across tests
      b0 <- rnorm(1, mu_b0[1], sqrt(s2_b0[1]))
      b1 <- rnorm(1, mu_b1[1], sqrt(s2_b1[1]))
    }
  } else {
    if (common_slopes) {
      # one slope per test, shared across classes
      b  <- rnorm(J, mu_b0, sqrt(s2_b0))  # (or average mu_b0/mu_b1)
    } else {
      # per-test, per-class
      b0 <- rnorm(J, mu_b0, sqrt(s2_b0))
      b1 <- rnorm(J, mu_b1, sqrt(s2_b1))
    }
  }
  I  <- rnorm(N, 0, 1)
  
  for (it in 1:(iterations + burnin)) {
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
    } else { 
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
    
    # 2) update a0, a1
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
    
    # 3) update I
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
      var_b <- 1 / (J * sumI2 + 1 / s2_b0[1])   # using b0 prior
      b <- rnorm(1, var_b * (sum_term + mu_b0[1]/s2_b0[1]), sqrt(var_b))
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
    } else if(model == "random" && common_slopes){
      II <- sum(I^2)
      diff_all <- Z - ifelse(matrix(D, N, J, byrow=FALSE)==1,
                             matrix(a1, N, J, byrow=TRUE),
                             matrix(a0, N, J, byrow=TRUE))
      IR <- colSums(matrix(I, N, J, byrow=FALSE) * diff_all)
      post_var_b  <- 1 / (II + 1 / s2_b0)
      post_mean_b <- post_var_b * (IR + mu_b0 / s2_b0)
      b <- rnorm(J, post_mean_b, sqrt(post_var_b))
    }
    
    # 5) update D
    if (model == "2LCR1" && common_slopes) {
      eta0 <- matrix(a0, N, J, byrow=TRUE) + matrix(b * I, N, J)
      eta1 <- matrix(a1, N, J, byrow=TRUE) + matrix(b * I, N, J)
    } else if (model == "2LCR1" && !common_slopes) {
      eta0 <- matrix(a0, N, J, byrow=TRUE) + matrix(b0 * I, N, J)
      eta1 <- matrix(a1, N, J, byrow=TRUE) + matrix(b1 * I, N, J)
    } else if (model != "2LCR1" && common_slopes) {
      eta0 <- matrix(a0, N, J, byrow=TRUE) + tcrossprod(I, b)
      eta1 <- matrix(a1, N, J, byrow=TRUE) + tcrossprod(I, b)
    } else { 
      eta0 <- matrix(a0, N, J, byrow=TRUE) + tcrossprod(I, b0)
      eta1 <- matrix(a1, N, J, byrow=TRUE) + tcrossprod(I, b1)
    }
    eta0 <- pmin(pmax(eta0, -cap), cap)
    eta1 <- pmin(pmax(eta1, -cap), cap)
    ll1 <- ifelse(Y == 1, pnorm(eta1, log.p=TRUE),
                  pnorm(eta1, lower.tail=FALSE, log.p=TRUE))
    ll0 <- ifelse(Y == 1, pnorm(eta0, log.p=TRUE),
                  pnorm(eta0, lower.tail=FALSE, log.p=TRUE))
    ll1[is.na(ll1)] <- 0; ll0[is.na(ll0)] <- 0
    loglik1 <- rowSums(ll1); loglik0 <- rowSums(ll0)
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
