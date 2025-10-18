


plot_prior_vs_posterior_2LCR <- function(
    fit,
    prior_input = NULL,                  # pass the same object used to fit
    model = NULL,                        # if NULL, uses fit$model
    common_slopes = FALSE,               # must match the fit
    nsim = 20000                         # number of prior predictive draws
){
  `%||%` <- function(x, y) if (is.null(x)) y else x
  .clamp01 <- function(x) pmin(pmax(x, 1e-6), 1-1e-6)
  
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
  
  .build_CI_tests_priors <- function(J, tests_beta = NULL, ranges = NULL,
                                     range_ci = 0.95, ranges_override = TRUE,
                                     default_beta = c(1,1)) {
    out <- vector("list", J)
    for (j in 1:J) out[[j]] <- list(sens = default_beta, spec = default_beta)
    if (!is.null(tests_beta)) {
      stopifnot(is.list(tests_beta), length(tests_beta)==J)
      for (j in 1:J) {
        tj <- tests_beta[[j]]
        if (!is.null(tj$sens)) out[[j]]$sens <- as.numeric(tj$sens)
        if (!is.null(tj$spec)) out[[j]]$spec <- as.numeric(tj$spec)
      }
    }
    if (!is.null(ranges)) {
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
  
  .eta_norm_from_ranges <- function(sens_range, spec_range, ci = 0.95) {
    z <- qnorm((1+ci)/2)
    Lse <- .clamp01(min(sens_range)); Use <- .clamp01(max(sens_range))
    m_se  <- qnorm((Lse + Use)/2)
    sd_se <- (qnorm(Use) - qnorm(Lse)) / (2*z); sd_se <- max(sd_se, 0.02)
    Lsp <- .clamp01(min(spec_range)); Usp <- .clamp01(max(spec_range))
    m_sp  <- qnorm((Lsp + Usp)/2)
    sd_sp <- (qnorm(Usp) - qnorm(Lsp)) / (2*z); sd_sp <- max(sd_sp, 0.02)
    list(m_se=m_se, sd_se=sd_se, m_sp=m_sp, sd_sp=sd_sp)
  }
  
  .build_RE_priors_from_ranges <- function(J, tests_norm = NULL, ranges = NULL,
                                           range_ci = 0.95, common_slopes = FALSE,
                                           default_norm = list(mean=0, sd=1),
                                           ndraw = 200000L) {
    mu_a0 <- rep(default_norm$mean, J); s2_a0 <- rep(default_norm$sd^2, J)
    mu_a1 <- rep(default_norm$mean, J); s2_a1 <- rep(default_norm$sd^2, J)
    mu_b0 <- rep(default_norm$mean, J); s2_b0 <- rep(default_norm$sd^2, J)
    mu_b1 <- rep(default_norm$mean, J); s2_b1 <- rep(default_norm$sd^2, J)
    
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
    
    if (!is.null(ranges)) {
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
  
  if (is.null(model)) model <- fit$model %||% stop("Provide 'model' or ensure fit$model exists.")
  Ypost_sens <- fit$sens
  Ypost_spec <- fit$spec
  Ypost_rho  <- fit$rho
  J <- if (!is.null(Ypost_spec)) ncol(Ypost_spec) else stop("fit$spec missing.")
  
  # If user didn't pass the priors, try to find them on the fit (if you stored them), else stop.
  if (is.null(prior_input)) {
    prior_input <- attr(fit, "prior_input")
    if (is.null(prior_input)) stop("Please pass 'prior_input' (the same object used in bayes_2LCR).")
  }
  
  # prevalence first
  prev_ab <- prior_input$prev %||% c(1,1)
  rho_prior <- rbeta(nsim, prev_ab[1], prev_ab[2])
  rho_df_prior <- tibble::tibble(
    value = rho_prior, quantity = "rho", test = "global", type = "Prior"
  )
  
  if (identical(tolower(model), "ci")) {
    tests_beta <- prior_input$tests
    tests_rng  <- prior_input$ranges
    range_ci   <- prior_input$range_ci %||% 0.95
    if (is.null(tests_beta) && is.null(tests_rng)) {
      tests_beta <- vector("list", J)
      for (j in 1:J) tests_beta[[j]] <- list(sens = c(1,1), spec = c(1,1))
    } else if (!is.null(tests_rng)) {
      tests_beta <- .build_CI_tests_priors(
        J, tests_beta = tests_beta, ranges = tests_rng, range_ci = range_ci,
        ranges_override = TRUE, default_beta = c(1,1)
      )
    } else {
      stopifnot(length(tests_beta)==J)
    }
    
    prior_tests_df <- purrr::map_dfr(seq_len(J), function(j){
      ab_se <- tests_beta[[j]]$sens; ab_sp <- tests_beta[[j]]$spec
      se <- rbeta(nsim, ab_se[1], ab_se[2])
      sp <- rbeta(nsim, ab_sp[1], ab_sp[2])
      tibble::tibble(test = paste0("Test ", j), Se = se, Sp = sp)
    }) |>
      tidyr::pivot_longer(c(Se, Sp), names_to = "quantity", values_to = "value") |>
      dplyr::mutate(type = "Prior")
    
  } else {
    # RANDOM / 2LCR1: build Normal priors for a0/a1/b0/b1 (allow ranges)
    # Reuse the same logic as in bayes_2LCR
    # If user provided ranges, they influence a0/a1 via slope priors.
    # Otherwise, explicit normals (or defaults) drive everything.
    # NOTE: we only need means/vars to *draw* from.
    .fill_default_tests_random <- function(J) {
      one <- list(
        a0 = list(mean = 0, sd = 1),
        a1 = list(mean = 0, sd = 1),
        b0 = list(mean = 0, sd = 1),
        b1 = list(mean = 0, sd = 1)
      )
      rep(list(one), J)
    }
    tests_list <- prior_input$tests %||% .fill_default_tests_random(J)
    
    REp <- .build_RE_priors_from_ranges(
      J = J,
      tests_norm   = tests_list,
      ranges       = prior_input$ranges,
      range_ci     = prior_input$range_ci %||% 0.95,
      common_slopes= common_slopes,
      default_norm = list(mean=0, sd=1),
      ndraw        = prior_input$range_ndraw %||% 200000L
    )
    mu_a0 <- REp$mu_a0; s2_a0 <- REp$s2_a0
    mu_a1 <- REp$mu_a1; s2_a1 <- REp$s2_a1
    mu_b0 <- REp$mu_b0; s2_b0 <- REp$s2_b0
    mu_b1 <- REp$mu_b1; s2_b1 <- REp$s2_b1
    
    # Draw and compute implied Se/Sp by model structure
    prior_tests_df <- purrr::map_dfr(seq_len(J), function(j){
      if (identical(tolower(model), "2lcr1")) {
        if (common_slopes) {
          # one slope b shared across tests AND classes (per draw)
          b_draw  <- rnorm(nsim, mean = mu_b0[1], sd = sqrt(s2_b0[1]))
          s_draw  <- sqrt(1 + b_draw^2)
          a0_draw <- rnorm(nsim, mu_a0[j], sqrt(s2_a0[j]))
          a1_draw <- rnorm(nsim, mu_a1[j], sqrt(s2_a1[j]))
          Se <- pnorm(a1_draw / s_draw)
          Sp <- pnorm(-a0_draw / s_draw)
        } else {
          # one slope per class, shared across tests (per draw)
          b0_draw <- rnorm(nsim, mu_b0[1], sqrt(s2_b0[1])); s0 <- sqrt(1 + b0_draw^2)
          b1_draw <- rnorm(nsim, mu_b1[1], sqrt(s2_b1[1])); s1 <- sqrt(1 + b1_draw^2)
          a0_draw <- rnorm(nsim, mu_a0[j], sqrt(s2_a0[j]))
          a1_draw <- rnorm(nsim, mu_a1[j], sqrt(s2_a1[j]))
          Se <- pnorm(a1_draw / s1)
          Sp <- pnorm(-a0_draw / s0)
        }
      } else { # "random"
        if (common_slopes) {
          # one slope per test, shared across classes
          b_draw  <- rnorm(nsim, mu_b0[j], sqrt(s2_b0[j]))  # use b0 prior
          s_draw  <- sqrt(1 + b_draw^2)
          a0_draw <- rnorm(nsim, mu_a0[j], sqrt(s2_a0[j]))
          a1_draw <- rnorm(nsim, mu_a1[j], sqrt(s2_a1[j]))
          Se <- pnorm(a1_draw / s_draw)
          Sp <- pnorm(-a0_draw / s_draw)
        } else {
          # per-test, per-class slopes
          b0_draw <- rnorm(nsim, mu_b0[j], sqrt(s2_b0[j])); s0 <- sqrt(1 + b0_draw^2)
          b1_draw <- rnorm(nsim, mu_b1[j], sqrt(s2_b1[j])); s1 <- sqrt(1 + b1_draw^2)
          a0_draw <- rnorm(nsim, mu_a0[j], sqrt(s2_a0[j]))
          a1_draw <- rnorm(nsim, mu_a1[j], sqrt(s2_a1[j]))
          Se <- pnorm(a1_draw / s1)
          Sp <- pnorm(-a0_draw / s0)
        }
      }
      tibble::tibble(test = paste0("Test ", j), Se = Se, Sp = Sp)
    }) |>
      tidyr::pivot_longer(c(Se, Sp), names_to = "quantity", values_to = "value") |>
      dplyr::mutate(type = "Prior")
  }
  
  # ---- posterior draws (already on Se/Sp scale in your fit) ----
  post_tests_df <- purrr::map_dfr(seq_len(J), function(j){
    tibble::tibble(
      test = paste0("Test ", j),
      Se   = Ypost_sens[, j],
      Sp   = Ypost_spec[, j]
    )
  }) |>
    tidyr::pivot_longer(c(Se, Sp), names_to = "quantity", values_to = "value") |>
    dplyr::mutate(type = "Posterior")
  
  rho_df_post <- tibble::tibble(
    value = as.numeric(Ypost_rho),
    quantity = "rho", test = "global", type = "Posterior"
  )
  
  dat <- dplyr::bind_rows(prior_tests_df, post_tests_df, rho_df_prior, rho_df_post)

  p_rho <- dat |>
    dplyr::filter(quantity == "rho") |>
    ggplot2::ggplot(ggplot2::aes(x = value, fill = type)) +
    ggplot2::geom_density(alpha = 0.35, adjust = 1.2) +
    ggplot2::labs(x = expression(rho), y = "Density",
                  title = "Prevalence (rho): Prior vs Posterior") +
    ggplot2::theme_minimal()
  
  p_tests <- dat |>
    dplyr::filter(quantity != "rho") |>
    ggplot2::ggplot(ggplot2::aes(x = value, fill = type)) +
    ggplot2::geom_density(alpha = 0.35, adjust = 1.2) +
    ggplot2::facet_grid(quantity ~ test, scales = "free") +
    ggplot2::labs(x = "Probability", y = "Density",
                  title = "Se/Sp: Prior vs Posterior") +
    ggplot2::theme_minimal()
  
  list(per_test = p_tests, rho = p_rho)
}
