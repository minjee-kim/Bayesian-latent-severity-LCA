

plot_prior_vs_posterior <- function(fit, nsim = 20000) {
  stopifnot(!is.null(attr(fit, "priors")))
    E_phi_betaS_minus_gamma <- function(beta, gamma, sev, rel.tol = 1e-6) {
    stype <- tolower(sev$type)
    if (stype == "ci") return(pnorm(beta - gamma))
    fs <- .severity_pdf(sev)
    upper <- if (stype == "gamma") Inf else sev$mu0 + 8*(sev$tau %||% 1.48495)
    mapply(function(b,g){
      integrate(function(s) pnorm(b*s - g) * fs(s),
                lower = 0, upper = upper, rel.tol = rel.tol)$value
    }, beta, gamma)
  }
  
  .calibrate_beta_from_range_quad <- function(
    Lse, Use, m_gamma, sd_gamma, sev,
    init_m = 1, init_sd = 0.6, nsim = 20000,
    z = qnorm(0.975),
    sd_bounds = c(1e-3, 5), m_bounds = c(1e-6, 15)
  ){
    stopifnot(Lse > 0, Use < 1, Lse < Use)
    target <- c(qL = Lse, qU = Use)
    obj <- function(par) {
      m  <- exp(par[1]); sd <- exp(par[2])
      beta  <- truncnorm::rtruncnorm(nsim, a=0, b=Inf, mean=m, sd=sd)
      gamma <- rnorm(nsim, mean = m_gamma, sd = sd_gamma)
      se_vals <- E_phi_betaS_minus_gamma(beta, gamma, sev)
      qs <- stats::quantile(se_vals, probs = c(0.025, 0.975), names = FALSE)
      sum((qs - target)^2)
    }
    par0 <- c(log(init_m), log(init_sd))
    obj_pen <- function(par){
      m <- exp(par[1]); sd <- exp(par[2]); pen <- 0
      if (m < m_bounds[1])  pen <- pen + (m_bounds[1]-m)^2
      if (m > m_bounds[2])  pen <- pen + (m-m_bounds[2])^2
      if (sd < sd_bounds[1]) pen <- pen + (sd_bounds[1]-sd)^2
      if (sd > sd_bounds[2]) pen <- pen + (sd-sd_bounds[2])^2
      obj(par) + 1e3*pen
    }
    o <- optim(par0, obj_pen, method="Nelder-Mead",
               control=list(maxit=200, reltol=1e-6))
    list(m_beta = exp(o$par[1]), sd_beta = exp(o$par[2]),
         value = o$value, converged = (o$convergence==0))
  }
  .severity_pdf <- function(sev) {
    stype <- tolower(sev$type)
    if (stype == "gamma") {
      a <- sev$aS; b <- sev$bS
      return(function(s) ifelse(s >= 0, dgamma(s, shape = a, rate = b), 0))
    } else if (stype %in% c("nm+","nm_plus","nmplus","normal moment")) {
      mu0 <- sev$mu0 %||% 0
      tau <- sev$tau %||% 1.48495
      g_unnorm <- function(s) ifelse(
        s >= 0,
        2 * (s^2) * dnorm((s - mu0)/tau) / (tau^4 * sqrt(2*pi)),
        0
      )
      U <- mu0 + 8*tau
      Z <- integrate(g_unnorm, lower = 0, upper = U, rel.tol = 1e-8)$value
      return(function(s) ifelse(s >= 0 & s <= U, g_unnorm(s)/Z, 0))
    } else if (stype == "ci") {
      return(function(s) 0) # unused
    }
    stop("Unknown severity type.")
  }
  prior_predictive_draws <- function(result, nsim = 50000){
    pr    <- attr(result, "priors")
    sev   <- pr$severity_prior
    stype <- tolower(sev$type)
    perBG <- pr$per_test
    J     <- length(perBG)
    
    # rho prior
    a_rho <- as.numeric(pr$rho_prior["a"])
    b_rho <- as.numeric(pr$rho_prior["b"])
    rho_vec <- rbeta(nsim, a_rho, b_rho)
    rho_draws <- tibble::tibble(
      value = rho_vec, quantity = "rho", test = "global", type = "Prior"
    )
    
    per_test_prior <- purrr::map_dfr(seq_len(J), function(j){
      m_b <- perBG[[j]]$beta$mean; sd_b <- perBG[[j]]$beta$sd
      m_g <- perBG[[j]]$gamma$mean; sd_g <- perBG[[j]]$gamma$sd
      beta  <- truncnorm::rtruncnorm(nsim, a=0, b=Inf, mean=m_b, sd=sd_b)
      gamma <- rnorm(nsim, m_g, sd_g)
      
      if (stype == "ci") {
        se <- pnorm(beta - gamma)
        sp <- pnorm(gamma)
      } else {
        se <- E_phi_betaS_minus_gamma(beta, gamma, sev) # <- quadrature
        sp <- pnorm(gamma)
      }
      
      tibble::tibble(
        test = paste0("Test ", j),
        Se   = se,
        Sp   = sp
      )
    }) |>
      tidyr::pivot_longer(c(Se, Sp), names_to = "quantity", values_to = "value") |>
      dplyr::mutate(type = "Prior")
    
    dplyr::bind_rows(per_test_prior, rho_draws)
  }
  
  posterior_draws_tidy <- function(result){
    sens <- result$sensitivity_Samples
    spec <- result$specificity_Samples
    rho  <- result$rho_Samples
    J    <- ncol(spec)
    
    per_test_post <- purrr::map_dfr(seq_len(J), function(j){
      tibble::tibble(
        test = paste0("Test ", j),
        Se   = sens[, j],
        Sp   = spec[, j]
      )
    }) |>
      tidyr::pivot_longer(c(Se, Sp), names_to = "quantity", values_to = "value") |>
      dplyr::mutate(type = "Posterior")
    
    rho_df <- tibble::tibble(value = rho, quantity = "rho", test = "global", type = "Posterior")
    dplyr::bind_rows(per_test_post, rho_df) |>
      dplyr::filter(is.finite(value))
  }
  
  prior <- prior_predictive_draws(fit, nsim = nsim)
  post  <- posterior_draws_tidy(fit)
  dat   <- dplyr::bind_rows(prior, post)
  
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
                  title = "Se/Sp: Model-aware Prior vs Posterior") +
    ggplot2::theme_minimal()
  
  list(per_test = p_tests, rho = p_rho)
}


