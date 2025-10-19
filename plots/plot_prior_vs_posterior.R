

plot_prior_vs_posterior <- function(fit, priors, nsim = 20000) {
  stopifnot(is.list(priors),
            !is.null(priors$severity_prior),
            !is.null(priors$rho_prior),
            !is.null(priors$per_test))
  
  `%||%` <- function(x, y) if (is.null(x)) y else x
  
  sev   <- priors$severity_prior
  stype <- tolower(sev$type %||% "ci")
  per_tests <- priors$per_test
  J <- length(per_tests)
  
  .severity_pdf <- function(sev) {
    s <- tolower(sev$type)
    if (s == "gamma") {
      a <- sev$aS; b <- sev$bS
      return(function(x) ifelse(x >= 0, stats::dgamma(x, shape = a, rate = b), 0))
    } else if (s == "nm+") {
      mu0 <- sev$mu0 %||% 0
      tau <- sev$tau %||% 1.48495
      g_unnorm <- function(x) ifelse(
        x >= 0,
        2 * (x^2) * stats::dnorm((x - mu0)/tau) / (tau^4 * sqrt(2*pi)),
        0
      )
      U <- mu0 + 8*tau
      Z <- stats::integrate(g_unnorm, lower = 0, upper = U, rel.tol = 1e-8)$value
      return(function(x) ifelse(x >= 0 & x <= U, g_unnorm(x)/Z, 0))
    } else if (s == "ci") {
      return(function(x) 0) # unused
    }
    stop("Unknown severity type.")
  }
  
  E_phi_betaS_minus_gamma <- function(beta, gamma, sev, rel.tol = 1e-6) {
    s <- tolower(sev$type)
    if (s == "ci") return(stats::pnorm(beta - gamma))
    fs <- .severity_pdf(sev)
    upper <- if (s == "gamma") Inf else sev$mu0 + 8*(sev$tau %||% 1.48495)
    mapply(function(b,g){
      stats::integrate(function(x) stats::pnorm(b*x - g) * fs(x),
                       lower = 0, upper = upper, rel.tol = rel.tol)$value
    }, beta, gamma)
  }
  
  # prior draws
  a_rho <- as.numeric(priors$rho_prior["a"])
  b_rho <- as.numeric(priors$rho_prior["b"])
  rho_prior <- stats::rbeta(nsim, a_rho, b_rho)
  prior_rho <- tibble::tibble(value = rho_prior, quantity = "rho",
                              test = "global", type = "Prior")
  
  prior_tb <- purrr::map_dfr(seq_len(J), function(j){
    mb <- per_tests[[j]]$m_beta; sdb <- per_tests[[j]]$sd_beta
    mg <- per_tests[[j]]$m_gamma; sdg <- per_tests[[j]]$sd_gamma
    beta  <- truncnorm::rtruncnorm(nsim, a = 0, b = Inf, mean = mb, sd = sdb)
    gamma <- stats::rnorm(nsim, mg, sdg)
    if (stype == "ci") {
      se <- stats::pnorm(beta - gamma)
      sp <- stats::pnorm(gamma)
    } else {
      se <- E_phi_betaS_minus_gamma(beta, gamma, sev)
      sp <- stats::pnorm(gamma)
    }
    tibble::tibble(test = paste0("Test ", j), Se = se, Sp = sp)
  }) |>
    tidyr::pivot_longer(c(Se, Sp), names_to = "quantity", values_to = "value") |>
    dplyr::mutate(type = "Prior")
  
  # posterior draws
  Se_post <- fit$sensitivity_Samples
  Sp_post <- fit$specificity_Samples
  rho_post <- fit$rho_Samples
  Jpost <- ncol(Sp_post)
  
  post_tb <- purrr::map_dfr(seq_len(Jpost), function(j){
    tibble::tibble(test = paste0("Test ", j),
                   Se = Se_post[, j],
                   Sp = Sp_post[, j])
  }) |>
    tidyr::pivot_longer(c(Se, Sp), names_to = "quantity", values_to = "value") |>
    dplyr::mutate(type = "Posterior")
  
  post_rho <- tibble::tibble(value = rho_post, quantity = "rho",
                             test = "global", type = "Posterior")
  
  dat <- dplyr::bind_rows(prior_tb, post_tb, prior_rho, post_rho) %>%
    dplyr::filter(is.finite(value)) %>%
    dplyr::mutate(
      type = stringr::str_trim(type),
      type = stringr::str_to_title(type),         
      type = factor(type, levels = c("Prior","Posterior"))
    )
  
  dat <- dat %>%
    dplyr::filter(is.finite(value)) %>%
    dplyr::mutate(
      type = stringr::str_to_title(trimws(type)),
      type = factor(type, levels = c("Prior","Posterior"))
    )
  
  p_rho <- dat %>%
    dplyr::filter(quantity == "rho") %>%
    ggplot2::ggplot(ggplot2::aes(x = value, fill = type, color = type)) +  # <- color mapped
    ggplot2::geom_density(alpha = 0.18, linewidth = 0.6, adjust = 1.2) +
    ggplot2::labs(x = expression(rho), y = "Density",
                  title = "Prevalence (\u03C1): Prior vs Posterior",
                  fill = NULL, color = NULL) +
    ggplot2::theme_minimal()
  
  p_tests <- dat %>%
    dplyr::filter(quantity != "rho") %>%
    ggplot2::ggplot(ggplot2::aes(x = value, fill = type, color = type)) +  # <- color mapped
    ggplot2::geom_density(alpha = 0.18, linewidth = 0.6, adjust = 1.2) +
    ggplot2::facet_grid(quantity ~ test, scales = "free_y") +
    ggplot2::labs(x = "Probability", y = "Density",
                  title = "Se/Sp: Prior vs Posterior",
                  fill = NULL, color = NULL) +
    ggplot2::theme_minimal()  
  list(per_test = p_tests, rho = p_rho)
}



