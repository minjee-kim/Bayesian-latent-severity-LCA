


plot_prior_vs_posterior_2LCR <- function(
    fit,
    nsim = 20000
){
  `%||%` <- function(x, y) if (is.null(x)) y else x
  theme_pub <- function(base = 12){
    ggplot2::theme_bw(base_size = base) +
      ggplot2::theme(
        legend.position = "bottom",
        strip.background = ggplot2::element_rect(fill = "grey95", color = NA),
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_line(linewidth = 0.25, colour = "grey90"),
        plot.title.position = "plot",
        plot.caption.position = "plot"
      )
  }
  scales_prior_post <- list(
    ggplot2::scale_colour_manual(NULL, values = c(Prior = "grey30", Posterior = "steelblue4")),
    ggplot2::scale_fill_manual  (NULL, values = c(Prior = "grey70", Posterior = "steelblue"))
  )
  
  pr <- fit$priors
  if (is.null(pr)) stop("No priors found on 'fit'. Re-run bayes_2LCR() from this version.")
  model <- tolower(pr$model %||% fit$model)
  J <- if (!is.null(fit$spec)) ncol(fit$spec) else stop("fit$spec missing.")
  a <- pr$prev[1]; b <- pr$prev[2]
  
  # prevalence
  xg <- seq(0, 1, length.out = 2000)
  prior_df <- data.frame(x = xg, y = dbeta(xg, a, b), grp = "Prior")
  rho_post <- data.frame(value = as.numeric(fit$rho), grp = "Posterior")
  
  # per-test Se/Sp priors
  if (identical(model, "ci")) {
    tests_beta <- pr$tests_beta
    prior_tests_df <- purrr::map_dfr(seq_len(J), function(j){
      ab_se <- tests_beta[[j]]$sens; ab_sp <- tests_beta[[j]]$spec
      tibble::tibble(
        test = paste0("Test ", j),
        Se = rbeta(nsim, ab_se[1], ab_se[2]),
        Sp = rbeta(nsim, ab_sp[1], ab_sp[2])
      )
    }) |>
      tidyr::pivot_longer(c(Se, Sp), names_to = "quantity", values_to = "value") |>
      dplyr::mutate(type = "Prior")
  } else {
    mu_a0 <- pr$norm$mu_a0; s2_a0 <- pr$norm$s2_a0
    mu_a1 <- pr$norm$mu_a1; s2_a1 <- pr$norm$s2_a1
    mu_b0 <- pr$norm$mu_b0; s2_b0 <- pr$norm$s2_b0
    mu_b1 <- pr$norm$mu_b1; s2_b1 <- pr$norm$s2_b1
    common_slopes <- isTRUE(pr$common_slopes)
    prior_tests_df <- purrr::map_dfr(seq_len(J), function(j){
      if (identical(model, "2lcr1")) {
        if (common_slopes) {
          b_draw  <- rnorm(nsim, mean = mu_b0[1], sd = sqrt(s2_b0[1]))
          s_draw  <- sqrt(1 + b_draw^2)
          a0_draw <- rnorm(nsim, mu_a0[j], sqrt(s2_a0[j]))
          a1_draw <- rnorm(nsim, mu_a1[j], sqrt(s2_a1[j]))
          Se <- pnorm(a1_draw / s_draw)
          Sp <- pnorm(-a0_draw / s_draw)
        } else {
          b0_draw <- rnorm(nsim, mu_b0[1], sqrt(s2_b0[1])); s0 <- sqrt(1 + b0_draw^2)
          b1_draw <- rnorm(nsim, mu_b1[1], sqrt(s2_b1[1])); s1 <- sqrt(1 + b1_draw^2)
          a0_draw <- rnorm(nsim, mu_a0[j], sqrt(s2_a0[j]))
          a1_draw <- rnorm(nsim, mu_a1[j], sqrt(s2_a1[j]))
          Se <- pnorm(a1_draw / s1)
          Sp <- pnorm(-a0_draw / s0)
        }
      } else {
        if (common_slopes) {
          b_draw  <- rnorm(nsim, mu_b0[j], sqrt(s2_b0[j]))
          s_draw  <- sqrt(1 + b_draw^2)
          a0_draw <- rnorm(nsim, mu_a0[j], sqrt(s2_a0[j]))
          a1_draw <- rnorm(nsim, mu_a1[j], sqrt(s2_a1[j]))
          Se <- pnorm(a1_draw / s_draw)
          Sp <- pnorm(-a0_draw / s_draw)
        } else {
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
  
  post_tests_df <- purrr::map_dfr(seq_len(J), function(j){
    tibble::tibble(
      test = paste0("Test ", j),
      Se   = fit$sens[, j],
      Sp   = fit$spec[, j]
    )
  }) |>
    tidyr::pivot_longer(c(Se, Sp), names_to = "quantity", values_to = "value") |>
    dplyr::mutate(type = "Posterior")
  
  p_rho <- ggplot2::ggplot() +
    ggplot2::geom_line(data = prior_df,
                       ggplot2::aes(x = x, y = y, colour = grp),
                       linewidth = 0.8, linetype = "dashed") +
    ggplot2::geom_density(data = rho_post,
                          ggplot2::aes(x = value, colour = grp, fill = grp),
                          alpha = 0.25, linewidth = 0.8, adjust = 1.1) +
    scales_prior_post +
    ggplot2::guides(fill = "none") +
    ggplot2::labs(x = expression(rho), y = "Density",
                  title = sprintf("Prevalence: Prior (Beta(%g,%g)) vs Posterior", a, b)) +
    theme_pub()
  
  dat <- dplyr::bind_rows(
    prior_tests_df,
    post_tests_df |> dplyr::rename(grp = type)
  )
  p_tests <- ggplot2::ggplot() +
    ggplot2::geom_density(
      data = dat |> dplyr::filter(grp == "Prior"),
      ggplot2::aes(x = value, colour = grp, fill = grp),
      alpha = 0.30, linewidth = 0.8, adjust = 1.1
    ) +
    ggplot2::geom_density(
      data = dat |> dplyr::filter(grp == "Posterior"),
      ggplot2::aes(x = value, colour = grp, fill = grp),
      alpha = 0.25, linewidth = 0.8, adjust = 1.1
    ) +
    ggplot2::facet_grid(quantity ~ test, scales = "free_y") +
    scales_prior_post +
    ggplot2::labs(x = "Probability", y = "Density",
                  title = "Sensitivity / Specificity: Prior vs Posterior") +
    theme_pub() +
    ggplot2::theme(panel.spacing = grid::unit(0.8, "lines"))
  
  list(per_test = p_tests, rho = p_rho)
}




