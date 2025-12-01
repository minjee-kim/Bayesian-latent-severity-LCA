


plot_prior_vs_posterior <- function(fit, nsim = 20000){
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(truncnorm)
  
  theme_pub <- function(base = 12){
    theme_bw(base_size = base) +
      theme(
        legend.position = "bottom",
        strip.background = element_rect(fill = "grey95", color = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.25, colour = "grey90"),
        plot.title.position = "plot",
        plot.caption.position = "plot"
      )
  }
  
  scales_prior_post <- list(
    scale_colour_manual(NULL, values = c(Prior = "grey30", Posterior = "steelblue4")),
    scale_fill_manual  (NULL, values = c(Prior = "grey70", Posterior = "steelblue"))
  )
  
  pr <- fit$priors
  severity <- tolower(pr$severity)
  J <- length(pr$per_test)
  a <- pr$rho_ab["a"]; b <- pr$rho_ab["b"]
  
  prior_tb <- do.call(rbind, lapply(seq_len(J), function(j){
    pt <- pr$per_test[[j]]
    gamma <- rnorm(nsim, pt$mu_gamma, pt$sd_gamma)
    mu_b <- if (!is.null(pt$mu_beta)) pt$mu_beta else pt$m_beta
    sd_b <- pt$sd_beta
    
    if (severity == "ci") {
      beta <- truncnorm::rtruncnorm(nsim, a = 0, b = Inf, mean = mu_b, sd = sd_b)
      S <- 1
    } else if (severity == "gamma") {
      beta <- truncnorm::rtruncnorm(nsim, a = 0, b = Inf, mean = mu_b, sd = sd_b)
      S <- rgamma(nsim, shape = pr$aS, rate = pr$bS)
    } else {
      stop("Unknown severity in priors: ", pr$severity)
    }
    
    data.frame(test = paste0("Test ", j),
               Se = pnorm(beta * S - gamma),
               Sp = pnorm(gamma),
               grp = "Prior")
  })) |>
    tidyr::pivot_longer(c(Se, Sp), names_to="quantity", values_to="value")
  
  Se_post  <- fit$sensitivity_Samples
  Sp_post  <- fit$specificity_Samples
  rho_post <- as.numeric(fit$rho_Samples)
  
  post_tb <- do.call(rbind, lapply(seq_len(ncol(Sp_post)), function(j)
    data.frame(test = paste0("Test ", j),
               Se   = Se_post[, j],
               Sp   = Sp_post[, j],
               grp  = "Posterior")
  )) |>
    tidyr::pivot_longer(c(Se, Sp), names_to="quantity", values_to="value")
  
  xg <- seq(0, 1, length.out = 2000)
  prior_df <- data.frame(x = xg, y = dbeta(xg, a, b), grp = "Prior")
  post_df  <- data.frame(value = rho_post, grp = "Posterior")
  
  p_rho <- ggplot() +
    geom_line(data = prior_df,
              aes(x = x, y = y, colour = grp),
              linewidth = 0.8, linetype = "dashed") +
    geom_density(data = post_df,
                 aes(x = value, colour = grp, fill = grp),
                 alpha = 0.25, linewidth = 0.8, adjust = 1.1) +
    scales_prior_post +
    guides(fill = "none") +
    labs(x = expression(rho), y = "Density",
         title = "Prevalence: Prior vs Posterior") +
    theme_pub()
  
  p_tests <- ggplot() +
    geom_density(data = prior_tb,
                 aes(x = value, colour = grp, fill = grp),
                 alpha = 0.30, linewidth = 0.8, adjust = 1.1) +
    geom_density(data = post_tb,
                 aes(x = value, colour = grp, fill = grp),
                 alpha = 0.25, linewidth = 0.8, adjust = 1.1) +
    facet_grid(quantity ~ test, scales = "free_y") +
    scales_prior_post +
    labs(x = "Probability", y = "Density",
         title = "Sensitivity / Specificity: Prior vs Posterior") +
    theme_pub() +
    theme(panel.spacing = grid::unit(0.8, "lines"))
  
  list(per_test = p_tests, rho = p_rho)
}
