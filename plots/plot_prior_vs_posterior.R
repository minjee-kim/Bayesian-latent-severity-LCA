


plot_prior_vs_posterior <- function(fit, nsim = 20000){
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(truncnorm)
  
  pr <- fit$priors
  J  <- length(pr$per_test)
  
  rho_prior <- rbeta(nsim, pr$rho_ab["a"], pr$rho_ab["b"])
  prior_rho <- data.frame(value = rho_prior, quantity = "rho", test = "global", type = "Prior")
  
  prior_tb <- do.call(rbind, lapply(seq_len(J), function(j){
    mb  <- pr$per_test[[j]]$m_beta; sdb <- pr$per_test[[j]]$sd_beta
    mg  <- pr$per_test[[j]]$m_gamma; sdg <- pr$per_test[[j]]$sd_gamma
    beta  <- rtruncnorm(nsim, 0, Inf, mb, sdb)
    gamma <- rnorm(nsim, mg, sdg)
    if (pr$severity == "ci"){
      se <- pnorm(beta - gamma)
      S  <- rep(1, nsim)
    } else if (pr$severity == "gamma"){
      S  <- rgamma(nsim, pr$aS, pr$bS)
      se <- pnorm(beta * S - gamma)
    } else {
      S  <- rep(1, nsim)
      se <- pnorm(beta * S - gamma)
    }
    sp <- pnorm(gamma)
    data.frame(test = paste0("Test ", j),
               Se = se, Sp = sp, type = "Prior")
  })) |> pivot_longer(c(Se, Sp), names_to="quantity", values_to="value")
  
  Se_post <- fit$sensitivity_Samples
  Sp_post <- fit$specificity_Samples
  rho_post <- fit$rho_Samples
  
  post_tb <- do.call(rbind, lapply(seq_len(ncol(Sp_post)), function(j)
    data.frame(test = paste0("Test ", j),
               Se = Se_post[, j],
               Sp = Sp_post[, j],
               type = "Posterior")
  )) |> pivot_longer(c(Se, Sp), names_to="quantity", values_to="value")
  
  post_rho <- data.frame(value = rho_post, quantity = "rho", test = "global", type = "Posterior")
  dat <- rbind(prior_tb, post_tb, prior_rho, post_rho)
  
  p_rho <- ggplot(dat[dat$quantity=="rho",], aes(x=value, fill=type)) +
    geom_density(alpha=0.25) + theme_minimal() +
    labs(x=expression(rho), y="Density", title="Prevalence: Prior vs Posterior")
  
  p_tests <- ggplot(dat[dat$quantity!="rho",], aes(x=value, fill=type)) +
    geom_density(alpha=0.25) +
    facet_grid(quantity ~ test, scales="free_y") +
    theme_minimal() + labs(x="Probability", y="Density", title="Se/Sp: Prior vs Posterior")
  
  list(per_test = p_tests, rho = p_rho)
}
