


####################################################################
################# Strongyloides data Prior vs Posterior ############
####################################################################

library(ggplot2)
library(dplyr)
library(purrr)
library(tidyr)
library(patchwork)
source("plot_prior_vs_posterior.R")
fit_RE_CI   <- readRDS("Strongyloides_CI.RDS")         # Bayesian Random Effect: CI
fit_RE_RAND <- readRDS("Strongyloides_random.RDS")     # Bayesian Random Effect: Random Effect
fit_BLS_CI    <- readRDS("Strongyloides_BLS_CI.RDS")   # Our Model: CI
fit_BLS_Gamma <- readRDS("Strongyloides_BLS_Gamma.RDS")# Our Model: Gamma
fit_BLS_NM    <- readRDS("Strongyloides_BLS_NM.RDS")   # Our Model: NM+

################# Dendukuri & Joseph Bayesian RE ############
theme_pub <- function() theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "grey95"),
        panel.grid.minor = element_blank())
fill_sc  <- scale_fill_manual(values = c(Prior = "grey70", Posterior = "steelblue"))
color_sc <- scale_color_manual(values = c(Prior = "grey45", Posterior = "steelblue4"))

get_rho  <- function(fit) fit$rho  %||% fit$rho_Samples
get_sens <- function(fit) fit$sens %||% fit$sensitivity_Samples
get_spec <- function(fit) fit$spec %||% fit$specificity_Samples

draw_prior_CI <- function(pr, nsim = 50000){
  J <- length(pr$tests)
  list(
    rho  = rbeta(nsim, pr$prev[1], pr$prev[2]),
    sens = lapply(seq_len(J), \(j) rbeta(nsim, pr$tests[[j]]$sens[1], pr$tests[[j]]$sens[2])),
    spec = lapply(seq_len(J), \(j) rbeta(nsim, pr$tests[[j]]$spec[1], pr$tests[[j]]$spec[2]))
  )
}
draw_prior_2LCR1 <- function(pr, nsim = 50000, common_slopes = FALSE){
  se_from <- function(a1,b1) pnorm(a1 / sqrt(1 + b1^2))
  sp_from <- function(a0,b0) pnorm(-a0 / sqrt(1 + b0^2))
  J <- length(pr$tests)
  
  if (common_slopes){
    b0 <- rnorm(nsim, pr$tests[[1]]$b0$mean, pr$tests[[1]]$b0$sd)
    b1 <- rnorm(nsim, pr$tests[[1]]$b1$mean, pr$tests[[1]]$b1$sd)
    s0 <- sqrt(1 + b0^2); s1 <- sqrt(1 + b1^2)
    sens <- spec <- vector("list", J)
    for (j in seq_len(J)){
      a0 <- rnorm(nsim, pr$tests[[j]]$a0$mean, pr$tests[[j]]$a0$sd)
      a1 <- rnorm(nsim, pr$tests[[j]]$a1$mean, pr$tests[[j]]$a1$sd)
      sens[[j]] <- pnorm(a1 / s1)
      spec[[j]] <- pnorm(-a0 / s0)
    }
  } else {
    # per-class slopes (shared across tests)
    b0 <- rnorm(nsim, pr$tests[[1]]$b0$mean, pr$tests[[1]]$b0$sd)
    b1 <- rnorm(nsim, pr$tests[[1]]$b1$mean, pr$tests[[1]]$b1$sd)
    s0 <- sqrt(1 + b0^2); s1 <- sqrt(1 + b1^2)
    sens <- spec <- vector("list", J)
    for (j in seq_len(J)){
      a0 <- rnorm(nsim, pr$tests[[j]]$a0$mean, pr$tests[[j]]$a0$sd)
      a1 <- rnorm(nsim, pr$tests[[j]]$a1$mean, pr$tests[[j]]$a1$sd)
      sens[[j]] <- pnorm(a1 / s1)
      spec[[j]] <- pnorm(-a0 / s0)
    }
  }
  list(rho = rbeta(nsim, pr$prev[1], pr$prev[2]), sens = sens, spec = spec)
}

tidy_prior <- function(draws, label){
  J <- length(draws$sens)
  bind_rows(
    tibble(quantity = "rho",  test = "global", value = draws$rho, type = "Prior", model = label),
    map_dfr(seq_len(J), \(j) tibble(quantity="Se", test=paste0("Test ",j), value=draws$sens[[j]], type="Prior", model=label)),
    map_dfr(seq_len(J), \(j) tibble(quantity="Sp", test=paste0("Test ",j), value=draws$spec[[j]], type="Prior", model=label))
  )
}
tidy_post <- function(fit, label){
  rho <- as.numeric(get_rho(fit)); Se <- get_sens(fit); Sp <- get_spec(fit); J <- ncol(Se)
  bind_rows(
    tibble(quantity = "rho", test = "global", value = rho, type = "Posterior", model = label),
    map_dfr(seq_len(J), \(j) tibble(quantity="Se", test=paste0("Test ",j), value=Se[,j], type="Posterior", model=label)),
    map_dfr(seq_len(J), \(j) tibble(quantity="Sp", test=paste0("Test ",j), value=Sp[,j], type="Posterior", model=label))
  )
}
make_row_plots <- function(prior_draws, fit, model_lab, width_ratio = c(1, 3.2)) {
  dd <- dplyr::bind_rows(tidy_prior(prior_draws, model_lab),
                         tidy_post (fit,         model_lab)) %>%
    dplyr::mutate(type = factor(type, levels = c("Prior","Posterior")))
  
  # LEFT: prevalence
  p_rho <- dd %>%
    dplyr::filter(quantity == "rho") %>%
    ggplot2::ggplot(ggplot2::aes(x = value, fill = type, color = type)) +
    ggplot2::geom_density(alpha = .18, linewidth = .6) +
    ggplot2::labs(title = paste0(model_lab, ": Prevalence (", "\u03C1", ")"),
                  x = NULL, y = "Density", fill = NULL, color = NULL) +
    fill_sc + color_sc + theme_pub()
  
  # RIGHT: 2×2 grid Se/Sp
  p_tests <- dd %>%
    dplyr::filter(quantity != "rho") %>%
    ggplot2::ggplot(ggplot2::aes(x = value, fill = type, color = type)) +
    ggplot2::geom_density(alpha = .18, linewidth = .6) +
    ggplot2::facet_grid(quantity ~ test, scales = "free_y") +
    ggplot2::labs(title = paste0(model_lab, ": Sensitivity (Se) and Specificity (Sp)"),
                  x = NULL, y = "Density", fill = NULL, color = NULL) +
    fill_sc + color_sc + theme_pub() +
    ggplot2::theme(panel.spacing = grid::unit(0.9, "lines"))
  
  # Widen the right side
  p_rho | p_tests + plot_layout(widths = width_ratio)
}

prior_CI <- list(
  prev = c(1,1),
  tests = list(
    list(sens=c(4.44, 13.31), spec=c(71.25, 3.75)),
    list(sens=c(21.96, 5.49), spec=c(4.1, 1.76))
  )
)
prior_rand <- list(
  prev = c(1,1),   # Beta(a,b) for prevalence
  tests = list(
    # Test 1 priors (means & SDs)
    list(a1 = list(mean = -0.811, sd = 0.380),   
         a0 = list(mean =  -2.171, sd = 0.261),  ## in the paper spec was pnorm(a0/ sqrt(1+b0^2))
         b1 = list(mean =  0.668, sd = 0.5),  
         b0 = list(mean =  0.861, sd = 0.5)), 
    # Test 2 priors
    list(a1 = list(mean =  1.012, sd = 0.268),
         a0 = list(mean =  -0.692, sd = 0.560),
         b1 = list(mean =  0.668, sd = 0.5),    
         b0 = list(mean =  0.861, sd = 0.5))
  )
)

nsim <- 500000
pr_draws_CI <- draw_prior_CI(prior_CI, nsim)
pr_draws_RE <- draw_prior_2LCR1(prior_rand, nsim, common_slopes = FALSE)

row_CI <- make_row_plots(prior_draws = pr_draws_CI,
                         fit         = fit_RE_CI,
                         model_lab   = "Conditional Independence",
                         width_ratio = c(2, 0.20))

row_RE <- make_row_plots(prior_draws = pr_draws_RE,
                         fit         = fit_RE_RAND,
                         model_lab   = "Random Effects",
                         width_ratio = c(2, 0.20))

DJ_plot <- row_CI / row_RE + plot_layout(guides = "collect")
DJ_plot <- DJ_plot & theme(legend.position = "bottom") 

print(DJ_plot)
ggsave("Strongyloides_DJ_Plot.png",
       DJ_plot, width = 15.5, height = 9.5, dpi = 300)


quantile(fit_RE_CI$rho, c(0.5, 0.025, 0.975))
quantile(fit_RE_CI$sens[,1], c(0.5, 0.025, 0.975))
quantile(fit_RE_CI$spec[,1], c(0.5, 0.025, 0.975))
quantile(fit_RE_CI$sens[,2], c(0.5, 0.025, 0.975))
quantile(fit_RE_CI$spec[,2], c(0.5, 0.025, 0.975))

quantile(fit_RE_RAND$rho, c(0.5, 0.025, 0.975))
quantile(fit_RE_RAND$sens[,1], c(0.5, 0.025, 0.975))
quantile(fit_RE_RAND$spec[,1], c(0.5, 0.025, 0.975))
quantile(fit_RE_RAND$sens[,2], c(0.5, 0.025, 0.975))
quantile(fit_RE_RAND$spec[,2], c(0.5, 0.025, 0.975))
quantile(fit_RE_RAND$a1[,1], c(0.5, 0.025, 0.975))
quantile(fit_RE_RAND$a1[,2], c(0.5, 0.025, 0.975))
quantile(fit_RE_RAND$a0[,1], c(0.5, 0.025, 0.975))
quantile(fit_RE_RAND$a0[,2], c(0.5, 0.025, 0.975))
quantile(fit_RE_RAND$b0[,1], c(0.5, 0.025, 0.975))
quantile(fit_RE_RAND$b1[,1], c(0.5, 0.025, 0.975))

################# Our Model ############

row_from_fit <- function(fit, title_left, title_right, nsim = 20000, widths = c(1, 3)) {
  pp <- plot_prior_vs_posterior(fit, nsim = nsim)
  p_rho   <- pp$rho      + ggplot2::labs(title = title_left)
  p_tests <- pp$per_test + ggplot2::labs(title = title_right)
  p_rho | p_tests + plot_layout(widths = widths)
}

# 2) Build rows for your three BLS models
row_CI    <- row_from_fit(fit_BLS_CI,
                          title_left  = "BLS–CI: Prevalence (ρ)",
                          title_right = "BLS–CI: Sensitivity (Se) and Specificity (Sp)",
                          nsim = 20000, widths = c(2, 0.2))

row_Gamma <- row_from_fit(fit_BLS_Gamma,
                          title_left  = "BLS–Gamma: Prevalence (ρ)",
                          title_right = "BLS–Gamma: Sensitivity (Se) and Specificity (Sp)",
                          nsim = 20000, widths = c(2, 0.2))

row_NM    <- row_from_fit(fit_BLS_NM,
                          title_left  = "BLS–NM+: Prevalence (ρ)",
                          title_right = "BLS–NM+: Sensitivity (Se) and Specificity (Sp)",
                          nsim = 20000, widths = c(2, 0.2))

# 3) Stack the three rows and collect a single legend
BLS_grid <- row_CI / row_Gamma / row_NM + plot_layout(guides = "collect")
BLS_grid <- BLS_grid & theme(legend.position = "bottom")

print(BLS_grid)










