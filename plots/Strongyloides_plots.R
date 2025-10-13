


####################################################################
################# Strongyloides data Prior vs Posterior ############
####################################################################

library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)

# Probit -> Se/Sp mapping for random-effect model
se_from <- function(a1, b1) pnorm(a1 / sqrt(1 + b1^2))
sp_from <- function(a0, b0) pnorm(-a0 / sqrt(1 + b0^2))

# Safe fetch of test-level prior Beta(a,b) lists for CI model
get_beta_ab <- function(pr, j, field){
  ab <- pr$tests[[j]][[field]]
  stopifnot(length(ab) == 2)
  ab
}
########################### Prior #########################################
# CI priors: Beta for prevalence, Se_j, Sp_j
CI_prior_input = list(
  prev = c(1,1),  # Beta(1,1) on prevalence
  tests = list(
    list(sens=c(4.44, 13.31), spec=c(71.25, 3.75)),  # Test 1
    list(sens=c(21.96, 5.49), spec=c(4.1, 1.76))   # Test 2
  )
)
set.seed(1)
n_prior <- 50000
prior_CI <- {
  pr <- CI_prior_input 
  J  <- length(pr$tests)
  
  rho_prior <- rbeta(n_prior, pr$prev[1], pr$prev[2])
  
  # per test Se/Sp priors
  se_list <- vector("list", J)
  sp_list <- vector("list", J)
  for(j in seq_len(J)){
    ab_se <- get_beta_ab(pr, j, "sens")
    ab_sp <- get_beta_ab(pr, j, "spec")
    se_list[[j]] <- rbeta(n_prior, ab_se[1], ab_se[2])
    sp_list[[j]] <- rbeta(n_prior, ab_sp[1], ab_sp[2])
  }
  
  tibble(
    draw = seq_len(n_prior),
    rho  = rho_prior
  ) %>%
    mutate(model = "CI prior") %>%
    bind_rows(
      map_dfr(seq_len(J), function(j){
        tibble(draw = seq_len(n_prior),
               param = paste0("Se (Test ", j, ")"),
               value = se_list[[j]])
      }) %>% mutate(model = "CI prior")
    ) %>%
    bind_rows(
      map_dfr(seq_len(J), function(j){
        tibble(draw = seq_len(n_prior),
               param = paste0("Sp (Test ", j, ")"),
               value = sp_list[[j]])
      }) %>% mutate(model = "CI prior")
    ) %>%
    # reshape rho into same (param, value) layout for faceting
    bind_rows(
      tibble(draw = seq_len(n_prior),
             param = "Prevalence",
             value = rho_prior,
             model = "CI prior")
    )
}

# 2LCR1 priors
prior_input_rand <- list(
  prev = c(1,1),   # Beta(a,b) for prevalence
  tests = list(
    # Test 1 priors (means & SDs)
    list(a1 = list(mean = -0.811, sd = 0.380),   # class 0
         a0 = list(mean =  2.171, sd = 0.261),   # class 1
         b1 = list(mean =  0.668, sd = 0.5),   # class 0 slope
         b0 = list(mean =  0.861, sd = 0.5)),   # class 1 slope
    # Test 2 priors
    list(a1 = list(mean =  1.012, sd = 0.268),
         a0 = list(mean =  0.692, sd = 0.560),
         b1 = list(mean =  0.668, sd = 0.5),           # class 0 slope
         b0 = list(mean =  0.861, sd = 0.5))
  )
)
Jr <- length(prior_input_rand$tests)

rand_prior_draws <- function(pr, n = 50000, common_slopes = TRUE){
  # draw a0,a1,b0,b1 from Normal priors per your list; compute implied Se/Sp
  a0 <- matrix(NA, n, Jr)
  a1 <- matrix(NA, n, Jr)
  
  if (common_slopes) {
    # class-specific slopes shared across tests: use test 1 for b0/b1 priors
    b0 <- rnorm(n, mean = pr$tests[[1]]$b0$mean, sd = pr$tests[[1]]$b0$sd)
    b1 <- rnorm(n, mean = pr$tests[[1]]$b1$mean, sd = pr$tests[[1]]$b1$sd)
    b0_mat <- matrix(b0, n, Jr)
    b1_mat <- matrix(b1, n, Jr)
  } else {
    # per-test slopes
    b0_mat <- do.call(cbind, lapply(pr$tests, function(tj) rnorm(n, tj$b0$mean, tj$b0$sd)))
    b1_mat <- do.call(cbind, lapply(pr$tests, function(tj) rnorm(n, tj$b1$mean, tj$b1$sd)))
  }
  
  for (j in seq_len(Jr)) {
    a0[, j] <- rnorm(n, pr$tests[[j]]$a0$mean, pr$tests[[j]]$a0$sd)
    a1[, j] <- rnorm(n, pr$tests[[j]]$a1$mean, pr$tests[[j]]$a1$sd)
  }
  
  # implied Se/Sp
  Se <- se_from(a1, b1_mat)
  Sp <- sp_from(a0, b0_mat)
  
  list(Se = Se, Sp = Sp)
}

set.seed(123)
rand_pr <- rand_prior_draws(prior_input_rand, n = n_prior, common_slopes = TRUE)

prior_rand <- bind_rows(
  map_dfr(seq_len(Jr), function(j){
    tibble(draw = seq_len(n_prior),
           param = paste0("Se (Test ", j, ")"),
           value = rand_pr$Se[, j],
           model = "2LCR1 prior")
  }),
  map_dfr(seq_len(Jr), function(j){
    tibble(draw = seq_len(n_prior),
           param = paste0("Sp (Test ", j, ")"),
           value = rand_pr$Sp[, j],
           model = "2LCR1 prior")
  })
)
# prevalence prior is still Beta(1,1) for the random model:
prior_rand <- bind_rows(
  prior_rand,
  tibble(draw = seq_len(n_prior),
         param = "Prevalence",
         value = rbeta(n_prior, prior_input_rand$prev[1], prior_input_rand$prev[2]),
         model = "2LCR1 prior")
)

########################### Posterior #########################################
# CI posterior
fit_CI <- readRDS("Strongyloides_CI.RDS")
post_CI <- {
  J <- ncol(fit_CI$sens)
  tibble(param = "Prevalence", value = fit_CI$rho, model = "CI posterior") %>%
    bind_rows(
      map_dfr(seq_len(J), function(j){
        tibble(param = paste0("Se (Test ", j, ")"),
               value = fit_CI$sens[, j],
               model = "CI posterior")
      }),
      map_dfr(seq_len(J), function(j){
        tibble(param = paste0("Sp (Test ", j, ")"),
               value = fit_CI$spec[, j],
               model = "CI posterior")
      })
    )
}

# 2LCR1 posterior
fit_rand <- readRDS("Strongyloides_random.RDS")
post_rand <- {
  J <- ncol(fit_rand$sens)
  tibble(param = "Prevalence", value = fit_rand$rho, model = "2LCR1 posterior") %>%
    bind_rows(
      map_dfr(seq_len(J), function(j){
        tibble(param = paste0("Se (Test ", j, ")"),
               value = fit_rand$sens[, j],
               model = "2LCR1 posterior")
      }),
      map_dfr(seq_len(J), function(j){
        tibble(param = paste0("Sp (Test ", j, ")"),
               value = fit_rand$spec[, j],
               model = "2LCR1 posterior")
      })
    )
}

########################### Plot #########################################

ci_df <- bind_rows(prior_CI, post_CI) %>%
  mutate(model = factor(model, levels = c("CI prior","CI posterior"))) %>%
  filter(!is.na(param), is.finite(value))

re_df <- bind_rows(prior_rand, post_rand) %>%
  mutate(model = factor(model, levels = c("2LCR1 prior","2LCR1 posterior"))) %>%
  filter(!is.na(param), is.finite(value))
plot_pp <- function(df, title, subtitle) {
  ggplot(df, aes(x = value, fill = model, color = model)) +
    geom_density(alpha = 0.18, linewidth = 0.6, adjust = 1.2) +
    facet_wrap(~ param, scales = "free", ncol = 3) +
    scale_fill_manual(values = c("grey70","steelblue")) +
    scale_color_manual(values = c("grey45","steelblue4")) +
    labs(x = NULL, y = "Density", title = title, subtitle = subtitle, fill = NULL, color = NULL) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom",
          strip.background = element_rect(fill = "grey95"),
          panel.grid.minor = element_blank())
}

gg_ci <- plot_pp(
  ci_df,
  title    = "Priors vs Posteriors (CI model)",
  subtitle = "Prevalence, Sensitivity, Specificity (Beta priors)"
)

gg_re <- plot_pp(
  re_df,
  title    = "Priors vs Posteriors (2LCR1 / Random-Effects)",
  subtitle = "Prevalence (Beta) + Se/Sp implied by Normal priors on probit parameters"
)

print(gg_ci)
print(gg_re)

ggsave("Strongyloides_CI_prior_posterior.png", gg_ci, width = 10, height = 6, dpi = 300)
ggsave("Strongyloides_2LCR1_prior_posterior.png", gg_re, width = 10, height = 6, dpi = 300)


