

theoretical_SeSp_CI <- function(beta, gamma){
  cbind(Se_true = pnorm(beta - gamma),
        Sp_true = pnorm(gamma))
}

theoretical_SeSp_RE <- function(beta, gamma, shapeS, rateS, nsim = 5e5){
  S <- rgamma(nsim, shape = shapeS, rate = rateS)
  J <- length(beta)
  Se <- vapply(seq_len(J), function(j) mean(pnorm(beta[j] * S - gamma[j])), numeric(1))
  Sp <- pnorm(gamma)
  cbind(Se_true = Se, Sp_true = Sp)
}

`%||%` <- function(x, y) if (is.null(x) || inherits(x, "try-error")) y else x

get_priors <- function(ranges, typ){
  build_priors_from_ranges(ranges, severity = typ) %||%
    build_priors_from_ranges(ranges, severity_type = typ)
}

safe_run <- function(expr, ok_path, err_path) {
  tryCatch({
    obj <- force(expr)
    saveRDS(obj, ok_path)
    invisible(ok_path)
  }, error = function(e) {
    saveRDS(list(error = conditionMessage(e)), err_path)
    stop(conditionMessage(e))
  })
}

run_job <- function(job) {
  Tij        <- job$Tij
  model_type <- job$model_type
  ranges     <- job$ranges
  out        <- job$out
  iters      <- job$iters
  burn       <- job$burn
  thin       <- job$thin
  
  if (model_type %in% c("2LCR_CI","2LCR1")) {
    prior_2LCR <- list(prev = c(1,1), ranges = ranges, range_ci = 0.95)
    mdl <- if (model_type == "2LCR_CI") "CI" else "2LCR1"
    safe_run(
      expr = bayes_2LCR(
        data = Tij, model = mdl, common_slopes = FALSE,
        iterations = iters, burnin = burn, thin = thin,
        prior_input = prior_2LCR
      ),
      ok_path  = sprintf("%s_bayes2LCR_%s.RDS", out, mdl),
      err_path = sprintf("%s_bayes2LCR_%s.ERR.RDS", out, mdl)
    )
  } else if (model_type == "BLS_CI") {
    pr <- get_priors(ranges, "CI")
    safe_run(
      expr = Bayesian_LCA_severity(
        data = Tij, iterations = iters, burnin = burn, thin = thin,
        severity = "CI",
        mu_beta = pr$mu_beta,  sd_beta = pr$sd_beta,
        mu_gamma = pr$mu_gamma, sd_gamma = pr$sd_gamma,
        rho_beta = c(1,1)
      ),
      ok_path  = sprintf("%s_fitBLS_CI.RDS", out),
      err_path = sprintf("%s_fitBLS_CI.ERR.RDS", out)
    )
  } else if (model_type == "BLS_gamma") {
    pr <- get_priors(ranges, "gamma")
    safe_run(
      expr = Bayesian_LCA_severity(
        data = Tij, iterations = iters, burnin = burn, thin = thin,
        severity = "gamma",
        mu_beta = pr$mu_beta,  sd_beta = pr$sd_beta,
        mu_gamma = pr$mu_gamma, sd_gamma = pr$sd_gamma,
        aS = pr$aS, bS = pr$bS,
        rho_beta = c(1,1)
      ),
      ok_path  = sprintf("%s_fitBLS_Gamma.RDS", out),
      err_path = sprintf("%s_fitBLS_Gamma.ERR.RDS", out)
    )
  } else {
    stop("Unknown model_type")
  }
}


source("~/Desktop/Bayesian-latent-severity-LCA/R/init.R", chdir = TRUE)

N <- 4000
J <- 4
rho_true <- 0.35
gamma_true <- qnorm(c(0.99, 0.80, 0.70, 0.98))

beta_true_CI <- c(1.0, 2.0, 2.8, 1.5)
D_CI <- rbinom(N, 1, rho_true)
Mu_CI <- outer(D_CI, beta_true_CI) - matrix(gamma_true, N, J, byrow = TRUE)
CI_Tij <- matrix(rbinom(N*J, 1, pnorm(Mu_CI)), nrow = N, ncol = J)
ci_thr  <- theoretical_SeSp_CI(beta_true_CI, gamma_true)


beta_true_RE <- c(0.5, 1.8, 2.95, 0.75)
D_RE <- rbinom(N, 1, rho_true)
S_RE <- ifelse(D_RE == 1, rgamma(N, shape = 4.5, rate = sqrt(4.5)), 0)
Mu_RE <- outer(D_RE * S_RE, beta_true_RE) - matrix(gamma_true, N, J, byrow = TRUE)
RE_Tij <- matrix(rbinom(N*J, 1, pnorm(Mu_RE)), nrow = N, ncol = J)
re_thr  <- theoretical_SeSp_RE(beta_true_RE, gamma_true, shapeS = 4.5, rateS = sqrt(4.5))


informative_ranges <- list(
  list(sens=c(0.01, 0.60), spec=c(0.90, 0.999)),
  list(sens=c(0.60, 0.99), spec=c(0.60, 0.99)),
  list(sens=c(0.60, 0.99), spec=c(0.60, 0.99)),
  list(sens=c(0.01, 0.60), spec=c(0.90, 0.999))
)

noninformative_ranges <- list(
  list(sens=c(0.01, 0.999), spec=c(0.01, 0.999)),
  list(sens=c(0.01, 0.999), spec=c(0.01, 0.999)),
  list(sens=c(0.01, 0.999), spec=c(0.01, 0.999)),
  list(sens=c(0.01, 0.999), spec=c(0.01, 0.999))
)

informative_ranges13 <- list(
  list(sens=c(0.01, 0.60), spec=c(0.90, 0.999)),
  list(sens=c(0.60, 0.99), spec=c(0.60, 0.99))
  )

inf_pr_CI <- build_priors_from_ranges(informative_ranges, severity="CI")
RE_inf_fitBLS_CI <- Bayesian_LCA_severity(
  data       = RE_Tij,
  iterations = 300000,
  burnin     = 100000,
  thin       = 10,
  severity   = "CI",   
  mu_beta    = inf_pr_CI$mu_beta,
  sd_beta    = inf_pr_CI$sd_beta,
  mu_gamma    = inf_pr_CI$mu_gamma,
  sd_gamma   = inf_pr_CI$sd_gamma,
  rho_beta   = c(1,1)
)

inf_pr_gamma <- build_priors_from_ranges(informative_ranges, severity="gamma", aS = 4.5, bS = sqrt(4.5))
RE_inf_fitBLS_Gamma <- Bayesian_LCA_severity(
  data       = RE_Tij,
  iterations = 300000,
  burnin     = 100000,
  thin       = 10,
  severity   = "gamma",   
  mu_gamma    = inf_pr_gamma$mu_gamma,
  sd_gamma   = inf_pr_gamma$sd_gamma,
  mu_beta = inf_pr_gamma$mu_beta, 
  sd_beta = inf_pr_gamma$sd_beta,
  aS = 4.5, bS = sqrt(4.5),
  rho_beta   = c(1,1)
)

inf_pr_gamma_13 <- build_priors_from_ranges(informative_ranges13, severity="gamma", aS = 4.5, bS = sqrt(4.5))
inf_fitBLS_Gamma13 <- Bayesian_LCA_severity(
  data       = RE_Tij,
  iterations = 300000,
  burnin     = 100000,
  thin       = 10,
  severity   = "gamma",   
  mu_gamma    = inf_pr_gamma_13$mu_gamma,
  sd_gamma   = inf_pr_gamma_13$sd_gamma,
  mu_beta = inf_pr_gamma_13$mu_beta, 
  sd_beta = inf_pr_gamma_13$sd_beta,
  aS = 4.5, bS = sqrt(4.5),
  rho_beta   = c(1,1)
)


