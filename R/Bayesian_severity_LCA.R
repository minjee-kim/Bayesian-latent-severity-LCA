

source("Gamma_LCA_severity.R")
source("NM_LCA_severity.R")
source("CI_LCA_probit.R")

Bayesian_LCA_severity <- function(
    data, iterations, burnin, thin = 1,
    severity_prior = list(type = "gamma", aS = 3, bS = sqrt(3)),  
    ranges    = NULL,   # list per test: list(sens=c(L,U), spec=c(L,U))
    re_priors = NULL,   # list per test: list(beta=list(mean,sd), gamma=list(mean,sd))
    rho_prior = list(beta = c(1, 1)),    # prevalence prior
    ci_level  = 0.95
){
  library(truncnorm)
  clamp_sd <- 1.5
  
  Tij <- as.matrix(data)
  if (!all(Tij %in% c(0,1))) stop("data must be 0/1.")
  J <- ncol(Tij)
  
  `%||%` <- function(x, y) if (is.null(x)) y else x
  .clamp01 <- function(x) pmin(pmax(x, 1e-6), 1 - 1e-6)
  
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
  
  stype <- tolower(severity_prior$type)
  if (!stype %in% c("gamma","normal moment","nm+","nmplus","nm_plus","ci"))
    stop("severity_prior$type must be 'gamma', 'normal moment' (nm+), or 'ci'.")
  
  # prevalence prior
  translate_rho_prior <- function(rp, ci_default = 0.95){
    if (!is.list(rp)) stop("rho_prior must be a list.")
    if (!is.null(rp$beta)) {
      ab <- rp$beta; if (length(ab)!=2 || any(ab<=0)) stop("rho_prior$beta must be c(a,b)")
      return(list(a=ab[1], b=ab[2]))
    }
    if (!is.null(rp$mean) && !is.null(rp$ess)) {
      m <- rp$mean; k <- rp$ess; if (!(m>0 && m<1 && k>0)) stop("rho_prior mean∈(0,1), ess>0")
      return(list(a=m*k, b=(1-m)*k))
    }
    if (!is.null(rp$range)) {
      L <- .clamp01(min(rp$range)); U <- .clamp01(max(rp$range))
      pL <- (1-(rp$ci %||% ci_default))/2; pU <- 1-pL
      o <- optim(c(log(2),log(2)),
                 function(par){ a<-exp(par[1]); b<-exp(par[2]); (qbeta(pL,a,b)-L)^2+(qbeta(pU,a,b)-U)^2 },
                 method="Nelder-Mead")
      return(list(a=exp(o$par[1]), b=exp(o$par[2])))
    }
    list(a=1,b=1)
  }
  ab <- translate_rho_prior(rho_prior); a_rho <- ab$a; b_rho <- ab$b
  z  <- qnorm((1+ci_level)/2)
  
  if (!is.null(ranges) && !is.null(re_priors)) stop("Provide exactly one of ranges or re_priors.")
  if (is.null(ranges) && is.null(re_priors)) {
    ranges <- replicate(J, list(sens=c(0.5,0.99), spec=c(0.5,0.99)), simplify=FALSE)
  }
  
  priors_beta_gamma <- vector("list", J)
  priors_probit     <- vector("list", J)
  
  if (!is.null(ranges)) {
    # build priors from ranges, model-aware (Se via quadrature)
    for (j in seq_len(J)) {
      rg <- ranges[[j]]
      Lsp <- .clamp01(min(rg$spec))
      Usp <- .clamp01(max(rg$spec))
      m_g  <- qnorm((Lsp + Usp)/2)                                 # NOTE: no minus
      sd_g <- (qnorm(Usp) - qnorm(Lsp)) / (2*z)
      sd_g <- min(max(sd_g, 0.02), clamp_sd)
      
      Lse <- .clamp01(min(rg$sens))
      Use <- .clamp01(max(rg$sens))
      m_se  <- qnorm((Lse + Use)/2)
      sd_se <- (qnorm(Use) - qnorm(Lse)) / (2*z)
      sd_se <- min(max(sd_se, 0.02), clamp_sd)
      
      if (stype == "ci") {
        ndraw <- 200000
        eta_se <- rnorm(ndraw, m_se, sd_se)
        gamma_d <- rnorm(ndraw, m_g, sd_g)
        beta_d  <- eta_se + gamma_d
        beta_d  <- pmax(beta_d, 0)  
        m_b  <- mean(beta_d)
        sd_b <- sd(beta_d)
        
      } else {
        cal <- .calibrate_beta_from_range_quad(Lse, Use, m_g, sd_g, severity_prior,
                                               nsim = 20000)
        m_b <- cal$m_beta; sd_b <- cal$sd_beta
      }
      
      priors_beta_gamma[[j]] <- list(
        beta  = list(mean = m_b,  sd = sd_b),
        gamma = list(mean = m_g,  sd = sd_g)
      )
      priors_probit[[j]] <- list( # for plotting back-translation if needed
        se_probit = list(mean = qnorm((Lse+Use)/2),
                         sd   = min(max((qnorm(Use)-qnorm(Lse))/(2*z),0.02),clamp_sd)),
        sp_probit = list(mean = qnorm((Lsp+Usp)/2), sd = sd_g)
      )
    }
  } else {
    # explicit priors supplied
    if (!is.list(re_priors) || length(re_priors)!=J)
      stop(sprintf("'re_priors' must be a list of length %d.", J))
    for (j in seq_len(J)) {
      stopifnot(!is.null(re_priors[[j]]$beta), !is.null(re_priors[[j]]$gamma))
      priors_beta_gamma[[j]] <- re_priors[[j]]
      priors_probit[[j]] <- NULL
    }
  }
  
  m_gamma <- sapply(priors_beta_gamma, function(p) p$gamma$mean)
  sd_gamma<- sapply(priors_beta_gamma, function(p) p$gamma$sd)
  mu_beta <- sapply(priors_beta_gamma, function(p) p$beta$mean)
  sd_beta <- sapply(priors_beta_gamma, function(p) p$beta$sd)
  
  
  # Fit according to severity type
  res <- switch(stype,
                "ci" = CI_LCA_probit(
                  data       = data,
                  iterations = iterations,
                  burnin     = burnin,
                  thin       = thin,
                  m_beta     = mu_beta,
                  sd_beta    = sd_beta,
                  m_gamma    = m_gamma,
                  sd_gamma   = sd_gamma,
                  a_rho      = a_rho,
                  b_rho      = b_rho
                ),
                "gamma" = Gamma_LCA_severity(
                  data       = data,
                  iterations = iterations,
                  burnin     = burnin,
                  thin       = thin,
                  aS         = severity_prior$aS,
                  bS         = severity_prior$bS,
                  mu_beta    = mu_beta,
                  sd_beta    = sd_beta,
                  m_gamma    = m_gamma,
                  sd_gamma   = sd_gamma,
                  a_rho      = a_rho,
                  b_rho      = b_rho
                ),
                "nm+" =,
                "normal moment" = ,
                "nm_plus" = NM_LCA_severity(
                  data       = data,
                  iterations = iterations,
                  burnin     = burnin,
                  thin       = thin,
                  mu0        = severity_prior$mu0 %||% 0,
                  tau        = severity_prior$tau %||% 1.48495,
                  mu_beta    = mu_beta,
                  sd_beta    = sd_beta,
                  m_gamma    = m_gamma,
                  sd_gamma   = sd_gamma,
                  a_rho      = a_rho,
                  b_rho      = b_rho
                )
  )
  
  attr(res, "priors") <- list(
    severity_prior = severity_prior,
    rho_prior      = c(a = a_rho, b = b_rho),
    per_test       = priors_beta_gamma,   # used by samplers
    per_test_probit= priors_probit,       # original probit priors (if provided/constructed)
    mode           = if (is.null(ranges)) "explicit" else "ranges"
  )
  res
}
