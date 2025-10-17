



set.seed(1)
N <- 300; J <- 2
rho_true <- 0.35
D <- rbinom(N, 1, rho_true)

gamma_true <- c(qnorm(0.95), qnorm(0.90)) # Sp ~ (0.95, 0.90)
beta_true  <- c(1.2, 0.8)
S <- ifelse(D==1, rgamma(N, shape=4.5, rate=sqrt(4.5)), 0)
Mu <- outer(D*S, beta_true) - matrix(gamma_true, N, J, byrow=TRUE)
Tij <- matrix( rbinom(N*J, 1, pnorm(Mu)), nrow=N, ncol=J)

# Per-test prior ranges for Se/Sp (interpreted as 95% intervals)
ranges <- list(
  list(sens=c(0.70, 0.95), spec=c(0.88, 0.98)),   # test 1
  list(sens=c(0.60, 0.90), spec=c(0.80, 0.95))    # test 2
)

# Run wrapper (short chain just to demo)
fitA <- Bayesian_LCA_severity(
  data       = data,
  iterations = 20000,
  burnin     = 5000,
  thin       = 2,
  severity_prior = list(type="gamma", aS=2, bS=sqrt(2)),
  ranges     = ranges,
  rho_prior  = list(range = c(0.55, 0.85), ci = 0.95), 
  ci_level   = 0.95,
  nsim       = 30000
)

# Peek at outputs
str(fitA, max.level = 1)
colMeans(fitA$sensitivity_Samples, na.rm=TRUE)   # posterior mean Se per test
pnorm(colMeans(fitA$gamma_Samples))              # posterior mean Sp per test
mean(fitA$rho_Samples)                           # posterior mean prevalence
