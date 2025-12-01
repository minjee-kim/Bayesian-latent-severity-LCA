

Bayesian_CI_Test <- function(test1, test2, D_posterior, use_yates = FALSE) {
  test1 <- as.integer(test1); test2 <- as.integer(test2)
  stopifnot(all(test1 %in% 0:1), all(test2 %in% 0:1))
  N <- length(test1)
  
  # rows are draws (S x N)
  if (ncol(D_posterior) == N) {
    S <- nrow(D_posterior)
  } else if (nrow(D_posterior) == N) {
    D_posterior <- t(D_posterior); S <- nrow(D_posterior)
  } else stop("D_posterior must be SxN or NxS (N = length(test1)).")
  
  chisq <- matrix(NA_real_, S, 2, dimnames = list(NULL, c("d=0","d=1")))
  for (s in seq_len(S)) {
    Ds <- D_posterior[s, ]
    for (d in 0:1) {
      idx <- which(Ds == d)
      if (length(idx) < 2) next
      tab <- table(factor(test1[idx], 0:1), factor(test2[idx], 0:1))
      n  <- sum(tab)
      rs <- rowSums(tab)
      cs <- colSums(tab)
      exp <- outer(rs, cs) / n
      if (any(exp == 0)) next
      if (use_yates) {
        diff <- pmax(abs(tab - exp) - 0.5, 0)
        x2 <- sum((diff^2) / exp)
      } else {
        x2 <- sum((tab - exp)^2 / exp)
      }
      chisq[s, d + 1L] <- x2
    }
  }
  list(chisq = chisq, summary = list(mean_chisq = colMeans(chisq, na.rm = TRUE)))
}


