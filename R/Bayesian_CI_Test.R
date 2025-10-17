

Bayesian_CI_Test <- function(test1, test2, D_posterior){
  
  S <- nrow(D_posterior)
  N <- length(test1)
  
  chisq <- matrix(NA_real_, nrow = S, ncol = 2, dimnames = list(NULL, c("d=0","d=1")))

  for (s in 1:S) {
    Ds <- D_posterior[s, ]
    for (d in 0:1) {
      idx <- which(Ds == d)
      if (length(idx) < 2) next
      
      t1 <- test1[idx] + 1L
      t2 <- test2[idx] + 1L
      tab <- matrix(0L, 2, 2)
      for (k in seq_along(t1)) tab[t1[k], t2[k]] <- tab[t1[k], t2[k]] + 1L
      
      n  <- sum(tab)
      rs <- rowSums(tab)
      cs <- colSums(tab)
      exp <- outer(rs, cs) / n
      if (any(exp == 0)) next
      
      # Yates correction
      diff <- pmax(abs(tab - exp) - 0.5, 0)
      x2 <- sum((diff^2) / exp)
      
      chisq[s, d + 1] <- x2
    }
  }
  
  list(
    chisq = chisq,
    summary = list(
      mean_chisq = colMeans(chisq, na.rm = TRUE)
    )
  )
}