#! Rank based Mahalanobis distance, from Paul Rosenbaum's Design of Observational Studies, p. 251
.smahal = function(z, X) {
  X = as.matrix(X)
  n = nrow(X)
  rownames(X) = seq_len(n)
  k = ncol(X)
  m = sum(z)

  for (j in seq_len(k)) {
    X[,j] <- rank(X[,j])
  }

  cv = cov(X)
  vuntied = var(seq_len(n))

  #! ***PENDING: correct this
  diag(cv)[diag(cv) == 0] <- .01
  rat = sqrt(vuntied / diag(cv))
  cv = diag(rat) %*% cv %*% diag(rat)
  out = matrix(NA_real_, nrow = m, ncol = n - m)
  Xc = X[z == 0, ]
  Xt = X[z == 1, ]
  rownames(out) = rownames(X)[z==1]
  colnames(out) = rownames(X)[z==0]

  icov <- MASS::ginv(cv)
  for (i in seq_len(m)) {
    out[i,] <- mahalanobis(Xc, Xt[i,], icov, inverted = TRUE)
  }
  out
}
