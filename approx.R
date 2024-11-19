approx <- function(s, w = 1 -diag(nrow(s)), p = 2, itmax = 10, eps = 1e-10, verbose = TRUE) {
  n <- nrow(s)
  xold <- matrix(rnorm(n * p), n, p)
  itel <- 1
  repeat {
    change <- 0
    xaux <- xold
    for (i in 1:n) {
      z <- drop((s[, i]  * w[, i]) %*% xold)
      v <- crossprod(xold, w[, i] * xold)
      xaux[i, ] <- solve(v, z)
    }
    xnew <- xaux
    for (i in 1:n) {
      z <- drop((s[, i]  * w[, i]) %*% xaux)
      v <- crossprod(xaux, w[, i] * xaux)
      xnew[i, ] <- solve(v, z)
    }
    change <-  max(abs(xold - xnew))
    ssq <- sum(w * (s - tcrossprod(xnew)) ^ 2)
    if (verbose) {
      cat("itel ", formatC(itel, format = "d", width = 4),
          "change ", formatC(change, digits = 10, format = "f"),
          "ssq ", formatC(ssq, digits = 10, format = "f"),
          "\n")
    }
    if ((change < eps) || (itel == itmax)) {
      break
    }
    itel <- itel + 1
    xold <- xnew
  }
  return(list(x = xnew, ssq = ssq, s = s, w = w))
}