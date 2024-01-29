w <- 1 - diag(4)
delta <- matrix(c(1, 2, 3, 4, 4, 3, 2, 1, 1, 2, 3, 4, 4, 3, 2, 1), 4, 4)
delta <- (delta + t(delta)) / 2
diag(delta) <- 0
x <- matrix(c(1, 1, -1, -1, 1, -1, 1, -1), 4, 2)
ei <- function (i, n) {
  return(ifelse(i == 1:n, 1, 0))
}

aij <- function (i, j, n) {
  df <- ei(i, n) - ei(j, n)
  return(outer(df, df))
}

flfMDS <- function (delta, w, x) {
  n <- nrow(x)
  p <- ncol(x)
  s <- matrix(0, n * p, n * p)
  t <- matrix(0, n * p, n * p)
  u <- matrix(0, n * p, n * p)
  xx <- as.vector(x)
  nn <- diag(p)
  d <- as.matrix(dist(x))
  for (i in 2:n) {
    for (j in 1:(i - 1)) {
      a <- kronecker(nn, aij(i, j, n))
      y <- drop(a %*% xx)
      v <- outer(y, y)
      if (w[i, j] == 0) {
        next
      }
      s <- s + (w[i, j] / delta[i, j] ^ 2) * a
      t <- t + (w[i, j] / d[i, j] ^ 2) * a
      u <- u + (w[i, j] / d[i, j] ^ 2) * v
    }
  }
  g <- drop((s - t) %*% xx)
  h <- s - t + 2 * u
  return(list(s = s, t = t, u = u, g = g, h = h))
}