JMCA <- function(cbrt, ncat, ndim = 2, itmax = 100, eps = 1e-10, verbose = TRUE) {
  cbrt <- as.matrix(cbrt)                 
  nvar <- length(ncat)
  msum <- sum(ncat)
  ntot <- sum(cbrt) / (nvar ^ 2)
  dbrt <- diag(cbrt)
  ebrt <- cbrt / sqrt(outer(dbrt, dbrt))
  sbrt <- eigen(ebrt)
  lbrt <- sbrt$values[-c(1, (msum - nvar + 2):msum)] / nvar
  vbrt <- sbrt$vectors[, -c(1, (msum - nvar + 2):msum)]
  qbrt <- lbrt / sum(lbrt)
  ubrt <- (nvar * lbrt - 1) ^ 2
  pbrt <- ubrt / sum(ubrt)
  tbrt <- ntot * sum(ubrt)
  ybrt <- sqrt(1 / dbrt) * vbrt
  return(list(
    lbrt = lbrt,
    vbrt = vbrt,
    ybrt = ybrt,
    tbrt = tbrt,
    pbrt = pbrt,
    qbrt = qbrt
  ))
}
