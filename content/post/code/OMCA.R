OMCA <-
  function(cbrt,
           ncat) {
    cbrt <- as.matrix(cbrt)                  ## handle dataframes
    nvar <- length(ncat)                     ## number of variables
    ntot <-
      sum(cbrt) / (nvar * nvar)        ## number of observations
    cbrt <- cbrt / ntot                      ## relative frequencies
    msum <- sum(ncat)                        ## total categories
    dbrt <- diag(cbrt)                       ## vector of marginals
    ll <- kk <- ww <- matrix(0, msum, msum)  ## reserve memory
    chi <- matrix(0, nvar, nvar)             ## reserve memory
    kup <- cumsum(ncat)
    klw <- kup - (ncat - 1)
    ind <- lapply(1:nvar, function(i)
      klw[i]:kup[i])
    ebrt <- cbrt / sqrt(outer(dbrt, dbrt))
    for (i in 1:nvar) {
      kk[ind[[i]], ind[[i]]] <- orthopol(1:ncat[i], dbrt[ind[[i]]])
    }
    kek <- crossprod (kk, cbrt %*% kk)
    for (i in 1:nvar) {
      for (j in 1:nvar) {
        ww[ind[[i]], ind[[j]]] <-
          ifelse(outer(1:ncat[i], 1:ncat[j], "=="), 1, 0)
        chi[i, j] <-
          ntot * sum (kek[ind[[i]], ind[[j]]][-1, -1] ^ 2)
      }
    }
    chitot <- sum(chi) - sum(diag(chi))
    ossq <-
      (sum(ww * kek ^ 2) - (msum + nvar * (nvar - 1))) / (chitot / ntot)
    kl <- unlist(sapply(ncat, function(i)
      1:i))
    pp <- ifelse(outer(1:msum, order(kl), "=="), 1, 0)
    pkekp <- t(pp) %*% kek %*% pp
    kp <- kk %*% pp
    km <- as.vector(table(kl))
    nm <- length(km)
    klw <- 1 + cumsum(c(0, km))[1:nm]
    kup <- cumsum(km)
    chidca <- matrix(0, nm, nm)
    for (i in 1:nm) {
      ind <- klw[i]:kup[i]
      ll[ind, ind] <- eigen(pkekp[ind, ind])$vectors
      for (j in 1:nm) {
        jnd <- klw[j]:kup[j]
        chidca[i, j] <- sum(pkekp[ind, jnd] ^ 2)
        if (i == j) {
          chidca[i, j] <- chidca[i, j] - length(ind)
        }
      }
    }
    chidca <- ntot * chidca[-1,-1]
    chiper <- chidca / chitot
    lpkekpl <- t(ll) %*% pkekp %*% ll
    kpl <- kp %*% ll
    return(
      list(
        kek = kek,
        pkekp = pkekp,
        lpkekpl = lpkekpl,
        k = kk,
        p = pp,
        l = ll,
        kp = kp,
        kpl = kpl,
        chisquares = chi,
        chipartition = chidca,
        chipercentages = chiper,
        func = ossq
      )
    )
  }

orthopol <- function (x, w) {
  n <- length(x)
  z <- outer(x, 0:(n - 1), "^")
  for (i in 1:n) {
    if (i > 1) {
      for (j in 1:(i - 1)) {
        z[, i] <- z[, i] - sum(w * z[, i] * z[, j]) * z[, j]
      }
    }
    z[, i] <- z[, i] / sqrt(sum(w * z[, i] * z[, i]))
  }
  return(z)
}
