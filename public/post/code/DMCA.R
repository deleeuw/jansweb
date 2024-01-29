DMCA <-
  function(cbrt,
           ncat,
           eps = 1e-8,
           itmax = 500,
           verbose = TRUE,
           vectors = TRUE) {
    cbrt <- as.matrix(cbrt)                  ## handle dataframes
    nvar <- length(ncat)                     ## number of variables
    ntot <- sum(cbrt) / (nvar * nvar)        ## number of observations
    cbrt <- cbrt / ntot                      ## relative frequencies
    msum <- sum(ncat)                        ## total categories
    dbrt <- diag(cbrt)                       ## vector of marginals
    ll <- kk <- ww <- matrix(0, msum, msum)  ## reserve memory
    chi <- matrix(0, nvar, nvar)             ## reserve memory
    itel <- 1
    kup <- cumsum(ncat)
    klw <- kup - (ncat - 1)
    ind <- lapply(1:nvar, function(i) 
      klw[i]:kup[i])
    ebrt <- cbrt / sqrt(outer(dbrt, dbrt))
    for (i in 1:nvar) {
      kk[ind[[i]], ind[[i]]] <- svd(ebrt[ind[[i]], ])$u
    }
    kek <- crossprod (kk, ebrt %*% kk)
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
    repeat {
      for (l in 1:nvar) {
        if (ncat[l] == 2)
          next()
        li <- ind[[l]]
        for (i in (klw[l] + 1):(kup[l] - 1))
          for (j in (i + 1):kup[l]) {
            bi <- kek[i, -li]
            bj <- kek[j, -li]
            wi <- ww[i, -li]
            wj <- ww[j, -li]
            acc <- sum(wi * bi ^ 2) + sum(wj * bj ^ 2)
            acs <- sum((wi - wj) * bi * bj)
            ass <- sum(wi * bj ^ 2) + sum(wj * bi ^ 2)
            u <-
              eigen(matrix(c(acc, acs, acs, ass), 2, 2))$vectors[, 1]
            c <- u[1]
            s <- u[2]
            kek[-li, i] <- kek[i, -li] <- c * bi + s * bj
            kek[-li, j] <- kek[j, -li] <- c * bj - s * bi
            if (vectors) {
              ki <- kk[li, i]
              kj <- kk[li, j]
              kk[li, i] <- c * ki + s * kj
              kk[li, j] <- c * kj - s * ki
            }
          }
      }
      nssq <-
        (sum(ww * kek ^ 2) - (msum + nvar * (nvar - 1))) / (chitot / ntot)
      if (verbose)
        cat(
          "Iteration ",
          formatC(itel, digits = 4),
          "ossq ",
          formatC(ossq, digits = 10, width = 12),
          "nssq ",
          formatC(nssq, digits = 10, width = 12),
          "\n"
        )
      if (((nssq - ossq) < eps) || (itel == itmax))
        break()
      itel <- itel + 1
      ossq <- nssq
    }
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
        itel = itel,
        func = nssq
      )
    )
  }
