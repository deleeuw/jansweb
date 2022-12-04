cburt <- as.matrix(burt) / 100
dburt <- diag(cburt)
eburt <- cburt / (4 * sqrt (outer(dburt, dburt)))
vburt <- eigen(eburt)
lburt <- vburt$values
xburt <- (1.0 / sqrt(dburt)) * vburt$vectors
hburt <- DMCA(cburt, c(3,3,2,2), verbose = TRUE, eps = 1e-15)
yburt <- (1.0 / sqrt(dburt)) * hburt$kpl
mburt <- diag(hburt$lpkbkpl) / 4
rburt <- crossprod (vburt$vectors, hburt$kpl)

