data(galo, package = "homals")
galo <- galo[, -(ncol(galo))]
ngalo <- nrow(galo)
mgalo <- ncol(galo)
kgalo <- sapply(galo, function(x) length(unique(x)))
ggalo <- NULL
for (j in 1:4) {
  ggalo <- cbind(ggalo, ifelse (outer (galo[[j]], unique(galo[[j]]), "=="), 1, 0))
}
cgalo <- crossprod(ggalo)
dgalo <- diag(cgalo)
egalo <- cgalo / sqrt(outer(dgalo, dgalo))
sgalo <- eigen(egalo)
vgalo <- sgalo$vectors
lgalo <- sgalo$values / mgalo
ygalo <- (1 / sqrt(dgalo)) * vgalo
hgalo <- DMCA(cgalo, kgalo, itmax = 100, verbose = TRUE)
wgalo <- diag(hgalo$lpkekpl) / mgalo
zgalo <- (1 / sqrt(dgalo)) * hgalo$kpl
rgalo <- crossprod (vgalo, hgalo$kpl)
pgalo1 <- hgalo$pkekp[1:4, 1:4]
pgalo2 <- hgalo$pkekp[5:8, 5:8]
pgalo3 <- hgalo$pkekp[9:11, 9:11]
pgalo4 <- hgalo$pkekp[12:14, 12:14]
pgalo5 <- hgalo$pkekp[15:17, 15:17]
pgalo6 <- hgalo$pkekp[18:20, 18:20]
pgalo7 <- hgalo$pkekp[21:22, 21:22]
pgalo8 <- hgalo$pkekp[23:24, 23:24]
