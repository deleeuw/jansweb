data(bfi, package = "psychTools")
bfi <- bfi[,1:25]
bfi <- bfi[which(rowSums(is.na(bfi)) == 0), ]
nbfi <- nrow(bfi)
mbfi <- ncol(bfi)
kbfi <- apply(bfi, 2, max)
gbfi <- NULL
for (j in 1:mbfi) {
  gbfi <- cbind(gbfi, ifelse(outer(bfi[, j], unique(bfi[, j]), "=="), 1, 0))
}
cbfi <- crossprod(gbfi)
dbfi <- diag(cbfi)
ebfi <- cbfi / sqrt(outer(dbfi, dbfi))
sbfi <- eigen(ebfi)
lbfi <- sbfi$values / mbfi
vbfi <- sbfi$vectors
ybfi <- sqrt(1 / dbfi) * vbfi
abfi <- lbfi[-c(1,127:150)]
ubfi <- (mbfi * abfi - 1) ^ 2
tbfi <- nbfi * sum(ubfi)
hbfi <- DCA(cbfi, rep(6, 25), verbose = TRUE)
wbfi <- diag(hbfi$lpkekpl)[-(1:25)] / mbfi
zbfi <- sqrt(1 / dbfi) * hbfi$kpl
rbfi <- crossprod(vbfi, hbfi$kpl)
ibfi <- outer(1:25, 1:25, ">")
bfir1 <- hbfi$pkekp[26:50,26:50][ibfi]
bfir2 <- hbfi$pkekp[51:75,51:75][ibfi]
bfir3 <- hbfi$pkekp[76:100,76:100][ibfi]
bfir4 <- hbfi$pkekp[101:125,101:125][ibfi]
bfir5 <- hbfi$pkekp[126:150,126:150][ibfi]


  
