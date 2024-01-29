grades<-c("A", "B", "C", "D", "F")
g1<-ifelse(outer(johnson[,1], grades, "=="),1,0)
g2<-ifelse(outer(johnson[,2], grades, "=="),1,0)
g3<-ifelse(outer(johnson[,3], grades, "=="),1,0)
g4<-ifelse(outer(johnson[,4], grades, "=="),1,0)
gjohnson <- cbind(g1, g2, g3, g4)
cjohnson <- crossprod(gjohnson)
djohnson <- diag(cjohnson)
ejohnson <- cjohnson / sqrt (outer(djohnson, djohnson))
vjohnson <- eigen(ejohnson)
ljohnson <- vjohnson$values / 4.0
xjohnson <- (1.0 / sqrt(djohnson)) * vjohnson$vectors
hjohnson <- DMCA(cjohnson, c(5,5,5,5), verbose = TRUE)
yjohnson <- (1.0 / sqrt(djohnson)) * hjohnson$kpl
mjohnson <- diag(hjohnson$lpkbkpl) / 4
rjohnson <- crossprod (vjohnson$vectors, hjohnson$kpl)


