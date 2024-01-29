matrixPrint <- function (x,
                         digits = 6,
                         width = 8,
                         format = "f",
                         flag = "+") {
  print (noquote (
    formatC (
      x,
      digits = digits,
      width = width,
      format = format,
      flag = flag
    )
  ))
}

superMatrixPrint <- function(x,
                             rows,
                             columns,
                             digits = 6,
                             width = 8,
                             format = "f",
                             flag = "+",
                             vsep = "|",
                             hsep = "",
                             space = 2) {
  nr <- nrow (x)
  nc <- ncol (x)
  char <- matrix("?", nr, nc)
  for (i in 1:nr) {
    for (j in 1:nc) {
      char[i, j] <- formatC (
        x[i, j],
        digits = digits,
        width = width,
        format = format,
        flag = flag
      )
    }
  }
  mr <- length(rows)
  mc <- length(columns)
  rup <- cumsum(rows)
  rlw <- rup - (rows - 1)
  rnd <- lapply(1:mr, function(i)
    rlw[i]:rup[i])
  cup <- cumsum(columns)
  clw <- cup - (columns - 1)
  cnd <- lapply(1:mc, function(i)
    clw[i]:cup[i])
  hhsep <- NULL
  kkk <- (width + 2) * nc + 2 * mc + 1
  for (k in 1:kkk) {
    hhsep <- paste(hhsep, hsep, sep = "")
  }
  for (i in 1:mr) {
    work <- matrix(vsep, rows[i], 1)
    for (j in 1:mc) {
      work <- cbind(work, char[rnd[[i]], cnd[[j]]])
      work <- cbind(work, matrix(vsep, rows[i], 1))
    }
    cat(rep(hsep, kkk), sep="")
    prmatrix(
      work,
      quote = FALSE,
      rowlab = rep("", nrow(work)),
      collab = rep("", ncol(work))
    )
  }
  cat(rep(hsep, kkk), sep="")
}
