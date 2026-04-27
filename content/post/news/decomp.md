---
cms_exclude: true
header:
  caption: ""
  image: ""
title: "Matrix Decomposition (2026-04-24)"
date: 2026-04-24
view: 1
---
<https://github.com/deleeuw/decomp> has R code and paper for matrix 
decomposition of X in the squared Kronecker norm
tr RXCX', with R and C positive semi-definite matrices. 
Specifically we look at decompositions of the form X = c + ae' + eb' + D. The four components are used to define four Kronecker-orthogonal subspaces of matrix space.
We project X on these four subspaces, using the Kronecker metric. We discuss

* applications to matrix approximation,
* generalizations to multidimensional arrays, 
* extensions to external or passive variables, and 
* ML estimation for the matrix variate normal distribution.