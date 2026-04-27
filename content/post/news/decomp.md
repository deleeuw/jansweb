---
cms_exclude: true
header:
  caption: ""
  image: ""
title: "Matrix Decomposition (2026-04-24)"
date: 2026-04-14
view: 1
---
<https://github.com/deleeuw/decomp> has R code and paper for matrix 
decomposition of Y by X in the Kronecker norm
tr R(Y-X)C(Y-X)', with R and C positive semi-definite matrices. 
Specifically we look at X of the form c + ae' + eb' + D. The four
components define four Kronecker-orthogonal subspaces of matrix space.
We project Y on these four subspaces, using the Kronecker metric.
Generalizations to multidimensional arrays, to external variables
different from vectors with all elements equal to one, and to
the matrix variate normal distribution are 
discussed.