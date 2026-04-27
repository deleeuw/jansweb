---
cms_exclude: true
header:
  caption: ""
  image: ""
title: "Matrix Approximation (2026-04-22)"
date: 2026-04-22
view: 1
---
"Matrix Approximation with Kronecker and Precision Weighting" by
De Leeuw and Graffelman is at 

<https://jansweb.netlify.app/publication/deleeuw-graffelman-e-26-d/>

and the files are at 

<https://github.com/deleeuw/triple>

There are two R programs (triple.R and nested.R) for minimizing the
loss function trace(R(W*(Y-X))C(W*(Y-X))') over a set in matrix space
(for example the matrices of rank ≤ p, or the Hankel matrices, or
whatever). Here Y is the matrix to be approximated, X is the
approximating matrix, R and C are psd weight matrices, W is a non-negative
matrix of weights, and * is Hadamard (elementwise) multiplication.
Both programs use majorization (MM), the function nested() uses
a nice nested majorization. The type of approximation wanted is
a parameter, the name of a function for unweighted least squares
approximation (defaults to eckart-young()). The file auxiliary.R
contains nine different approximation routines that can be used
in both nested() and triple().

As always, everything is public domain, and both paper and code may
change over time.