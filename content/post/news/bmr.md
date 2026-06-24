---
cms_exclude: true
header:
  caption: ""
  image: ""
title: "Binary Monotone Regression (2026-06-24)"
date: 2026-06-24
view: 1
---

In the file monotone.R in the github repository

<https://github.com/deleeuw/voronoi>

there is a new routine binaryMonotoneRegression()
for the special case of least squares monotone
regression in which the partial order is defined
by a partition of the target into two sets.
All elements in  the first set are required to be
smaller than or equal to all elements in the second set.
Previously I used Kruskal's primary approach
to ties to handle this case, which requires
a lot of sorting and index manipulation in each
iteration.

The new algorithm reduces the monotone regression
problem to minimizing a differentiable and piecewise
quadratic convex function of a single variable. We
sort the target and use the fact that the function is quadratic
between each interval of successive values in the sorted
vector. We start at one end of the scale, evaluate the
derivative at successive values until it changes sign, 
and then interpolate linearly to find the zero.

The new function is used in the voronoiHomogeneityAnalysis(),
voronoiCentroidAnalysis(), and voronoiSpericalAnalysis()
functions for unfolding categorical data that can be
found in the same github repository.



