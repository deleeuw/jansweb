---
cms_exclude: true
header:
  caption: ""
  image: ""
title: "Voronoi Homogeneity Analysis (2026-05-09)"
date: 2026-05-09
view: 1
---
In Homogeneity Analysis (HA), better known as Multiple Correpondence
Analysis (MCA), there are n *objects* and m *variables* that measure
the objects. By "measure" we mean that each variable j maps
the n objects into a set of k(j) possible values, the
*categories* of the variable. The number of categories can be
finite or infinite and the categories can be nominal, ordinal, or
numerical. In actual data analysis the number of categories
that are actually present in the data will necessarily be 
finite.

HA maps objects and categories into a p-dimensional Euclidean
space. The numerical output is an n x p matrix X of *object scores*
and m matrices Y(j), of dimension k(j) x p, of *category 
quantifications* for nominal and ordinal variables or *category transformations*
for numerical variables. In HA, in the
form in which it is implemented in Gifi (1990), we require
that X is column-centered with X'X = nI. The category 
quatifications are the centroids of the scores of the
objects in the category.

The most interesting graphical outputs of HA are the *star plots*, which
are joint plots 
of objects scores and category quantifications/transformations. 
There is one star plot for each variable. In the star plot we
draw stars by connecting each object with the categories the
object is in. The resulting k(j) graphs are star shapes, because
the category point is the centroid of "its" object points.