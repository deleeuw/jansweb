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
the n objects into a set of k<sub>j</sub> possible values, the
*categories* of the variable. The number of categories can be
finite or infinite and the categories can be nominal, ordinal, or
numerical. In actual data analysis the number of categories
that are actually present in the data will necessarily be 
finite.

HA maps objects and categories into a p-dimensional Euclidean
space. The numerical output is an n x p matrix X of *object scores*
and m matrices Y<sub>j</sub>, of dimension k<sub>j</sub> x p, called *category 
quantifications* for nominal and ordinal variables and *category transformations*
for numerical variables. In HA, in the
form in which it is implemented in Gifi (1990), we require
that X is column-centered with X'X = nI. The category 
quantifications Y<sub>j</sub>
are the centroids of the scores of the
objects in the category (we also say that objects "belong to" a category
or are "in" a category)

The most interesting graphical outputs of HA are the *star plots*, which
are joint plots 
of objects scores and category quantifications/transformations. 
There is one star plot for each variable. In the star plot we
draw stars by connecting each object with the categories the
object belongs to. The resulting k<sub>j</sub> graphs are star-shaped, because
the category point is the centroid of "its" object points.

In the paper we first reformulate HA as a form of non-metric multidimensional
scaling, using Kruskal's stress and requiring that the distances between
object-points and the category-points the objects belong to are zero.
We then relax this to the requirement that the distance of an object-point
to the point of the category the object belongs to must be less than or equal to
the distance to the other category points of the same variable. We adapt the
monotone regression routines and the normalization of the solution to
this requirement. The classical HA solution is used as an initial
estimate of the smacof iterations.

Geometrically the category-points of a variable define a partition of the space into Voronoi regions. Our new requirement is that all object-points 
are in the Voronoi region of the category they belong to. For each variable the star 
graphs are replaced by regions partitioning the space. There is no centroid
relation between category-points and object-points are more -- the object points
can be anywhere in the Voronoi region of the category they belong to.

The paper is at

<https://jansweb.netlify.app/publication/deleeuw-e-26-f/>

and the files are at 

<https://github.com/deleeuw/voronoi>


