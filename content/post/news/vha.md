---
cms_exclude: true
header:
  caption: ""
  image: ""
title: "More Voronoi Homogeneity Analysis (2026-06-04)"
date: 2026-06-04
view: 1
---

My brief introductory remarks for the April 10 workshop
 
<https://deleeuwworkshop.stat.ucla.edu/#schedule>

had the title "A Tale of Two Loss Functions". The
two loss functions were, of course, the Gifi meet-loss
of homogeneity analysis (MCA) and the stress of smacof.
Patrick's presentation at the workshop uses the same
two chapters to discuss my work in "dimension reduction".

I am currently working on the paper "Voronoi Homogeneity Analysis",
latest version at 

<https://github.com/deleeuw/voronoi>

In the
introductory section there is a new derivation of the fact
that homogeneity analysis is a special case of non-metric
smacof (and that consequently the Gifi system can be embedded
in smacof with restrictions). More specifically for m
variables MCA consists of m linked nonmetric unfolding problems,
which are linked because they have the same X but
different Y_j. Thus it is now "A Tale of one Loss Function".

Voronoi MCA relaxes the ordinal constraints on the disparities.
As a consequence, as in unfolding, there are the usual
degenerate solutions with perfect fit. But because
of the weak ordinal constraints other perfect fit solutions exist,
which are not necessarily degenerate (or only partly degenerate).

The paper has software and examples to illustrate the 
workings of the algorithm and the influence of various
constraints (some of which emulate MCA properties).
The appendix introduces MALS, which is Majorized
Alternating Least Squares.
