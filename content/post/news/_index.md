---
cms_exclude: true
header:
  caption: ""
  image: ""
title: "News"
date: 2026-04-14
view: 1
---

2026-04-14

Over the last couple of years I have written a large number of R codes for
matrix approximation, multidimensional scaling, and multivariate analysis.
The general approach is to write a stand-alone R program and in addition an R driver that
uses .C() to call a shared library, which is by itself a stand-alone
C program. The main recent projects are C versions of the smacof MDS
majorization (MM) method with extensions to Sammon mapping, Elastic Scaling, 
strain, fStress, rStress -- all written in both R and C, all with both
metric and non-metric versions.

The papers/manuals/vignettes, largely unpublished and often unfinished, are at
https://jansweb.netlify.app/publication and the corresponding C, h, R, pdf,
qmd, tex, bib files are at https://github.com/deleeuw?tab=repositories. I am not
interested in transforming these papers/codes to full-fledged publications,
but in the current form they maybe useful to some. Everything in the repositories
is open source with a CC0 license, so the files can be used in whatever way you see fit. 
Attribution is appreciated, but not required.
On the other hand, all suggestions for improvement of paper and code are welcome.
