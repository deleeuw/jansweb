---
cms_exclude: true
header:
  caption: ""
  image: ""
title: 'On Computational Statistics'
date: 2021-11-24
view: 1
---
&#9658; **Computational Statistics** is the discipline that proposes **techniques**, provides descriptions and implementations of these techniques, and compares and analyzes their properties using real and artificial data. 

&#9658; A **technique** is a map of **carrier space** into **target space**.
Both spaces are often, but not necessarily, real linear spaces.

&#9658;  Computational Statistics is not **data theory**. Design, collection, coding, and cleaning of data prepare the input fot the technique. Our job has not started yet.

&#9658; Computational Statistics is not **data analysis**. Inference, induction, probabilistic reasoning, and causal analysis use the output of the technique. Our job is already done.

&#9658; **Data space** is a finite subset of carrier space. It is the set of codings of all possible outcomes of a given experiment or other data collection effort. For example: coded versions of all possible results of a multiple choice (attitude/intelligence/personality) test on a number of subjects, outcomes of a clinical trial on a number of subjects, outcomes of a survey of a number of voters, essays in an essay test of a number of students, event histories in a longitudinal study of a number of individuals, heights and weights of a number of soldiers, dissimilarity matrices of a number of colors.

&#9658; **Representation space** is a finite subset of target space. It is the set of all outputs of the technique applied to the data space. For example: a number of correlation matrices, cross tables, graphs and plots, regression coefficients.

&#9658; Continuity and smoothness, defined using suitable topologies on carrier and target space, are among the most interesting aspects of statistical techniques. They define the **stability** of the technique under various types of perturbations.

&#9658; A **package** is a technique, including a description of carrier/data and target/representation space, plus an algorithm, plus an implementation, plus a manual, plus examples, plus a stability analysis, plus a license. A package is **incomplete** if it does not all these components, and **inaccessible** if its components are not open source and open access.

&#9658; For more details, see Albert Gifi, *Nonlinear Multivariate Analysis*,
John Wiley & Sons, Chichester, 1990, Chapter 1, Section 5.


