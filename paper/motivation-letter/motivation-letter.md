---
output:
  pdf_document:
    latex_engine: pdflatex
fontsize: 12pt
geometry: margin=1in
---

\thispagestyle{empty}
\today

Editors\
The R Journal

Dear Editors,

Please consider our article, "statnnet: Statistical Selection and Interpretation
for Neural Networks Fitted with nnet", for publication in *The R Journal*.

The article presents the `statnnet` package, which supplies a coherent workflow
for selecting and interpreting a deliberately restricted class of shallow
neural networks fitted with `nnet`. The package selects input terms and
hidden-layer size using an information criterion, retains the evaluated
candidate path and convergence diagnostics, reconstructs the fitted objective
and curvature, and provides covariate-level Wald summaries, predictive
finite-difference effects, partial-dependence contrasts, diagnostic plots, and
prediction methods through a familiar R modelling interface.

We believe the article is relevant to R Journal readers because `nnet` remains
a widely used and stable R implementation of shallow neural networks, while
its standard output is primarily directed toward fitting and prediction. The
contribution of `statnnet` is intentionally narrower than a general neural-
network framework: selection and interpretation use `nnet` as the sole fitting
engine and preserve its exact parameter ordering, objective, Hessian,
convergence information, and fitted object. The article demonstrates the
complete workflow using the open Auto MPG data, places the contribution in
context with model-agnostic tools and neural-network packages, and states
clearly where their estimands and inferential mechanisms differ.

The manuscript and its examples are fully reproducible from the supplied R
Markdown source, scripts, and package data. The package source, tests,
documentation, and development history are maintained publicly at
<https://github.com/andrew-mcinerney/statnnet>. We understand that an R Journal
package submission requires CRAN availability and will submit this manuscript
only after the corresponding `statnnet` release has been accepted on CRAN.

This manuscript is original, has not been published previously, and is not
under consideration by another journal. Both authors have approved its
submission.

Yours sincerely,

Andrew McInerney\
Department of Mathematics and Statistics\
University of Limerick\
Limerick, Ireland\
andrew.mcinerney@ul.ie

Kevin Burke\
Department of Mathematics and Statistics\
University of Limerick\
Limerick, Ireland\
kevin.burke@ul.ie
