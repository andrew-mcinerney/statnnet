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

Please consider our article, "statnnet: Regression-Style Statistical Summaries
for Neural Networks Fitted with nnet", for publication in *The R Journal*.

The article presents the `statnnet` package, which augments a deliberately
restricted class of fitted `nnet` models with statistical information that is
not supplied by the fitting package itself. The package reconstructs the fitted
objective and curvature, reports individual- and covariate-level Wald
summaries, distinguishes observed-row predictive finite differences from
partial-dependence contrasts, supplies pointwise uncertainty intervals, and
provides diagnostic, plotting, and prediction methods through a familiar R
modelling interface.

We believe the article is relevant to R Journal readers because `nnet` remains
a widely used and stable R implementation of shallow neural networks, while
its standard output is primarily directed toward fitting and prediction. The
contribution of `statnnet` is intentionally narrower than a general neural-
network interpretation framework: it uses the exact parameter ordering,
objective, Hessian, convergence information, and fitted object produced by
`nnet`. The article places this contribution in context with model-agnostic
tools and with `NeuralNetTools`, `NeuralSens`, and `cito`, and states clearly
where their estimands and inferential mechanisms differ.

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
