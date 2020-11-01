\name{RXshrink-package}
\alias{RXshrink-package}
\docType{package}
\title{Maximum Likelihood (ML) Shrinkage using Generalized Ridge or Least Angle Regression Methods}
\description{The functions in this package augment basic calculations of Generalized Ridge and
  Least Angle Regression with important visualization tools.  Specifically, they display TRACEs
  of (biased) estimates for 5 KEY quantities [fitted coefficients, relative MSE, excess eigenvalues,
  inferior direction cosines and shrinkage factors] that completely characterize the effects of
  shrinkage along Paths with 2 or more parameters: (m-Extent plus q-Shape or Unrestricted). All
  Paths start at the Ordinary Least-Squares estimate [m = 0] and end at the origin, (0, 0, ..., 0),
  where all coefficient estimates have been shrunken to zero [m = p = rank(X).] Three different
  types of Likelihood of minimal MSE risk (Classical Normal-Theory, Empirical Bayes, and Random
  Coefficients) can be monitored to suggest an optimal m-Extent of shrinkage for the given
  matrix of X-variables and the observed y-outcome vector.  
}
\details{
  \tabular{ll}{
  Package: \tab RXshrink\cr
  Type: \tab Package\cr
  Version: \tab 1.4.3\cr
  Date: \tab 2020-11-01\cr
  License: \tab GNU GENERAL PUBLIC LICENSE, Version 2, June 1991\cr
  }
  
  unr.ridge() calculates and displays TRACEs for an Unrestricted (p-parameter) shrinkage Path
  that always passes through the Beta coefficient point-estimate that is most likely to achieve
  optimal MSE risk reductions under Normal distribution-theory.
  
  MLboot(), MLcalc() and MLhist() support use of Bootstrap resampling to study the bias and MSE
  risk characeristics of non-linear (unrestricted) Generalized Ridge Regression (GRR) estimators.
  
  When true regression parameters have user-specified (KNOWN) numerical values, MLtrue() uses this
  information and generates a new data.frame that contains a y-outcome vector of the expected form
  with "disturbance" terms that are IID Normal errors-in-measurement. Arguments to MLtrue() must
  include the "formula" for a desired linear model and a data.frame containing the specified
  X-variables.
  
  qm.ridge() calculates and displays TRACEs for traditional Paths defined by 2-parameters:
  \bold{q}-Shape and \bold{m}-Extent of Shrinkage. The restricted Path of most likely
  \bold{q}-Shape is found via search on a lattice of 21 values within [-5,+5]. Default lattices for
  both q and m searches are easy to modify using the qmax, qmin, nq and steps arguments to qm.ridge().
  The "ordinary" ridge Path of Hoerl and Kennard has \bold{q}-Shape = 0, while "uniform" shrinkage
  corresponds to \bold{q}-Shape = +1. None of these qm-Paths generally achieve overall minimum MSE
  risk when p > 2 because they restrict attention to Shrinkage-factors that are "monotome" (increasing
  or decreasing).  
  
  aug.lars() augments the Efron-Hastie lars() R-function to perform Least Angle Regression with
  MSE risk calculations and Maximum Likelihood TRACE displays ...like those of qm.ridge() and
  unr.ridge().
  
  uc.lars() applies Least Angle Regression methods to the \bold{Uncorrelated Components} of a possibly
  ill-conditioned set of x-variables. Calculations use a closed-form expression for lars/lasso
  shrinkage delta-factors that apply because NO Ill-Conditioning is present in these "uc" cases.
  
  correct.signs() displays the Normal-theory maximum likelihood estimate of the regression
  coefficient vector that minimizes MSE Risk in the UNKNOWN direction of p-space PARALLEL to
  the true Beta vector. This estimate corrects "wrong-sign" problems in the sense that its
  coefficients have the same relative magnitudes and numerical signs as those of the "Correlation
  Form" of the X'y vector.
  
  RXpredict() makes predictions (i.e. computes "fitted.values") for all 6 types of RXshrink
  estimation ...either at a user-specified m-Extent of Shrinkage or at the Normal-theory "minMSE"
  m-Extent.
}
\author{Bob Obenchain <wizbob@att.net>}
\references{
  Efron B, Hastie T, Johnstone I, Tibshirani R. (2003) Least angle regression.
  \emph{Annals of Statistics} \bold{32}, 407-499.

  Goldstein M, Smith AFM. (1974) Ridge-type estimators for regression analysis.
  \emph{J. Roy. Stat. Soc. B} \bold{36}, 284-291. (The 2-parameter shrinkage family.)

  Obenchain RL. (1975) Ridge Analysis Following a Preliminary Test of the Shrunken
  Hypothesis. \emph{Technometrics} \bold{17}, 431-441. <doi:10.1080/00401706.1975.10489369>

  Obenchain RL. (1978) Good and Optimal Ridge Estimators. \emph{Annals of Statistics}
  \bold{6}, 1111-1121. <doi:10.1214/aos/1176344314>
  
  Obenchain RL. (2020) Ridge TRACE Diagnostics. \url{https://arxiv.org/abs/2005.14291}
  
  Obenchain RL. (2020) \emph{RXshrink_in_R.PDF}  RXshrink package vignette-like document.
  \url{http://localcontrolstatistics.org}
}
\keyword{package}
\examples{
  demo(longley2)
}
