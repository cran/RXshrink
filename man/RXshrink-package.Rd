\name{RXshrink-package}
\alias{RXshrink-package}
\docType{package}
\title{Maximum Likelihood Shrinkage using Generalized Ridge or Least Angle Regression Methods}
\description{The functions in this package augment basic calculations of Generalized Ridge and
  Least Angle Regression with important visualization tools.  Specifically, they display TRACEs
  of (biased) estimates for 5 KEY quantities [fitted coefficients, MSE risks, excess eigenvalues,
  inferior direction cosines and shrinkage factors] that completely characterize the effects of
  shrinkage along PATHS with 2 or more parameters (M-extent plus Q-shape or Unrestricted). All
  PATHS start at the Ordinary Least-Squares estimate [M = 0] and end at the origin, (0, 0, ..., 0),
  where all coefficient estimates have been shrunken to zero [M = rank(X).] Three different types
  of Likelihood of minimal MSE risk (Classical Normal-Theory, Empirical Bayes, and Random
  Coefficients) can be monitored to suggest an optimal M-extent of shrinkage.  
}
\details{
  \tabular{ll}{
  Package: \tab RXshrink\cr
  Type: \tab Package\cr
  Version: \tab 1.3.1\cr
  Date: \tab 2020-01-20\cr
  License: \tab GNU GENERAL PUBLIC LICENSE, Version 2, June 1991\cr
  }
  
  qm.ridge() calculates and displays TRACEs for the 2-parameter PATH (i.e. \bold{M}-extent of
  shrinkage and \bold{Q}-shape) that is most likely under Normal-distribution theory to yield
  optimal reductions in MSE Risk.
  
  unr.ridge() calculates and displays TRACEs for a "new" 1-parameter PATH that passes through the
  Unrestricted (p-parameter) Coefficient Estimates most likely to yield optimal MSE reductions.
  
  When regression parameters have specified (KNOWN) numerical values, true.risk() calculates and
  displays the corresponding \bold{True} MSE Risk profiles. Similarly, true.simu() simulates
  Y-outcome data and calculates True Squared Error LOSSES associated with \bold{Q}-shape shrinkage.  
  
  aug.lars() augments the Efron/Hastie lars() R-function to perform Least Angle Regression with
  MSE risk canculations and Maximum Likelihood TRACE displays ...like those of qm.ridge() and
  unr.ridge().
  
  uc.lars() applies Least Angle Regression methods to the Uncorrelated Components of a possibly
  ill-conditioned set of X-variables. Calculations use a closed-form expression for lars/lasso
  shrinkage delta-factors that apply because NO Ill-Conditioning is present in these "uc" cases.
  
  correct.signs() displays the Normal-theory maximum likelihood estimate of the regression
  coefficient vector that minimizes MSE Risk in the UNKNOWN direction of p-space PARALLEL to
  the true Beta vector. This estimate corrects "wrong-sign" problems in the sense that its
  coefficients have the same relative magnitudes and numerical signs as those of the "Correlation
  Form" of the X'y vector.
}
\author{Bob Obenchain <wizbob@att.net>}
\references{
  Efron B, Hastie T, Johnstone I, Tibshirani R. (2004) Least angle regression.
  \emph{Annals of Statistics} \bold{32}, 407-499.

  Goldstein M, Smith AFM. (1974) Ridge-type estimators for regression analysis.
  \emph{J. Roy. Stat. Soc. B} \bold{36}, 284-291. (The 2-parameter shrinkage family.)

  Obenchain RL. (2005) \emph{Shrinkage Regression: ridge, BLUP, Bayes, spline and
  Stein}.  Electronic book-in-progress (185+ pages.) \url{http://localcontrolstatistics.org}

  Obenchain RL. (2018) \emph{RXshrink_in_R.PDF}  RXshrink package vignette-like file.
  \url{http://localcontrolstatistics.org}
  
  Obenchain RL. (2019) Maximum Likelihood Shrinkage in Regression. Submitted to
  \emph{Technometrics}. [The "new" unr.ridge() PATH is defined in equation (7.1).]
}
\keyword{package}
\examples{
  demo(longley2)
}