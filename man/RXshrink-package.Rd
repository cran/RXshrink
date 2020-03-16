\name{RXshrink-package}
\alias{RXshrink-package}
\docType{package}
\title{Maximum Likelihood Shrinkage using Generalized Ridge or Least Angle Regression Methods}
\description{The functions in this package augment basic calculations of Generalized Ridge and
  Least Angle Regression with important visualization tools.  Specifically, they display TRACEs
  of (biased) estimates for 5 KEY quantities [fitted coefficients, MSE risks, excess eigenvalues,
  inferior direction cosines and shrinkage factors] that completely characterize the effects of
  shrinkage along PATHS with 2 or more parameters (m-Extent plus q-Shape or Unrestricted). All
  PATHS start at the Ordinary Least-Squares estimate [m = 0] and end at the origin, (0, 0, ..., 0),
  where all coefficient estimates have been shrunken to zero [m = rank(X).] Three different types
  of Likelihood of minimal MSE risk (Classical Normal-Theory, Empirical Bayes, and Random
  Coefficients) can be monitored to suggest an optimal m-Extent of shrinkage.  
}
\details{
  \tabular{ll}{
  Package: \tab RXshrink\cr
  Type: \tab Package\cr
  Version: \tab 1.3.2\cr
  Date: \tab 2020-03-20\cr
  License: \tab GNU GENERAL PUBLIC LICENSE, Version 2, June 1991\cr
  }
  
  unr.ridge() calculates and displays TRACEs for an Unrestricted (p-parameter) PATH that passes
  through the Unrestricted Maximum Likelihood Coefficient point-Estimate that yields optimal
  MSE risk reductions under Normal distribution-theory.
  
  qm.ridge() calculates and displays TRACEs for PATHs defined by 2-parameters: \bold{q}-Shape
  and \bold{m}-Extent of shrinkage. The restricted PATH of most likely \bold{q}-Shape is found
  via search on a lattice of 21 values within [-5,+5]. Default lattices for q and m searches are
  easy to modify using the qmax, qmin, nq and steps arguments to qm.ridge(). The "ordinary" ridge
  PATH of Hoerl and Kennard has \bold{q}-Shape = 0, while "uniform" shrinkage corresponds to
  \bold{q}-Shape = +1.  
  
  When true regression parameters have user-specified (KNOWN) numerical values, true.risk()
  calculates and displays \bold{True} MSE Risk profiles. Similarly, true.simu() simulates
  Y-outcome data and calculates True Squared Error LOSSES associated with \bold{q}-Shape shrinkage.  
  
  aug.lars() augments the Efron-Hastie lars() R-function to perform Least Angle Regression with
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

  Obenchain RL. (2020) \emph{RXshrink_in_R.PDF}  RXshrink package vignette-like document.
  \url{http://localcontrolstatistics.org}
}
\keyword{package}
\examples{
  demo(longley2)
}