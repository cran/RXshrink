\name{RXshrink-package}
\alias{RXshrink-package}
\docType{package}
\title{Maximum Likelihood Shrinkage via Generalized Ridge or Least Angle Regression}
\description{The functions in this package augment the basic calculations of Generalized
  Ridge and Least Angle Regression with important visualization tools.  Specifically, they display
  TRACEs of normal-distribution-theory Maximum Likelihood estimates of the key quantities that
  completely characterize the effects of shrinkage on the MSE Risk of fitted coefficients.
}
\details{
  \tabular{ll}{
  Package: \tab RXshrink\cr
  Type: \tab Package\cr
  Version: \tab 1.0-5\cr
  Date: \tab 2010-09-10\cr
  License: \tab GNU GENERAL PUBLIC LICENSE, Version 2, June 1991\cr
  }
  
  RXridge() calculates and displays TRACEs for the \bold{Q}-shaped shrinkage path, including
  the \bold{M}-extent of shrinkage along that path, that are most likely under normal
  distribution theory to yield optimal reducions in MSE Risk.
  
  When regression parameters have specified, KNOWN numerical values, RXtrisk() calculates and
  displays the corresponding True MSE Risk profiles and RXtsimu() first simulates Y-outcome data 
  then calculates true Squared Error Losss associated with \bold{Q}-shape shrinkage.  
  
  RXlarlso() calls the Efron/Hastie lars() R-function to perform Least Angle Regression then
  augments these calculations with Maximum Likelihood TRACE displays like those of RXridge().
  
  RXuclars() applies Least Angle Regression to the uncorrelated components of a possibly
  ill-conditioned set of X-variables using a closed-form expression for the lars/lasso
  shrinkage delta factors that exits in this special case.
}
\author{Bob Obenchain <wizbob@att.net>}
\references{
  Efron B, Hastie T, Johnstone I, Tibshirani R. (2004) Least angle regression.
  \emph{Ann. Statis.} \bold{32}, 407-499.

  Goldstein M, Smith AFM. (1974) Ridge-type estimators for regression analysis.
  \emph{J. Roy. Stat. Soc. B} \bold{36}, 284-291. (2-parameter shrinkage family.)

  Obenchain RL. (2005)
  \emph{Shrinkage Regression: ridge, BLUP, Bayes, spline and Stein}.  Electronic
  book-in-progress (200+ pages.) \url{http://members.iquest.net/~softrx/}.

  Obenchain RL. (2010) \emph{RXshrink-R.PDF}  ../R\_HOME/library/RXshrink/doc
}
\keyword{package}
\examples{
  demo(longley2)
}