\name{RXshrink-internal}
\title{Internal RXshrink functions}
\alias{RXshrink-internal}
\alias{print.aug.lars}
\alias{print.correct.signs}
\alias{print.MLboot}
\alias{print.MLcalc}
\alias{print.MLhist}
\alias{print.MLtrue}
\alias{print.qm.ridge}
\alias{print.RXpredict}
\alias{print.uc.lars}
\alias{print.unr.ridge}
\alias{mstep}
\alias{kofm}
\alias{mofk}
\description{These are internal functions for printing and calculating quantities such
 as additive K-factors and multiplicative Delta-Factors that determine a specified
 m-Extent of Shrinkage. mstep() is called by qm.ridge(), true.risk(), true.simu(), aug.lars()
 and uc.lars(). kofm() and mofk() are called by unr.ridge(),
}
\usage{
  mstep(mobj, kinc, p, qp, eqm1)
  kofm(muobj, p, dMSE, delmax)
  mofk(p, k, dMSE)
}
\arguments{
  \item{mobj}{Objective value for the M-extent-of-shrinkage along a path of shape qp.}
  \item{kinc}{Previous additive eigenvalue inflation constant yielding a smaller mobj.}
  \item{p}{Number of elements in the regression coefficient vector.}
  \item{qp}{Shrinkage path \bold{Q}-shape parameter.}
  \item{eqm1}{Vector of X'X eigenvalues raised to the power (qp-1).}
  \item{muobj}{Objective value for the M-extent-of-shrinkage along the Unrestricted path.}
  \item{delmax}{Positive maximum value for Shrinkage delta-factors that are strictly less
    than 1 whenever m-Extent exceeds ZERO (default = 0.999999.)}
  \item{dMSE}{Maximum Likelihood Shrinkage Factors.}
  \item{k}{Unrestricted Ridge k-factor = non-negative multiplier.}
}
\value{An output list consisting of either two objects (k-factor & d-matrix) or one (m):
  \item{kinc or k}{New k-factor producing an m-Extent-of-shrinkage of mobj or muobj.}
  \item{d}{A pxp diagonal matrix of multiplicative shrinkage (delta) factors for the
     uncorrelated components of the regression coefficient vector.}
  \item{m}{New m-Extent for Shrinkage.}
}
\author{Bob Obenchain <wizbob@att.net>}
\keyword{internal}