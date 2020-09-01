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
\alias{print.unr.biv}
\alias{mstep}
\description{These are internal functions for printing and calculating quantities such
 as additive K-factors and multiplicative Delta-Factors that determine a specified
 m-Extent of Shrinkage. mstep() is called by qm.ridge(), true.risk(), true.simu(), aug.lars()
 and uc.lars().
}
\usage{
  mstep(mobj, kinc, p, qp, eqm1)
}
\arguments{
  \item{mobj}{Objective value for the M-extent-of-shrinkage along a path of shape qp.}
  \item{kinc}{Previous additive eigenvalue inflation constant yielding a smaller mobj.}
  \item{p}{Number of elements in the regression coefficient vector.}
  \item{qp}{Shrinkage path \bold{Q}-shape parameter.}
  \item{eqm1}{Vector of X'X eigenvalues raised to the power (qp-1).}
}
\value{An output list consisting of two objects (k-factor & d-matrix):
  \item{kinc}{New k-factor producing an m-Extent-of-shrinkage of mobj.}
  \item{d}{A pxp diagonal matrix of multiplicative shrinkage (delta) factors for the
   uncorrelated components of the regression coefficient vector.}
}
\author{Bob Obenchain <wizbob@att.net>}
\keyword{internal}