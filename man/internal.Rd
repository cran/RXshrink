\name{RXshrink-internal}
\title{Internal RXshrink functions}
\alias{RXshrink-internal}
\alias{print.RXlarlso}
\alias{print.RXridge}
\alias{print.RXtrisk}
\alias{print.RXtsimu}
\alias{print.RXuclars}
\alias{RXmstep}
\description{These are internal functions for printing and calculating quantities such
 as additive K-factors and multiplicative Delta-Factors that determine a specified M
 extent-of-shrinkage.  RXmstep() is called by RXridge(), RXtrisk() and RXtsimu().
}
\usage{
  RXmstep(mobj, kinc, p, qp, eqm1)
}
\arguments{
  \item{mobj}{Newest objective value for the M-extent-of-shrinkage along a path of shape qp.}
  \item{kinc}{Previous additive eigenvalue inflation constant yielding a smaller mobj.}
  \item{p}{Number of elements in the regression coefficient vector.}
  \item{qp}{Shrinkage path \bold{Q}-shape parameter.}
  \item{eqm1}{Vector of X'X eigenvalues raised to the power (qp-1).}
}
\value{An output list consisting of two objects:
  \item{kinc}{New, increased value of kinc that produces an M-extent-of-shrinkage = mobj.}
  \item{d}{A pxp diagonal matrix of multiplicative shrinkage (delta) factors for the
     uncorrelated components of the regression coefficient vector.}
}
\author{Bob Obenchain <wizbob@att.net>}
\keyword{internal}