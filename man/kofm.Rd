\name{kofm}
\alias{kofm}
\title{k-Multipliers and delta-Factors for unr.ridge() Shrinkage.}
\description{
  The kofm() function computes the "kStar" multipliers and the diagonal "d" Matrices
  used by the unr.ridge() function. While these calculations are extensive, they are
  greatly simplified by using interpolation between "knots" on "piecewise
  linear-spline" functions.
}
\usage{
kofm(muobj, p, dMSE, delmax = 0.999999) 
}
\arguments{
  \item{muobj}{A "numeric" value that is >= 0 and <= p.}
  \item{p}{The integer number of non-constant x-variables used in defining the linear
   model being fitted to ill-conditioned (intercorrelated, confounded) data. Note that
   p must also be rank of the given X-matrix.}
  \item{dMSE}{Maximum Likelihood estimates of minimum MSE risk Shrinkage Delta-Factors.}
  \item{delmax}{Maximum allowed value for Shrinkage delta-factors that is strictly less
   than 1. (default = 0.999999, which prints as 1 when rounded to fewer than 6 decimal
   places.)}   
}
\value{An output list of two objects:
  \item{kStar}{The k-factor needed to multiply each given dMSE-factor.}
  \item{d}{The (p x p) diagonal matrix of requested shrinkage-factors.}
}
\references{
  Obenchain RL. (2020) \emph{Ridge TRACE Diagnostics.} \url{https://arxiv.org/abs/2005.14291}
  
  Obenchain RL. (2020) \emph{The Unrestricted Shrinkage Path: Technical Details.} "unrPath.pdf" 
  \url{http://localcontrolstatistics.org}
}
\author{Bob Obenchain <wizbob@att.net>}
\seealso{\code{\link{unr.ridge}} and \code{\link{mofk}}.}
\keyword{regression}
