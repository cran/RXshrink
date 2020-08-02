\name{mofk}
\alias{mofk}
\title{m-Extents of Shrinkage used in unr.ridge() Calculations.}
\description{
  The mofk() function computes the numerical m-Extent of Shrinkage corresponding to a given
  value for "k", the scalar that multiplies observed ML estimates of "dMSE"-factors for
  each of the "p" uncorrelated components in GRR estimation using unr.ridge(). These
  calculations are quite simple; the function is many-to-one in the sense that all truly
  large values of k imply that m = 0.0, and k values that are either 0.0 or negative imply
  that m = p.
}
\usage{
  mofk(p, k, dMSE)
}
\arguments{
  \item{p}{The integer number of non-constant x-variables used in defining the linear
   model being fitted to ill-conditioned (intercorrelated, confounded) data. Note that
   p must also be rank of the given X-matrix.}
  \item{k}{The given k-factor that multiplies each given dMSE-factor.}
  \item{dMSE}{Maximum Likelihood estimates of Shrinkage Delta-Factors with minimum MSE risk.}
}
\value{The appropriate scalar value for m:
  \item{m}{The computed m-Extent of Shrinkage ...a scalar within [0, p].}
}
\references{
  Obenchain RL. (2020) \emph{Ridge TRACE Diagnostics.} \url{http://arxiv.org/abs/2005.14291}
  
  Obenchain RL. (2020) \emph{The Unrestricted Shrinkage Path: Technical Details.} "unrPath.pdf" 
  \url{http://localcontrolstatistics.org}
}
\author{Bob Obenchain <wizbob@att.net>}
\seealso{\code{\link{unr.ridge}} and \code{\link{kofm}}.}
\keyword{regression}
