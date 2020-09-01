\name{unr.aug}
\alias{unr.aug}
\title{Augment calculations performed by unr.ridge() to prepare for display of eliptical confidence
 regions for pairs of biased coefficient estimates using plot.unr.biv()}
\description{
  This function makes classical (rather than Bayesian) Normal distribution-theory calculations of the
  form proposed in Obenchain(1977). Instead of providing "new" confidence regions for estimable linear 
  functions, Generalized Ridge Regression (GRR) can focus interest on estimates that are within traditional
  confidence intervals and regions but which deviate reasonably from the interval or region centroid.
}
\usage{
  unr.aug(urobj)
}
\arguments{
  \item{urobj}{An output object of class "unr.ridge".}
}
\value{ An output list object of class "unr.aug"...
  \item{p}{Number of regression predictor variables.}
  \item{LMobj}{The lm() output object for the model fitted using unr.ridge().}
  \item{bstar}{The p by p+1 matrix of shrunken GRR coefficients. The first p correspond to "knots"
   on piecewise linear splines, and the last column contains minimum MSE risk estimates.}
  \item{mcal}{p+1 increasing measures of shrinkage "Extent". The first is 0 for the OLS (BLUE)
   estimate, while the last denotes the "Extent" most likely to yield minimum MSE risk.}
  \item{vnams}{All p+1 names for variables used in the GRR model.} 
}
\references{
  
  Obenchain RL. (1977) Classical F-Tests and Confidence Regions for Ridge Regression.
  \emph{Technometrics} 19, 429-439. http://doi.org/10.1080/00401706.1977.10489582
  
  Obenchain RL. (2020) \emph{Ridge TRACE Diagnostics.} \url{https://arxiv.org/abs/2005.14291}
  
  Obenchain RL. (2020) \emph{RXshrink_in_R.PDF}  RXshrink package vignette-like file.
  \url{http://localcontrolstatistics.org}
}
\author{Bob Obenchain <wizbob@att.net>}
\seealso{\code{\link{unr.ridge}}, \code{\link{mofk}} and \code{\link{kofm}}}
\keyword{regression}
