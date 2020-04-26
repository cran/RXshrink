\name{unr.ridge}
\alias{unr.ridge}
\title{Unrestricted Maximum Likelihood Shrinkage using Linear-Spline Patterns}
\description{
  Compute and display TRACEs for the p-paramater Shrinkage PATH passing through the
  (classical) Normal-theory point-estimate of Regression Coefficients in p-Space that
  is most likely to yield optimal reducions in MSE Risk.  
}
\usage{
  unr.ridge(form, data, rscale = 1, steps = 8, delmax = 0.999999)
}
\arguments{
  \item{form}{A regression formula [y~x1+x2+...] suitable for use with lm().}
  \item{data}{Data frame containing observations on all variables in the formula.}
  \item{rscale}{One of three possible choices (0, 1 or 2) for rescaling of variables
   as they are being "centered" to remove non-essential ill-conditioning:  0 implies no
   rescaling; 1 implies divide each variable by its standard error; 2 implies rescale as
   in option 1 but re-express answers as in option 0.}
  \item{steps}{Number of equally spaced values per unit change along the horizontal
   M-extent-of-shrinkage axis for estimates to be calculated and displayed in TRACES
   (default = 8.)}
  \item{delmax}{Maximum allowed value for Shrinkage delta-factors that is strictly less
   than 1. (default = 0.999999, which prints as 1 when rounded to fewer than 6 decimal
   places.)}
}
\details{Ill-conditioned and/or nearly multi-collinear regression models are unlikely to
  produce Ordinary Least Squares (OLS) regression coefficient estimates that are very
  close, numerically, to their unknown true values.  Specifically, OLS estimates can then
  tend to have "wrong" numerical signs and/or unreasonable relative magnitudes, while
  shrunken (generalized ridge) estimates chosen to maximize their likelihood of reducing
  Mean Squared Error (MSE) Risk (expected loss) can be much more stable and reasonable,
  numerically.  On the other hand, because only OLS estimates are guaranteed to be minimax
  when risk is matrix valued (truly multivariate), no guarantee of an actual reduction in
  MSE Risk is necessarily associated with shrinkage. 
}
\value{ An output list object of class unr.ridge:
  \item{data}{Name of the data.frame object specified as the second argument.}
  \item{form}{The regression formula specified as the first argument.}
  \item{p}{Number of regression predictor variables.}
  \item{n}{Number of complete observations after removal of all missing values.}
  \item{r2}{Numerical value of R-square goodness-of-fit statistic.}
  \item{s2}{Numerical value of the residual mean square estimate of error.}
  \item{prinstat}{Listing of principal statistics.}
  \item{gmat}{Orthogonal Matrix of Direction Cosines for Principal Axes.}
  \item{coef}{Matrix of shrinkage-ridge regression coefficient estimates.}
  \item{rmse}{Matrix of MSE risk estimates for fitted coefficients.}
  \item{exev}{Matrix of excess MSE eigenvalues (ordinary least squares minus ridge.)}
  \item{infd}{Matrix of direction cosines for the estimated inferior direction, if any.}
  \item{spat}{Matrix of shrinkage pattern multiplicative delta factors.}
  \item{mlik}{Listing of criteria for maximum likelihood selection of M-extent-of-shrinkage.}
  \item{sext}{Listing of summary statistics for all M-extents-of-shrinkage.}
  \item{mUnr}{Unrestricted m-Extent of Shrinkage corresponding to k == 1 on TRACE displays.}  
  \item{mClk}{Most Likely Extent of Shrinkage Observed: best multiple of (1/steps) <= p.} 
  \item{minC}{Minimum Observed Value of Normal-theory -2*log(Likelihood).}    
}
\references{

  Thompson JR. (1968) Some shrinkage techniques for estimating the mean. \emph{Journal
  of the American Statistical Association} \bold{63}, 113-122. (The ``cubic'' estimator.)
  
  Obenchain RL. (1978) Good and Optimal Ridge Estimators. \emph{Annals of Statistics}
  \bold{6}, 1111-1121.

  Obenchain RL. (2005) \emph{Shrinkage Regression: ridge, BLUP, Bayes, spline and Stein}.
  Electronic book-in-progress; 185+ pages. \url{http://localcontrolstatistics.org}

  Obenchain RL. (2018) \emph{RXshrink_in_R.PDF}  RXshrink package vignette-like file.
  \url{http://localcontrolstatistics.org}
}
\author{Bob Obenchain <wizbob@att.net>}
\seealso{\code{\link{qm.ridge}}, \code{\link{correct.signs}}, \code{\link{true.risk}} 
and \code{\link{true.simu}}.}
\examples{
  data(longley2)
  form <- GNP~GNP.deflator+Unemployed+Armed.Forces+Population+Year+Employed
  rxuobj <- unr.ridge(form, data=longley2)
  rxuobj
  plot(rxuobj)
  str(rxuobj)
}
\keyword{regression}
