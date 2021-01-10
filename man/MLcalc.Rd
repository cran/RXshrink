\name{MLcalc}
\alias{MLcalc}
\title{Calculate Unrestricted Maximum Likelihood (ML) point-estimates for a Linear Model that
  are either Unbiased (OLS) or nearly Optimally Biased under Normal-distribution theory.}
\description{
  Compute MSE risk-optimal point-estimates of Beta-Coefficients and their Relative MSE risks.
  Much of the code for this function is identical to that of unr.ridge(), which computes
  multiple points along the "Unrestricted" shrinkage-path. MLcalc() restricts attention to
  only 2 points: Unbiased OLS (BLUE) and the best (nearly optimally) Biased estimate with
  nearly minimum MSE.}
\usage{
  MLcalc(form, data, rscale = 1, delmax = 0.999999)}
\arguments{
  \item{form}{A regression formula [y~x1+x2+...] suitable for use with lm().}
  \item{data}{Data frame containing observations on all variables in the formula.}
  \item{rscale}{One of three possible choices (0, 1 or 2) for "rescaling" of variables
   (after being "centered") to remove all "non-essential" ill-conditioning:  0 implies no
   rescaling; 1 implies divide each variable by its standard error; 2 implies rescale as
   in option 1 but re-express answers as in option 0.}
  \item{delmax}{Maximum allowed value for Shrinkage delta-factors that is strictly less
   than 1. (default = 0.999999, which prints as 1 when rounded to fewer than 6 decimal
   places.)}
}
\details{Ill-conditioned and/or nearly multi-collinear regression models are unlikely to
  produce Ordinary Least Squares (OLS) regression coefficient estimates that are very
  close, numerically, to their unknown true values.  Specifically, OLS estimates can then
  tend to have "wrong" numerical signs and/or unreasonable relative magnitudes, while
  shrunken (generalized ridge) estimates chosen to maximize their likelihood of reducing
  Mean Squared Error (MSE) Risk (expected squared-error loss) can be more stable
  numerically.  On the other hand, because only OLS estimates are guaranteed to be minimax
  when risk is matrix valued (truly multivariate), no guarantee of an expected reduction in
  MSE Risk is necessarily associated with Generalized Ridge Regression shrinkage. 
}
\value{ An output list object of class MLcalc:
  \item{data}{Name of the data.frame object specified as the second argument.}
  \item{form}{The regression formula specified as the first argument.}
  \item{p}{Number of regression predictor variables.}
  \item{n}{Number of complete observations after removal of all missing values.}
  \item{r2}{Numerical value of R-square goodness-of-fit statistic.}
  \item{s2}{Numerical value of the residual mean square estimate of error.}
  \item{prinstat}{Listing of principal statistics.}
  \item{gmat}{Orthogonal Matrix of Direction Cosines for Principal Axes [1:p, 1:p].}
  \item{beta}{Numerical shrinkage-ridge regression coefficient estimates [1:2, 1:p].}
  \item{rmse}{Numerical MSE risk estimates for fitted coefficients [1:2, 1:p].}
  \item{dMSE}{Numerical delta-factors for shrinking OLS components [1:p].}
  \item{ys}{Numerical rescaling factor for y-outcome variable [1, 1].}
  \item{xs}{Numerical rescaling factors for given x-variables [1:p].}  
}
\references{

  Thompson JR. (1968) Some shrinkage techniques for estimating the mean. \emph{Journal
  of the American Statistical Association} \bold{63}, 113-122. (The ``cubic'' estimator.)
  
  Obenchain RL. (2020) Ridge TRACE Diagnostics. \url{https://arxiv.org/abs/2005.14291}

  Obenchain RL. (2021) \emph{RXshrink_in_R.PDF} \url{http://localcontrolstatistics.org}
  RXshrink package vignette-like file.
}
\author{Bob Obenchain <wizbob@att.net>}
\seealso{\code{\link{unr.ridge}}, \code{\link{MLboot}}}
\keyword{regression}
