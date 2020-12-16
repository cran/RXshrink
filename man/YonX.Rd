\name{YonX}
\alias{YonX}
\title{Maximum Likelihood (ML) Shrinkage in Simple Linear Regression}
\description{
  Compute and display Normal-theory ML Shrinkage statistics when a Y-variagle is regressed
  upon a Single X-variable. Since the OLS (BLUE) estimate is a scalar in this case, the
  MSE optimal Shrinkage factor, dMSE, is also a scalar (less than +1 and greater than 0).  
  The m-Extent of overall Optimal Shrinkage corresponding to this solution is marked by
  the "purple" vertical dashed-line on all TRACE displays for YonX() Shrinkage. Note that
  YonX() has no "steps" argument; YonX() calculations always use 101 equally spaced m-Values
  between 0.00 and 1.00 (in steps of 0.01).
}
\usage{
  YonX(form, data, delmax = 0.999999)
}
\arguments{
  \item{form}{A regression formula [y ~ x] suitable for use with lm().}
  \item{data}{Data frame containing observations on both variables in the formula.}
  \item{delmax}{Maximum allowed value for Shrinkage delta-factor that is strictly less
   than 1. (default = 0.999999, which prints as 1 when rounded to fewer than 6 decimal
   places.)}
}
\details{Since only a single X-variable is being used, these "simple" models are
  (technically) never "Ill-conditioned". Of course, the Y-variable may be nearly
  multi-collinear with the given X-variable, but this simply means that the model
  then has low "lack-of-fit". In fact, the OLS estimate can never have the "wrong"
  numerical sign in these simple p = 1 models!  Furthermore, since "risk" estimates
  are scalar-valued, no "emse" TRACE is displayed because its content would duplicate
  information in the "rmse" TRACE. Similarly, no "infd" TRACE is displayed because
  any "inferior direction" COSINE would be either +1 ("upwards") or -1 ("downwards").
  Finally, see the documentation of plot.YonX() for information on some curious
  differences between estimates "Most Likely" to minimize risk under Normal-theory
  and the amount of shrinkage that numerically minimizes traditional (nonlinear)
  estimates of MSE risk.
}
\value{ An output list object of class YonX:
  \item{data}{Name of the data.frame object specified as the second argument.}
  \item{form}{The regression formula specified as the first argument to YonX() must
    have only one right-hand-side X-variable.}
  \item{p}{Number of X-variables must be p = 1.}
  \item{n}{Number of complete observations after removal of all missing values.}
  \item{r2}{Numerical value of R-square goodness-of-fit statistic.}
  \item{s2}{Numerical value of the residual mean square estimate of error.}
  \item{prinstat}{Vector of five Principal Statistics: eigval, sv, b0, rho & tstat.} 
  \item{yxnam}{Character Names of "Y" and "X" data vectors.}
  \item{yvec}{"Y" vector of data values.}
  \item{xvec}{"X" vector of data values.}
  \item{coef}{Vector of Shrinkage regression beta-coefficient estimates: delta * B0.}
  \item{rmse}{Vector of MSE Risk unbiased estimates for shrunken coefficients.}
  \item{spat}{Vector of Shrinkage (multiplicative) delta-factors: 1.00 to 0.00 by -0.01.}
  \item{mlik}{Normal-theory Likelihood ...for Maximum Likelihood estimation of Shrinkage m-Extent.}
  \item{sext}{Listing of summary statistics for all M-extents-of-shrinkage.}
  \item{mUnr}{Unrestricted optimal m-Extent of Shrinkage from the dMSE estimate; mUnr = 1 - dMSE.}
  \item{Good}{"GOOD" m-Extent of Shrinkage corresponding to max(0, 2*dMSE-1).}
  \item{mMSE}{Minimum MSE risk estimate.}   
  \item{mClk}{Most Likely Observed m-Extent of Shrinkage: best multiple of (1/steps) <= 1.} 
  \item{minC}{Minimum Observed Value of CLIK Normal-theory -2*log(Likelihood-Ratio).}
  \item{minE}{Minimum Observed Value of EBAY (Empirical Bayes) criterion.}
  \item{minR}{Minimum Observed Value of RCOF (Random Coefficients) criterion.}  
  \item{mRmin}{m-Extent where the numerical "rmse" estimate is minimized.}
  \item{mReql}{m-Extent where the numerical "rmse" estimate is first >= its observed OLS value
    in rmse[1] (i.e. the rescaled variance)at m = 0.}
  \item{Phi2}{Phi-Squared noncentrality parameter of the F-ratio for testing H: true
    beta-coefficient = zero.}  
  \item{dMSE}{Most Likely to be Optimal-value for the Shrinkage Delta-factor.}  
}
\references{

  Obenchain RL. (1978) Good and Optimal Ridge Estimators. \emph{Annals of Statistics}
  \bold{6}, 1111-1121. <doi:10.1214/aos/1176344314>
  
  Obenchain RL. (2020) \emph{RXshrink_in_R.PDF}  RXshrink package vignette-like file.
  \url{http://localcontrolstatistics.org}
}
\author{Bob Obenchain <wizbob@att.net>}
\seealso{\code{\link{correct.signs}} and \code{\link{MLtrue}}}
\examples{
  data(haldport)
  form <- heat ~ p4caf
  YXobj <- YonX(form, data=haldport)
  YXobj
  plot(YXobj)
}
\keyword{regression}
