\name{RXsigns}
\alias{RXsigns}
\title{Normal-theory Maximum Likelihood Estimation of Beta Coefficients with "Correct" Signs}
\description{
  RXsigns displays the Beta Coefficient estimate, denoted by \bold{B(=)}, that is most likely to have
  minimum MSE risk in the one, unknown direction PARALLEL to the true Beta in p-dimensional
  likelihood space. Shrinkage to ZERO of any components ORTHOGONAL to the true Beta is MSE optimal.
  Obenchain(1978) shows that \bold{B(=)} is of the form k * X'y where the scalar k is given by
  equation (4.2) on page 1118, the optimal shrinkage factors are proportional to known eigenvalues,
  and the formula for the maximum likelihood estimate of k given on page 1119 is corrected.  
}
\usage{
  RXsigns(form, data)
}
\arguments{
  \item{form}{A regression formula [y~x1+x2+...] suitable for use with lm().}
  \item{data}{Data frame containing observations on all variables in the formula.}
}
\details{Ill-conditioned (nearly multi-collinear) regression models can produce Ordinary Least
  Squares (OLS) regression coefficient estimates with numerical signs that differ from those
  of the X'y vector. This is disturbing because X'y contains the sample \bold{correlations} between
  the X-predictor variables and y-response variable if these variables have first been "centered"
  by subtracting off their mean values and then rescaled to vectors of length one.  Besides
  displaying the B(=) estimate, the RXsigns() function also displays the OLS vector, the
  "correlation form" of X'y, the estimated Delta shrinkage-factors and the rescaled coefficients,
  \bold{Bfit = f * B(=)}, where f is the positive scalar that minimizes the Residual Sum-of-Squares;
  RSS(Bfit) >= RSS(OLS).
}
\value{ An output list object of class RXsigns:
  \item{data}{Name of the data.frame object specified as the second argument.}
  \item{form}{The regression formula specified as the first argument.}
  \item{p}{Number of regression predictor variables.}
  \item{n}{Number of complete observations after removal of all missing values.}
  \item{r2}{Numerical value of R-square goodness-of-fit statistic.}
  \item{s2}{Numerical value of the residual mean square estimate of error.}
  \item{prinstat}{Listing of principal statistics (p by 5) from RXridge().}
  \item{kpb}{Maximum likelihood estimate of k-factor in equation (4.2) of Obenchain(1978).}
  \item{bmf}{Rescaling factor for B(=) to minimize the Residual Sum-of-Squares.}
  \item{signs}{Listing of five Beta coefficient statistics (p by 5): OLS, X'y, Delta, B(=) and Bfit.}
  \item{loff}{Lack-of-Fit statistics: Residual Sum-of-Squares for OLS, X'y, B(=) and Bfit.}
  \item{mcor}{Squared Correlation between the y-vector and its predicted values. The two values
    displayed are for OLS predictions or for predictions using Bfit, X'y or B(=). These two values
    are the familiar R^2 goodness-of-fit statistics for OLS and Bfit.}
}
\references{
  Obenchain RL. (1978) Good and Optimal Ridge Estimators. \emph{Annals of Statistics}
  \bold{6}, 1111-1121.

  Obenchain RL. (2005)
  \emph{Shrinkage Regression: ridge, BLUP, Bayes, spline and Stein}.  Electronic
  book-in-progress (185+ pages.) \url{http://localcontrolstatistics.org}

  Obenchain RL. (2018) \emph{RXshrink_in_R.PDF}  RXshrink package vignette-like file.
  \url{http://localcontrolstatistics.org}
}
\author{Bob Obenchain <wizbob@att.net>}
\seealso{\code{\link{RXridge}}, \code{\link{RXtrisk}} and \code{\link{RXtsimu}}.}
\examples{
  data(longley2)
  form <- GNP~GNP.deflator+Unemployed+Armed.Forces+Population+Year+Employed
  rxsobj <- RXsigns(form, data=longley2)
  rxsobj
  str(rxsobj)
}
\keyword{regression}
