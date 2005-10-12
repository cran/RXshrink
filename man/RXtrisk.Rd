\name{RXtrisk}
\alias{RXtrisk}
\title{True MSE Risk of Shrinkage Resulting from Known Regression Parameters}
\description{
  By specifying numerical values for regression parameters (uncorrelated components and error
  sigma) that usually are unknown, these functions can calculate and display True MSE Risk
  statistics associated with shrinkage along a given \bold{Q}-shaped path.  
}
\usage{
RXtrisk(form, data, trug, trus, Q = 0, rscale = 1,
          steps = 8, qmax = 5, qmin = -5)
}
\arguments{
  \item{form}{A regression formula [y~x1+x2+...] suitable for use with lm().}
  \item{data}{Data frame containing observations on all variables in the formula.}
  \item{trug}{Column vector of numerical values for the true uncorrelated components
   of the regression coefficient vector.}
  \item{trus}{Numerical value for the true error standard deviation, Sigma.}
  \item{Q}{Numerical value for the shape parameter controlling shrinkage path curvature.
   Default shape is Q = 0 for Hoerl-Kennard "ordinary" ridge regression.}
  \item{rscale}{One of three possible choices (0, 1 or 2) for rescaling of variables
   as they are being "centered" to remove non-essential ill-conditioning:  0 implies no
   rescaling; 1 implies divide each variable by its standard error; 2 implies rescale as
   in option 1 but re-express answers as in option 0.}
  \item{steps}{Number of equally spaced values per unit change along the horizontal
   M-extent-of-shrinkage axis where estimates are calculated and displayed in TRACES
   (default = 8.)}
  \item{qmax}{Maximum allowed Q-shape (default = +5.)}
  \item{qmin}{Minimum allowed Q-shape (default = -5.)}
}
\details{The RXridge() functions calculate maximum likelihood estimates (corrected, if
  necessary, so as to have correct range) for typical statistical inference situations where
  regression parameters are unknowns.  In sharp contrast with this usual situation, the
  RXtrisk() functions show exactly how expected regression coefficients and true Mean
  Squared Error Risk actually do change with shrinkage when regression parameters take on
  specified, KNOWN numerical values.
}
\value{ An output list object of class RXtrisk:
  \item{form}{The regression formula specified as the first argument.}
  \item{data}{Name of the data.frame object specified as the second argument.}
  \item{trug}{Vector of numerical values for the true uncorrelated gamma components.}
  \item{trus}{Numerical value for the true error standard deviation, Sigma.}
  \item{qp}{Numerical value of the Q-shape actually used for shrinkage.}
  \item{p}{Number of regression predictor variables.}
  \item{n}{Number of complete observations after removal of all missing values.}
  \item{prinstat}{Listing of principal statistics.}
  \item{coef}{Matrix of expected shrinkage-ridge regression coefficients.}
  \item{rmse}{Matrix of true MSE risk values for shrunken coefficients.}
  \item{exev}{Matrix of true excess eigenvalues (ordinary least squares minus ridge.)}
  \item{infd}{Matrix of direction cosines for the true inferior direction, if any.}
  \item{spat}{Matrix of shrinkage pattern delta factors.}
  \item{sext}{Listing of summary statistics for all M-extents-of-shrinkage.}
}
\author{Bob Obenchain <softrx@iquest.net>}
\seealso{\code{\link{RXridge}} and \code{\link{RXtsimu}}.}
\examples{
    data(mpg)  
    form <- mpg~cylnds+cubins+hpower+weight
    rxrobj <- RXridge(form, data=mpg)
    # define true parameter values.
    trugam <- matrix(c(-.5,-.1,.1,-.6),4,1)
    trusig <- 0.4
    # create true shrinkage MSE risk scenario.
    trumse <- RXtrisk(form, data=mpg, trugam, trusig, Q=-1, steps=4)
    plot(trumse)
}
\keyword{regression}
\keyword{hplot}