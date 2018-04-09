\name{RXlarlso}
\alias{RXlarlso}
\title{Maximum Likelihood Estimation of Effects in Least Angle Regression}
\description{
  Identify whether least angle regression estimates are generalized
  ridge shrinkage estimates and generate TRACE displays for estimates
  that do correspond to ridge shrinkage factors between 0.00 and 0.99.
}
\usage{
RXlarlso(form, data, rscale = 1, type = "lar", trace = FALSE, 
    eps = .Machine$double.eps, omdmin = 9.9e-13, \dots) 
}
\arguments{
  \item{form}{A regression formula [y~x1+x2+...] suitable for use with lm().}
  \item{data}{Data frame containing observations on all variables in the formula.}
  \item{rscale}{One of three possible choices (0, 1 or 2) for rescaling of variables
   as they are being "centered" to remove non-essential ill-conditioning:  0 implies no
   rescaling; 1 implies divide each variable by its standard error; 2 implies rescale as
   in option 1 but re-express answers as in option 0.}
  \item{type}{One of "lasso", "lar" or "forward.stagewise" for function lars(). Names can be
   abbreviated to any unique substring. Default in RXlarlso() is "lar".}
  \item{trace}{If TRUE, lars() function prints out its progress.}
  \item{eps}{The effective zero for lars().}
  \item{omdmin}{Strictly positive minimum allowed value for one-minus-delta (default = 9.9e-013.)}
  \item{\dots}{Optional argument(s) passed on to the lars() function from the lars R-package.}
}
\details{RXlarlso() calls the Efron/Hastie lars() function to perform Least Angle Regression on
 X-variables that have been centered and possibly rescaled but which may be (highly) correlated.
 Maximum likelihood TRACE displays paralleling those of RXridge are also computed and (optionally)
 plotted. 
}
\value{ An output list object of class RXlarlso:
  \item{form}{The regression formula specified as the first argument.}
  \item{data}{Name of the data.frame object specified as the second argument.}
  \item{p}{Number of regression predictor variables.}
  \item{n}{Number of complete observations after removal of all missing values.}
  \item{r2}{Numerical value of R-square goodness-of-fit statistic.}
  \item{s2}{Numerical value of the residual mean square estimate of error.}
  \item{prinstat}{Listing of principal statistics.}
  \item{gmat}{Orthogonal matrix of direction cosines for regressor principal axes.}
  \item{lars}{An object of class lars.}
  \item{coef}{Matrix of shrinkage-ridge regression coefficient estimates.}
  \item{risk}{Matrix of MSE risk estimates for fitted coefficients.}
  \item{exev}{Matrix of excess MSE eigenvalues (ordinary least squares minus ridge.)}
  \item{infd}{Matrix of direction cosines for the estimated inferior direction, if any.}
  \item{spat}{Matrix of shrinkage pattern multiplicative delta factors.}
  \item{mlik}{Listing of criteria for maximum likelihood selection of M-extent-of-shrinkage.}
  \item{sext}{Listing of summary statistics for all M-extents-of-shrinkage.}
}
\references{
  Breiman L. (1995) Better subset regression using the non-negative garrote.
  \emph{Technometrics} \bold{37}, 373-384.

  Efron B, Hastie T, Johnstone I, Tibshirani R. (2004)
  Least angle regression. \emph{Annals of Statistics} \bold{32}, 407-499.

  Obenchain RL. (2005)
  \emph{Shrinkage Regression: ridge, BLUP, Bayes, spline and Stein}.  Electronic
  book-in-progress (200+ pages.) \url{http://localcontrolstatistics.org}

  Obenchain RL. (2018) \emph{RXshrink_in_R.PDF}  RXshrink package vignette-like file.
  \url{http://localcontrolstatistics.org}

  Tibshirani R. (1996) Regression shrinkage and selection via the lasso.
  \emph{J. Roy. Stat. Soc. B}  \bold{58}, 267-288.
}
\author{Bob Obenchain <wizbob@att.net>}
\seealso{\code{\link{RXuclars}}.}
\examples{
  data(longley2)
  form <- GNP~GNP.deflator+Unemployed+Armed.Forces+Population+Year+Employed
  rxlobj <- RXlarlso(form, data=longley2)
  rxlobj
  names(rxlobj)
  plot(rxlobj)
}
\keyword{regression}
\keyword{hplot}