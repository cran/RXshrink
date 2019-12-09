\name{qm.ridge}
\alias{qm.ridge}
\title{Maximum Likelihood Shrinkage in Regression}
\description{
  Compute and display TRACEs for the \bold{Q}-shaped shrinkage path, including the \bold{M}-extent
  of shrinkage along that path, that are most likely under Normal-distribution theory to
  yield optimal reducions in MSE Risk.  
}
\usage{
  qm.ridge(form, data, rscale = 1, Q = "qmse", steps = 8, nq = 21,
              qmax = 5, qmin = -5, omdmin = 9.9e-13)
}
\arguments{
  \item{form}{A regression formula [y~x1+x2+...] suitable for use with lm().}
  \item{data}{Data frame containing observations on all variables in the formula.}
  \item{rscale}{One of three possible choices (0, 1 or 2) for rescaling of variables
   as they are being "centered" to remove non-essential ill-conditioning:  0 implies no
   rescaling; 1 implies divide each variable by its standard error; 2 implies rescale as
   in option 1 but re-express answers as in option 0.}
  \item{Q}{Shape parameter that controls the curvature of the shrinkage path through
   regression-coefficient likelihood space (default = "qmse" implies use the value found 
   most likely to be optimal.)  Use Q = 0 to specify Hoerl-Kennard "ordinary" ridge regression.}
  \item{steps}{Number of equally spaced values per unit change along the horizontal
   M-extent-of-shrinkage axis for estimates to be calculated and displayed in TRACES
   (default = 8.)}
  \item{nq}{Number of equally spaced values on the lattice of all possible values for
   shrinkage Q-shape between the "qmin" and "qmax" parameter settings (default = 21.)}
  \item{qmax}{Maximum allowed Q-shape (default = +5.)}
  \item{qmin}{Minimum allowed Q-shape (default = -5.)}
  \item{omdmin}{Strictly positive minimum allowed value for one-minus-delta (default = 9.9e-013.)}
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
\value{ An output list object of class qm.ridge:
  \item{form}{The regression formula specified as the first argument.}
  \item{data}{Name of the data.frame object specified as the second argument.}
  \item{p}{Number of regression predictor variables.}
  \item{n}{Number of complete observations after removal of all missing values.}
  \item{r2}{Numerical value of R-square goodness-of-fit statistic.}
  \item{s2}{Numerical value of the residual mean square estimate of error.}
  \item{prinstat}{Listing of principal statistics.}
  \item{mx}{Matrix containing mean values of X-predictors.}
  \item{crlqstat}{Listing of criteria for maximum likelihood selection of path Q-shape.}
  \item{qmse}{Numerical value of Q-shape most likely to be optimal.}
  \item{qp}{Numerical value of the Q-shape actually used for shrinkage.}
  \item{coef}{Matrix of shrinkage-ridge regression coefficient estimates.}
  \item{risk}{Matrix of MSE risk estimates for fitted coefficients.}
  \item{exev}{Matrix of excess MSE eigenvalues (ordinary least squares minus ridge.)}
  \item{infd}{Matrix of direction cosines for the estimated inferior direction, if any.}
  \item{spat}{Matrix of shrinkage pattern multiplicative delta factors.}
  \item{mlik}{Listing of criteria for maximum likelihood selection of M-extent-of-shrinkage.}
  \item{sext}{Listing of summary statistics for all M-extents-of-shrinkage.}
  \item{mClk}{Most Likely Extent of Shrinkage Observed: best multiple of (1/steps) <= p.} 
  \item{minC}{Minimum Observed Value of Normal-theory -2*log(Likelihood).}
}
\references{
  Burr TL, Fry HA. (2005)
  Biased Regression: The Case for Cautious Application. \emph{Technometrics}
  \bold{47}, 284-296.

  Goldstein M, Smith AFM. (1974)
  Ridge-type estimators for regression analysis. \emph{J. Roy. Stat. Soc. B}
  \bold{36}, 284-291. (The 2-parameter shrinkage family.)
  
  Obenchain RL. (1978) Good and Optimal Ridge Estimators. \emph{Annals of Statistics}
  \bold{6}, 1111-1121.

  Obenchain RL. (2005)
  \emph{Shrinkage Regression: ridge, BLUP, Bayes, spline and Stein}.  Electronic
  book-in-progress (200+ pages.) \url{http://localcontrolstatistics.org}

  Obenchain RL. (2018) \emph{RXshrink_in_R.PDF}  RXshrink package vignette-like file.
  \url{http://localcontrolstatistics.org}
}
\author{Bob Obenchain <wizbob@att.net>}
\seealso{\code{\link{unr.ridge}}, \code{\link{correct.signs}}, \code{\link{true.risk}} 
and \code{\link{true.simu}}.}
\examples{
  data(longley2)
  form <- GNP~GNP.deflator+Unemployed+Armed.Forces+Population+Year+Employed
  rxrobj <- qm.ridge(form, data=longley2)
  rxrobj
  plot(rxrobj)
  str(rxrobj)
}
\keyword{regression}