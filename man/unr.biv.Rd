\name{unr.biv}
\alias{unr.biv}
\title{Specify pairs of GRR Coefficient Estimates for display in Bivariate Confidence Regions}
\description{
  This function specifies which Pair of GRR estimates to display and the single (or dual) Confidence
  Level(s) of the Ellipse(s) displayed. Reguested confidence levels must both be equal to or greater
  than 0.05 and less than or equal to 0.95.
}
\usage{
  unr.biv(uraug, x1 = 1, x2 = 2, conf1 = 0.95, conf2 = 0.50) 
}
\arguments{
  \item{uraug}{An output list object of class "unr.aug" for a GRR model with p >= 2 x-Variables.}
  \item{x1}{Integer index value >= 1 and <= p for the x-Coefficient displayed on the horizontal axis.}
  \item{x2}{Integer index value >= 1 and <= p for the x-Coefficient displayed on the vertical axis.
   Index x2 must differ from x1 to display a plot.}
  \item{conf1}{This first Confidence level must be >= 0.05 and <= 0.95 to display an Ellipse.}
  \item{conf2}{When the second Confidence level is >= 0.05 and <= 0.95, its Ellipse is displayed. No
   plot is displayed when both conf1 and conf2 are outside of the [0.05, 0.95] range.}   
}
\value{ An output list object of class "unr.biv"...
  \item{p}{Number of regression predictor variables.}
  \item{LMobj}{The lm() output object for the model fitted using unr.ridge().}
  \item{bstar}{The p by p+2 matrix of shrunken GRR coefficients. The first p correspond to "knots"
   on piecewise linear splines, the next to last contains minimum MSE risk estimates, and the last
   represents the shrinkge terminus.}
  \item{mcal}{p+2 increasing measures of shrinkage "Extent". The first value is 0 for the OLS = BLU
   Estimate, plotted as a Point of color "blue". The next to last corresponds to the "Extent" most
   likely to yield minimum MSE risk is colored "purple", and the last may be outside the plot frame.}
  \item{ellip1}{matrix[100, 2] of points on confidence ellipse 1.}
  \item{conf1}{confidence level of ellipse 1 within [0.05, 0.95].}
  \item{ecor1}{Pearson correlation between x1 and x2 in ellipse 1.}
  \item{ellip2}{matrix[100, 2] of points on confidence ellipse 2.}
  \item{conf2}{confidence level of ellipse 2 within [0.05, 0.95].}
  \item{ecor2}{Pearson correlation between x1 and x2 in ellipse 2.}
}
\references{
  
  Obenchain RL. (1977) Classical F-Tests and Confidence Regions for Ridge Regression.
  \emph{Technometrics} 19, 429-439. http://doi.org/10.1080/00401706.1977.10489582
  
  Murdoch DJ. and Chow ED. (1996). A graphical display of large correlation matrices.
  \emph{The American Statistician} 50, 178-180.
  
  Murdoch DJ. ellipse: Functions for Drawing Ellipses and Ellipse-Like Confidence Regions.
  \url{https://CRAN.R-project.org/package=ellipse}
}
\author{Bob Obenchain <wizbob@att.net>}
\seealso{\code{\link{unr.aug}}, \code{\link{ellipse.lm}}, \code{\link{unr.ridge}}}
\keyword{regression}