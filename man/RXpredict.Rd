\name{RXpredict}
\alias{RXpredict}
\title{Predictions from Models Fit using RXshrink Generalized Ridge Estimation Methods.}
\description{
  RXpredict() makes in-sample predictions (i.e. computes "fitted.values") for all 5 forms of
  RXshrink estimation either at user specified m-Extents of shrinkage or at the normal-theory
  Maximum Lilelihood (ML) m-Extent.
}
\usage{
  RXpredict(x, data, m="ML", rscale=1)
}
\arguments{
  \item{x}{An object output by one of the 5 RXshrink estimation functions. Thus class(x) must
    be "qm.ridge", "unr.ridge", "aug.lars", "uc.lars" or "correct.signs".}
  \item{data}{Data frame containing observations on all variables used by the RXshrink function
    used in estimation.}
  \item{m}{Either [i] a "numeric" value >= 0 and <= rank(X) or [ii] the (default) character string
    "ML" to request use of the m-Extent of shrinkage most likely to be MSE optimal under Normal
    distribution-theory.}
  \item{rscale}{One of three possible choices (0, 1 or 2) for "rescaling" of variables
   (after being "centered") to remove all "non-essential" ill-conditioning:  0 implies no
   rescaling; 1 implies divide each variable by its standard error; 2 implies rescale as
   in option 1 but re-express answers as in option 0.}
}
\value{ An output list object of class RXpredict:
  \item{cryprd}{Predicted values for the "centered" and possibly "rescaled" outcome y-vector, cry.
    These values correspond, for example, to the default "predicted.values" from lm().}
  \item{cry}{This the "centered" and possibly "rescaled" outcome y-vector from the input data.frame.}
  \item{yvecprd}{Predicted values for the Y-outcome variable, yvec.}
  \item{yvec}{The Y-outcome vector from the input data.frame specified by the "data" argument.}
  \item{m}{"numeric" Value of m-Extent requested in the call to RXpredict(), possibly via a default
    call with m="ML". Restriction: 0 <= m <= rank(X).} 
  \item{mobs}{Observed m-Extent most close to the above m, yet is on the lattice of m-Extents stored
    within the x-object.}
}
\references{
  Obenchain RL. (1978) Good and Optimal Ridge Estimators. \emph{Annals of Statistics}
  \bold{6}, 1111-1121. <doi:10.1214/aos/1176344314>

  Obenchain RL. (2005)
  \emph{Shrinkage Regression: ridge, BLUP, Bayes, spline and Stein}.  Electronic
  book-in-progress (185+ pages.) \url{http://localcontrolstatistics.org}

  Obenchain RL. (2020) \emph{RXshrink_in_R.PDF}  RXshrink package vignette-like file.
  \url{http://localcontrolstatistics.org}
}
\author{Bob Obenchain <wizbob@att.net>}
\seealso{\code{\link{qm.ridge}}, \code{\link{unr.ridge}} and \code{\link{MLtrue}}.}
\examples{
  data(tycobb)
  form <- batavg~atbats+seasons+CMspl
  rxuobj <- unr.ridge(form, data=tycobb)
  tycfit <- RXpredict(rxuobj, tycobb, m="ML")
  plot(tycfit)
  abline(h = 0.400, col = "red", lty = 2, lwd = 2)
}
\keyword{regression}
