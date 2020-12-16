\name{plot.YonX}
\alias{plot.YonX}
\title{Plot method for YonX objects}
\description{Graphics for Shrinkage in "Simple" Linear Regression: Models with only p=1
 X-variable. The default is to first [1] display three TRACES in one plot and then to [2]
 display the "Y on X" scatter-plot plus Three Fitted Lines: the OLS fit is BLUE, the
 optimally Shrunken fit is "purple", and the most shrinkage with estimated relative risk
 not exceeding OLS is "red". Note: these 3 fitted lines can be quite close to each other
 whenever model lack-of-fit is small.
}
\usage{
  \method{plot}{YonX}(x, trace = "all", \dots )
}
\arguments{
  \item{x}{Output list object of class "YonX".}
  \item{trace}{One of six possible options: "all" to display 3 traces on the first plot,
    then the "YonX" scatter plot; "seq" to display 4 full-sized plots in sequence (in
    response to user prompts); "coef" to display only the estimated beta coefficient trace
    (a straight line); "rmse" to display only the (quadratic) estimated relative MSE risk
    trace; "spat" to display only the delta-factor trace (a straight line); or "YonX" to
    display only the Y-vs-X scatter plot with 3 fitted-lines: the OLS fitted line (BLUE),
    the "purple" Maximum-Likelihood (optimally biased) line, and the "red" line marking
    the most shrinkage with estimated relative MSE risk less than or equal to that of OLS.}
  \item{\dots}{Optional argument(s) passed on to plot().}
}
\details{The effects of Shrinkage on Simple Linear Regression models (p = 1) are, in reality,
 NOT easier to illustrate than the corresponding effects on Multiple Linear Regressiom models
 (p > 1). In both situations, alternative estimates of effects and risks abound. For example,
 the estimate chosen can be [1] Maximum Likelihood under Normal-theory, [2] Unbiased under
 Normal-theory or [3] have "Correct Range". [See Obenchain (1978), equations (3.3) to (3.5),
 and corresponding text.] When a graphic contains only a single curve, a "reality" is that
 the general "shape" of the curve (plus any highlighted "points" on that curve) should "look
 right" or, at the very least, "reasonable".
 
 My choices among alternative estimates of (nonlinear) MSE risk were initially made roughly
 30 years ago ...and have remained unchanged (primarily) for consistency with earlier versions
 of RXshrink.

 The NEW "rmse" TRACE for class "YonX" shows the ML "purple" point and dotted-line, but attemps 
 to focus user attention upon the "green" point (and dotted-line) that numerically minimizes the
 traditional quadratic MSE risk estimate. Similarly, the "red" point (and dotted-line) marks the
 shrunken estimate with MSE relative risk that equals the relative Variance of (unbiased) OLS.
 Note that this "red" upper limit suggests a lesser m-Extent of Shrinkage than the ML "Good"
 upper limit: max(0, 2*dMSE-1).
}
\value{NULL}
\references{

  Obenchain RL. (1978) Good and Optimal Ridge Estimators. \emph{Annals of Statistics}
  \bold{6}, 1111-1121. <doi:10.1214/aos/1176344314>
  
  Obenchain RL. (2020) \emph{RXshrink_in_R.PDF}  RXshrink package vignette-like file.
  \url{http://localcontrolstatistics.org}
}
\author{Bob Obenchain <wizbob@att.net>}
\examples{
  data(haldport)
  form <- heat ~ p4caf
  YXobj <- YonX(form, data=haldport)
  plot(YXobj)
}
\keyword{trace}
