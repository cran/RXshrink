\name{plot.eff.ridge}
\alias{plot.eff.ridge}
\title{Plot method for eff.ridge objects}
\description{Plot TRACE displays for the "efficient" Path that is [1]
 as short as possible and [2] consists of Piecewise Linear-Splines with a
 emph{Single Interior Knot} at the MSE Optimal m-Extent of Shrinkage. The
 default is to display all five Traces on one page with no legends.
}
\usage{
  \method{plot}{eff.ridge}(x, trace = "all", trkey = FALSE, \dots )
}
\arguments{
  \item{x}{Output list object of class eff.ridge.}
  \item{trace}{One of seven possible options: "all" to display 5 traces in one graph,
   "seq" to display 5 full-sized traces in sequence in response to user prompts,
   "coef" to display only the estimated beta coefficient trace,
   "rmse" to display only the estimated relative MSE risk trace,
   "exev" to display only the estimated excess MSE eigenvalue (OLS minus larlso),
   "infd" to display only the estimated inferior direction cosine trace, or
   "spat" to display only the delta-factor pattern trace.}
  \item{trkey}{If TRUE, display a crude legend at the bottom of each trace plot.}
  \item{\dots}{Optional argument(s) passed on to plot().}
}
\value{NULL}
\author{Bob Obenchain <wizbob@att.net>}
\examples{
  data(longley2)
  form <- GNP~GNP.deflator+Unemployed+Armed.Forces+Population+Year+Employed
  rxefobj <- eff.ridge(form, data=longley2)
  plot(rxefobj)
}
\keyword{trace}
