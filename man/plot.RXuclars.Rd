\name{plot.RXuclars}
\alias{plot.RXuclars}
\title{Plot method for RXuclars objects}
\description{
 Plot trace displays for RXuclars regression coefficients.  The
 default is to display all five traces on one page with no legends.
}
\usage{
  \method{plot}{RXuclars}(x, trace = "all", trkey = FALSE, \dots)
}
\arguments{
  \item{x}{Output list object of class RXuclars.}
  \item{trace}{One of seven possible options: "all" to display 5 traces in one graph,
   "seq" to display 5 full-sized traces in sequence in response to user prompts,
   "coef" to display only the estimated shrunken beta coefficient trace,
   "rmse" to display only the estimated relative MSE risk trace,
   "exev" to display only the estimated excess MSE eigenvalue (OLS minus ridge) trace,
   "infd" to display only the estimated inferior direction cosine trace, or
   "spat" to display only the shrinkage (delta) factor pattern trace.}
  \item{trkey}{If TRUE, display a crude legend at the bottom of each trace plot.}
  \item{\dots}{Optional argument(s) passed on to plot().}
}
\value{NULL}
\author{Bob Obenchain <wizbob@att.net>}
\examples{
  data(longley2)
  form <- GNP~GNP.deflator+Unemployed+Armed.Forces+Population+Year+Employed
  rxuobj <- RXuclars(form, data=longley2)
  plot(rxuobj)
}
\keyword{trace}
