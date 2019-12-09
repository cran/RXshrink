\name{plot.true.risk}
\alias{plot.true.risk}
\title{Plot method for true.risk objects}
\description{Plot trace displays for the "true" (known components) risk in coefficient
 estimates.  The default is to display all five traces on one page with no legends.
}
\usage{
  \method{plot}{true.risk}(x, trace = "all", trkey = FALSE, \dots )
}
\arguments{
  \item{x}{Output list object of class true.risk.}
  \item{trace}{One of seven possible options: "all" to display 5 traces in one graph,
   "seq" to display 5 full-sized traces in sequence in response to user prompts,
   "coef" to display only the estimated beta coefficient trace,
   "rmse" to display only the estimated relative MSE risk trace,
   "exev" to display only the estimated excess MSE eigenvalue (OLS minus larlso),
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
  # set true regression parameter values not unlike their OLS estimates from RXridge().
  trugam <- matrix(c(.5,-.1,.2,.2,.2,-.2),6,1)
  trusig <- 0.04
  rxtobj <- true.risk(form, data=longley2, trugam, trusig, Q=-1.5)
  plot(rxtobj)
}
\keyword{trace}
