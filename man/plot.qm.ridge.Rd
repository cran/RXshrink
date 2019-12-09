\name{plot.qm.ridge}
\alias{plot.qm.ridge}
\title{Plot method for qm.ridge objects}
\description{Plot TRACE displays for 2-parameter (q, m) generalized ridge coefficients.  The
 default is to display all five traces on one page with no legends.
}
\usage{
  \method{plot}{qm.ridge}(x, trace = "all", trkey = FALSE, \dots )
}
\arguments{
  \item{x}{Output list object of class qm.ridge.}
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
  rxrobj <- qm.ridge(form, data=longley2)
  plot(rxrobj)
}
\keyword{trace}
