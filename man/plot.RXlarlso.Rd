\name{plot.RXlarlso}
\alias{plot.RXlarlso}
\title{Plot method for RXlarlso objects}
\description{Plot trace displays for RXlarlso regression coefficients.  The default
 is to display all five traces on one page with no legends.
}
\usage{
  \method{plot}{RXlarlso}(x, trace = "all", trkey = FALSE, \dots)
}
\arguments{
  \item{x}{Output list object of class RXlarlso.}
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
\author{Bob Obenchain <softrx@iquest.net>}
\examples{
  data(longley2)
  form <- GNP~GNP.deflator+Unemployed+Armed.Forces+Population+Year+Employed
  rxlobj <- RXlarlso(form, data=longley2)
  plot(rxlobj)
}
\keyword{hplot}
\keyword{methods}