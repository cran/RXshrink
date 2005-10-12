\name{plot.RXtrisk}
\alias{plot.RXtrisk}
\title{Plot method for RXtrisk objects}
\description{
Plot trace displays of RXtrisk shrinkage statistics.  The default is to display all
 five traces on one page with no legends.
}
\usage{plot.RXtrisk(x, trace = "all", trkey = FALSE, \dots )}
\arguments{
  \item{x}{Output list object of class RXtrisk.}
  \item{trace}{One of seven possible options: "all" to display 5 traces in one graph,
   "seq" to display 5 full-sized traces in sequence in response to user prompts,
   "coef" to display only the estimated beta coefficient trace,
   "rmse" to display only the estimated relative MSE risk trace,
   "exev" to display only the estimated excess MSE eigenvalue (OLS minus larlso),
   "infd" to display only the estimated inferior direction cosine trace, or
   "spat" to display only the shrinkage (delta) factor pattern trace.}
  \item{trkey}{If TRUE, a crude legend is displayed at the bottom of each trace plot.}
  \item{\dots}{Optional argument(s) passed on to plot().}
}
\value{NULL}
\author{Bob Obenchain <softrx@iquest.net>}
\examples{
    data(longley2)
    form <- GNP~GNP.deflator+Unemployed+Armed.Forces+Population+Year+Employed
    # set true regression parameter values not unlike their OLS estimates from RXridge().
    trugam <- matrix(c(.5,-.1,.2,.2,.2,-.2),6,1)
    trusig <- 0.04
    rxtobj <- RXtrisk(form, data=longley2, trugam, trusig, Q=-1.5)
    plot(rxtobj)
}
\keyword{hplot}
\keyword{methods}