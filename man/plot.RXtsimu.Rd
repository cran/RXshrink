\name{plot.RXtsimu}
\alias{plot.RXtsimu}
\title{Plot method for RXtsimu objects}
\description{Plot trace displays for RXtsimu simulated regression coefficients.  The
 default is to display all three traces on one page with no legends.
}
\usage{
  \method{plot}{RXtsimu}(x, trace = "all", trkey = FALSE, \dots )
}
\arguments{
  \item{x}{Output list object of class RXtsimu.}
  \item{trace}{One of five possible options: "all" to display 3 traces in one graph,
   "seq" to display 3 full-sized traces in sequence in response to user prompts,
   "coef" to display only the estimated shrunken beta coefficient trace,
   "rsel" to display only the relative squared error loss trace, or
   "spat" to display only the shrinkage (delta) factor pattern.}
  \item{trkey}{If TRUE, display a crude legend at the bottom of each trace plot.}
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
  rxsobj <- RXtsimu(form, data=longley2, trugam, trusig, Q=-1.5)
  plot(rxsobj)
}
\keyword{hplot}
\keyword{methods}