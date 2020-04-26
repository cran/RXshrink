\name{plot.RXpredict}
\alias{plot.RXpredict}
\title{Plot method for RXpredict objects}
\description{Plot Predicted and/or Fitted.Values for all 5 RXshrink regression estimation
methods. The default is to plot Predictions for the y-Outcome variable.
}
\usage{
  \method{plot}{RXpredict}(x, fit = "yvecprd", \dots )
}
\arguments{
  \item{x}{Output list object of class RXpredict.}
  \item{fit}{One of three possible options:,
   "yvecprd" to display Predictions of the Observed y-Outcomes in a single plot,
   "cryprd" to display Fitted.Values for the Centered and Rescaled y-Outcomes,
   "both" to display "yvecprd" and "cryprd" plots as rows on one page.}
  \item{\dots}{Optional argument(s) passed on to plot().}
}
\value{NULL}
\author{Bob Obenchain <wizbob@att.net>}
\examples{
  data(longley2)
  form <- Employed~GNP+GNP.deflator+Unemployed+Armed.Forces+Population+Year
  rxuobj <- unr.ridge(form, longley2)
  rxuprd <- RXpredict(rxuobj, longley2)
  plot(rxuprd)
  # Clearly Biased predictions can still represent "Optimal" Variance-Bias Trade-Offs...
}
\keyword{predict}