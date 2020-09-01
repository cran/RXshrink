\name{plot.unr.biv}
\alias{plot.unr.biv}
\title{Plot method for unr.biv objects}
\description{Plot Bivariate Confidence Ellipse(s) or a TRACE display for unrestricted
 ridge shrinkage coefficients.  The default is to display two ellipses at confidence
 levels 0.95 and 0.50 on a single plot for the first 2 of p >= 2 coefficients. The
 projection of the unr.ridge() Path from the OLS solution (m == 0) to the Minimum
 MSE Risk solution (m == mUnr) onto this bivariate-plane is also shown in "red".
}
\usage{
  \method{plot}{unr.biv}(x, type = "ellip", ...)
}
\arguments{
  \item{x}{Output list object of class unr.ridge.}
  \item{type}{One of 2 options: "ellip" or "trace". The default option of "ellip"
   displays the Confidence Ellipse(s) specified by arguments to unr.biv().}
  \item{\dots}{Optional argument(s) passed on to plot().}
}
\value{NULL}
\author{Bob Obenchain <wizbob@att.net>}
\examples{
  # Cost-Effectiveness inferences using Linear Models and GRR...
  \dontrun{
    # library(ICEinfer)
    data(sepsis)
    ndr <- ICEpref(sepsis$qalypres, sepsis$totcost, lambda=50000, beta=0.1)
    sndr <- data.frame(cbind(ndr, sepsis)) # ndr: non-linear diminishing returns...
    form4 <- ndr ~ icu + age + orgfails + apache
    usndra <- unr.aug(unr.ridge(form4, sndr)) # compare ndr of 2 Intensive Care Units...
    plot(urobj <- unr.biv(usndra, 2, 4))
    urobj               # implicit print...
    # plot(urobj, type = "tr")
  }
}
\keyword{trace}
