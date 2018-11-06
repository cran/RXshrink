\name{RXtsimu}
\alias{RXtsimu}
\title{True Squared Error LOSS of Shrinkage for a Simulated Response Y-vector}
\description{
  By specifying numerical values for regression parameters (uncorrelated components and error
  sigma) that usually are unknown, this function allows the user to simulate response data and
  display the True Squared Error Loss associated with shrinkage along a given \bold{Q}-shaped path.  
}
\usage{
  RXtsimu(form, data, trug, trus, Q = 0, rscale = 1, steps = 8,
    seed = 0123, qmax = 5, qmin = -5)
}
\arguments{
  \item{form}{A regression formula [y~x1+x2+...] suitable for use with lm().}
  \item{data}{Data frame containing observations on all variables in the formula.}
  \item{trug}{Column vector of numerical values for the true uncorrelated components
   of the regression coefficient vector.}
  \item{trus}{Numerical value for the true error standard deviation, Sigma.}
  \item{Q}{Numerical value for the shape parameter controling the shrinkage path curvature.
   Default shape is Q = 0 for Hoerl-Kennard "ordinary" ridge regression.}
  \item{rscale}{One of three possible choices (0, 1 or 2) for rescaling of variables
   as they are being "centered" to remove non-essential ill-conditioning:  0 implies no
   rescaling; 1 implies divide each variable by its standard error; 2 implies rescale as
   in option 1 but re-express answers as in option 0.}
  \item{steps}{Number of equally spaced values per unit change along the horizontal
   M-extent-of-shrinkage axis where estimates are calculated and displayed in TRACEs
   (default = 8.)}
  \item{seed}{Seed for random number generator.  To get two different simulated response
   vectors and different sets of coefficients and SE losses, invoke the RXtsimu() function
   twice with different seed values. (default = 0123.)}
  \item{qmax}{Maximum allowed Q-shape (default = +5.)}
  \item{qmin}{Minimum allowed Q-shape (default = -5.)}
}
\details{The RXridge() function calculates maximum likelihood estimates (corrected, if
  necessary, so as to have correct range) for typical statistical inference situations where
  regression parameters are unknowns.  In sharp contrast with this usual situation, the
  RXtsimu() function shows exactly how regression coefficient estimates and their true
  Squared Error Losses change with shrinkage for a simulated response Y-vector generated
  using given numerical values for regression parameters.  In fact, it is interesting to
  compare the output from RXtrisk() and RXtsimu() for given regression parameters to the
  corresponding output from RXridge() on the data.frame in which any original response
  Y-vector has been replaced by the ydat object from the RXtsimu() output list. 
}
\value{ An output list object of class RXtsimu:
  \item{form}{The regression formula specified as the first argument.}
  \item{data}{Name of the data.frame object specified as the second argument.}
  \item{trug}{Vector of numerical values for the true uncorrelated gamma components.}
  \item{trus}{Numerical value for the true error standard deviation, Sigma.}
  \item{qp}{Numerical value of the Q-shape actually used for shrinkage.}
  \item{p}{Number of regression predictor variables.}
  \item{n}{Number of complete observations after removal of all missing values.}
  \item{prinstat}{Listing of principal statistics.}
  \item{ydat}{Matrix with simulated Y-responses in its first column and the true expected
   values of those responses in its second column.}
  \item{coef}{Matrix of shrinkage-ridge regression coefficient estimates.}
  \item{rsel}{Matrix of true relative SE losses in regression coefficient estimates.}
  \item{spat}{Matrix of shrinkage pattern multiplicative delta factors.}
  \item{sext}{Listing of summary statistics for all M-extents-of-shrinkage.}
}
\author{Bob Obenchain <wizbob@att.net>}
\seealso{\code{\link{RXridge}} and \code{\link{RXtrisk}}.}
\examples{
  data(haldport)
  form <- heat~p3ca+p3cs+p4caf+p2cs
  rxrobj <- RXridge(form, data=haldport)
  plot(rxrobj)
  # define true parameter values.
  trugam <- matrix(c(.8,.0,.3,.5),4,1)
  trusig <- 0.2
  # create true shrinkage MSE risk scenario.
  trumse <- RXtrisk(form, data=haldport, trugam, trusig, Q=-5)
  # calculate true shrinkage squared error losses.
  trusim <- RXtsimu(form, data=haldport, trugam, trusig, Q=-5)
  haldpsim <- haldport
  haldpsim[,5] <- trusim$ydat[,1]
  rxsobj <- RXridge(form, data=haldpsim) # analysis as if parameters unknown
  plot(rxsobj)
}
\keyword{regression}
