\name{haldport}
\alias{haldport}
\docType{data}
\title{Portland Cement benchmark of Hald(1952).}
\description{Heat evolved during setting of 13 cement mixtures of four basic ingredients.  Each
 ingredient percentage appears to be rounded down to a full integer.  The sum of the four mixture
 percentages varies from a maximum of 99\% to a minimum of 95\%.  If all four regressor
 X-variables always summed to 100\%, the centered X-matrix would then be of rank only 3.  Thus, the
 regression of heat on four X-percentages is ill-conditioned, with an approximate rank deficiency
 of MCAL = 1. 
}
\usage{data(haldport)}
\details{The (RXshrink) haldport data are identical to the (MASS) cement data except for
 variable names.
}
\format{
  A data frame with 13 observations on the following 5 variables.
  \describe{
    \item{\code{p3ca}}{Integer percentage of 3CaO.Al2O3 in the mixture.}
    \item{\code{p3cs}}{Integer percentage of 3CaO.SiO2 in the mixture.}
    \item{\code{p4caf}}{Integer percentage of 4CaO.Al2O3.Fe2O3 in the mixture.}
    \item{\code{p2cs}}{Integer percentage of 2CaO.SiO2 in the mixture.} 
    \item{\code{heat}}{Heat (cals/gm) evolved in setting, recorded to nearest tenth.}
  }
}
\source{
Woods H, Steinour HH, Starke HR. "Effect of composition of Portland cement on heat
 evolved during hardening. Industrial Engineering and Chemistry 1932; 24: 1207-1214.
}
\references{
Hald A. Statistical Theory with Engineering Applications. 1952 (page 647.) New York; Wiley.
}
\keyword{datasets}