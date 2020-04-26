\name{tycobb}
\alias{tycobb}
\docType{data}
\title{Ty Cobb batting statistics for 1905--1928 and Carl Morris' Linear Spline.}
\description{Linear Regression models can be used to predict Ty Cobb's Expected true
yearly batting averages from his observed yearly averages, batavg, from the 5 other
variables stored in this "tycobb" data.frame. The predictions from such "models" can
address the question: "Was Ty Cobb ever a TRUE .400 hitter?" Since a player's
seasonal batavg is a "random variable," the fact that Cobb's batavg was 0.4196 in
1911 and 0.4105 in 1912 does NOT imply that his Expected Seasonal batavg was truly
over .400 in either of these consecutive seasons. For example, his batavg was
0.4011 in 1922, which was 10 years later. But his seasonal batavg had dipped to
"only" 0.3341 in 1920. Revealing true "time-series" features using the tycobb
data.frame requires more than just using "Weighted" Least-Squares estimation with a
24x24 Variance - Covariance matrix ASSUMED to be "non-diagonal."  Major League
Baseball Hall of Fame executive, Branch Rickey: "Cobb lived off the field as though
he wished to live forever. He lived on the field as though it was his last day." 
}
\usage{data(tycobb)}
\format{
  A data frame with 24 observations (years) on the following 6 variables.
  \describe{
    \item{\code{year}}{Ty Cobb's 24 American League Seasons: 1905 through 1928.}
    \item{\code{hits}}{Total number of Hits that season.}
    \item{\code{atbats}}{Total number of times at Bat that season.}
    \item{\code{CMspl}}{Carl Morris' Linear Spline term with "knot" at 1910, season 6.}
    \item{\code{seasons}}{Linear "Trend" term: 1, 2, ..., 24.}
    \item{\code{batavg}}{Cobb's Seasonal Batting Average ...accurate to 7 decimal places.} 	
  }
}
\references{
Carl Morris (1982) "Was Ty Cobb ever a TRUE .400 hitter?" One-page Handout for his JSM
    Lecture on August 18 in Cincinnati, Ohio.
}
\keyword{datasets}