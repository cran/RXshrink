###################################################
#  Reference:
#
#  Carl Morris. "Was Ty Cobb ever a TRUE .400 hitter?"  ASA JSM Lecture, August 18, 1982. Cincinnati.
#
#  CMspl = Carl Morris' Piecewise Linear Spline term with "knot" at year = 1910 (Cobb's 6th Season.)
#  seasons = 1, 2, ..., 24 will be a pure linear trend term in the fitted model.
#
#  This RXshrink demo() illustrates use of the RXpredict() function.
#
require(RXshrink)
data(tycobb)
str(tycobb)
#
# formula suitable for use with lm() for p=3 x-variables...
form <- batavg~atbats+seasons+CMspl
#
# Fit a generalized linear regression (GRR) model using unr.ridge()...
tycuobj <- unr.ridge(form, data=tycobb, steps=100)
# tycuobj
#     Only the 5 final lines of printed output are listed below...
# 
#     Most Likely UNRestricted Shrinkage Extent, mUnr = 1.039595
#     Corresponding -2*log(LikelihoodRatio) statistic = 0.0
#     Most Likely m-Value on the Lattice,        mClk = 1.04 <note: steps=100 above>
#     Smallest Observed -2*log(LikelihoodRatio), minC = 0.0003832262
#     dMSE Estimates = 0.01700621 0.9698978 0.973501 
#
plot(tycuobj)    # Show all 5 unr.ridge() TRACE Diagnostic plots...
#
# Display the first 3 "k-star" values for "knots"; final "knot" at m=4 & "k-star"=0...
rep(1,3)/tycuobj$dMSE
#
testLS <- RXpredict(tycuobj, data=tycobb, m=0)     # OLS fit at m == 0
testML <- RXpredict(tycuobj, data=tycobb)          # ML  fit for tycuobj...
testSZ <- RXpredict(tycuobj, data=tycobb, m=3)     # All Coefficients Shrunken to Zero...
#
ym <- mean(tycobb$batavg)              # 0.3610738
ys <- sqrt(var(tycobb$batavg))         # 0.03848413
# Calculate value of batavg == 0.400 on the "cry" scale...
crx400 <- ( 0.4 - ym )/ys              # 1.011486  ...slightly greater than 1.00...
#
plot( tycobb$year, testLS$cry, ann = FALSE, type = "b") 
lines(tycobb$year, testLS$cryprd, lty=2, lwd=2, col="blue")
title(main="Ty Cobb's Batting Averages and Fitted Values",
  xlab="Year", ylab="Centered and Rescaled Batting Averages")
lines(tycobb$year, testML$cryprd, lty=3, lwd=2, col="green")
abline(h = crx400, col = "red", lty = 2, lwd = 2)
abline(h = 0, col = "gray", lty = 1, lwd = 1)
# Next, some additional plot() annotation...
text(x=1907, y=+1.2, labels="400 Batting", col="red")
text(x=1918, y=-1.5, labels="_o_ Ty Cobb's Season Averages", col="black")
text(x=1918, y=-1.8, labels="_ _ OLS fitted Averages", col="blue")
text(x=1918, y=-2.1, labels="... ML optimally biased Averages", col="darkgreen") 
#
maxLS <- max(testLS$yvecprd)  # 0.4000915 >> Yes, Ty Cobb was a "True 400-hitter ...in 1911 !!! 
maxML <- max(testML$yvecprd)  # 0.3994661 >> NO, not even in 1911. Sorry to any Ty Cobb fans out there !!!
#
################## End of "tycobb" DEMO...
