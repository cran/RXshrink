  library(RXshrink)
  # Input the "revised" Longley dataset of Art Hoerl(1979)...
  #
  data(longley2)
  #
  # Specify the "formula" for a regression model...
  #
  form <- GNP~GNP.deflator+Unemployed+Armed.Forces+Population+Year+Employed
  #
  # Fit this model using 2-parameter Generalized Ridge Regression...
  #
  rxqmobj <- qm.ridge(form, data=longley2)
  #
  # Now [1] Print all qm.ridge() Summary Statistics
  # and [2] Display all 5-types of generalized ridge TRACE plots...
  #
  rxqmobj
  # SCROLL ^^^ UP ^^^ to see the PRINTED output from qm.ridge()...
  #
  plot(rxqmobj)
  #
  # Next, show the corresponding results for the Piecewise-Linear
  # Spline PATH passing through the Unrestricted Maximum Likelihood
  # point-estimate of "MSE Optimal" Regression Coefficients...
  #
  rxunobj <- unr.ridge(form, data=longley2)
  rxunobj
  # SCROLL ^^^ UP ^^^ to see the PRINTED output from unr.ridge()...
  #
  plot(rxunobj)
  #
  # Finally, print information about Beta-coefficient estimates with
  # Guaranteed "Correct" SIGNS...
  #
  # These Coefficients have Minimum MSE Risk in the "Unknown"
  # direction PARALLEL to true Beta vector...
  rxcsobj <- correct.signs(form, data=longley2)
  rxcsobj
  #
  # The PLOT didn't change, only the above TEXT is NEW...
  #
  ################## End of "longley2" DEMO...
  