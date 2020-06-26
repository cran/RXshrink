  require(RXshrink)
  # Input the "revised" Longley dataset of Hoerl(2000)...
  #
  data(longley2)
  #
  # Specify the "formula" for a regression model...
  #
  form <- GNP~GNP.deflator+Unemployed+Armed.Forces+Population+Year+Employed
  #
  # Fit this model using 2-parameter Generalized Ridge Regression...
  #
  rxrobj <- qm.ridge(form, data=longley2)
  #
  # Now [1] Print all qm.ridge() Summary Statistics
  # and [2] Display all 5-types of generalized ridge TRACE plots...
  #
  rxrobj
  # SCROLL ^^^ UP ^^^ to see PRINTED output from qm.ridge()...
  #
  plot(rxrobj)
  #
  # Next, see the corresponding results for the Piecewise-Linear
  # Spline PATH passing through the Unrestricted Maximum Likelihood
  # point-estimate of "MSE Optimal" Regression Coefficients...
  #
  rxuobj <- unr.ridge(form, data=longley2)
  rxuobj
  # SCROLL ^^^ UP ^^^ to see PRINTED output from unr.ridge()...
  #
  plot(rxuobj)
  #
  # Finally, print Beta-coefficient estimates with
  # Guaranteed "Correct" SIGNS...
  #
  # These Coefficients have Minimum MSE Risk in the "Unknown"
  # direction PARALLEL to true Beta vector...
  rxcsobj <- correct.signs(form, data=longley2)
  rxcsobj
  #
  ################## End of "longley2" DEMO...
  