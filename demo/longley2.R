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
  # Press ENTER to [1] Print all qm.ridge() Summary Statistics and
  # [2] Display all 5-types of generalized ridge TRACE plots...
  #
  scan()
  rxrobj
  plot(rxrobj)
  #
  # Press ENTER again to see the corresponding results for the
  # Linear Spline PATH passing the Unrestricted Maximum Likelihood
  # point-estimate of MSE Optimal Regression Coefficients...
  #
  scan()
  rxuobj <- unr.ridge(form, data=longley2)
  rxuobj
  plot(rxuobj)
  #
  # Finally, Press ENTER to display fitted Coefficients with
  # Guaranteed "Correct" SIGNS...
  # 
  scan()
  # These Coefficients have Minimum MSE Risk in the unknown 
  # direction PARALLEL to true Beta vector...
  rxcsobj <- correct.signs(form, data=longley2)
  rxcsobj
  #
  # END of demo(longley2)...
  #