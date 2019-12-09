  require(RXshrink)
  # Input revised Longley dataset of Hoerl(2000).
  data(longley2)
  # Specify form of (linear) regression model...
  form <- GNP~GNP.deflator+Unemployed+Armed.Forces+Population+Year+Employed
  # Fit model of this form using 2-parameter Generalized Ridge Regression...
  rxrobj <- qm.ridge(form, data=longley2)
  rxrobj
  names(rxrobj)
  plot(rxrobj)
  cat("\n Press ENTER for Unrestricted Shrinkage model fitting...")
  scan()
  # Fit model of same form using the Shrinkage Path that passes through
  # the Unrestricted Maximum Likelihood estimate... 
  rxuobj <- unr.ridge(form, data=longley2)
  rxuobj
  names(rxuobj)
  plot(rxuobj)
  cat("\n Press ENTER for Coefficients with guarenteed Correct SIGNS...")
  scan()
  # Show Coefficients with Minimum MSE Risk Parallel to unknown, true Beta...
  rxcsobj <- correct.signs(form, data=longley2)
  rxcsobj