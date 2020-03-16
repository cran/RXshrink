  require(RXshrink)
  # Read in the Portland Cement data.frame...
  data(haldport)
  # Next, define a regression model formula...
  form <- heat~p3ca+p3cs+p4caf+p2cs
  # Then perform QM (2-parameter) generalized Ridge Shrinkage...
  rxrobj <- qm.ridge(form, data=haldport)
  # Print Summary Statistics...
  rxrobj
  # Display generalized Ridge TRACE plots...
  plot(rxrobj)
  #
  # Press ENTER to see corresponding results for
  # the Unrestricted Path...
  #
  scan()
  rxuobj <- unr.ridge(form, data=haldport)
  rxuobj
  plot(rxuobj)
  #
  # Press ENTER to see results for Least Angle Regression that
  # are Augmented with 5 more generalized TRACE displays...
  #
  scan()
  rxaobj <- aug.lars(form, data=haldport)
  # Print Summary and Plot TRACES...
  rxaobj
  plot(rxaobj)
  #
  # Press ENTER to see the lars() fit to all 4 UNCORRELATED
  # X-COMPONENTS plus their generalized ridge TRACE displays...
  #  
  scan()
  rxcobj <- uc.lars(form, data=haldport)
  # Print Summary and Plot TRACES...
  rxcobj
  plot(rxcobj)
  #
  # END of demo(haldport)...
  #
  
  