  require(RXshrink)
  #
  # Input the Gasoline Mileage data of Hocking(1976)...
  data(mpg)
  #
  # Specify a linear regression model with four predictors of gas mileage...
  form <- mpg~cylnds+cubins+hpower+weight
  #
  # Fit this model using 2-parameter Generalized Ridge Regression...
  #
  rxrobj <- qm.ridge(form, data=mpg)
  #
  # Press ENTER to [1] Print all qm.ridge() Summary Statistics and
  # [2] Display all 5-types of generalized ridge TRACE plots...
  #
  scan()
  # Print...
  rxrobj
  # Plot...
  plot(rxrobj)
  #  
  # Press ENTER again to [1] Declare TRUE vales for Beta Coefficients and
  # the Standard Deviation of Errors, then [2] Calculate and Plot TRUE
  # MSE Risks associated with QM (2-parameter) Shrinkage...
  scan()
  #
  # Define true parameter values...
  trugam <- matrix(c(-.5,-.1,.1,-.6),4,1)
  trusig <- 0.4
  #
  # Create true shrinkage MSE risk scenario.
  trumse <- true.risk(form, data=mpg, trugam, trusig, Q=-1, steps=4)
  #
  # Print Summary Statistics...
  trumse
  #
  plot(trumse)
  #
  # Finally, Press ENTER to Simulate, Print and Plot TRUE Squared Error LOSS...
  scan()
  trusel <- true.simu(form, data=mpg, trugam, trusig, Q=-1, steps=4)
  # Print...
  trusel
  # Plot...
  plot(trusel)
  #
  # END of demo(mpg)...
  #
  