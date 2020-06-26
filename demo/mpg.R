  require(RXshrink)
  #
  # Input the Gasoline Mileage data of Hocking(1976)...
  data(mpg)
  #
  # Specify a linear regression model with four predictors of gas mileage...
  form <- mpg~cylnds+cubins+hpower+weight
  #
  # Fit this model using unrestricted generalized ridge regression (GRR)...
  rxuobj <- unr.ridge(form, data=mpg)
  #
  rxuobj
  # SCROLL ^^^ UP ^^^ to see PRINTED output from qm.ridge()...
  #
  plot(rxuobj)
  #
  # Define true parameter values...
  trugam <- matrix(c(-.5,-.1,.1,-.6),4,1)
  truvar <- 0.16
  #
  # Simulate a Correct Linear Model with IID errors from a Normal-distribution...
  trudat <- MLtrue(form, data=mpg, truc=trugam, truv=truvar)
  #
  form2 <- Yvec~cylnds+cubins+hpower+weight
  #  
  ideal <- unr.ridge(form2, data=trudat$new)
  #
  plot(ideal)
  #
  ################## End of "mpg" DEMO...
  