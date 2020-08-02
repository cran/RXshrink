"MLboot" <-  
function (form, data, reps=100, seed, rscale=1) 
{
    if (missing(form) || class(form) != "formula") 
        stop("First argument to MLboot must be a valid linear regression formula.")
    yvar <- deparse(form[[2]])
    if (missing(data) || !inherits(data, "data.frame")) 
        stop("Second argument to MLboot must be an existing Data Frame.")
    dfname <- deparse(substitute(data))
    if (!is.element(yvar, dimnames(data)[[2]])) 
        stop("Response variable in the MLboot formula must be an existing variable.")
    if (reps < 10) reps <- 10
    if (reps > 10000)
        cat("\nMore than 10,000 MLboot replications are NOT recommended...")
    if (missing(seed)) seed <- 1 + floor(1000 * runif(1))
    set.seed(seed)
    if (rscale != 1) rscale <- 2
    dmat <- as.matrix(data)
    samp <- dmat
    n <- length(data[,1])
    for (j in 1:reps) {
      if (j == 1) dfsam <- data        # use the given data rows...
      else {
        sam <- sample(1:n, replace=T)
        for (i in 1:n ) {
          samp[i,] <- dmat[sam[i],]    # sample from the given rows with replacement...
          } 
        dfsam <- data.frame(samp)	
      }
      out <- MLcalc(form, dfsam, rscale)
      if (out[1] == "Rsq1") {
        cat("Rsq1 exception at step: j =", j,"\n")
        sam <- sample(1:n, replace=T)
        for (i in 1:n ) {
          samp[i,] <- dmat[sam[i],]    # sample again from the given rows with replacement...
          } 
        dfsam <- data.frame(samp)	
        out <- MLcalc(form, dfsam, rscale)       # Second attempt at j-th sample...
        if (out[1] == "Rsq1") break              # Second consecutive attempt also failed...
        }	  
      if (j == 1) {
        ols.beta <- out$beta[1,]
        ols.rmse <- out$rmse[1,]
        opt.dmse <- out$dMSE
        opt.beta <- out$beta[2,]
        opt.rmse <- out$rmse[2,]
      }
      else {
        ols.beta <- rbind(ols.beta, out$beta[1,])
        ols.rmse <- rbind(ols.rmse, out$rmse[1,])
        opt.dmse <- rbind(opt.dmse, out$dMSE)
        opt.beta <- rbind(opt.beta, out$beta[2,])
        opt.rmse <- rbind(opt.rmse, out$rmse[2,])
      }
    }
    p <- out$p
    Totreps <- length(opt.beta[,1])
    RXolist <- list(data = dfname, form = form, reps = Totreps, seed = seed,
        p = p, n = n)
    RXolist <- c(RXolist, list(ols.beta = ols.beta, ols.rmse = ols.rmse, 
        opt.dmse = opt.dmse, opt.beta = opt.beta, opt.rmse = opt.rmse))
    class(RXolist) <- "MLboot"
    RXolist
}

"print.MLboot" <-
function (x, ...) 
{
    cat("\nMLboot Object: Resampling of data.frame Observations WITH Replacement...\n")
    cat("Data Frame:", x$data, "\n")
    cat("Regression Equation:\n")
    print(x$form)
    cat("\n    Number of Replications,       reps =", x$reps, "\n")
    cat("    Random Number Seed Value,     seed =", x$seed, "\n")
    cat("    Number of Predictors,            p =", x$p, "\n")
    cat("    Number of Observations,          n =", x$n, "\n\n")
    cat("OLS Beta Coefficient matrix            = ols.beta\n")
    cat("ML Optimally Biased Coefficient matrix = opt.beta\n")
    cat("OLS Relative MSE Risk matrix           = ols.rmse\n")
    cat("ML Optimally Biased Coefficient matrix = opt.rmse\n") 		
    cat("ML Shrinkage Delta-factor matrix       = opt.dmse\n\n")	
}

"MLhist" <- 
function(x, comp="opt.beta", xvar=1, npct=95, bins=50)
{
  if (missing(x) || !inherits(x, "MLboot"))
    stop("\nFirst argument to MLhist must be an existing MLboot object.\n\n")
  if (comp != "ols.beta" && comp != "ols.rmse" && comp != "opt.rmse" &&
    comp != "opt.dmse") comp <- "opt.beta"
  if (xvar <= 1) xvar <- 1
  if (xvar > x$p) xvar <- x$p
  if (npct < 66) npct <- 66
  if (npct > 100) npct <- 100
  npct <- round(npct)
  xout <- paste0("x$", comp, "[,", xvar, "]")
  xv <- eval(str2lang(xout))
  xlng <- length(xv)
  xost <- xv[order(xv)]                # xv order-statistics...
  noin <- round(xlng * npct / 100)     # chosen number of order-statistics
  nlo <- 1 + round((xlng-noin)/2)
  nup <- xlng - nlo + 1
  hst <- hist(xost[nlo:nup], breaks=bins, main = paste("Histogram of a MLboot Distribution"),
    xlab = paste(npct, "% of MLboot order-statistics"))
  abline(v=xv[1], col="blue", lwd=2, lty=2)  # Add location of observed sample estimate...
  RXolist <- list(x = xout, npct = npct, rbins = bins, dbins=hst$breaks, ntot = xlng, p = x$p,
    nlo = nlo, nup = nup, noin = noin, xmn = xv[1])
  class(RXolist) <- "MLhist"
  RXolist
}

"print.MLhist" <-
function (x, ...) 
{
    cat("\nMLhist Object: Characteristics of the displayed MLboot() Distribution...\n")
    cat("\n    MLboot() Distribution x$comp[,xvar]        : ", x$x, "\n")
    cat("\n    Percentage of Distribution Displayed, npct =", x$npct, "\n")
    cat("    Number of Histogram Bins requested,  rbins =", x$rbins, "\n")
    cat("    Number of Histogram Bins displayed,  dbins =", x$dbins, "\n")
    cat("    Total Number of Estimates available,  ntot =", x$ntot, "\n")
    cat("    First Order-Statistic Included,        nlo =", x$nlo, "\n")
    cat("    Last Order-Statistic Included,         nup =", x$nup, "\n")
    cat("    Number of >>Middle<< Estimates shown, noin =", x$noin, "\n")
    cat("    Observed Mean-Value of All Estimates,  xmn =", x$xmn, "\n\n")
}
