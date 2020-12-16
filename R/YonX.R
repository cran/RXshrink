"YonX" <- function (form, data, delmax = 0.999999) 
{
    if (missing(form) || class(form) != "formula") 
        stop("First argument to YonX() must be a valid linear regression formula.")
    yvar <- deparse(form[[2]])
    xvar <- deparse(form[[3]])
    if (missing(data) || !inherits(data, "data.frame")) 
        stop("Second argument to YonX() must be an existing Data Frame.")
    dfname <- deparse(substitute(data))
    if (!is.element(yvar, dimnames(data)[[2]])) 
        stop("Response Y-variable in the YonX() formula must be an existing variable.")
    if (!is.element(xvar, dimnames(data)[[2]])) 
        stop("Regressor X-variable in the YonX() formula must be an existing variable.")
    lmobj <- lm(form, data)
    yvec <- as.matrix(lmobj$model[,1])
    xvec <- as.matrix(lmobj$model[,2])
    p <- ncol(xvec)
    if (p > 1) 
        stop("The formula must contain only 1 (non-constant) Regressor X-variable.")
    n <- length(xvec)
    if (n < 5) 
        stop("Number of observations must be at least 5.")
    mx <- mean(xvec)
    crx <- xvec - mx * matrix(1, n, 1)	
    sx <- svd(crx)                               # Singular Value Decomposition of Centered xvec
    if (sx$v == -1) sx$u <- -sx$u         
    eigval <- sx$d^2
    cry <- yvec - mean(yvec) * matrix(1, n, 1)   # n x 1
    bstar <- as.numeric(t(sx$u) %*% ( cry / sx$d ))
    risk <- 1/eigval
    exev <- 0
    delta <- 1
    d <- 1
    cold <- 1
    sv <- sx$d
    ssy <- as.numeric(t(cry) %*% cry)
    rho <- (sv * bstar)/sqrt(ssy)
    arho <- abs(rho)
    r2 <- arho^2             # OLS "R^2" statistic....
    if (r2 >= 1) 
        stop(" Maximum Likelihood Shrinkage is not applicable when RSQUARE=1.")
    res <- cry - crx %*% bstar
    s2 <- as.numeric(t(res) %*% res/(n - 2))
    varrho <- s2/ssy
    tstat <- rho/sqrt(varrho)
    frat <- tstat^2
    stat <- cbind(eigval, sv, bstar, rho, tstat)
    dimnames(stat) <- list(1, c("LAMBDA", "SV", "COMP", "RHO", "TRAT"))
    yxnam <- c(yvar,xvar)
    RXolist <- list(dfname = dfname, form = form, p = 1, n = n, 
        r2 = r2, s2 = s2, prinstat = stat, yxnam = yxnam, yvec = yvec, xvec = xvec)
    OmR2dN <- (1 - r2)/n
    dMSE <- rho^2 / (rho^2 + OmR2dN )    # Maximum-Likelihood estimate...
    mMSE <- 1 - dMSE         # is positive...
    mcal <- 0
    kinc <- 1/dMSE           # innitially > 1
    kstar <- kinc
    const <- (n - 4)/(n - 2)
    srat <- tstat / sv
    MCAL <- 0
    KSTAR <- kinc
    C <- Inf
    E <- Inf
    R <- Inf
    steps <- 100
    leng <- steps + 1
    MCAL <-  rep(MCAL, leng)
    KSTAR <- rep(KSTAR, leng)
    C <- rep(C, leng)
    E <- rep(E, leng)
    R <- rep(R, leng)
    bstar <- rep(bstar, leng)
    risk <- rep(risk, leng)
    exev <- rep(exev, leng)
    delta <- rep(delta, leng)
    kstar <- rep(kstar,leng)
    mcal <- rep(mcal, leng)
    term <- n * log((1 - r2)/n)
    for (inc in 1:steps) {
        ip <- inc + 1
        mobj <- inc/steps
        iter <- kofm1(mobj, dMSE)  # Function defined below...
        kinc <- iter$kStar
        d <- iter$d
        omd <- max(1 - d, 1 - delmax)   # Strictly Positive values...
        ddomd <- d/omd
        rxi <- arho * sqrt(ddomd)
        slik <- 2/(rxi + sqrt(4 * n + rxi^2))
        clik <- 2 * n * log(slik) + ddomd - (rxi/slik) - term
        if (clik < 1 - delmax) clik <- 1 - delmax
        ebay <- frat * omd - log(omd)
        sr2d <- d * rho^2
        if( sr2d >= 1 ) {rcof <- Inf} else {
            rcof <- abs(log(omd)) + n * log((1 - sr2d)/(1 - r2))
        }
        minc <- 1 - d
        MCAL[ip] <- minc
        KSTAR[ip] <- kinc
        C[ip] <- clik
        E[ip] <- ebay
        R[ip] <- rcof
        binc <- d * bstar[1]
        vecr <- (1 - d) * srat
        risk[ip] <- const * vecr^2 + (2*d - 1)/eigval 
        rinc <- max(risk[ip], d^2/eigval)
        emse <- rinc - risk[ip]
        cinc <- 1
        if (risk[ip] > rinc) cinc <- as.numeric(NA)
        sfac <- abs(emse)/100
        if (sfac < 1e-05) sfac <- 1e-05
        eign <- emse/sfac
        einc <- eign*sfac
        if (is.na(einc)) 
            einc <- 0
        if (einc + 1 - delmax < 0) {
            cinc <- eign
            cold <- cinc
        }
        bstar[ip] <- binc
        exev[ip] <- einc
        delta[ip] <- d
        kstar[ip] <- kinc
        mcal[ip] <- minc
    }
    mlik <- cbind(MCAL, KSTAR, C, E, R)
    dimnames(mlik) <- list(1:leng, c("M", "KSTAR", "CLIK", "EBAY", "RCOF"))
    sext <- cbind(risk, kstar, mcal)
    dimnames(sext) <- list(1:leng, c("TSMSE", "KSTAR", "MCAL"))
    mUnr <- mofk1(k=1, dMSE)
    minC <- min(mlik[,3])
    minE <- min(mlik[,4])
    minR <- min(mlik[,5])
    Good <- max(0,2*dMSE-1)            # Potentially WRONG or "misleading" ???
    Phi2 <- dMSE/(1-dMSE)
    minRMSE <- min(risk)               # Minimum estimated MSE risk...
    eqlRMSE <- risk[1]                 # estimated risk = relative variance of OLS estimate... 
    mClk <- 0; mRmin <- 0; mReql <- 0	
    for (inc in 2:leng) {
        if (mlik[inc,3] <= minC) {mClk <- (inc-1)/steps}
        if (risk[inc] <= minRMSE) {mRmin <- (inc-1)/steps} 
        if (risk[inc] <= eqlRMSE) {mReql <- (inc-1)/steps} 		
    }	
    RXolist <- c(RXolist, list(coef = bstar, rmse = risk, spat = delta,
        mlik = mlik, sext = sext, mUnr = mUnr, Good = Good, mMSE = mMSE,
        mClk = mClk, minC = minC, minE = minE, minR = minR, mRmin = mRmin,
        mReql = mReql, Phi2 = Phi2, dMSE = dMSE))
    class(RXolist) <- "YonX"
    RXolist
}

"kofm1" <- function (muobj, dMSE, delmax = 0.999999) 
{
    kM <- 1/dMSE 
    if (muobj <= 0) {
        d <- delmax         # d strictly < 1...
        kStar <- kM
        return(list(kStar = kStar, d = d))
    }
    if (muobj >= 1) {
        d <- 0
        kStar <- 0                     # Terminus of the shrinkage path
        return(list(kStar = kStar, d = d))
    }
    mVar <- mofk1(k=kM, dMSE)
    kMin <- (1-muobj)/dMSE
    if( muobj > mVar ) {               # Tail; Large m, small k Cases...
        d <- kMin*dMSE
        kStar <- kMin
        return(list(kStar = kStar, d = d))
    } else {                           # Locally Linear cases... 
        if( mVar < muobj ) {
            B <- muobj - mVar
            D <- 1 - muobj 
            kStar <- (B*D/(B+D))*kM/B
        }
        dp <- min(delmax,kStar*dMSE)
        return(list(kStar = kStar, d = dp))		
    }
}

"mofk1" <- function(k, dMSE)  # Many-to-One Function for large k-values...
{
  if( k < 0 ) k <- 0
  kM <- 1/dMSE
  if( k > kM ) k <- kM
  m <- 1 - min(1,k*dMSE)
  m
}

"mapp" <- function(mVal, YXobj)  # Approximate index on YonX() lattice of m = mVal...
{
  steps <- length(YXobj$spat)-1
  m <- 1  
  if (mVal <= 0) return(m)
  m <- steps+1
  if (mVal >= 1) return(m)
  m <- round(mVal*steps,0)+1   
}

"print.YonX" <- function (x, ...) 
{
    cat("\nYonX Object: Shrinkage in Simple Regression [One X-variable]\n")
    cat("Data Frame:", x$dfname, "\n")
    cat("Regression Equation [with Implicit Intercept]:\n")
    print(x$form)
    cat("\n    Number of Observations, n =", x$n, "\n")
    cat("\nPrincipal Axis Summary Statistics of Ill-Conditioning...\n")
    print.default(x$prinstat, quote = FALSE)
    cat("\n    Residual Mean Square for Error =", x$s2, "\n")
    cat("    Estimate of Residual Std. Error =", sqrt(x$s2), "\n")
    cat("\nThe extent of shrinkage (M value) most likely to be optimal\n")
    cat("depends upon whether one uses the Classical, Empirical Bayes, or\n")
    cat("Random Coefficient criterion.  In each case, the objective is to\n")
    cat("minimize the minus-two-log-likelihood statistics listed below:\n")
    print.default(x$mlik, quote = FALSE)
    cat("\nExtent of Shrinkage Statistics...\n")
    print.default(x$sext, quote = FALSE)
    if (x$r2 > 0.999) {
        cat("\n\nWARNING! R-squared exceeds 0.999;")
        cat("  Details of calculations are possibly Misleading...\n")
    }
    b0 <- x$prinstat[3]
    bm <- b0 * x$dMSE	
    cat("\n    OLS Beta Coefficient (BLUE),                bHat =", b0)
    cat("\n    Minimum MSE Risk (Optimally Biased) Beta,   bStar=", bm)
    cat("\n    Most Likely UNRestricted Shrinkage Extent,  mUnr =", x$mUnr)
    cat("\n    Corresponding Expected -2*log(Likelihood Ratio)  = 0.0")  	
    cat("\n    Most Likely m-Value on Observed Lattice,    mClk =", x$mClk)
    cat("\n    Smallest Observed -2*log(Likelihood Ratio), minC =", x$minC)
    cat("\n    Smallest Observed EBAY criterion,           minE =", x$minE) 
    cat("\n    Smallest Observed RCOF criterion,           minR =", x$minR) 
    cat("\n    Estimated NonCentrality of F-test,          Phi2 =", x$Phi2) 	
    cat("\n    dMSE Estimate (Optimal Shrinkage Factor)    dMSE =", x$dMSE, "\n\n")	
}

"plot.YonX" <- function (x, trace = "all", ...)
{ 
    mcal <- x$sext[,3]    # MCAL values are in the 3rd column... 
    if (trace != "coef" && trace != "rmse" && trace != "spat" &&
        trace != "seq"  && trace != "YonX") trace <- "all"
    steps <- length(x$spat)-1		
    mO <- x$mUnr    # ML Optimal m-extent ...with k*==1 in YonX()...
    mM <- x$mRmin   # m-extent of Minimum Estimated Risk...
    mE <- x$mReql   # m-extent with Estimated Risk equal to OLS Relative Variance...	
    b0 <- x$prinstat[3]
    bm <- b0 * x$dMSE
    be <- b0 * (1 - x$mReql)	
    opar <- par(no.readonly = TRUE)  
    on.exit(par(opar))
    if (trace == "all") par(mfrow=c(2,2)) else par(mfrow=c(1,1))  
    if (trace == "all" || trace == "seq" || trace == "coef") {
        plot(mcal, x$coef, ann = FALSE, type = "n")
        lines(mcal, x$coef, col = 1, lty = 1, lwd = 2)
        abline(v = 0,  col = "blue", lty = 2, lwd = 1)
        points(0, b0, col = "blue", lwd = 2)  		
        abline(v = mO, col = "purple", lty = 2, lwd = 1)
        points(mO, bm, col = "purple", lwd = 2)
        abline(v = mM, col = "limegreen", lty = 2, lwd = 1)		  		 		 		
        abline(h = 0, col = gray(0.9), lwd = 2)  
        title(main = paste("COEFFICIENT TRACE"),
            xlab = "m = Multicollinearity Allowance", ylab = "Fitted Coefficient")   
    }
    if (trace == "seq") {
        cat("\nPress the Enter key to view the RMSE trace...")  
        scan()  
    }  
    if (trace == "all" || trace == "seq" || trace == "rmse") {
        rmax <- max(x$rmse)
        plot(mcal, x$rmse, ann = FALSE, type = "n", ylim=c(0,rmax))
        abline(v = 0,  col = "blue", lty = 2, lwd = 1)
        abline(v = mO, col = "purple", lty = 2, lwd = 1)
        abline(v = mM, col = "limegreen", lty = 2, lwd = 1)		  		
        abline(v = mE, col = "red", lty = 2, lwd = 1)		
        lines(mcal, x$rmse, col = 1, lty = 1, lwd = 2)	 
        points(0, x$rmse[1], col = "blue", lwd = 2)  			
        points(mM, x$rmse[mapp(mM,x)], col = "green", lwd = 2)
        points(mE, x$rmse[mapp(mE,x)], col = "red", lwd = 2)
        abline(h = 0, col = gray(0.9), lwd = 2) 		
        title(main = paste("RELATIVE MSE"),
            xlab = "m = Multicollinearity Allowance", ylab = "Scaled MSE Risk")   
    }
    if (trace == "seq") {
        cat("\nPress the Enter key to view the SPAT trace...")  
        scan()  
    }	
    if (trace == "all" || trace == "seq" || trace == "spat") {
        plot(mcal, x$spat, ann = FALSE, type = "n")  
        abline(v = 0,  col = "blue", lty = 2, lwd = 1)  			
        abline(v = mO, col = "purple", lty = 2, lwd = 1)
        abline(v = mM, col = "limegreen", lty = 2, lwd = 1)		  				
        lines(mcal, x$spat, col = 1, lty = 1, lwd = 2)
        points(0, x$spat[1], col = "blue", lwd = 2)  		
        points(mO, x$spat[mapp(mO,x)], col = "purple", lwd = 2) 	
        abline(h = 0, col = gray(0.9), lwd = 2)  
        title(main = paste("SHRINKAGE PATTERN"), 
            xlab = "m = Multicollinearity Allowance", ylab = "Delta Factor") 
    }
    if (trace == "all" || trace == "seq") {
        cat("\nPress the Enter key to view the YonX Scatter plot...")  
        scan()  
    }
    if (trace == "YonX" || trace == "all" || trace == "seq") {
        par(mfrow=c(1,1))
        yvec <- x$yvec
        xvec <- x$xvec
        v1 <- matrix(1, x$n, 1)
        ym <- mean(yvec)
        xm <- mean(xvec)
        b0fit <- (ym - b0*xm)*v1 + b0*xvec
        bmfit <- (ym - bm*xm)*v1 + bm*xvec
        befit <- (ym - be*xm)*v1 + be*xvec		
        plot(xvec, yvec, xlab = paste("Xvar =", x$yxnam[2]),
            ylab = paste("Yvar =", x$yxnam[1]), main = "YonX Shrinkage Plot")
        lines(xvec, b0fit, lwd=2, col="blue")
        lines(xvec, bmfit, lwd=2, col="purple")
        lines(xvec, befit, lwd=2, col="red")
        abline(h = ym, v = xm, lty = 2)
    }
}
