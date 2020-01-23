"unr.ridge" <-     # for RXswhrink version 1.3.1.
function (form, data, rscale = 1, steps = 8, omdmin = 9.9e-13) 
{
    if (missing(form) || class(form) != "formula") 
        stop("First argument to unr.ridge must be a valid linear regression formula.")
    yvar <- deparse(form[[2]])
    if (missing(data) || !inherits(data, "data.frame")) 
        stop("Second argument to unr.ridge must be an existing Data Frame.")
    dfname <- deparse(substitute(data))
    if (!is.element(yvar, dimnames(data)[[2]])) 
        stop("Response variable in the unr.ridge formula must be an existing variable.")
    lmobj <- lm(form, data)
    yvec <- as.matrix(lmobj$model[, 1])
    xmat <- as.matrix(lmobj$model[, 2:length(lmobj$model)])
    if (rscale < 1) 
        rscale <- 0
    if (rscale > 1) 
        rscale <- 2
    p <- ncol(xmat)
    if (p < 2) 
        stop("Number on non-constant regressor variables must be at least 2.")
    n <- nrow(xmat)
    if (n != nrow(yvec)) 
        stop("Numbers of observations in XMAT and YVEC must match.")
    if (n < p + 4) 
        stop("Number of observations must exceed number of regressors by at least 4.")
    mx <- matrix(apply(xmat, 2, "mean"), nrow = 1)
    crx <- xmat - matrix(1, n, 1) %*% mx
    xscale <- diag(p)
    if (rscale >= 1) {
        xscale <- diag(sqrt(diag(var(crx))))
        crx <- crx %*% solve(xscale)
    }	
    sx <- svd(crx)
    eigval <- matrix(sx$d^2, ncol = 1)         # p x 1
    eiginv <- solve(diag(sx$d^2, ncol = p))    # p x p
    cry <- matrix(yvec - mean(yvec), ncol = 1) # n x 1
    yscale <- 1
    if (rscale >= 1) {
        yscale <- sqrt(var(cry))
        cry <- cry/yscale[1, 1]
    }
    smse <- sx$v %*% eiginv %*% t(sx$v)
    risk <- diag(smse)
    tsmse <- sum(risk)
    comp <- solve(diag(sx$d, ncol = p)) %*% t(sx$u) %*% cry
    bstar <- sx$v %*% comp
    exev <- matrix(0, p, 1)
    infd <- matrix(0, p, 1)
    delta <- matrix(1, p, 1)
    idty <- diag(p)
    d <- idty
    cold <- delta
    sv <- matrix(sx$d, p, 1)
    ssy <- t(cry) %*% cry
    rho <- (sv * comp)/sqrt(ssy[1, 1])
    arho <- matrix(abs(rho), nrow = 1)
    r2 <- sum(arho^2)     # OLS "R^2" statistic....
    if (r2 >= 1) 
        stop(" Maximum Likelihood Shrinkage is not applicable when RSQUARE=1.")
    res <- cry - crx %*% bstar
    s2 <- t(res) %*% res/(n - p - 1)
    varrho <- s2[1, 1]/ssy[1, 1]
    tstat <- rho/sqrt(varrho)
    frat <- rho^2/varrho
    stat <- cbind(eigval, sv, comp, rho, tstat)
    dimnames(stat) <- list(1:p, c("LAMBDA", "SV", "COMP", "RHO", "TRAT"))
    RXolist <- list(data = dfname, form = form, p = p, n = n, 
        r2 = r2, s2 = s2, prinstat = stat, gmat = sx$v)
    OmR2dN <- (1 - r2 )/n
    dMSE <- rep(1,p)      # meaningless initial values...
    for( i in 1:p ) {
        dMSE[i] <- (rho[i])^2 / ( (rho[i])^2 + OmR2dN )    # Maximum-Likelihood estimate...
    }
    mcal <- 0
    kinc <- 1/min(dMSE)  # innitially > 1 but decreasing...
    konst <- kinc
    const <- (n - p - 3)/(n - p - 1)
    srat <- solve(diag(as.vector(sv), ncol = p)) %*% tstat
    MCAL <- 0
    KONST <- kinc
    C <- Inf
    E <- Inf
    R <- Inf
    maxinc <- p * steps   # i.e. total number of steps with "m" strictly > ZERO.
    for (inc in 1:maxinc) {
        mobj <- inc/steps
        iter <- m.ukd(mobj, p, dMSE)  # ...new function
        kinc <- iter$kinc
        d <- matrix(iter$d, p, p)
        dinc <- diag(d)
        omd <- pmax(1 - dinc, omdmin)
        ddomd <- dinc/omd
        rxi <- sum(arho * sqrt(ddomd))
        slik <- 2/(rxi + sqrt(4 * n + rxi^2))
        clik <- 2 * n * log(slik) + sum(ddomd) - (rxi/slik) - 
            n * log((1 - r2)/n)
        ebay <- sum(frat * omd - log(omd))
        sr2d <- sum(dinc * rho^2)
        if( sr2d >= 1 ) rcof <- Inf else {
            rcof <- -sum(log(omd)) + n * log((1 - sr2d)/(1 - r2)) }
        minc <- p - sum(dinc)
        C <- rbind(C, clik)
        E <- rbind(E, ebay)
        R <- rbind(R, rcof)
        KONST <- rbind(KONST, kinc)
        MCAL <- rbind(MCAL, minc)
        binc <- sx$v %*% d %*% comp
        vecr <- (idty - d) %*% srat
        compr <- const * vecr %*% t(vecr) + (2 * d - idty) * 
            eiginv
        diagc <- diag(diag(compr), ncol = p)
        lowr <- eiginv * d^2
        maxd <- matrix(pmax(diagc, lowr), p, p)
        compr <- compr - diagc + maxd
        smse <- sx$v %*% compr %*% t(sx$v)
        rinc <- diag(smse)
        lowb <- sx$v %*% lowr %*% t(sx$v)
        rinc <- pmax(rinc, diag(lowb))
        tinc <- sum(rinc)
        emse <- eiginv - compr
        sfac <- min(abs(diag(emse)))/100
        if (sfac < 1e-05) 
            sfac <- 1e-05
        eign <- eigen(emse/sfac)
        einc <- rev(eign$values) * sfac
        cinc <- matrix(0, p, 1)
        if (is.na(einc[1])) 
            einc[1] <- 0
        if (einc[1] < 0) {
            eign$vectors <- sx$v %*% eign$vectors
            cinc <- eign$vectors[, p]
            if (rscale == 2) {
                cinc <- cinc %*% xscale
                cinc <- cinc/sqrt(sum(cinc^2))
            }
            if (t(cold) %*% cinc < 0) 
                cinc <- -1 * cinc
            cold <- cinc
        }
        bstar <- cbind(bstar, binc)
        risk <- cbind(risk, rinc)
        exev <- cbind(exev, einc)
        infd <- cbind(infd, cinc)
        delta <- cbind(delta, dinc)
        tsmse <- rbind(tsmse, tinc)
        konst <- rbind(konst, kinc)
        mcal <- rbind(mcal, minc)
    }
    if (rscale == 2) {
        bstar <- yscale * solve(xscale) %*% bstar
        risk <- yscale^2 * solve(xscale^2) %*% risk
    }
    mlik <- cbind(MCAL, KONST, C, E, R)
    dimnames(mlik) <- list(0:maxinc, c("M", "K", "CLIK", "EBAY", "RCOF"))
    bstar <- t(bstar)
    risk <- t(risk)
    exev <- t(exev)
    infd <- t(infd)
    delta <- t(delta)
    sext <- cbind(tsmse, konst, mcal)
    dimnames(sext) <- list(0:maxinc, c("TSMSE", "KONST", "MCAL"))    mUnr <- mofk(p, k=1.0, dMSE)    minC <- min(mlik[,3])
    for( i in 1:maxinc ) {
        if( mlik[i,3] <= minC ) {       
            iUnr <- i - 1            break        }    }
    mClk <- iUnr/steps
    RXolist <- c(RXolist, list(coef = bstar, rmse = risk, 
        exev = exev, infd = infd, spat = delta, mlik = mlik, 
        sext = sext, mUnr = mUnr, mClk = mClk, minC = minC,
        dMSE = dMSE))
    class(RXolist) <- "unr.ridge"
    RXolist
}

"plot.unr.ridge" <-
function (x, trace = "all", trkey = FALSE, ...)
{ 
    mcal <- x$sext[, 3]    # MCAL values are in the 3rd column...  
    mcalp <- rep(mcal, times = x$p)  
    if (trace != "coef" && trace != "rmse" && trace != "exev" && 
        trace != "infd" && trace != "spat" && trace != "seq") 
        trace <- "all" 
    mV <- x$mClk    # m-extent with min Classical -2*log(Like)  		
    opar <- par(no.readonly = TRUE)  
    on.exit(par(opar))  
    if (trace == "all") par(mfrow=c(3,2)) else
        par(mfrow=c(1,1))  
    if (trace == "all" || trace == "seq" || trace == "coef") {
        plot(mcalp, x$coef, ann = FALSE, type = "n") 
        abline(v = mV, col = "gray", lty = 2)  		
        abline(h = 0, col = gray(0.9))  
        for (i in 1:x$p) lines(mcal, x$coef[, i], col = i, lty = i, 
            lwd = 2)  
        title(main = paste("COEFFICIENT TRACE"),
            xlab = "m = Multicollinearity Allowance", ylab = "Fitted Coefficients")  
        if( trkey ) 
            legend("bottom",all.vars(x$form)[2:(x$p+1)], col=1:(x$p), lty=1:(x$p), lwd=2)  
    }  
    if (trace == "seq") {
        cat("\nPress the Enter key to view the RMSE trace...")  
        scan()  
    }  
    if (trace == "all" || trace == "seq" || trace == "rmse") {
        plot(mcalp, x$rmse, ann = FALSE, type = "n") 
        abline(v = mV, col = "gray", lty = 2)  		
        abline(h = 0, col = gray(0.9))  
        for (i in 1:x$p) lines(mcal, x$rmse[, i], col = i, lty = i, 
            lwd = 2)  
        title(main = paste("RELATIVE MEAN SQ. ERROR"),
            xlab = "m = Multicollinearity Allowance", ylab = "Scaled MSE Risk")  
        if( trkey )
            legend("bottom",all.vars(x$form)[2:(x$p+1)], col=1:(x$p), lty=1:(x$p), lwd=2)  
    }  
    if (trace == "seq") {
        cat("\nPress the Enter key to view the EXEV trace...")  
        scan()  
    }  
    if (trace == "all" || trace == "seq" || trace == "exev") {
        plot(mcalp, x$exev, ann = FALSE, type = "n")  
        abline(v = mV, col = "gray", lty = 2)  
        abline(h = 0, col = gray(0.9)) 
        for (i in 1:x$p) lines(mcal, x$exev[, i], col = i, lty = i, 
            lwd = 2)  
        title(main = paste("EXCESS EIGENVALUES"),
            xlab = "m = Multicollinearity Allowance", ylab = "Least Squares minus UNRes")  
        if( trkey )
            legend("bottom",paste("Component",1:(x$p)), col=1:(x$p), lty=1:(x$p), lwd=2)  
    }  
    if (trace == "seq") {
        cat("\nPress the Enter key to view the INFD trace...")  
        scan()  
    }  
    if (trace == "all" || trace == "seq" || trace == "infd") {
        plot(mcalp, x$infd, ann = FALSE, type = "n")  
        abline(v = mV, col = "gray", lty = 2)  		
        abline(h = 0, col = gray(0.9))  
        for (i in 1:x$p) lines(mcal, x$infd[, i], col = i, lty = i, 
            lwd = 2)  
        title(main = paste("INFERIOR DIRECTION"), 
            xlab = "m = Multicollinearity Allowance", ylab = "Direction Cosines")  
        if( trkey )
            legend("bottom",all.vars(x$form)[2:(x$p+1)], col=1:(x$p), lty=1:(x$p), lwd=2)  
    }  
    if (trace == "seq") {
        cat("\nPress the Enter key to view the SPAT trace...")  
        scan()  
    }  
    if (trace == "all" || trace == "seq" || trace == "spat") {
        plot(mcalp, x$spat, ann = FALSE, type = "n")  
        abline(v = mV, col = "gray", lty = 2)  		
        abline(h = 0, col = gray(0.9))  
        for (i in 1:x$p) lines(mcal, x$spat[, i], col = i, lty = i, 
            lwd = 2)  
        title(main = paste("SHRINKAGE PATTERN"), 
            xlab = "m = Multicollinearity Allowance", ylab = "UNRes Delta Factors")  
        if( trkey )
            legend("bottom",paste("Component",1:(x$p)), col=1:(x$p), lty=1:(x$p), lwd=2)  
    }  
}  

"print.unr.ridge" <-
function (x, ...) 
{
    cat("\nunr.ridge Object: Shrinkage via a PATH of UNRestricted Shape\n")
    cat("Data Frame:", x$data, "\n")
    cat("Regression Equation:\n")
    print(x$form)
    cat("\n    Number of Regressor Variables, p =", x$p, "\n")
    cat("    Number of Observations, n =", x$n, "\n")
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
    cat("\n    Most Likely UNRestricted Shrinkage Extent, mUnr =", x$mUnr)
    cat("\n    Corresponding -2*log(LikelihoodRatio) statistic = 0.0")  	
    cat("\n    Most Likely m-Value on the Lattice, mClk =       ", x$mClk)
    cat("\n    Smallest Observed -2*log(LikelihoodRatio), minC =", x$minC)    	
    cat("\n    dMSE Estimates =", x$dMSE, "\n\n")	
} 
"m.ukd" <-function (muobj, p, dMSE) {    kM <- matrix(1,1,p)/dMSE[order(dMSE)]    if (muobj <= 0) {        d <- diag(p)        kinc <- kM[1]        return(list(kinc = kinc, d = d))    }    if (muobj >= p) {        d <- matrix(0, p, p)        kinc <- 0        return(list(kinc = kinc, d = d))    }    mVar <- matrix(1,1,p)    # initial values...    for( j in 1:p ) {        mVar[j] <- mofk(p, k=kM[j], dMSE)    }    kMin <- (p-muobj)/sum(dMSE)    if( muobj > mVar[p] ) {  # large m (truly small k) cases...        d <- kMin*diag(dMSE)        kinc <- kMin        return(list(kinc = kinc, d = d))    } else {        for( j in 2:p ) {             if( mVar[j-1] < muobj && muobj <= mVar[j] ) {	            if( muobj == mVar[j] ) {                    kinc <- kM[j]                } else {                    B <- muobj - mVar[j-1]                    D <- mVar[j] - muobj                     kinc <- (B*D/(B+D))*( kM[j-1]/B + kM[j]/D )                }            }        }        for( j in 1:p ) {            dd <- min(1,kinc*dMSE[j])            if( j == 1 ) dp <- dd else dp <- c(dp,dd)        }        d <- diag(dp)        return(list(kinc = kinc, d = d))		    }}"mofk" <- function(p, k, dMSE)  # Many-to-One Function for large k-values...{  p <- as.integer(p)  if( k < 0 ) k <- 0  kM <- 1/min(dMSE)  if( k > kM ) k <- kM  m <- as.double(0)  for( j in 1:p ) {    m <- m + 1 - min(1,k*dMSE[j])  }  m}
