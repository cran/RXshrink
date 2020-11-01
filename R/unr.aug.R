"unr.aug" <- function (urobj) 
{
    if (missing(urobj) || class(urobj) != "unr.ridge") 
        stop("First argument to unr.aug() must be a valid unr.ridge() output object.")
    p <- urobj$p
    form <- urobj$form
    vnams <- all.vars(form)            # p+1 variable names in called order...
    data <- urobj$data                 # data.frame containing all p+1 "vnams"
    lmraw <- lm(form, data)            # preliminary calculations...
    yvec <- as.matrix(lmraw$model[, 1])
    xmat <- as.matrix(lmraw$model[, 2:length(lmraw$model)])
    rscale <- urobj$rscale             # rscale = 0, 1 or 2 setting from unr.ridge()
    n <- nrow(xmat)
    mx <- matrix(apply(xmat, 2, "mean"), nrow = 1)
    crx <- xmat - matrix(1, n, 1) %*% mx
    xscale <- diag(p)
    if (rscale >= 1) {
        xscale <- diag(sqrt(diag(var(crx))))
        crx <- crx %*% solve(xscale)
    }	
    cry <- matrix(yvec - mean(yvec), ncol = 1) # n x 1
    yscale <- 1
    if (rscale >= 1) {
        yscale <- sqrt(var(cry))
        cry <- cry/yscale[1, 1]
    }
    crdata <- data.frame(cbind(cry, crx))
    names(crdata) <- vnams
    LMobj <- lm(form, crdata) # centered & rescaled calculations; large object
    # Next, calculate p+1 estimates at the p-unr "knots" plus the optimal shrinkage extent...
    dMSE <- urobj$dMSE
    gmat <- urobj$gmat
    beta0 <- LMobj$coef[2:(p+1)]       # p OLS estimates..
    cvec <- t(gmat) %*% beta0          # uncorrelated components
    dINC <- dMSE[order(dMSE)]          # dMSE order statistics...
    kM <- matrix(1,1,p+2)              # dummy starting values...
    for (j in 1:p) { kM[1,j] <- 1.0/dINC[j] }    # knots at ordered reciprocals
    kM[1,p+1] <- 1.0                   # optimal shrinkage at k-star == 1
    kM[1,p+2] <- 0.0                   # complete shrinkage at k-star == 0
    mM <- matrix(1,1,p+2)              # dummy mcal starting values...
    for (j in 1:(p+2)) { mM[1,j] <- mofk(p, kM[1,j], dMSE) }
    bstar <- beta0
    mcal <- mM[1,1]
    for (inc in 2:(p+2)) {
        mobj <- mM[1,inc]
        iter <- kofm(mobj, p, dMSE)
        d <- iter$d          # Diagonal Matrix, p x p...
        binc <- gmat %*% d %*% cvec
        bstar <- cbind(bstar, binc)
        mcal <- rbind(mcal, mobj)
    }
    RXolist <- list(p = p, LMobj = LMobj, bstar = bstar, mcal = mcal, vnams = vnams)
    class(RXolist) <- "unr.aug"
    RXolist
}
	
"unr.biv" <- function (uraug, x1 = 1, x2 = 2, conf1 = 0.95, conf2 = 0.50) 
{
    if (missing(uraug) || class(uraug) != "unr.aug") {
        cat("\nNOTE: Call unr.aug() and save its output list before making")
        cat("\nAny calls to unr.biv().\n")	
        stop("First argument to unr.biv() must be a valid unr.aug() output object.")
    }
    p <- uraug$p
    confOK <- 0
    if (conf1 >= 0.05 && conf1 <= 0.95) confOK <- 1
    if (conf2 >= 0.05 && conf2 <= 0.95) confOK <- confOK + 1
    if (confOK == 1) cat("\nWARNING: Only 1 Confidence Ellipse will be computed.\n\n")
    if (confOK == 0) stop("Neither Confidence Ellipse can be computed.\n\n")		
    x1 <- as.integer(x1)
    if (x1 >= 1 && x1 <= p) { x1var <- uraug$vnams[1+x1] } else { 
        stop("Variable x1 in unr.biv() Call is Out-of-Range.")
    }
    x2 <- as.integer(x2)
    if (x2 >= 1 && x2 <= p && x2 != x1 ) { x2var <- uraug$vnams[1+x2] } else { 
        stop("Variable x2 in unr.biv() Call is Out-of-Range or Equal to x1.")
    }
    LMobj <- uraug$LMobj
    bstar <- uraug$bstar
    mcal <- uraug$mcal
    RXolist <- list(p = p, LMobj = LMobj, bstar = bstar, mcal = mcal, x1 = x1, x2 = x2)
    if (conf1 >= 0.05 && conf1 <= 0.95) {
        ellip1 <- ellipse(LMobj, which = c((x1 + 1), (x2 + 1)), conf1)
        corr1 <- cor(ellip1)
        RXolist <- c(RXolist, list(ellip1 = ellip1, conf1 = conf1, ecor1 = corr1[1,2]))
    }
    if (conf2 >= 0.05 && conf2 <= 0.95) {
        ellip2 <- ellipse(LMobj, which = c((x1 + 1), (x2 + 1)), conf2)
        corr2 <- cor(ellip2)
        RXolist <- c(RXolist, list(ellip2 = ellip2, conf2 = conf2, ecor2 = corr2[1,2]))
    }
    class(RXolist) <- "unr.biv"
    RXolist
}

"plot.unr.biv" <- function (x, type = "ellip", ...)
{
    if (type != "ellip") type <- "trace"
    if (type == "ellip") {
        OK <- 0 
        if (length("x$ellip1") > 0) {
            plot(x$ellip1, type = 'l')
            sub1 <- paste("corr = ", round(x$ecor1, 4))
            OK <- 1 
        }
        if (length("x$ellip2") > 0) {
            if (OK == 1) { lines(x$ellip2, type = 'l') }
            else {
                plot(x$ellip2, type = 'l')
            	sub1 <- paste("corr = ", round(x$ecor2, 4))
            }
        }
        if (OK == 0) stop("No Confidence Ellipse is Available.")
        abline(h=0, lty=2, col="lightgray" )
        abline(v=0, lty=2, col="lightgray" )
        lines(x$bstar[x$x1,], x$bstar[x$x2,], type = "l", lwd=2, col="red")
        points(x$bstar[x$x1,1], x$bstar[x$x2,1], type = "p", lwd=2, col="blue")
        pp1 <- 1 + x$p; pp2 <- 1 + pp1
        points(x$bstar[x$x1,pp1], x$bstar[x$x2,pp1], type = "p", lwd=2, col="purple")
        points(x$bstar[x$x1,pp2], x$bstar[x$x2,pp2], type = "p", lwd=2, col="black")
        title(sub = sub1, col.sub="darkgreen", cex.sub=1.0) 
        title(main = "Confidence Ellipses & BLUE_minMSE_ZERO Path")
    }
    else {
        mcal <- t(x$mcal)
        yup <- max(x$bstar)
        ylo <- min(x$bstar)	
        plot(rep(mcal,x$p), x$bstar, ann = FALSE, type = "n", ylim = c(ylo, yup))
        abline(h = 0, col = gray(0.9), lwd = 2)  
        for (i in 1:x$p) lines(mcal, x$bstar[i,], col = i, lty = i, lwd = 2)
        for (i in 1:x$p) points(mcal, x$bstar[i,], col = i, lwd = 2)  		
        title(main = paste("Piecewise-Linear Splines"),
            xlab = "m = Multicollinearity Allowance", ylab = "Shrunken Coefficients") 
    }	
}  

"print.unr.biv" <- function (x, ...) 
{
    cat("\nunr.biv Object: Bivariate displays of UNRestricted Shrinkage\n")
    cat("\n    Current Horizontal Coefficient Number =", x$x1)
    cat("\n    Current Vertical   Coefficient Number =", x$x2)
    cat("\n    Matrix of Fitted Coefficients and their mcal-Extents:\n\n")
    pp1 <- 1 + x$p; pp2 <- 1 + pp1 
    mat <- as.matrix(cbind(t(x$bstar), x$mcal))
    rownames(mat) <- c(1:pp2)
    colnames(mat)[pp1] <- c("mcal")
    print(mat)
}
