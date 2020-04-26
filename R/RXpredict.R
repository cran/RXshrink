"RXpredict" <- function (x, data, m="ML", rscale=1)
{
    cond <- class(x) %in% c("qm.ridge", "unr.ridge", "aug.lars", "uc.lars", "correct.signs" )
    if (missing(x) || cond == FALSE) {
        cat("\nFirst argument to RXpredict() should be an output object from one of Five\n")
        stop("RXshrink functions: qm.ridge, unr.ridge, aug.lars, uc.lars or correct.signs")
    }		
    mMax <- as.double(x$p)
    if (class(x) == "correct.signs") {tsteps <- 1} else {tsteps <- nrow(x$coef)-1}
    if (missing(m) || is.logical(m) || is.character(m)) {
        if (class(x) == "correct.signs") {m <- 0} else {m <- x$mClk}
    }
    if (m < 0) m <- 0
    if (m > mMax) m <- mMax
    if (class(x) == "correct.signs") {bstar <- x$signs[,5]} else {bstar <- x$coef}
    stsize <- mMax/tsteps              # Step-Size == 1/steps in unr.ridge()
    midx <- round(m/stsize, 0)         # m-Index: "row" number == midx + 1
    mobs <- stsize*midx                # Nearest "observed" value of of m...
    lmobj <- lm(x$form, data)
    yvec <- as.matrix(lmobj$model[, 1])
    n <- nrow(yvec)
    xmat <- as.matrix(lmobj$model[,-1])
    if (class(x) == "correct.signs") {coef <- bstar} else {coef <- bstar[midx+1,]}
    # bstar coefficients for specified m-extent
    mx <- matrix(apply(xmat, 2, "mean"), nrow = 1)
    cx <- xmat - matrix(1, n, 1) %*% mx          # Centering is always applied (implicit intercept)...
    if (rscale >= 1) {
        xscale <- diag(sqrt(diag(var(cx))))      # Rescaling to equal variances & std. errors...
        crx <- cx %*% solve(xscale)              # Centered and Rescaled X-matrix...
    }
    cry <- matrix(yvec - mean(yvec), ncol = 1)   # n x 1 ; more Centering...
    yscale <- 1
    if (rscale >= 1) {
        yscale <- sqrt(var(cry))
        cry <- cry/yscale[1, 1]
    }
    cryprd <- crx %*% coef
    if (class(x) == "correct.signs" && 0 < m) {
        cryprd <- cryprd * (mMax - m) / mMax
        mobs <- m
    }	
    yvecprd <- matrix( mean(yvec), n, 1) + cryprd %*% sqrt(var(yvec))
    if (class(x) == "aug.lars" || class(x) == "uc.lars") mobs <- x$mlik[mobs+1,1]	
    RXolist <- list(class=class(x), cryprd = as.vector(cryprd), cry = as.vector(cry),
                    yvecprd = as.vector(yvecprd), yvec = as.vector(yvec), m = m,
                    mobs = mobs)
    class(RXolist) <- "RXpredict"
    RXolist
}

"print.RXpredict" <-
function (x, ...) 
{
    cat("\nRXpredict Object: In-Sample Predictions (fitted.values) for RXshrink functions...\n")
    cat("\nClass of Estimation Method:", x$class, "\n")
    cat("\nCentered and Rescaled y-Outcome Vector: cry =\n")
    print.default(x$cry, quote = FALSE)
    cat("\nFitted.Values = Predictions of cry: cryprd =\n")
    print.default(x$cryprd, quote = FALSE)
    cat("\nObserved y-Outcome Vector: yvec =\n")
    print.default(x$yvec, quote = FALSE)
    cat("\nPredictions of yvec: yvecprd =\n")
    print.default(x$yvecprd, quote = FALSE)
    cat("\nShrinkage m-Extent requested: m =", x$m, "\n")
    cat("\nObserved m-Extent most close to the requested m is: mobs =", x$mobs, "\n")
}

"plot.RXpredict" <-
function (x, fit = "yvecprd", ...) 
{
    if (fit != "yvecprd" && fit != "cryprd") fit <- "both"
    p <- length(x$cryprd)
    obs <- 1:p
    m <- x$mobs    # m-extent actually found and saved... 	
    op <- par(no.readonly = TRUE)
    if (fit == "both") {par(mfrow=c(2,1))} else {par(mfrow=c(1,1))}
	
    if (fit == "both" || fit == "yvecprd") {
        plot(obs, x$yvec, ann = FALSE, type = "n")
        points(obs, x$yvec)
        lines(obs, x$yvecprd, lwd=2, col="blue")
        title(main = paste(x$class, "Predictions for m-Extent =", round(m,2)), 
            xlab = "Observation Numbers", ylab = "Observed and Predicted y-Outcomes")
    }
    
    if (fit == "both" || fit == "cryprd") {
        plot(obs, x$cry, ann = FALSE, type = "n")
        abline(h = 0, col = gray(0.9), lty=2, lwd = 2)
        points(obs, x$cry)
        lines(obs, x$cryprd, lwd=2, col="blue")
        title(main = paste(x$class, "Fitted.Values for m-Extent =", round(m,2)), 
            xlab = "Observation Numbers", ylab = "Centered and Rescaled y-Outcomes")
    }
    par(op)
}
