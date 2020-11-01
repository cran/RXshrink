"MLcalc" <-  
function (form, data, rscale = 1, delmax = 0.999999) 
{
    if (missing(form) || class(form) != "formula") 
        stop("First argument to MLcalc must be a valid linear regression formula.")
    yvar <- deparse(form[[2]])
    if (missing(data) || !inherits(data, "data.frame")) 
        stop("Second argument to MLcalc must be an existing Data Frame.")
    dfname <- deparse(substitute(data))
    if (!is.element(yvar, dimnames(data)[[2]])) 
        stop("Response variable in the MLcalc formula must be an existing variable.")
    if (rscale != 1) rscale <- 2  # just two rscale options from unr.ridge()
    lmobj <- lm(form, data)
    yvec <- as.matrix(lmobj$model[, 1])
    xmat <- as.matrix(lmobj$model[, 2:length(lmobj$model)])
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
    xscale <- diag(sqrt(diag(var(crx))))
    crx <- crx %*% solve(xscale)
    sx <- svd(crx) # sx object is Singular Value Decomposition of Centered X-matrix
    eigval <- matrix(sx$d^2, ncol = 1) # p x 1
    eiginv <- rep(0, p)                # p x 1
    svainv <- rep(0, p)                # p x 1
    for (i in 1:p) {
      if (eigval[i] > .Machine$double.eps) {
        eiginv[i] <- 1/eigval[i]
        svainv[i] <- sqrt(eiginv[i])
        }
    }
    eiginv <- diag(eiginv, ncol = p)   # p x p
    svainv <- diag(svainv, ncol = p)   # p x p
    cry <- matrix(yvec - mean(yvec), ncol = 1) # n x 1
    yscale <- sqrt(var(cry))
    cry <- cry/yscale[1, 1]
    smse <- sx$v %*% eiginv %*% t(sx$v)
    risk <- diag(smse)
    tsmse <- sum(risk)
    comp <- svainv %*% t(sx$u) %*% cry
    betaols <- sx$v %*% comp
    exev <- matrix(0, p, 1)
    infd <- matrix(0, p, 1)
    delta <- matrix(1, p, 1)
    idty <- diag(p)
    d <- idty
    sv <- matrix(sx$d, p, 1)
    ssy <- t(cry) %*% cry
    rho <- (sv * comp)/sqrt(ssy[1, 1])
    arho <- matrix(abs(rho), nrow = 1)
    r2 <- sum(arho^2)     # OLS "R^2" statistic....
    if (r2 >= 1) {
        cat("\nMaximum Likelihood Shrinkage is not applicable when RSQUARE=1.\n")
        return("Rsq1")
        }
    res <- cry - crx %*% betaols
    s2 <- t(res) %*% res/(n - p - 1)
    varrho <- s2[1, 1]/ssy[1, 1]
    tstat <- rho/sqrt(varrho)
    frat <- rho^2/varrho
    stat <- cbind(eigval, sv, comp, rho, tstat)
    dimnames(stat) <- list(1:p, c("LAMBDA", "SV", "COMP", "RHO", "TRAT"))
    RXolist <- list(data = dfname, form = form, p = p, n = n, 
        r2 = r2, s2 = s2, prinstat = stat, gmat = sx$v)
    OmR2dN <- (1 - r2 )/n
    dMSE <- rep(1,p)         # meaningless initial values...
    for( i in 1:p ) {
        dMSE[i] <- (rho[i])^2 / ( (rho[i])^2 + OmR2dN )    # Maximum-Likelihood shrinkage factors...
    }
    kinc <- 1                          # kinc*dMSE values then have minimum MSE risk...
    dFact <- (n - p - 3)/(n - p - 1)
    srat <- solve(diag(as.vector(sv), ncol = p)) %*% tstat
    mobj <- mofk(p, kinc, dMSE)
    dinc <- dMSE	# The p Optimal Shrinkage-Factors...
    diag( d ) <- dMSE
    omd <- pmax(1 - dinc, 1 - delmax)  # Strictly Positive values...
    ddomd <- dinc/omd                  # Diagonal Elements...
    rxi <- sum(arho * sqrt(ddomd))
    slik <- 2/(rxi + sqrt(4 * n + rxi^2))
    clik <- 2 * n * log(slik) + sum(ddomd) - (rxi/slik) - 
        n * log((1 - r2)/n)
    if (clik < 1 - delmax) clik <- 1 - delmax
    minc <- p - sum(dinc)
    bstar <- sx$v %*% d %*% comp       # MLbeta
    vecr <- (idty - d) %*% srat
    compr <- dFact * vecr %*% t(vecr) + (2 * d - idty) * eiginv
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
    if (sfac < 1e-05) sfac <- 1e-05
    eign <- eigen(emse/sfac)
    einc <- sort(eign$values) * sfac   # Increasing order; negative values first...
    cinc <- matrix(0, p, 1)
    if (is.na(einc[1])) einc[1] <- 0
    if (einc[1] + 1 - delmax < 0) {    # New Correction...	
        eign$vectors <- sx$v %*% eign$vectors
        cinc <- eign$vectors[, p]
        if (rscale == 2) {
            cinc <- cinc %*% xscale
            cinc <- cinc/sqrt(sum(cinc^2))
        }
    }
    bstar <- t(cbind(betaols, bstar))
    risk <- t(cbind(risk, rinc))
    if (rscale == 2) {
        bstar <- as.double(yscale) * (bstar %*% solve(xscale))
        risk <- as.double(yscale^2) * (risk %*% solve(xscale^2)) 
    }                  
    RXolist <- c(RXolist, list(beta = bstar, rmse = risk, dMSE = dMSE, ys = yscale, xs = diag(xscale)))  
    class(RXolist) <- "MLcalc"         
    RXolist                            
}

"print.MLcalc" <-
function (x, ...) 
{
    cat("\nMLcalc Object: Most Likely Coefficient Estimates under Normal-Theory\n")
    cat("Data Frame:", x$data, "\n")
    cat("Regression Equation:\n")
    print(x$form)
    cat("\n    Number of Regressor Variables, p =", x$p, "\n")
    cat("    Number of Observations, n =", x$n, "\n")
    cat("\nPrincipal Axis Summary Statistics of Ill-Conditioning...\n")
    print.default(x$prinstat, quote = FALSE)
    cat("\n    Residual Mean Square for Error =", x$s2, "\n")
    cat("    Estimate of Residual Std. Error =", sqrt(x$s2), "\n")
    cat("\nOLS Beta Coefficients =\n", x$beta[1,], "\n")
    cat("\nML Optimally Biased Coefficients =\n", x$beta[2,], "\n") 
    cat("\nOLS Relative MSE Risks =\n", x$rmse[1,], "\n")
    cat("\nML Minimum Relative Risks =\n", x$rmse[2,], "\n")  	
    cat("\ndMSE Estimates =\n", x$dMSE, "\n\n")	
} 