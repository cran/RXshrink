"MLtrue" <-  
function (form, data, seed, go=TRUE, truv, trub, truc) 
{
    if (missing(form) || class(form) != "formula") 
        stop("First argument to MLtrue must be a valid linear model formula.")
    yvar <- deparse(form[[2]])
    if (missing(data) || !inherits(data, "data.frame")) 
        stop("Second argument to MLtrue must be an existing Data Frame.")
    dfname <- deparse(substitute(data))
    if (!is.element(yvar, dimnames(data)[[2]])) 
        stop("Response y-variable in the MLtrue() <formula> must be within the input data.frame.")
    lmobj <- lm(form, data)
    yvec <- as.matrix(lmobj$model[, 1])
    mnam <- dimnames(lmobj$model)[[2]]
    ynam <- mnam[1]
    xmat <- as.matrix(lmobj$model[, 2:length(lmobj$model)])
    xnam <- mnam[2:length(mnam)]
    p <- ncol(xmat)
    if (p < 2) {
        cat("\nNumber on non-constant regressor variables must be at least 2.\n")
        return("p_less_than_2")        
    }
    n <- nrow(xmat)
    if (n != nrow(yvec)) 
        stop("Numbers of observations in XMAT and YVEC must match.")
    if (n < p + 4) {
        cat("\nNumber of observations must exceed number of regressors by at least 4.\n")
        return("n_less_than_p+4")   
    }
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
    sv <- matrix(sx$d, p, 1)
    ssy <- t(cry) %*% cry
    rho <- (sv * comp)/sqrt(ssy[1, 1])
    r2 <- sum(rho^2)     # OLS "R^2" statistic....
    if (r2 >= 1) {
        cat("\nMaximum Likelihood Shrinkage is not applicable when RSQUARE=1.\n")
        return("Rsq_is_1")
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
    RXolist <- c(RXolist, list(beta = betaols, comp = comp, rmse = risk, ys = yscale, xs = diag(xscale)))  
    class(RXolist) <- "MLtrue"	
    if (!go) {
        cat("\nPrint the MLtrue output.list to View (default) parameter values...\n\n")
        return(RXolist)
    }
    else {
        tvar  <- s2          # Set defaults...
        tbeta <- betaols
        tcomp <- comp
        useb <- FALSE
        if (!missing(truv) && truv >= 0) tvar <- truv
        if (!missing(trub)) {
            useb <- TRUE
            if(length(trub) == p) tbeta <- as.matrix(trub, p, 1)
        }
        if (!useb && !missing(truc)) {      # ignore truc when useb==TRUE
            if(length(truc) == p) tcomp <- as.matrix(truc, p, 1)
        }			
    }
    tsig <- sqrt(tvar)
    if (useb) {
        Yhat <- crx %*% tbeta
    }
    else {
        Yhat <- crx %*% sx$v %*% tcomp   
    }
    if (missing(seed)) 
        seed <- sample(1001, 1) - 1
    set.seed(seed)
    Yvec <- Yhat + as.matrix(rnorm(n, sd=tsig), n, 1)
    new <- data.frame(cbind(Yvec, Yhat, crx, cry))
    names(new) <- c("Yvec", "Yhat", xnam, ynam)
    RXolist <- c(list(new = new, Yvec = Yvec, Yhat = Yhat, seed = seed, tvar = tvar,
        tbeta = tbeta, tcomp = tcomp, useb = useb), RXolist)
    class(RXolist) <- "MLtrue"
    RXolist		
}

"print.MLtrue" <-
function (x, ...) 
{
    cat("\nMLtrue Object: OLS Beta-Coefficient Estimates and Related Statistics...\n")
    cat("Data Frame:", x$data, "\n")
    cat("Regression Equation:\n")
    print(x$form)
    cat("\n    Number of Regressor Variables, p =", x$p, "\n")
    cat("    Number of Observations, n =", x$n, "\n")
    cat("\nPrincipal Axis Summary Statistics of Ill-Conditioning...\n")
    print.default(x$prinstat, quote = FALSE)
    cat("\nOLS Residual Mean Square for Error = Estimate of truv =\n", x$s2, "\n")
    cat("\nOLS Beta Coefficients   =\n", x$beta, "\n")
    cat("\nUncorrelated Components =  (COMP-column Above)\n", x$comp, "\n")
    if (length(x$new) > 0) {
        cat("\nRandom Number SEED value  =", x$seed, "\n")
        cat("\nTrue error Variance, tvar =", x$tvar, "\n")
        if (x$useb) {
            cat("\nTrue OLS Beta Coefficients, tbeta =\n", x$tbeta, "\n")
            }
        else cat("\nTrue Uncorrelated Components, tcomp =\n", x$tcomp, "\n")
    }
    cat("\n")
} 