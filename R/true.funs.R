"true.risk" <-
function (form, data, trug, trus, Q = 0, rscale = 1, steps = 8, 
    qmax = 5, qmin = -5) 
{
    if (missing(form) || class(form) != "formula") 
        stop("First argument to true.risk must be a valid linear regression formula.")
    yvar <- deparse(form[[2]])
    if (missing(data) || !inherits(data, "data.frame")) 
        stop("Second argument to true.risk must be an existing Data Frame.")
    dfname <- deparse(substitute(data))
    lmobj <- lm(form, data)
    xmat <- as.matrix(lmobj$model[, 2:length(lmobj$model)])
    if (rscale < 1) 
        rscale <- 0
    if (rscale > 1) 
        rscale <- 2
    qp <- min(max(Q, qmin), qmax)
    n <- nrow(xmat)
    p <- ncol(xmat)
    if (p < 2) 
        stop("Number on non-constant regressor variables must be at least 2.")
    if (p != nrow(trug)) 
        stop("Number of true components must equal number of columns in XMAT.")
    mx <- matrix(apply(xmat, 2, "mean"), nrow = 1)
    crx <- xmat - matrix(1, n, 1) %*% mx
    xscale <- diag(p)
    if (rscale >= 1) {
        xscale <- diag(sqrt(diag(var(crx))))
        crx <- crx %*% solve(xscale)
    }
    sx <- svd(crx)
    eigval <- matrix(sx$d^2, ncol = 1)
    if (eigval[p] <= 0) {
        cat("\n\ntrue.risk requires XMAT to have Full Column Rank.")
        print(sx$d)
        print(sx$v)
        stop()
    }
    eiginv <- solve(diag(sx$d^2, ncol = p))
    smse <- sx$v %*% eiginv %*% t(sx$v)
    risk <- diag(smse)
    tsmse <- sum(risk)
    fact <- trus^2
    for (i in 1:p) fact <- fact + trug[i]^2 * eigval[i]/(n - 
        1)
    trug <- trug/sqrt(fact)
    trus <- trus/sqrt(fact)
    beta <- sx$v %*% trug
    maxinc <- p * steps
    qm1 <- qp - 1
    exev <- matrix(0, p, 1)
    infd <- matrix(0, p, 1)
    delta <- matrix(1, p, 1)
    d <- diag(delta)
    cold <- delta
    sv <- matrix(sx$d, p, 1)
    idty <- diag(p)
    if (qp == 1) 
        eqm1 <- matrix(1, p, 1)
    else eqm1 <- exp(qm1 * log(eigval))
    stat <- cbind(eigval, sv, trug, beta)
    dimnames(stat) <- list(1:p, c("LAMBDA", "SV", "GAMMA", "BETA"))
    RXolist <- list(data = dfname, form = form, trug = trug, 
        trus = trus, qp = qp, p = p, n = n, prinstat = stat)
    mcal <- 0
    konst <- 0
    kinc <- 0
    for (inc in 1:maxinc) {
        mobj <- inc/steps
        iter <- mstep(mobj, kinc, p, qp, eqm1)
        kinc <- iter$kinc
        d <- matrix(iter$d, p, p)
        binc <- sx$v %*% d %*% trug
        rbias <- (idty - d) %*% trug/trus
        compr <- rbias %*% t(rbias) + eiginv * d^2
        smse <- sx$v %*% compr %*% t(sx$v)
        rinc <- diag(smse)
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
        dinc <- diag(d)
        minc <- p - sum(dinc)
        beta <- cbind(beta, binc)
        risk <- cbind(risk, rinc)
        exev <- cbind(exev, einc)
        infd <- cbind(infd, cinc)
        delta <- cbind(delta, dinc)
        tsmse <- rbind(tsmse, tinc)
        konst <- rbind(konst, kinc)
        mcal <- rbind(mcal, minc)
    }
    beta <- t(beta)
    risk <- t(risk)
    exev <- t(exev)
    infd <- t(infd)
    delta <- t(delta)
    sext <- cbind(tsmse, konst, mcal)
    dimnames(sext) <- list(0:maxinc, c("TSMSE", "KONST", "MCAL"))
    RXolist <- c(RXolist, list(coef = beta, rmse = risk, exev = exev, 
        infd = infd, spat = delta, sext = sext))
    class(RXolist) <- "true.risk"
    RXolist
}

"true.simu" <-
function (form, data, trug, trus, Q = 0, rscale = 1, steps = 8, 
    seed = 123, qmax = 5, qmin = -5) 
{
    if (missing(form) || class(form) != "formula") 
        stop("First argument to true.simu must be a valid linear regression formula.")
    yvar <- deparse(form[[2]])
    if (missing(data) || !inherits(data, "data.frame")) 
        stop("Second argument to true.simu must be an existing Data Frame.")
    dfname <- deparse(substitute(data))
    lmobj <- lm(form, data)
    xmat <- as.matrix(lmobj$model[, 2:length(lmobj$model)])
    if (rscale < 1) 
        rscale <- 0
    if (rscale > 1) 
        rscale <- 2
    qp <- min(max(Q, qmin), qmax)
    n <- nrow(xmat)
    p <- ncol(xmat)
    if (p < 2) 
        stop("Number on non-constant regressor variables must be at least 2.")
    if (p != nrow(trug)) 
        stop("Number of true components must equal number of columns in XMAT.")
    mx <- matrix(apply(xmat, 2, "mean"), nrow = 1)
    crx <- xmat - matrix(1, n, 1) %*% mx
    xscale <- diag(p)
    if (rscale >= 1) {
        xscale <- diag(sqrt(diag(var(crx))))
        crx <- crx %*% solve(xscale)
    }
    sx <- svd(crx)
    eigval <- matrix(sx$d^2, ncol = 1)
    if (eigval[p] <= 0) {
        cat("\n\ntrue.simu requires XMAT to have Full Column Rank.")
        print(sx$d)
        print(sx$v)
        stop()
    }
    eiginv <- solve(diag(sx$d^2, ncol = p))
    fact <- trus^2
    for (i in 1:p) fact <- fact + trug[i]^2 * eigval[i]/(n - 
        1)
    trug <- trug/sqrt(fact)
    trus <- trus/sqrt(fact)
    beta <- sx$v %*% trug
    expy <- crx %*% beta
    if (missing(seed)) 
        seed <- sample(1001, 1) - 1
    set.seed(seed)
    error <- matrix(rnorm(n, 0, trus), ncol = 1)
    yvec <- matrix(expy + error - mean(error), ncol = 1)
    ydat <- cbind(yvec, expy)
    dimnames(ydat) <- list(1:n, c("YVEC", "EXPY"))
    ssy <- t(yvec) %*% yvec
    comp <- solve(diag(sx$d)) %*% t(sx$u) %*% yvec
    bstar <- sx$v %*% comp
    loss <- (bstar - beta)^2/trus^2
    totlos <- sum(loss)
    maxinc <- p * steps
    qm1 <- qp - 1
    delta <- matrix(1, p, 1)
    d <- diag(p)
    sv <- matrix(sx$d, p, 1)
    idty <- diag(p)
    if (qp == 1) 
        eqm1 <- matrix(1, p, 1)
    else eqm1 <- exp(qm1 * log(eigval))
    stat <- cbind(eigval, sv, trug, beta)
    dimnames(stat) <- list(1:p, c("LAMBDA", "SV", "GAMMA", "BETA"))
    RXolist <- list(data = dfname, form = form, trug = trug, 
        trus = trus, qp = qp, p = p, n = n, ydat = ydat, prinstat = stat)
    rho <- (sv * comp)/sqrt(ssy[1, 1])
    res <- yvec - crx %*% bstar
    s2 <- t(res) %*% res/(n - p - 1)
    varrho <- s2[1, 1]/ssy[1, 1]
    trat <- rho/sqrt(varrho)
    mcal <- 0
    konst <- 0
    kinc <- 0
    for (inc in 1:maxinc) {
        mobj <- inc/steps
        iter <- mstep(mobj, kinc, p, qp, eqm1)
        kinc <- iter$kinc
        d <- matrix(iter$d, p, p)
        binc <- sx$v %*% d %*% comp
        linc <- (binc - beta)^2/trus^2
        tinc <- sum(linc)
        dinc <- diag(d)
        minc <- p - sum(dinc)
        bstar <- cbind(bstar, binc)
        loss <- cbind(loss, linc)
        delta <- cbind(delta, dinc)
        totlos <- rbind(totlos, tinc)
        konst <- rbind(konst, kinc)
        mcal <- rbind(mcal, minc)
    }
    bstar <- t(bstar)
    loss <- t(loss)
    delta <- t(delta)
    sext <- cbind(totlos, konst, mcal)
    dimnames(sext) <- list(0:maxinc, c("TLOSS", "KONST", "MCAL"))
    RXolist <- c(RXolist, list(coef = bstar, rsel = loss, spat = delta, 
        sext = sext))
    class(RXolist) <- "true.simu"
    RXolist
}

"plot.true.risk" <-
function (x, trace = "all", trkey = FALSE, ...) 
{
    mcal <- x$sext[, 3]
    mcalp <- rep(mcal, times = x$p)
    if (trace != "coef" && trace != "rmse" && trace != "exev" && 
        trace != "infd" && trace != "spat" && trace != "seq")
        trace <- "all"
    opar <- par(no.readonly = TRUE)
    on.exit(par(opar))
    if (trace == "all")
        par(mfrow=c(3,2))
    else
        par(mfrow=c(1,1))
    if (trace == "all" || trace == "seq" || trace == "coef") {
        plot(mcalp, x$coef, ann = FALSE, type = "n")
        abline(h = 0, col = gray(0.9), lwd = 2)
        for (i in 1:x$p) lines(mcal, x$coef[, i], col = i, lty = i, 
            lwd = 2)
        title(main = paste("EXPECTED COEFFICIENTS: Q-shape =", 
            x$qp), xlab = "m = Multicollinearity Allowance", 
            ylab = "Expected Coefficients")
        if( trkey )
            legend("bottom",all.vars(x$form)[2:(x$p+1)], col=1:(x$p), lty=1:(x$p), lwd=2)
    }
    if (trace == "seq") {
        cat("\nPress the Enter key to view the RMSE trace...")
        scan()
    }
    if (trace == "all" || trace == "seq" || trace == "rmse") {
        plot(mcalp, x$rmse, ann = FALSE, type = "n")
        abline(h = 0, col = gray(0.9), lwd = 2)
        for (i in 1:x$p) lines(mcal, x$rmse[, i], col = i, lty = i, 
            lwd = 2)
        title(main = paste("RELATIVE TRUE MeanSqError: Q-shape =", 
            x$qp), xlab = "m = Multicollinearity Allowance", 
            ylab = "Scaled True MSE Risk")
        if( trkey )
            legend("bottom",all.vars(x$form)[2:(x$p+1)], col=1:(x$p), lty=1:(x$p), lwd=2)
    }
    if (trace == "seq") {
        cat("\nPress the Enter key to view the EXEV trace...")
        scan()
    }
    if (trace == "all" || trace == "seq" || trace == "exev") {
        plot(mcalp, x$exev, ann = FALSE, type = "n")
        abline(h = 0, col = gray(0.9), lwd = 2)
        for (i in 1:x$p) lines(mcal, x$exev[, i], col = i, lty = i, 
            lwd = 2)
        title(main = paste("TRUE EXCESS EIGENVALUES: Q-shape =", 
            x$qp), xlab = "m = Multicollinearity Allowance", 
            ylab = "Least Squares minus Ridge")
        if( trkey )
            legend("bottom",paste("Component",1:(x$p)), col=1:(x$p), lty=1:(x$p), lwd=2)
    }
    if (trace == "seq") {
        cat("\nPress the Enter key to view the INFD trace...")
        scan()
    }
    if (trace == "all" || trace == "seq" || trace == "infd") {
        plot(mcalp, x$infd, ann = FALSE, type = "n")
        abline(h = 0, col = gray(0.9), lwd = 2)
        for (i in 1:x$p) lines(mcal, x$infd[, i], col = i, lty = i, 
            lwd = 2)
        title(main = paste("TRUE INFERIOR DIRECTION: Q-shape =", 
            x$qp), xlab = "m = Multicollinearity Allowance", 
            ylab = "Inferior Direction Cosines")
        if( trkey )
            legend("bottom",all.vars(x$form)[2:(x$p+1)], col=1:(x$p), lty=1:(x$p), lwd=2)
    }
    if (trace == "seq") {
        cat("\nPress the Enter key to view the SPAT trace...")
        scan()
    }
    if (trace == "all" || trace == "seq" || trace == "spat") {
        plot(mcalp, x$spat, ann = FALSE, type = "n")
        abline(h = 0, col = gray(0.9), lwd = 2)
        for (i in 1:x$p) lines(mcal, x$spat[, i], col = i, lty = i, 
            lwd = 2)
        title(main = paste("SHRINKAGE PATTERN: Q-shape =", x$qp), 
            xlab = "m = Multicollinearity Allowance", ylab = "Ridge Delta Factors")
        if( trkey )
            legend("bottom",paste("Component",1:(x$p)), col=1:(x$p), lty=1:(x$p), lwd=2)
    }
}

"plot.true.simu" <-
function (x, trace = "all", trkey = FALSE, ...) 
{
    mcal <- x$sext[, 3]
    mcalp <- rep(mcal, times = x$p)
    if (trace != "coef" && trace != "rsel" &&
        trace != "spat" && trace != "seq")
        trace <- "all"
    opar <- par(no.readonly = TRUE)
    on.exit(par(opar))
    if (trace == "all")
        par(mfrow=c(2,2))
    else
        par(mfrow=c(1,1))
    if (trace == "all" || trace == "seq" || trace == "coef") {
        plot(mcalp, x$coef, ann = FALSE, type = "n")
        abline(h = 0, col = gray(0.9), lwd = 2)
        for (i in 1:x$p) lines(mcal, x$coef[, i], col = i, lty = i, 
            lwd = 2)
        title(main = paste("SIMULATED COEFFICIENTS: Q-shape =", 
            x$qp), xlab = "m = Multicollinearity Allowance", 
            ylab = "Simulated Coefficients")
        if( trkey )
            legend("bottom",all.vars(x$form)[2:(x$p+1)], col=1:(x$p), lty=1:(x$p), lwd=2)
    }
    if (trace == "seq") {
        cat("\nPress the Enter key to view the RSE Loss trace...")
        scan()
    }
    if (trace == "all" || trace == "seq" || trace == "rsel") {
        plot(mcalp, x$rsel, ann = FALSE, type = "n")
        abline(h = 0, col = gray(0.9), lwd = 2)
        for (i in 1:x$p) lines(mcal, x$rsel[, i], col = i, lty = i, 
            lwd = 2)
        title(main = paste("RELATIVE SqError LOSS: Q-shape =", 
            x$qp), xlab = "m = Multicollinearity Allowance", 
            ylab = "Scaled True SE Loss")
        if( trkey )
            legend("bottom",all.vars(x$form)[2:(x$p+1)], col=1:(x$p), lty=1:(x$p), lwd=2)
    }
    if (trace == "seq") {
        cat("\nPress the Enter key to view the SPAT trace...")
        scan()
    }
    if (trace == "all" || trace == "seq" || trace == "spat") {
        plot(mcalp, x$spat, ann = FALSE, type = "n")
        abline(h = 0, col = gray(0.9), lwd = 2)
        for (i in 1:x$p) lines(mcal, x$spat[, i], col = i, lty = i, 
            lwd = 2)
        title(main = paste("SHRINKAGE PATTERN: Q-shape =", x$qp), 
            xlab = "m = Multicollinearity Allowance", ylab = "Ridge Delta Factors")
        if( trkey )
            legend("bottom",paste("Component",1:(x$p)), col=1:(x$p), lty=1:(x$p), lwd=2)
    }
}

"print.true.risk" <-
function (x, ...) 
{
    cat("\ntrue.risk Object: True Risk of Shrinkage in Regression\n")
    cat("Data Frame:", x$data, "\n")
    cat("Regression Equation:\n")
    print(x$form)
    cat("\n    Number of Regressor Variables, p =", x$p, "\n")
    cat("    Number of Observations, n =", x$n, "\n")
    cat("\nPrincipal Axis Summary Statistics of Ill-Conditioning...\n")
    print.default(x$prinstat, quote = FALSE)
    cat("\n\ntrue.risk: Shrinkage PATH Shape =", x$qp, "\n")
    cat("\nExtent of shrinkage statistics...\n")
    print.default(x$sext, quote = FALSE)
}

"print.true.simu" <-
function (x, ...) 
{
    cat("\ntrue.simu Object: Simulated Loss of Shrinkage in Regression\n")
    cat("Data Frame:", x$data, "\n")
    cat("Regression Equation:\n")
    print(x$form)
    cat("\n    Number of Regressor Variables, p =", x$p, "\n")
    cat("    Number of Observations, n =", x$n, "\n")
    cat("\nPrincipal Axis Summary Statistics of Ill-Conditioning...\n")
    print.default(x$prinstat, quote = FALSE)
    cat("\n\ntrue.simu: Shrinkage PATH Shape =", x$qp, "\n")
    cat("\nExtent of shrinkage statistics...\n")
    print.default(x$sext, quote = FALSE)
    cat("\nSimulated Response and Expected Values...\n")
    print.default(x$ydat, quote = FALSE)
}
