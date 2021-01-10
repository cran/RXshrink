"kofm" <- function (muobj, p, dMSE, delmax = 0.999999) 
{ 
    kM <- matrix(1,1,p)/dMSE[order(dMSE)]   # Vector of decreasing (1/dMSE) values...
    if (muobj <= 0) { 
        d <- diag(delmax, p)         # Matrix with Diagonal Values strictly < 1...
        kStar <- kM[1] 
        return(list(kStar = kStar, d = d)) 
    } 
    if (muobj >= p) {
        d <- matrix(0, p, p) 
        kStar <- 0                     # Terminus of the shrinkage path
        return(list(kStar = kStar, d = d))
    } 
    mVar <- matrix(1,1,p)              # Meaningless initial values...
    for( j in 1:p ) {                  # 0 < muobj < 1 below here...
        mVar[j] <- mofk(p, k=kM[j], dMSE) 
    } 
    kMin <- (p-muobj)/sum(dMSE )
    if( muobj > mVar[p] ) {            # Tail; Large m, small k Cases...
        d <- kMin*diag(dMSE) 
        kStar <- kMin 
        return(list(kStar = kStar, d = d)) 
    } else {                           # Locally Linear cases...
        for( j in 2:p ) { 
            if( mVar[j-1] < muobj && muobj <= mVar[j] ) { 
	            if( muobj == mVar[j] ) { 
                    kStar <- kM[j] 
                } else { 
                    B <- muobj - mVar[j-1] 
                    D <- mVar[j] - muobj 
                    kStar <- (B*D/(B+D))*( kM[j-1]/B + kM[j]/D ) 
                } 
            } 
        } 
        for( j in 1:p ) { 
            dd <- min(delmax,kStar*dMSE[j]) 
            if( j == 1 ) dp <- dd else dp <- c(dp,dd) 
        } 
        d <- diag(dp) 
        return(list(kStar = kStar, d = d)) 
    } 
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
  
"mofk" <- function(p, k, dMSE)  # Many-to-One Function for large k-values...
{ 
    p <- as.integer(p)
    if( k < 0 ) k <- 0
    kM <- 1/min(dMSE)
    if( k > kM ) k <- kM
    m <- as.double(0)
    for( j in 1:p ) {
      m <- m + 1 - min(1,k*dMSE[j])
    } 
    m 
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
  m 
} 
