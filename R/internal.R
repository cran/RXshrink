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
  m 
} 
