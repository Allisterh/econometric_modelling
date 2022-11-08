#============================================================================
#
#     Program to estimate a symmetric BEKK MGARCH model of US yields
#       with conditionsl multivariate Student t disturbance
#
#============================================================================

rm(list = ls(all=T))
graphics.off()

#
# ------------------------ Helper Functiona ---------------------------------
#
# Load required functions - trimr
source("EMTSUtil.R")

#----------------------------------------------------------------------------
# Log-likelihood function for symmetric BEKK model
#----------------------------------------------------------------------------
neglog <- function( b,y ) {
    t <- nrow(y)
    n <- ncol(y)
    
    f <- rep(0, t)
    
    h <- cov(y)     # Initialise  conditional covariance matrix
    u <- apply(y, 2, sd)     # Initialise the disturbance vector
    #u <- t(u)
    
    c  <- matrix(c(b[1], 0.0 ,
                   b[2],   b[3]), nrow=2, byrow=T)
    
    a <-  matrix(c(b[4],b[5],
                   b[5],  b[6]), nrow=2, byrow=T)
    
    d <-  matrix(c(b[7],  b[8],
                   b[8], b[9]), nrow=2, byrow=T)
    nu <- b[12]
    for (i in seq(t)) {
      m    <- b[c(10, 11)]                   # Update conditional mean
      u    <- y[i,]- m                   # Update residuals 
      const <- gamma( (nu+1)/2 ) / ( sqrt(pi*(nu-2)) * gamma( nu/2 ) )                   
      f[i]  <- log(const) - 0.5*log(det(h)) - 0.5*(nu+1)*log( 1 + u %*% inv(h) %*% cbind(u/(nu-2)) )   
            
      h    <- c %*% t(c) + a %*% (u %*% t(u)) %*% t(a) + d %*% h %*% t(d)  # Update conditional covariance matrix
    }
    lf <- -mean( f )
    return(lf)
}
#
#-------------------- Mgarch Student t Model of US Yields  ------------------
#
mgarch_student <- function() {
  # Load data
  load('yields_us.Rdata')
  
  # Choose variables
  r  <- rdata[,c(1,2)]
  y  <- 100*(trimr( r,1,0 ) - trimr( r,0,1 )) 
  t  <- length( y )
  
  # Estimate the BEKK Aymmetric MGARCH(1,1) model             
  start <- c(0.8062958419864796,
             0.2816730578640683, 
             0.4297482870507287, 
             0.2685150842768220, 
             -0.0075466373201676, 
             0.1615508008596290, 
             0.9269147664364472, 
             0.0083884990798515, 
             0.9724763092352805, 
             -0.1336772926059021,
             -0.0906508454299734, 
             8.3673242655198212)
  
  estResults <- optim(start, neglog, y=y, method="BFGS")
  theta <- estResults$par
  
  cat('\nBEKK parameters')
  print(cbind(theta))
  
}
