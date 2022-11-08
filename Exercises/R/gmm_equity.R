#============================================================================
#
#   Decompose equity returns based on GMM estimates
#   Data are daily, July 29th 2004 to March 3rd 2009, SP500,FTSE100,EURO50
#
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()

# 
#---------------------------------- Helper functions ------------------------
#
# Load require functions - vech
source("EMTSUtil.R")
#----------------------------------------------------------------------------
# GMM unconstrained objective function  
#----------------------------------------------------------------------------
gmmcrit <- function(b,e2) {
  lam  <- b[1:3]
  phi  <- b[4:6]
  
  a        <- cbind(lam, diag(phi))
  
  m <- matrix(as.vector(e2) - as.vector(vech(a %*% t(a))), nrow(e2), ncol(e2), byrow=T)
  
  w <- t(m) %*% m/nrow(m)
  q <- 0.5*colMeans(m) %*% inv(w) %*%  cbind(colMeans(m))
  return(q)
}


gmm_equity <- function( ) {
  # Read in data - 1199x3 array named 'p'
  p <- read.table("equity_decomposition.dat")
  
  # Compute centered equity returns
  rt <- 100*(log(trimr(p,1,0)) - log(trimr(p,0,1)))                          
  et <- rt - colMeans(rt)
  t  <- nrow(et)
  
  # Generate squares and cross products of columns of et
  rows <- nrow(et)
  cols <- ncol(et)
  e2          <- array(0, c(rows,cols*cols))
  k<-1
  
  for (j in seq(cols)) {  
    for (i in seq(cols)) {    
      e2[,k] <- et[,j]*et[,i]    
      k<-k+1        
    }
  }
  
  # Drop repeated columns
  e2 <- e2[,vech( reshapeg(seq(cols*cols),cols,cols) ) ]
  
  # Estimate model
  start <- runif(6)
  estResults <- optim(start, gmmcrit, e2=e2, method="BFGS", hessian=T)
  theta <- estResults$par
  qu <- estResults$val
  hess <- estResults$hess
  
  vc <- inv(hess)/t
  cat('\n')
  cat('\nValue of objective function = ', qu)
  cat('\n')
  print(cbind("Params"=theta,   "Std Errors"=sqrt(diag(vc)) ))
  
  
  cat('\nJ statistic             = ', 2*t*qu)
  nu <- size(e2,2)-length(theta)
  if  (nu > 0.0)
    cat('\np-value             = ', 1-pchisq(2*t*qu,nu))
  
  
  cat('\n')
  cat('\n Covariance matrix of data \n')
  et <- as.matrix(et)
  print( t(et) %*% et/t )
  
  a <- cbind(theta[1:3], diag(theta[4:6]))
  cat('\n ')
  cat('\nCovariance matrix based on decomposition\n')
  print( a %*% t(a) )
  
  
  total <- theta[1:3]^2 + theta[4:6]^2
  cat('\n ')
  cat('\n   Common   Idiosyncratic\n ') 
  print( 100*( cbind((theta[1:3]^2)/total, (theta[4:6]^2)/total) ) )
}
