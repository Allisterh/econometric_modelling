#============================================================================
#
#   Compute GMM estimates of a contagion model
#
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()
set.seed(123, kind="Mersenne-Twister")

#
#------------------------- Helper Functions----------------------------------
#

# Load required functions - inv, trimr, reshapeg
source("EMTSUtil.R")

# Load require Library - repmat
library("matlab")

#----------------------------------------------------------------------------
# GMM uncontrained objective function  
#----------------------------------------------------------------------------
gmmcrit <- function(b,e2) {
  lam  <- b[1:7]
  phi  <- b[8:14]
  gam  <- b[15:20]
  a       <- cbind(lam, diag(phi))
  a[1:6,8] <- gam
  
  m <- t( apply(e2, 1, '-', t( vech(a %*% t(a)))) )   
  w <- t(m) %*% m/nrow(m)  
  q <- 0.5*colMeans(m) %*% inv(w) %*% cbind( colMeans(m))
  
  return(q)
}

#----------------------------------------------------------------------------
# GMM uncontrained objective function  
#----------------------------------------------------------------------------
gmmcritc <- function(b,e2) {
  lam  <- b[1:7]
  phi  <- b[8:14]
  gam  <- rep(0, 6)
  
  a        <- cbind(lam, diag(phi))
  a[1:6,8] <- gam
  m <- t( apply(e2, 1, '-', t( vech(a %*% t(a)) )))
  w <- t(m) %*% m/nrow(m)
  q <- 0.5*colMeans(m) %*% inv(w) %*% cbind( colMeans(m))
}


#
#--------------- Modelling Contagion in the Asian Crisis --------------------
#
gmm_contagion <- function() {
  # Load daily data 2 June 1997 - 31 August 1998
  # South Korea, Indonesia, Thailand, Malaysia, Australia, NZ, Japan
  st <- as.matrix(read.table("contagion.dat"))
  rt <- 100*(log(trimr(st,1,0)) - log(trimr(st,0,1)))  

  rt <- t( apply(rt, 1, '-', colMeans(rt)))
  
  et <- rt[,c(1, 2, 4, 7, 5, 6, 3)]
  t  <- nrow(et)

  # Generate squares and cross products of columns of et
  rows <- nrow(et)
  cols <- ncol(et)
  e2          <- array(0, c(rows,cols*cols))
  k <- 1
  for (j in seq(cols)) {
    for (i in seq(cols)) {
      e2[,k] <- et[,j]*et[,i]
      k <- k+1      
    }    
  }
  # Drop repeated columns
  e2 <- e2[,vech( reshapeg(seq(cols*cols),cols,cols) ) ]

  # Estimate the uncontrained model  
  theta0 <- runif(20)
  estResults <- optim(theta0, gmmcrit, e2=e2, method="BFGS")
  theta <- estResults$par
  qu <- estResults$value
 
  cat('\n ')
  cat('\nValue of objective function = ', qu)
  cat('\n ')
  cat('\nJ statistic                 = ', 2*t*qu)
  nu <- ncol(e2)-length(theta)
  if(nu > 0.0)
    cat('\np-value                     = ', 1-pchisq(2*t*qu,nu))
  
  cat('\n')
  total <- theta[1:7]^2 + theta[8:14]^2 + c(theta[15:20]^2, 0)  
  common <- (theta[1:7]^2)/total
  idiosync <- (theta[8:14]^2)/total
  contagion <- c(theta[15:20]^2, 0)/total
  print( 100*( cbind( common, idiosync,  contagion )))
  
  # Estimate the contrained model  
  theta0 <- runif(14)
  estResults <- optim(theta0, gmmcritc, e2=e2, method="BFGS")
  theta <- estResults$par
  qc <- estResults$value
  
  stat <- 2*t*(qc - qu)
  cat('\n ')
  cat('\nValue of objective function = ', qc)
  cat('\nTest of contagion = ', stat)
  cat('\nP-value           = ', 1-pchisq(stat,6))  
}

