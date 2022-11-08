#============================================================================
#
#   Recursive structural model estimated in four ways
#
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()

#
#------------------------- Helper Functions----------------------------------
#

#load required functions -  trimr, inv
source("EMTSUtil.R")


#--------------------------------------------------------------------------
# Log-likelihood function for the VAR
#--------------------------------------------------------------------------
neglog <- function(b,y,v) {
  t   <- nrow(v)
  n   <- ncol(y)
  lnl <- rep(0, t)
  b0 <-   matrix(c(1, 0,  0,  0,
                   -b[1], 1, 0, 0,
                   -b[2], -b[3],  1,  0,
                   -b[4], -b[5], -b[6], 1), nrow=4, byrow=T)
  # Structural residual variances 
  d  <- diag(n)
  d  <- diag(abs(b[7:10]))
	vc <- inv(b0) %*% d %*% t(inv(b0))
  
  for (i in seq(t)) {
    lnl[i] <- -0.5*n*log(2*pi) - 0.5*log(det(vc)) - 0.5*v[i,] %*% inv(vc) %*% cbind(v[i,])
  }
  f <- -mean(lnl)
  return(f) 
}


#
#-------------------- Impulse Responses and Variance Decompositions ---------
#

stsm_recursive <- function( ) {
    # Read the data: quarterly US data from Jan-1959 to Dec-1998
  ytdata <- as.matrix(read.table("sims_data.dat"))
 
    
  # Define variables
  r    <- ytdata[,1]        
  #lex  <- log( ytdata[,2] )
  #lcp  <- log( ytdata[,3] )
  lm   <- log( ytdata[,4] )
  lp   <- log( ytdata[,5] )
  lo   <- log( ytdata[,6] )
  #sdum <- ytdata[,7:17]

  t <- nrow(ytdata)

  # Construct variables for use in VAR
  # interest rate and the annual percentage growth rates in money, price and output
  yvar <- cbind(r,   lm,   lp,   lo)    
  tmp  <- 100*(trimr(yvar[,2:4],12,0) - trimr(yvar[,2:4],0,12))
  yvar <- cbind(trimr(yvar[,1],12,0), tmp) 
  lags <- cbind(trimr(yvar,1,1),   trimr(yvar,0,2) )

  # Structural equation approach  
  y  <- trimr(yvar[,1],2,0)
  x  <- cbind(rep(1, length(y)),   lags)
  a1 <- lm(y ~ x - 1)$coef
  u1 <- y - x %*% a1

  y  <- trimr(yvar[,2],2,0)
  x  <- cbind(rep(1, length(y)),   trimr(yvar[,1],2,0),   lags)
  a2 <- lm(y ~ x - 1)$coef
  u2 <- y - x %*% a2

  y  <- trimr(yvar[,3],2,0)
  x  <- cbind(rep(1, length(y)),   trimr(yvar[,c(1, 2)],2,0),   lags)
  a3 <- lm(y ~ x - 1)$coef
  u3 <- y - x %*% a3

  y  <- trimr(yvar[,4],2,0)
  x  <- cbind(rep(1, length(y)),   trimr(yvar[,c(1, 2, 3)],2,0),   lags)
  a4 <- lm(y ~ x - 1)$coef
  u4 <- y - x %*% a4

  # Structural residuals   and covariance 
  b01 <- matrix(c(1,  0,  0,  0,
                  -a2[2], 1,  0, 0,
                  -a3[2], -a3[3], 1, 0,
                  -a4[2],  -a4[3], -a4[4], 1), nrow=4, byrow=T)
  u  <- cbind(u1, u2, u3, u4)
  d1 <- t(u) %*% u/nrow(u)   	 

  # Reduced form approach	 
  y <- trimr(yvar,2,0)
  x <- cbind(rep(1, nrow(yvar)-2),   lags)

  bar <- lm(y ~ x - 1)$coef
  v   <- y - x %*% bar
  u1 <- v[,1]
    
  a2 <- lm(v[,2] ~ v[,1] - 1)$coef
  u2 <- v[,2] - v[,1] * a2
    
  a3 <- lm(v[,3] ~ v[,c(1,2)] - 1)$coef
  u3 <- v[,3] - v[,c(1,2)] %*% a3
  
  a4 <- lm(v[,4] ~ v[,c(1, 2, 3)] - 1)$coef
  u4 <- v[,4] - v[,c(1, 2, 3)] %*% a4

  # Structural residuals 	and covariance 
  b02 <- matrix(c(1,  0,  0,  0,
                  -a2[1], 1,  0,  0,
                  -a3[1], -a3[2],  1, 0,
                  -a4[1], -a4[2], -a4[3], 1), nrow=4, byrow=T)
  u  <- cbind(u1, u2, u3, u4)
  d2 <- t(u) %*% u/nrow(u) 		 

  # Choleski decomposition approach 
  y <- trimr(yvar,2,0)
  x <- cbind(rep(1, nrow(y)),   lags)
    
  bar    <- lm(y ~ x - 1)$coef
  v      <- y - x %*% bar
  vc     <- t(v) %*% v/nrow(v)
  s      <- t( chol(vc) )
  tmp    <- diag(s)
  b0inv  <- t( apply(s, 1, '/', tmp) )
  b03    <- inv(b0inv)
  d3    <- rep(0, 4)
  d3     <- diag(tmp^2) 

  # MLE approach  
  theta0 <- 0.1*rep(1, 10)
  # Converge in one iteration
  estResults <- optim(theta0, neglog, y=y, v=v, method="BFGS")
  theta <- estResults$par
  fval <- estResults$value

  lf <- -fval
  cat('\nLog-likelihood value     = ',lf)
  cat('\nT x Log-likelihood value = ', nrow(v)*lf)

  b04 <- matrix(c(1,  0,  0,  0,
                  -theta[1], 1, 0,  0,
                  -theta[2], -theta[3], 1, 0,
                  -theta[4], -theta[5], -theta[6], 1), nrow=4, byrow=T)
  d4 <- rep(0, 4)
  d4  <- diag(abs(theta[7:10]))
  cat('\n ')
  cat('\nB0: structural equation approach\n')	 
  print( b01 )
  cat('\nB0: reduced form approach\n' )	 
  print( b02 )
  cat('\nB0: Choleski decomposition approach\n' )	 
  print( b03 )
  cat('\nB0: maximum likelihood approach\n' )	 
  print( b04 )


  cat('\n ')
  cat('\nD: structural equation approach\n')	 
  print( d1 )
  cat('\nD: reduced form approach\n' )	 
  print( d2 )
  cat('\nD: Choleski decomposition approach\n')	 
  print( d3 )
  cat('\nD: maximum likelihood approach\n')	 
  print( d4 )
}

