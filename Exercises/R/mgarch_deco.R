#=============================================================================
#
#     Program to estimate a DECO model of US yields
#
#=============================================================================
rm(list = ls(all=T))
graphics.off()

#
# ------------------------ Helper Functions ----------------------------------
#
# Load required functions - trimr, inv
source("EMTSUtil.R")

#
#----------------------------------------------------------------------------
# Likelihood function for a GARCH(1,1) model
#----------------------------------------------------------------------------
negloglgarch <- function( b,y ) {
  u <- y       
  h <- recserar(cbind(b[1] + pnorm(b[2])*trimr(c(0.0,u^2),0,1)),cbind(sd(u)^2),cbind(pnorm(b[3])))
  z <- u/sqrt(h)
  f <- - 0.5*log( 2*pi ) - 0.5*log( h ) - 0.5*z^2
  
  lf <- -mean( f )
  return(lf)
}

#----------------------------------------------------------------------------
# Likelihood wrapper function
# This function returns the average log-likelihood.
#----------------------------------------------------------------------------
neglogldeco <- function( b,y,hv ) {
  b       <- pnorm( b )    
  t <- nrow(y)
  n <- ncol(y)
  f       <- rep(0, t)
  
  u    <- y              # Residuals                                                        
  z    <- u/sqrt(hv)    # Standardised residuals                                          
  qbar <- cov( z )
  q    <- qbar                     
  
  for (i in seq(t)) {
    # Diagonal matrix of conditional standard deviations
    s  <- diag( sqrt(hv[i,]) ) 
    
    # Conditional correlation matrix 
    tmp <- inv( diag( sqrt(diag(q)) ) )
    r   <- tmp %*% q %*% tmp     
    
    # Compute mean of covariances
    r[upper.tri(r, diag=T)] <- 0.0  # get the lower triangle
    ind  <- as.vector(r) != 0 
    rbar <- mean( r[ind] )
    
    # Redefine r for DECO    
    r <- (1 - rbar) * diag(n) + rbar * array(1, c(n,n))
    
    # Determinant and inverse of r for DECO
    rdet <- (1 - rbar)^(n-1) * ( 1 + (n - 1)*rbar)                                              
    rinv <- (1/(1 - rbar)) * diag(n) - (rbar/( (1-rbar)*(1 + (n - 1)*rbar) ) )*array(1, c(n,n))
    
    # Log of the likelihood at observation i   
    f[i] <- - 0.5*log(rdet) - 0.5*z[i,] %*% rinv %*% cbind(z[i,])

    # Update Q
    q  <- abs(1 - b[1] - b[2])*qbar + b[1]*cbind(z[i,]) %*% z[i,] + b[2]*q
  }
  lf <- -mean(f)
  return(lf)
}


#----------------------------------------------------------------------------
# Likelihood function for FULL DECO model
# This function returns the average log-likelihood.
#----------------------------------------------------------------------------
neglogfull <- function( b,y ) {
  ind <- c(2, 3, 5, 6, 8, 9, 11, 12, 14, 15, 16, 17)
  b[ind] <- pnorm(b[ind])
  
  t <- nrow(y)
  n <- ncol(y)
  f <- rep(0, t)
  
  u <- y                                       
  hv <- array(0, c(t,n))
  
  # Construct conditional variances
  hv[,1] <- recserar(cbind(b[1]+b[2]*trimr(c(0.0,u[,1]^2),0,1)),cbind(sd(u[,1])^2),cbind(b[3]))
  hv[,2] <- recserar(cbind(b[4]+b[5]*trimr(c(0.0,u[,2]^2),0,1)),cbind(sd(u[,1])^2),cbind(b[6]))      
  hv[,3] <- recserar(cbind(b[7]+b[8]*trimr(c(0.0,u[,3]^2),0,1)),cbind(sd(u[,1])^2),cbind(b[9]))
  hv[,4] <- recserar(cbind(b[10]+b[11]*trimr(c(0.0,u[,4]^2),0,1)),cbind(sd(u[,1])^2),cbind(b[12]))     
  hv[,5] <- recserar(cbind(b[13]+b[14]*trimr(c(0.0,u[,5]^2),0,1)),cbind(sd(u[,1])^2),cbind(b[15]))    
  
  # Unconditional covariance matrix  standardized residuals
  z    <- u/sqrt( hv )                                                
  qbar <- cov( z )
  q    <- qbar                                                                        
  
  for (i in seq(t)) {
    # Diagonal matrix of conditional standard deviations
    s  <- diag( sqrt(hv[i,]) ) 
    
    # Conditional correlation matrix 
    tmp <- inv( diag( sqrt(diag(q)) ) )
    r   <- tmp %*% q %*% tmp     
    
    # Compute mean of covariances
    r[upper.tri(r, diag=T)] <- 0.0  # get the lower triangle
    ind  <- as.vector(r) != 0 
    rbar <- mean( r[ind] )
    
    # Redefine r for DECO    
    r <- (1 - rbar) * diag(n) + rbar * array(1, c(n,n))
    
    # Update conditional variance-covariance matrix  
    h  <-  s %*% r %*% s                                                                                    
    
    # Likelihood function
    f[i] <- - 0.5*n*log(2*pi) - 0.5*log(det(h)) - 0.5*u[i,] %*% inv(h) %*% cbind(u[i,])
       
    # Update Q
    q  <- abs(1 - b[16] - b[17])*qbar + b[16]*cbind(z[i,]) %*% z[i,] + b[17]*q 
  }
  lf <- -mean( f )
  return(lf)
}

#
# ------------------------ MGARCH DECO Model --------------------------------
#
mgarch_deco <- function( ) { 

  # Load data: US daily yields 3-Jan-2000 to 21-Aug-2006
  data <- as.matrix(read.table("daily_finance.dat"))
  
  # Data manipulation
  data <- data[,c(1, 5, 10, 15, 20)]
  y    <- 100*(trimr(data,1,0) - trimr(data,0,1))
  y    <- y - colMeans(y)
  
  t <- nrow(y)
  n <- ncol(y)
  
  # Estimate univariate GARCH models for each variable   
  b1 <- array(0, c( 3,n ))
  h1 <- array(0, c( t,n ))
  h  <- rep(0, t)
  
  for (k in seq(n)) {
    bstart      <- c(0.1,  -2,  2)    
    estResults  <- optim(bstart, negloglgarch, y=y[,k], method="BFGS")
    b1[,k]     <- estResults$par 
    u <- y[,k]
    h1[,k] <- recserar(cbind(b1[1,k] + pnorm(b1[2,k])*trimr(c(0.0,u^2),0,1)),cbind(sd(u)^2),cbind(pnorm(b1[3,k])))
  }
  b1[c(2, 3),] <- pnorm( b1[c(2, 3),] )
  cat('\nGarch parameters\n' )
  print(b1)
  
  # Estimate the DECO model
  bstart <- c(-2, 2 )
  estResults <- optim(bstart, neglogldeco, y=y, hv=h1, method="BFGS")
  bc <- estResults$par  
  bc <- pnorm( bc )
  cat('\n', bc, '\n')
  
  start <-          c(0.999314709903659,
                           0.152115338735339,
                           0.843167761478777,
                           0.789746885859107,
                           0.064814542202766,
                           0.912043327137570,
                           0.630618257788049,
                           0.054182853876372,
                           0.920133291548674,
                           0.815035172227143,
                           0.055163278439432,
                           0.912538859743870,
                           0.796164122428437,
                           0.055205181793798,
                           0.911644837026472,
                           0.180319393683112,
                           0.953053140443055)
  
  estResults <- optim(start, neglogfull, y=y, method="BFGS")
  bf <- estResults$par
  ind <- c(2, 3, 5, 6, 8, 9, 11, 12, 14, 15, 16, 17)
  bf[ind] <- normcdf(bf(ind))
  print(cbind(cbind(c(b1), c(bc)),bf))
}

