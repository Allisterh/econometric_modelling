#============================================================================
#
#   Simulate and estimate a VARMA model with multivariate bilinearity
#
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()
set.seed(1234567, kind="Mersenne-Twister")

#
#------------------------- Helper Functions----------------------------------
#

#load required functions -  inv, trimr
source("EMTSUtil.R")

#--------------------------------------------------------------------------
# Log-likelihood function for an unrestricted VARMA model
#--------------------------------------------------------------------------
neglogl <- function( b,y ){
  t <- nrow(y)
  n <- ncol(y)
  
  e1      <- rep(0,t)
  e2      <- rep(0, t)  
  lf      <- rep(0,  t-1)
  
  for (i in 2:t) {
    e1[i] <- y[i,1]-b[1]-b[2]*y[i-1,1]-b[3]*e1[i-1]-b[4]*e2[i-1] - b[9]*y[i-1,1]*e1[i-1]
    e2[i] <- y[i,2]-b[5]-b[6]*y[i-1,2]-b[7]*e1[i-1]-b[8]*e2[i-1] - b[10]*y[i-1,2]*e2[i-1]    
  }
  e  <- cbind(trimr( e1,1,0 ), trimr( e2,1,0 ))
  vc <- t(e) %*% e/(t-1)
  
  for (i in seq(t-1)) {
    lf[i] <- -0.5*n*log(2*pi) - 0.5*log(det(vc))- 0.5*e[i,] %*% inv(vc) %*% cbind(e[i,])    
  }
  f <- -mean( lf )    
  return(f)
}

#--------------------------------------------------------------------------
# Log-likelihood function for an restricted VARMA model
#--------------------------------------------------------------------------
negloglr <- function( b,y ) {
  t <- nrow(y)
  n <- ncol(y)
  
  e1      <- rep(0,t)
  e2      <- rep(0, t)  
  lf      <- rep(0,  t-1)
  
  # First loop over MA part  
  for (i in 2:t) {
     e1[i] <- y[i,1]-b[1]-b[2]*y[i-1,1]-b[3]*e1[i-1]-b[4]*e2[i-1]
     e2[i] <- y[i,2]-b[5]-b[6]*y[i-1,2]-b[7]*e1[i-1]-b[8]*e2[i-1]    
  }
  e  <- cbind(trimr( e1,1,0 ), trimr( e2,1,0 ))
  vc <- t(e) %*% e/(t-1)
  
  for (i in seq(t-1)) {
    lf[i] <- -0.5*n*log(2*pi) - 0.5*log(det(vc)) - 0.5*e[i,] %*% inv(vc) %*% cbind(e[i,])   
  }
  f <- -mean( lf )
  return(f) 
}

#
#------------------------  VARMA Model simulation ----------------------------
#

stsm_varmab <- function() {
    # Generate the data
  # Equation 1
  nobs <- 1000                    
  mu1  <- 0.0 
  phi111 <- 0.6 
  si111 <- 0.2 
  si112<- -0.5 
  
  # Equation 2    
  mu2  <- 0.0
  phi122 <- 0.4 
  si121 <- 0.2 
  si122 <- 0.6 

  nobs <- nobs + 100          # Allow for startup  
  y1   <- rep(0, nobs)
  y2   <- rep(0, nobs)
  v1   <- rnorm(nobs)
  v2   <- rnorm(nobs)

  for (t in 2:nobs) {
    y1[t] <- mu1 + phi111*y1[t-1] + v1[t] + si111*v1[t-1] + si112*v2[t-1]
    y2[t] <- mu2 + phi122*y2[t-1] + v2[t] + si121*v1[t-1] + si122*v2[t-1]    
  }
  y <- cbind(trimr(y1,100,0), trimr(y2,100,0))

  t <- nrow(y)
  # Estimate the unrestricted model
  bstart <- c(mu1, phi111,si111,si112,mu2,phi122,si121,si122,0.01,0.01)
  estResults <- optim(bstart, neglogl, y=y, method="BFGS", hessian=T)              
  theta <- estResults$par
  fvalu <- estResults$value
  hess <- estResults$hessian

  # Wald test 
  vc <- (1/t)*inv( hess )
  r <- matrix(c(0,   0,   0,   0,   0,   0,   0,   0,  1,  0,
                0,   0,   0,   0,   0,   0,   0,   0,  1,  1), nrow=2, byrow=T)
  q <- rbind(0, 0)
  w <- t( (r %*% theta - q) ) %*% inv(r %*% vc %*% t(r)) %*% (r %*% theta - q)

  cat('\nWald test  = ', w) 
  cat('\np-value    = ', 1-pchisq(w,4))
     
  # Estimate the restricted model
  bstart <- c(mu1, phi111, si111, si112, mu2, phi122, si121, si122)
  estResults <- optim(bstart, negloglr, y=y, method="BFGS")                
  fvalr <- estResults$value

  # Likelihood Ratio test 
  fvalr <- -fvalr
  fvalu <- -fvalu
  lr <- -2*(t-1)*(fvalr - fvalu)
  
  cat('\n ')
  cat('\nUnconstrained log-likelihood  = ', fvalu)
  cat('\nConstrained log-likelihood    = ', fvalr)
  cat('\n ')

  cat('\nLR test  = ', lr) 
  cat('\np-value    = ', 1-pchisq(lr,4))
}
