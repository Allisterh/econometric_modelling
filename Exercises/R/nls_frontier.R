#================================================================================
#
#   Program to estimate a stochastic frontier model 
#
#================================================================================

rm (list = ls(all=TRUE))
graphics.off()
set.seed(1, kind="Mersenne-Twister")

#
#--------------------------- Helper Functions -----------------------------------
#
# load required functions - figure
source("EMTSUtil.R")

#--------------------------------------------------------------------------------
#   Log-likelihood function of the stochastic frontier model
#--------------------------------------------------------------------------------
neglog <- function(b,y,x) {
  u    <- y - b[1] - b[2]*x
  sig1 <- abs(b[3])
  sig2 <- abs(b[4])
        
  # Normal - exponential likelihood
  lf <- - log(sig2) + 0.5*sig1^2/sig2^2 + u/sig2 + log(pnorm( -(u + sig1^2/sig2 )/sig1 ) )
    
  lnf <- -mean( lf )
  return(lnf)  
}


#
#--------------------------- Stochastic Frontier Model  -------------------------
#

nls_frontier <- function() {
  # Set parameter values
  beta0 <- 1         
  beta1 <- 0.5        
  sig1  <- 1          
  sig2  <- 1.5          

  t <- 1000     
  x <- rnorm(t) 
    
  # Compute density and cumulative distribution function    
  h <- 0.01        
  u <- seq(-10, 5, h)

  f <- (1/sig2)*exp( sig1^2/(2*sig2^2) + u/sig2 )*pnorm( -(u + sig1^2/sig2)/sig1 )


  #********************************************************************
  #***
  #***     Generate graph
  #***
  #********************************************************************        
  figure()
  par(xaxs="i", yaxs="i")
  plot(u,f, type = "l",
       xlab = expression(u),
       ylab = expression(f(u)),
       bty = "l")
  
  # Monte Carlo simulation
  ndraws <- 5000
  theta  <- array(0, c(ndraws,4))
    
  # cdf to be used as a look-up table in inverse cdf method 
  cdf    <- h*cumsum(f) 
  
  pb <- txtProgressBar(min=0, max=ndraws, style=3)
  y <- array(0, c(t,1))
  for (i in seq(ndraws)) {
    # Generate random numbers for y  
    for (j in seq(t)) {
      ind <- which.min(abs(cdf - runif(1)))
      y[j]       <- beta0 + beta1*x[j] + u[ind]
    }
    # Estimate the model
    theta0     <- c(beta0, beta1, sig1,sig2)
    estResults <- optim(theta0, neglog, y=y, x=x, method="Nelder-Mead")
    bhat       <- estResults$par
        
    theta[i,] <- abs(bhat)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  cat('\n                           beta0     beta1     sig1     sig2')
  cat('\nPopulation parameter  =    ', theta0[1], '      ', theta0[2], '       ', theta0[3], '    ', theta0[4])
  cat('\nMean                  = ', colMeans(theta))    
}

