#==========================================================================
#
#    Program to estimate an exponential model
#    using the GAUSS-NEWTON algorithm 
#     
#==========================================================================

rm (list = ls(all=TRUE))
graphics.off()
set.seed(123457, kind="Mersenne-Twister")

#
#--------------------------- Helper Functions -----------------------------------
#
# load required functions - inv
source("EMTSUtil.R")


#--------------------------------------------------------------------------------
# Simulate data an return the results as a list of x and y
#--------------------------------------------------------------------------------
simulatedata <- function(t) {
  b0  <- 1.0
  b1  <- 0.05
  sig <- 0.5
  
  u <- sig*rnorm(t)
  x <- 1:t  
  
  y <- b0*exp( b1*x ) + u  
  return(list(y=y, x=x))
}

#
#--------------------------- Exponential Model by Gauss-Newton -----------------------------------
#
nls_exponential_gn <- function() {
  t <- 50
  simResults <- simulatedata(t)
  x <- simResults$x
  y <- simResults$y
  
  # GAUSS-NEWTON algorithm  

  b <- c(0.1, 0.1)     # Choose staring values       
  crit  <- 0.00001     # Convergence criterion         
  maxit <- 20          # Maximum number of iterations
  for (i in seq(maxit)) {
   e  <- y - b[1]*exp(b[2]*x)
   z1 <- exp(b[2]*x)
   z2 <- b[1]*x*exp(b[2]*x)
   
   z  <- cbind(z1, z2)
     
   cat('\nIteration        = ', i)
   bchange <- lm(e ~ z - 1)$coef
     
   if ( (t(bchange) %*% bchange) < crit) {
     break
   } else {
     b <- b + bchange       # Update new parameters     
   }        
  
   if (i == maxit)
      cat('\nFailed to converge after iteration number =', maxit, '\n')
  } 
  s2 <- c(t(e) %*% e/t)                       # Residual variance       
  
  omega <- s2 * inv( t(z) %*% z)
  
  sterr <- sqrt(diag(omega))
  tstat <- b/sterr
  
  cat('\n\nParameter estimates, std errors and t-stats\n')
  print(cbind(b, sterr, tstat))
     
  cat('\n\nEstimated asymptotic covariance matrix\n') 
  print(omega)  
}






