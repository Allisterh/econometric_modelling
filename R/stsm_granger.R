#============================================================================
#
#   Program to perform Granger causality tests on a VAR
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
# Log-likelihood function for an unrestricted VARMA model
#--------------------------------------------------------------------------
neglog <- function( b,y ) {
  t <- nrow(y)
  n <- ncol(y)
  v      <- array(0, c(t,4))
  lf      <- array(0, c(t-2,1))
  
  # First loop over MA part  
  for (i in 3:t) {
      v[i,1] <- y[i,1] -  b[1] -  b[2]*y[i-1,1] -  b[3]*y[i-1,2] - b[4]*y[i-1,3] -  b[5]*y[i-1,4] -  b[6]*y[i-2,1] - b[7]*y[i-2,2]  - b[8]*y[i-2,3]  - b[9]*y[i-2,4]

      v[i,2] <- y[i,2] - b[10] - b[11]*y[i-1,1] - b[12]*y[i-1,2] - b[13]*y[i-1,3] - b[14]*y[i-1,4] - b[15]*y[i-2,1] - b[16]*y[i-2,2] - b[17]*y[i-2,3] - b[18]*y[i-2,4]

      v[i,3] <- y[i,3] - b[19] - b[20]*y[i-1,1] - b[21]*y[i-1,2] - b[22]*y[i-1,3] - b[23]*y[i-1,4] - b[24]*y[i-2,1] - b[25]*y[i-2,2] - b[26]*y[i-2,3] - b[27]*y[i-2,4]

      v[i,4] <- y[i,4] - b[28] - b[29]*y[i-1,1] - b[30]*y[i-1,2] - b[31]*y[i-1,3] - b[32]*y[i-1,4] - b[33]*y[i-2,1] - b[34]*y[i-2,2] - b[35]*y[i-2,3] - b[36]*y[i-2,4]      
  }
  v  <- trimr( v,2,0 )    
  vc <- t(v) %*% v/nrow(v)

  for (i in seq(t-2)) {
     lf[i] <- -0.5*n*log(2*pi) - 0.5*log(det(vc)) - 0.5*v[i,] %*% inv(vc) %*% cbind(v[i,])      
  }
  f <- -mean(lf)
  return(f)
}


#
#--------------------------------- Granger Causality ------------------------
#

stsm_granger <- function( ) {

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

  

  # Construct variables for use in VAR
    # interest rate and the annual percentage growth rates in money, price and output
  yvar <- cbind(r,lm, lp, lo)
  tmp  <- 100*(trimr(yvar[,2:4],12,0) - trimr(yvar[,2:4],0,12))
  y    <- cbind(trimr(yvar[,1],12,0), tmp)

  # Granger causality tests (based on the Wald test)
  t   <- nrow(y)

  # Obtain estimates of the VAR(2) by OLS to each equation  
  bols   <- lm(trimr(y,2,0) ~ cbind(rep(1, t-2), trimr(y,1,1), trimr(y,0,2)) - 1)$coef               
  theta0 <- c(bols)   
  
  cat('Estimating', length(theta0), 'parameters. This may take several minutes...')
           
  # Estimate model: will converge in one step
  estResults <- optim(theta0, neglog, y=y, method="BFGS", hessian=T)
  theta1 <- estResults$par
  fval <- estResults$value
  H <- estResults$hessian
  
  cov1 <- (1/t)*inv(H)
  
  
  # Wald test of no causality from money to interest rates
  r <- array(0, c(2,36))  
  r[1,3] <- 1.0
  r[2,7] <- 1.0
  q      <- c(rep(0, nrow(r)))
  wd     <- t( (r %*% theta1 - q) ) %*% inv(r %*% cov1 %*% t(r)) %*% (r %*% theta1 - q)
  dof    <- nrow(r)
  
  cat('\n ')
  cat('\nWald test (money to interest rates)  = ',wd) 
  cat('\nNumber of degrees of freedom         = ',dof) 
  cat('\np-value                              = ',1-pchisq(wd,dof)) 
  cat('\n ')

  # Wald test of no causality from prices to interest rates
  r <- array(0, c(2,36))
  r[1,4] <- 1.0
  r[2,8] <- 1.0
  q      <- c(rep(0, nrow(r)))
  wd     <- t( (r %*% theta1 - q) ) %*% inv(r %*% cov1 %*% t(r)) %*% (r %*% theta1 - q)
  dof    <- nrow(r)
  cat('\n ')
  cat('\nWald test (prices to interest rates) = ',wd) 
  cat('\nNumber of degrees of freedom         = ',dof) 
  cat('\np-value                              = ',1-pchisq(wd,dof)) 
  cat('\n ')

  # Wald test of no causality from output to interest rates
  r <- array(0, c(2,36))
  r[1,5] <- 1.0
  r[2,9] <- 1.0
  q      <- c(rep(0, nrow(r)))
  wd     <- t( (r %*% theta1 - q) ) %*% inv(r %*% cov1 %*% t(r)) %*% (r %*% theta1 - q)
  
  cat('\n ')
  cat('\nWald test (output to interest rates) = ',wd) 
  cat('\nNumber of degrees of freedom         = ',dof) 
  cat('\np-value                              = ',1-pchisq(wd,dof)) 
  cat('\n ')

    
  # Wald test of no causality from all to interest rates  
  r      <- r <- array(0, c(6,36))
  r[1,3] <- 1.0
  r[2,7] <- 1.0
  r[3,4] <- 1.0
  r[4,8] <- 1.0
  r[5,5] <- 1.0
  r[6,9] <- 1.0
  q      <- c(rep(0, nrow(r)))
  wd     <- t( (r %*% theta1 - q) ) %*% inv(r %*% cov1 %*% t(r)) %*% (r %*% theta1 - q)
  dof    <- nrow(r)
  cat('\n ')
  cat('\nWald test (all to interest rates)    = ',wd) 
  cat('\nNumber of degrees of freedom         = ',dof) 
  cat('\np-value                              = ',1-pchisq(wd,dof)) 
  cat('\n ')

  # Wald test of no causality from interest rates to money
  r      <- array(0, c(2,36))
  r[1,11] <- 1.0
  r[2,15] <- 1.0
  q      <- c(rep(0, nrow(r)))
  wd     <- t( (r %*% theta1 - q) ) %*% inv(r %*% cov1 %*% t(r)) %*% (r %*% theta1 - q)
  dof    <- ncol(r)
  cat('\n ')
  cat('\nWald test (interest rates to money)  = ',wd) 
  cat('\nNumber of degrees of freedom         = ',dof) 
  cat('\np-value                              = ',1-pchisq(wd,dof)) 
  cat('\n ')

  # Wald test of no causality from prices to money
  r      <- array(0, c(2,36))
  r[1,13] <- 1.0
  r[2,17] <- 1.0
  q      <- c(rep(0, nrow(r)))
  wd     <- t( (r %*% theta1 - q) ) %*% inv(r %*% cov1 %*% t(r)) %*% (r %*% theta1 - q)
  dof    <- nrow(r)
  dof <- ncol(r)

  cat('\n ')
  cat('\nWald test (prices to money)          = ',wd) 
  cat('\nNumber of degrees of freedom         = ',dof) 
  cat('\np-value                              = ',1-pchisq(wd,dof)) 
  cat('\n ')
 
  # Wald test of no causality from output to money    
  r      <- array(0, c(2,36))
  r[1,14] <- 1.0
  r[2,18] <- 1.0
  q      <- c(rep(0, nrow(r)))
  wd     <- t( (r %*% theta1 - q) ) %*% inv(r %*% cov1 %*% t(r)) %*% (r %*% theta1 - q)
  dof    <- nrow(r)
  dof <- ncol(r)

  cat('\n ')
  cat('\nWald test (output to money)          = ',wd) 
  cat('\nNumber of degrees of freedom         = ',dof) 
  cat('\np-value                              = ',1-pchisq(wd,dof)) 
  cat('\n ')

    
  # Wald test of no causality from all to money    
  r      <- array(0, c(6,36))
  r[1,11] <- 1.0
  r[2,15] <- 1.0
  r[3,13] <- 1.0
  r[4,17] <- 1.0
  r[5,14] <- 1.0  
  r[6,18] <- 1.0
  q      <- c(rep(0, nrow(r)))
  wd     <- t( (r %*% theta1 - q) ) %*% inv(r %*% cov1 %*% t(r)) %*% (r %*% theta1 - q)
  dof    <- nrow(r)
  cat('\n ')
  cat('\nWald test (all to money)             = ',wd)
  cat('\nNumber of degrees of freedom         = ',dof)
  cat('\np-value                              = ',1-pchisq(wd,dof))
  cat('\n ')
    
  # Wald test of no causality from interest rates to prices  
  r      <- array(0, c(2,36))
  r[1,20] <- 1.0
  r[2,24] <- 1.0
  q      <- c(rep(0, nrow(r)))
  wd     <- t( (r %*% theta1 - q) ) %*% inv(r %*% cov1 %*% t(r)) %*% (r %*% theta1 - q)
  dof    <- ncol(r)
  
    
  cat('\n ')
  cat('\nWald test (interest rates to prices) = ',wd)
  cat('\nNumber of degrees of freedom         = ',dof)
  cat('\np-value                              = ',1-pchisq(wd,dof))
  cat('\n ')
    
  # Wald test of no causality from money to prices  
  r      <- array(0, c(2,36))
  r[1,21] <- 1.0
  r[2,25] <- 1.0
  q      <- c(rep(0, nrow(r)))
  wd     <- t( (r %*% theta1 - q) ) %*% inv(r %*% cov1 %*% t(r)) %*% (r %*% theta1 - q)
  dof    <- ncol(r)

  cat('\n ')
  cat('\nWald test (money to prices)          = ',wd)
  cat('\nNumber of degrees of freedom         = ',dof)
  cat('\np-value                              = ',1-pchisq(wd,dof))
  cat('\n ')

  # Wald test of no causality from output to prices  
  r      <- array(0, c(2,36))
  r[1,23] <- 1.0
  r[2,27] <- 1.0
  q      <- c(rep(0, nrow(r)))
  wd     <- t( (r %*% theta1 - q) ) %*% inv(r %*% cov1 %*% t(r)) %*% (r %*% theta1 - q)
  dof    <- ncol(r)
    
  cat('\n ')
  cat('\nWald test (output to prices)         = ',wd)
  cat('\nNumber of degrees of freedom         = ',dof)
  cat('\np-value                              = ',1-pchisq(wd,dof))
  cat('\n ')

  # Wald test of no causality from all to prices   
  r      <- array(0, c(6,36))
  r[1,20] <- 1.0
  r[2,24] <- 1.0
  r[3,21] <- 1.0
  r[4,25] <- 1.0
  r[5,23] <- 1.0  
  r[6,27] <- 1.0
  q      <- c(rep(0, nrow(r)))
  wd     <- t( (r %*% theta1 - q) ) %*% inv(r %*% cov1 %*% t(r)) %*% (r %*% theta1 - q)
  dof    <- nrow(r)
  
  cat('\n ')
  cat('\nWald test (all to prices)            = ',wd)
  cat('\nNumber of degrees of freedom         = ',dof)
  cat('\np-value                              = ',1-pchisq(wd,dof))
  cat('\n ')

  # Wald test of no causality from interest rates to output  
  r      <- array(0, c(2,36))
  r[1,29] <- 1.0
  r[2,33] <- 1.0
  q      <- c(rep(0, nrow(r)))
  wd     <- t( (r %*% theta1 - q) ) %*% inv(r %*% cov1 %*% t(r)) %*% (r %*% theta1 - q)
  dof    <- ncol(r)
    
  cat('\n ')
  cat('\nWald test (interest rates to output) = ',wd)
  cat('\nNumber of degrees of freedom         = ',dof)
  cat('\np-value                              = ',1-pchisq(wd,dof))
  cat('\n ')

  # Wald test of no causality from money to output
  r       <- array(0, c(2,36))
  r[1,30] <- 1.0
  r[2,34] <- 1.0
  q      <- c(rep(0, nrow(r)))
  wd     <- t( (r %*% theta1 - q) ) %*% inv(r %*% cov1 %*% t(r)) %*% (r %*% theta1 - q)
  dof    <- ncol(r)
  
  cat('\n ')
  cat('\nWald test (money to output)          = ',wd)
  cat('\nNumber of degrees of freedom         = ',dof)
  cat('\np-value                              = ',1-pchisq(wd,dof))
  cat('\n ')
     
  # Wald test of no causality from prices to output   **/
  r      <- array(0, c(2,36))
  r[1,31] <- 1.0
  r[2,35] <- 1.0
  q      <- c(rep(0, nrow(r)))
  wd     <- t( (r %*% theta1 - q) ) %*% inv(r %*% cov1 %*% t(r)) %*% (r %*% theta1 - q)
  dof    <- ncol(r)
   
  cat('\n ')
  cat('\nWald test (prices to output)         = ',wd)
  cat('\nNumber of degrees of freedom         = ',dof)
  cat('\np-value                              = ',1-pchisq(wd,dof))
  cat('\n ')
   # Wald test of no causality from all to output
  r       <- array(0, c(6,36))
  r[1,29] <- 1.0
  r[2,33] <- 1.0
  r[3,30] <- 1.0
  r[4,34] <- 1.0
  r[5,31] <- 1.0  
  r[6,35] <- 1.0
  q      <- c(rep(0, nrow(r)))
  wd     <- t( (r %*% theta1 - q) ) %*% inv(r %*% cov1 %*% t(r)) %*% (r %*% theta1 - q)
  dof    <- nrow(r)
  cat('\n ')
  cat('\nWald test (all to output)            = ',wd)
  cat('\nNumber of degrees of freedom         = ',dof)
  cat('\np-value                              = ',1-pchisq(wd,dof))

}

