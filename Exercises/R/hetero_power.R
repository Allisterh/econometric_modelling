#=============================================================================
#
#    Simulation example to generate the power function of the Wald test
#    of heteroskedasticity.
#                                                         
#=============================================================================

rm (list = ls(all=TRUE))
graphics.off()
set.seed(123456, kind="Mersenne-Twister")

#
#------------------------- Helper Functions------------------------------------
#
#load required functions - inv, figure
source("EMTSUtil.R")

#-----------------------------------------------------------------------------
#      Negative log-likelihood function   
#-----------------------------------------------------------------------------
neglog <- function(theta,y,x) {
  lf <- -mean( loglt(theta,y,x) )
  return(lf)
}

#-----------------------------------------------------------------------------
#      Log-likelihood function at each observation
#-----------------------------------------------------------------------------
loglt <- function(p,y,x) {
  mu   <- p[1] + p[2]*x
  sig2 <- exp(p[3] + p[4]*x)
  lft  <- -(0.5)*log(2*pi*sig2) - (y - mu)^2/(2*sig2)
  return(lft)
}

#------------------------- Power Properties ----------------------------------

hetero_power <- function() {
  t <- 50  
  n <- 10000            
  x <- runif(t)              

  # Population parameters                                    
  beta0  <- 1.0
  beta1  <- 2.0
  gam0   <- 1.0
  gam1   <- seq(0.0, 5.0, 0.5)                      

  # Pre-allocate memory for estimates
  wd        <- array(0, c(n,1))
  wd_power  <- array(0, c(length(gam1),1))

  # Restriction matrices
  R <- rbind(c(0,  0,  0,  1))
  Q <- 0

  # Loop over elements of gamma
  for (j in seq(gam1)) {
    start <- c(beta0, beta1, gam0, gam1[j])
    # Simulation loop
    pb <- txtProgressBar(min=0, max=n, style=3)
    for (k in seq(n)) {
      # Generate data    
      u <- sqrt( exp(gam0 + gam1[j]*x) )*rnorm(t)   
      y <- beta0 + beta1*x + u                    

      # Estimate parameters
      estResults <- optim(start, neglog, y=y, x=x, method="BFGS", hessian=T)
      theta <- estResults$par
      H <- estResults$hessian       
        
      # Wald test
      wd[k] <- t* t( (R %*% theta - Q) ) %*% inv(R %*% inv(H) %*% t(R)) %*% (R %*% theta - Q)
      setTxtProgressBar(pb, k)
    }    
    if (j==1) {
      c_5         <- quantile(wd,probs=0.95)
      wd_h0       <- wd 
      wd_power[j] <- 100*mean(wd > c_5 )
      cat('\n\nSize of test')
      cat('\n------------------------------')
      cat('\n5%  critical value (empirical)')
      cat('\n', c_5)
      cat('\nPower ( 5% nominal sig. level)')
      cat('\n', 100*mean( wd > qchisq(p=0.95, df=1) ) )
      cat('\nPower ( 5% nominal sig. level)')
      cat('\n', wd_power[j], '\n' )       
    } else {
      wd_power[j] <- 100*mean(wd > c_5 )
      cat('\n\nPower of test for gamma value')
      cat('\n', gam1[j])
      cat('\n------------------------------')
      cat('\nUnadjusted power')
      cat('\n', 100*mean( wd > qchisq(p=0.95, df=1) ) )
      cat('\nSize adjusted power')
      cat('\n', wd_power[j], '\n')       
    }
    close(pb)
  }
    
  #**************************************************************************
  #**
  #**     Generate graphs
  #**
  #**************************************************************************
  figure()

  # Plot the sampling distribution of the Wald statistic under the null hypothesis
  hist(wd_h0,21, xaxs="i", yaxs="i",
      main = "",
      xlab = "Values of the Wald test statistic",
      ylab = "Empirical Distribution",
      xlim = c(0, 15),
      bty = "l")
  # Plot the power function of the Wald statistic
  figure() 
  plot(gam1,wd_power, xaxs="i", yaxs="i", type="l",
      xlab = expression(paste("Values of the parameter ", gamma[1])),
      ylab = "Power (%)", 
      bty = "l")  
}
