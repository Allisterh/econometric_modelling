#==============================================================================
#
#     Simulation example to reproduce the Beach and MacKinnon (1978)
#     Econometrica, pp.51-58 study which derives the sampling
#     distributions MLE estimators of regression models with
#     autocorrelation.
#
#==============================================================================
rm (list = ls(all=TRUE))
graphics.off()
set.seed(12345, kind="Mersenne-Twister")

#
#------------------------- Helper Functions------------------------------------
#
#load required functions - recserar, trimr
source("EMTSUtil.R")

#-----------------------------------------------------------------------
#      Negative log-likelihood function   
#-----------------------------------------------------------------------
neglog <- function(b,y,x,flag) {
  beta0 <- b[1]                       
  beta1 <- b[2]  
  rho1  <- tanh(b[3])   #  rho1 stays in the unit circle
  sig2  <- abs(b[4])    #  variance is positive     

  u     <- y - beta0 - beta1*x
  v     <- trimr(u,1,0) - rho1*trimr(u,0,1)
  
  lnl_0 <- -0.5*log(2*pi) - 0.5*log(sig2) + 0.5*log(1 - rho1^2) - 0.5*(u[1] - 0)^2/(sig2/(1 - rho1^2))                
  lnl_1 <- -0.5*log(2*pi) - 0.5*log(sig2) - 0.5*v^2/sig2                                                           

  if (flag) {
    lnl <- -mean(lnl_1)
  }    
  else {    
    lnl <- -mean(c(lnl_0, lnl_1))
  }
  return(lnl)  
}
  
#
#--------------------- Beach-MacKinnon MC Study ----------------------------
#
  
auto_beachmack <-function() {
  # Simulate a regression model with an AR(1) disturbance term     
  t      <- 20                                                                                       
  beta0  <- 1.0
  beta1  <- 1.0
  rho1   <- 0.6
  sigma2 <- 0.0036
  theta0 <- c(beta0, beta1, rho1, sigma2)
    
  # xt is fixed in repeated samples
  x <- exp(0.04*seq(t)) + sqrt(0.0009)*rnorm(t)   
                                                     

  ndraws <- 200          # used by Beach and MacKinnon  
  zeros <- array(0, c(ndraws,4))
  theta_exact <- zeros
  theta_cond  <- zeros
  theta_ols   <- zeros
  pb <- txtProgressBar(min = 0, max=ndraws, style=3)
  for (k in seq(ndraws)) {
    #     Generate data      
    v <- sqrt(sigma2)*rnorm(t) 
    u <- recserar( cbind(v) , cbind(sqrt(1/(1-rho1^2))*v[1]) , cbind(rho1) )
    y <- beta0 + beta1*x + u

    # Exact MLE
    flag <- FALSE    
    estResults <- optim(theta0, neglog, y=y, x=x, flag=flag, method="BFGS")
    theta <- estResults$par      
    theta_exact[k,] <- c(theta[1:2], tanh(theta[3]), abs(theta[4]))
    

    # Conditional MLE       
    flag <-  TRUE
    estResults <- optim(theta0, neglog, y=y, x=x, flag=flag, method="BFGS")
    theta <- estResults$par       
    theta_cond[k,] <- c(theta[1:2], tanh(theta[3]), abs(theta[4]))      
  
    # OLS
    xx <- cbind(rep(1,t),x) 
    b_ols <- lm(y ~ xx - 1)$coef
    e    <- y - xx %*% b_ols
    sig2_ols <- mean(e^2)
    rho_ols <- lm(trimr(e,1,0) ~ trimr(e, 0,1) - 1)$coef
    theta_ols[k,] <- c(b_ols , sig2_ols , rho_ols) 
      
    setTxtProgressBar(pb, k)
  }
  close(pb)
  
  #  Compute statistics of sampling distributions and print results    
  rmse_exact <- sqrt( mean( (colMeans(theta_exact) - theta0)^2 ) )
  rmse_cond  <- sqrt( mean( (colMeans(theta_cond) - theta0)^2 ) )
  rmse_ols   <- sqrt( mean( (colMeans(theta_ols) - theta0)^2 ) )

  cat('\n                                 Beta0   Beta1      Rho1        Sigma^2')
  cat('\nPopulation parameter       =      ', theta0[1], '     ', theta0[2], '       ', theta0[3], '      ', theta0[4])
  cat('\n------------------------------------------------------------------------------')  
  cat('\nMean    (exact MLE)        = ', colMeans(theta_exact))
  cat('\nBias (exact MLE)           = ', theta0-colMeans(theta_exact))
  
  cat('\n')
  cat('\nMean (cond. MLE)           = ', colMeans(theta_cond))
  cat('\nBias (cond. MLE)           = ', theta0-colMeans(theta_cond))
  
  cat('\n')
  cat('\nMean (OLS)                 = ', colMeans(theta_ols))
  cat('\nBias (OLS)                 = ', theta0-colMeans(theta_ols))
        
  cat('\n')
  cat('\n\nRMSE (exact MLE)         = ', rmse_exact)
  cat('\n\nRMSE (cond. MLE)         = ', rmse_cond)
  cat('\n\nRMSE (OLS)               =', rmse_ols)
        
  cat('\n\n')
  cat('\nEfficiency (cond/exact)    = ', rmse_cond/rmse_exact)
  cat('\nEfficiency (ols/exact)     = ', rmse_ols/rmse_exact)  
}    


