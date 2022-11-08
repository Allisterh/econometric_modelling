#=======================================================================
#
#  Simulation example to compare the asymptotic distribution of the MLE 
#  and the OLS estimators of the autocorrelated regression model for 
#  alternative assumptions about the explanatory variable xt
#
#=======================================================================

rm (list = ls(all=TRUE))
graphics.off()
set.seed(1, kind="Mersenne-Twister")

#
#------------------------- Helper Functions-----------------------------
#
#load required functions - inv, recserar, trimr
source("EMTSUtil.R")
# load library functions - zeros, ones. repmat
library("matlab")

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

#--------------------- Efficiency of MLE and LS----------------------

auto_efficiency <- function( ) {
  #  Simulate a regression model with an AR(1) disturbance term      
  t <- 200                                                    
  beta0  <- 1.0
  beta1  <- 1.0
  rho1   <- 0.6
  sigma2 <- 0.0036
  theta0 <- c(beta0, beta1, rho1, sigma2)
  phi1   <- 0.6
  sigma2_w <- 0.0036

  ndraws <- 500             
  theta_exact  <- zeros(ndraws,4)
  theta_cond   <- zeros(ndraws,4)
  theta_ols    <- zeros(ndraws,4)

  pb <- txtProgressBar(min = 0, max=ndraws, style=3)
  for (k in seq(ndraws)) {
    #   Generate data       
    v <- sqrt(sigma2)*rnorm(t)
    u <- recserar( cbind(v) , cbind(sqrt(1/(1-rho1^2))*v[1]) , cbind(rho1) )      
    w  <- sqrt(sigma2_w)*rnorm(t)      
      
    # For the case where xt is not fixed
    x <- recserar( cbind(w) , cbind(sqrt(1/(1-phi1^2))*w[1]) , cbind(phi1) )      
    
    y <- beta0 + beta1 * x + u

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
    xx <- cbind(ones(t,1),x) 
    b_ols <- lm(y ~ xx - 1)$coef
    e    <- y - xx %*% b_ols
    sig2_ols <- mean(e^2)
    rho_ols <- lm(trimr(e,1,0) ~ trimr(e, 0,1) - 1)$coef
    theta_ols[k,] <- c(b_ols , sig2_ols , rho_ols) 
    setTxtProgressBar(pb, k)
      
  }
  close(pb)
  #  Compute statistics of sampling distributions and print results    
  mse_exact <- colMeans(((theta_exact - repmat(theta0,ndraws,1))^2))
  mse_cond  <- colMeans(((theta_cond  - repmat(theta0,ndraws,1))^2))
  mse_ols   <- colMeans(((theta_ols   - repmat(theta0,ndraws,1))^2))
    


  cat('\n                                 Beta0   Beta1      Rho1        Sigma^2')
  cat('\nPopulation parameter       =      ', theta0[1], '     ', theta0[2], '       ', theta0[3], '      ', theta0[4])
  cat('\n------------------------------------------------------------------------------')
  cat('\nMean        (exact MLE)           = ', colMeans(theta_exact))
  cat('\nBias (x100) (exact MLE)           = ' , 100*(colMeans(theta_exact)-theta0))
  cat('\nMSE  (x100) (exact MLE)           = ', 100*mse_exact)
  cat('\nRMSE        (exact MLE)           = ', sqrt(mse_exact))
  cat('\n ')
  cat('\nMean        (cond. MLE)           = ', colMeans(theta_cond))
  cat('\nBias (x100) (cond. MLE)           = ', 100*(colMeans(theta_cond)-theta0))
  cat('\nMSE  (x100) (cond. MLE)           = ', 100*mse_cond)
  cat('\nRMSE        (cond. MLE)           = ', sqrt(mse_cond))
  cat('\n ')
  cat('\nMean        (OLS)                 = ', colMeans(theta_ols))
  cat('\nBias (x100) (OLS)                 = ', 100*(colMeans(theta_ols)-theta0))
  cat('\nMSE  (x100) (OLS)                 = ', 100*mse_ols)
  cat('\nRMSE        (OLS)                 = ', sqrt(mse_ols))
  cat('\n ')
  cat('\nEfficiency (cond/exact)           = ', mse_cond/mse_exact)
  cat('\nEfficiency (ols/exact)            = ', mse_ols/mse_exact)
  
}



