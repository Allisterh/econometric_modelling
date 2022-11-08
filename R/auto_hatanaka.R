#================================================================================
#  
#  Simulation example to compure the asymptotic distribution of the MLE
#  estimator and the Hatanaka estimator for the regression model with 
#  autocorrelation.
#
#================================================================================

rm (list = ls(all=TRUE))
graphics.off()
set.seed(1, kind="Mersenne-Twister")

#
#--------------------------- Helper Functions -----------------------------------
#

#load required functions - inv, recserar
source("EMTSUtil.R")
library("matlab")

#--------------------------------------------------------------------------------
#   Log-likelihood function at each observation            
#--------------------------------------------------------------------------------
neglog <- function(b,y,x) {
  beta0  <- b[1]                        
  beta1  <- b[2]
  alpha1 <- tanh(b[3])     # rho1 stays in the unit circle    
  rho1   <- tanh(b[4])     # rho2 stays in the unit circle    
  sig2   <- abs(b[5])      # variance is positive     
  u      <- trimr(y,1,0) - beta0 - beta1*trimr(x,1,0) - alpha1*trimr(y,0,1)
  v      <- trimr(u,1,0) - rho1*trimr(u,0,1)
  lnl    <- - 0.5*log(2*pi) - 0.5*log(sig2) - 0.5*v^2/sig2
  lf     <- -mean( lnl )
  return(lf)
}

#--------------------------------------------------------------------------------
# Hatanaka 2-step efficient estimator      
#--------------------------------------------------------------------------------
hatanaka <- function(y,x,t) {
  yvar <- trimr(y,1,0)
  xvar <- cbind(ones(t-1,1), trimr(x,1,0), trimr(y,0,1))
  zvar <- cbind(ones(t-1,1), trimr(x,1,0), trimr(x,0,1))
 
    
  #   IV initial estimates of mean parameters
  biv  <- inv(t(zvar) %*% xvar) %*% (t(zvar) %*% yvar)         
  u    <- yvar - (xvar %*% biv)
    
  # Estimate rho
  rho  <- lm(trimr(u,1,0) ~ trimr(u,0,1) - 1)$coef    
  yvar <- trimr(y,2,0) - rho*trimr(y,1,1)             
  xvar <- cbind(ones(t-2,1), (trimr(x,2,0) - rho*trimr(x,1,1)), (trimr(y,1,1) - rho*trimr(y,0,2)), trimr(u,0,1))
   
  
  # Regression on transformed variables and update rho
  b    <- lm(yvar ~ xvar - 1)$coef
  rho  <- rho + b[4]                           
  #    Compute residual variance 
  v    <- trimr(y,2,0)-cbind(cbind(ones(t-2,1), trimr(x,2,0), trimr(y,1,1)) %*% b[1:3], -rho*(trimr(y,1,1)-cbind(ones(t-2,1), trimr(x,1,1), trimr(y,0,2)) %*% b[1:3]))
  sig2 <- mean(v^2)                                    
  h1   <- c(b[1:3], rho, sig2)  
  return(h1)
}

#
#-------------------------- Hatanaka Estimator ----------------------------------
#
            
auto_hatanaka <- function(){
  t      <- 1000000
  beta0  <- 1.0
  beta1  <- 1.0
  alpha1 <- 0.5
  rho1   <- 0.6
  sigma2 <- 0.1
  theta0 <- c(beta0, beta1, alpha1, rho1, sigma2)
    
  x <- runif(t) - 0.5  
  v <- sqrt(sigma2)*rnorm(t)

  u <- recserar(cbind(v), cbind(sqrt(1/(1-rho1^2))*v[1]), cbind(rho1))  # Beach and MacKinnon, fn3     
  y <- recserar(cbind(beta0 + beta1*x + u) , cbind(0.0), cbind(alpha1))                                                                      
    
  cov_theoretical <- rho1*sigma2/((1-alpha1*rho1)*(1-rho1^2))
  cov_simulated   <- cov(trimr(y,0,1), trimr(u,1,0)) 

  cat('\nTheoretical and simulated covariance\n')
  print( cbind(cov_theoretical, cov_simulated) )

  #   Monte Carlo replications  
  t      <- 1000
  ndraws <- 1000  
  x      <- runif(t) - 0.5
  theta_cond     <- zeros(ndraws,5)
  theta_hatanaka <- zeros(ndraws,5)
  theta_ols      <- zeros(ndraws,5)

  pb <- txtProgressBar(min=0, max=ndraws, style=3)
  for (k in seq(ndraws)) {
    #   Generate data       
    v <- sqrt(sigma2)*rnorm(t)
    u <- recserar( cbind(v) , cbind(sqrt(1/(1-rho1^2))*v[1]) , cbind(rho1) )   
    y <- recserar( cbind(beta0 + beta1*x + u) , cbind(0.0), cbind(alpha1) )    
        
    #   Conditional MLE
    estResults <- optim(theta0, neglog, y=y, x=x, method="Nelder-Mead")
    theta <- estResults$par
    
    theta_cond[k,] <- c( theta[1], theta[2], tanh(theta[3]), tanh(theta[4]),  abs(theta[5]) )                         

    #   Hatanaka 2-step efficient estimator       
    theta <- hatanaka(y,x,t)
    theta_hatanaka[k,] <- theta        

    #   OLS         
    b_ols <- lm(trimr(y,1,0) ~ cbind(ones(t-1,1), trimr(x,1,0), trimr(y,0,1)) - 1)$coef        
    u_hat <- trimr(y,1,0) - cbind(ones(t-1,1), trimr(x,1,0), trimr(y,0,1)) %*% b_ols
    rho_ols <- lm(trimr(u_hat,1,0) ~ trimr(u_hat,0,1) - 1)$coef
    sig2_ols <- mean(u_hat^2)
    theta_ols[k,] <- c(b_ols, rho_ols, sig2_ols)
        
    setTxtProgressBar(pb, k)      
  }
  close(pb)

  # Compute statistics of sampling distributions and print results     
  mse_cond     <- colMeans((theta_cond - repmat(theta0,ndraws,1))^2)
  mse_hatanaka <- colMeans((theta_hatanaka - repmat(theta0,ndraws,1))^2)
  mse_ols      <- colMeans((theta_ols - repmat(theta0,ndraws,1))^2)

  cat('\n                                        Beta0   Beta1      Alpha1    Rho1    Sigma^2')
  cat('\nPopulation parameter              =      ', theta0[1], '     ', theta0[2], '       ', theta0[3], '     ', theta0[4],  '    ', theta0[5])
  cat('\n----------------------------------------------------------------------------------------------------')
  cat('\nMean        (cond. MLE)           = ', colMeans(theta_cond))
  cat('\nBias (x100) (cond. MLE)           = ', 100*(colMeans(theta_cond)-theta0))
  cat('\nMSE  (x100) (cond. MLE)           = ', 100*mse_cond)
  cat('\nRMSE (x100) (cond. MLE)           = ', 100*sqrt(mse_cond))

  cat('\n ')    
  cat('\nMean        (Hatanaka)            = ', colMeans(theta_hatanaka))
  cat('\nBias (x100) (Hatanaka)            = ', 100*(colMeans(theta_hatanaka)-theta0))
  cat('\nMSE  (x100) (Hatanaka)            = ', 100*mse_hatanaka)
  cat('\nRMSE (x100) (Hatanaka)            = ', 100*sqrt(mse_hatanaka))

  cat('\n ')
  cat('\nMean        (OLS)                 = ', colMeans(theta_ols))
  cat('\nBias (x100) (OLS)                 = ', 100*(colMeans(theta_ols)-theta0))
  cat('\nMSE  (x100) (OLS)                 = ', 100*mse_ols)
  cat('\nRMSE (x100) (OLS)                 = ', 100*sqrt(mse_ols))
  cat('\n ')
  cat('\nEfficiency (hatanaka/cond)        = ', mse_hatanaka/mse_cond)
  cat('\nEfficiency (ols/cond)             = ', mse_ols/mse_cond)            
}


