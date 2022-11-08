#=======================================================================
#
#    Simulation example to generate the sampling distributions the MLE 
#    and OLS estimator of a regression model with heteroskedasticity.
#
#=======================================================================

rm (list = ls(all=TRUE))
graphics.off()
set.seed(1234567, kind="Mersenne-Twister")

#
#------------------------- Helper Functions-----------------------------
#
#load require lib - repmat
library("matlab")
#-----------------------------------------------------------------------
#      Negative log-likelihood function   
#-----------------------------------------------------------------------
neglog <- function(theta,y,x ) {
   lf <- -mean( lnlt(theta,y,x) )
   return(lf)  
}

#-----------------------------------------------------------------------
#      Log-likelihood function at each observation
#-----------------------------------------------------------------------
lnlt <- function(p,y,x) {
  mu   <- p[1] + p[2]*x
  sig2 <- exp(p[3] + p[4]*x)
  lft  <- -(0.5)*log(2*pi*sig2) - (y - mu)^2/(2*sig2)
  return(lft)
}


#
#----------------------- Sampling Properites ---------------------------
#

hetero_sampling <- function() {
  # Simulate a regression model with multiplicative heteroskedasticity
  t      <- 500                                                                                      
  beta0  <- 1.0
  beta1  <- 2.0
  gam0   <- 1.0
  gam1   <- 5.0
  theta0 <- c(beta0, beta1, gam0, gam1)

  # Fixed exogenous variable
  x <- runif(t)       
                     
  # Monte Carlo replications
  ndraws   <- 10000   
  zeros    <- array(0, c(ndraws,4) )
  theta_mle <- zeros
  theta_ols <- zeros
    
  pb <- txtProgressBar(min=0, max=ndraws, style=3)
  for (k in seq(ndraws)) {
    # Generate data       
    u <- sqrt( exp(gam0 + gam1*x) )*rnorm(t)                         
    y <- beta0 + beta1*x + u                                                          
        
    # ML Estimation
    estResults <- optim(theta0, neglog, y=y, x=x, method="BFGS")
    theta_hat  <- estResults$par     
    theta_mle[k,] <- theta_hat
        
    # OLS
    xx             <- cbind(ones(t,1),x)
    b_ols          <- lm(y ~ xx - 1)$coef
    e              <- y - cbind(ones(t,1),x) %*% b_ols
    gam0_ols       <- log(mean(e^2))
    theta_ols[k,] <- c(b_ols , gam0_ols , 0)
    setTxtProgressBar(pb, k)
  }
  close(pb)
        
  
  # Compute statistics of sampling distributions and print results
  rmse_mle   <- sqrt( colMeans( (theta_mle - repmat(theta0,ndraws,1))^2 ) )
  rmse_ols   <- sqrt( colMeans( (theta_ols - repmat(theta0,ndraws,1))^2 ) )

  cat('\n\nML results\n')
  cnames <- c("True", "Estimate", "Bias", "RMSE")
  results <- matrix(c(theta0, colMeans(theta_mle), colMeans(theta_mle)-theta0, rmse_mle), 
                nrow=4, dimnames = list(rep("", 4), cnames))
  print(results)
  
  cat('\n\n OLS results\n')
  results <- matrix(c(theta0, colMeans(theta_ols), colMeans(theta_ols)-theta0, rmse_ols),
         nrow=4, dimnames=list(rep("", 4), cnames))
  print(results)
    
  cat('\n\nEfficiency (OLS/ML)\n')
  print(unname(cbind((rmse_ols/rmse_mle))))
}
