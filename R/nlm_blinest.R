#============================================================================
#
#   Estimate and test a bilinear time series model by maximum likelihood
#
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()
set.seed(123456, kind="Mersenne-Twister")
#
#--------------------------- Helper Functions -------------------------------
#
# load required functions - inv
source("EMTSUtil.R") 

#----------------------------------------------------------------------------
#  Unrestricted log-likelihood function
#----------------------------------------------------------------------------
neglog <- function(b,y) 
{
  nobs <- length(y)     
  v    <- rep(0, nobs)
  
  for (t in 2:nobs) {
    v[t] <- y[t] - b[1] - b[2]*y[t-1] - b[3]*y[t-1]*v[t-1]
  }  
  v <- trimr(v,1,0)     
  s2 <- mean(v^2) 
  z  <- v/sqrt(s2)
  
  lf <- -mean( -0.5*log(2*pi)-0.5*log(s2)-0.5*z^2 )
  return(lf)
}

#----------------------------------------------------------------------------
#  Restricted log-likelihood function
#----------------------------------------------------------------------------
neglogr <- function(b,y) {
  nobs <- length(y)     
  v    <- rep(0, nobs)
  
  for (t in 2:nobs) {
    v[t] <- y[t] - b[1] - b[2]*y[t-1]
  } 
  
  v <- trimr(v,1,0)   	
  s2 <- mean(v^2) 
  z  <- v/sqrt(s2)
  
  lf <- -mean( -0.5*log(2*pi)-0.5*log(s2)-0.5*z^2 )
  return(lf)
}


#
#------------------------------- Bilinearity  -------------------------------
#
nlm_blinest <- function( ) {
  # Simulate data
  t <- 1000
  y <- rep(0, t+100)
  v <- sqrt(0.1)*rnorm(t+100)
  
  for (k in 2:(t+100)) {
    y[k] <- 0.1 + 0.4*y[k-1] + 0.4*y[k-1]*v[k-1] + v[k]
  }
  
  y <- trimr(y,100,0)
  
  #  Estimate the unrestricted model 
  theta_0 <- c(0.1,  0.1,  0.1)
  estResults <- optim(theta_0, neglog, y=y, method="BFGS", hessian=T)
  theta1 <- estResults$par
  lf1 <- estResults$val
  hess <- estResults$hessian
  
  vc  <- (1/t)*inv(hess)
  se  <- sqrt(diag(vc))
  lf1 <- -lf1
  
  cat('\nUnrestricted log-likelihood     = ', lf1)
  cat('\n')
  print(cbind(Estimates=theta1, StdErrors= se))
  cat('\n')
  
  # Restricted model
  theta_0 <- c(0.1,  0.1,  0)
  
  estResults <- optim(theta_0, neglogr, y=y, method="BFGS")
  lf0 <- estResults$val
  
  lf0 <- -lf0
  
  # LR test       
  lr <- -2*t*(lf0 - lf1)
  
  cat('\nLikelihood ratio test  =  ', lr)
  cat('\np-value                =  ', 1 - pchisq(lr, 1))
  cat('\n')
  
  # Wald test       
  w <- (theta1[3] - 0)^2/vc[3,3]
  
  cat('\nWald test              =  ', w)
  cat('\np-value                =  ', 1 - pchisq(w, 1))
  cat('\n')
  
  # LM test       
  # First stage regression      
  x <- cbind(rep(1, length(y)-1), trimr(y,0,1))
  y <- trimr(y,1,0)           
  k <- ncol(x)
  v <- lm(y ~ x - 1)$residuals                            
  
  # Second stage regression      
  x <- cbind(trimr(x,1,0),  trimr(x[,2],1,0)*trimr(v,0,1))              
  v <- trimr(v,1,0)
  e <- lm(v ~ x - 1)$residuals
  
  t  <- length(y)
  r2 <- 1 - sum(e^2)/sum( (v-mean(v))^2 )
  lm <- t*r2
  
  cat('\nLM test              =  ', lm)
  cat('\np-value                =  ', 1 - pchisq(lm, 1))  
}

