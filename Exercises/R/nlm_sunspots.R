#============================================================================
#
#   Sunspot data
#
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()
#
#--------------------------- Helper Functions -------------------------------
#
# load required functions - inv, trimr
source("EMTSUtil.R")


#----------------------------------------------------------------------------
#  Compute the LM statistic to test a threshold autoregressive model 
#  assuming one lag in the auxiliary model
#----------------------------------------------------------------------------
tar_test <- function(yvar,p) {
  # First stage regression      
  y <- trimr(yvar,1,0)
  x <- cbind(rep(1, length(y)) , trimr(yvar,0,1))
  k <- ncol(x)
  u <- lm(y ~ x - 1)$residuals
  
  # Second stage regression      
  if (p == 1)
  x <- cbind(x,  x[,2]*(x[,2])) 
  else if (p == 2)
    x <- cbind(x,  x[,2]*(x[,2])^1,  x[,2]*(x[,2]^2)) 
  else if (p == 3)
    x <- cbind(x, x[,2]*(x[,2]^1) , x[,2]*x[,2]^2 , x[,2]*x[,2]^3)
  else if (p == 4)
    x <- cbind(x, x[,2]*(x[,2]^1), x[,2]*(x[,2]^3))       
  e <- lm(u ~ x - 1)$residuals
  
  # Compute LM statistic        
  r2  <- 1 - sum(e^2)/sum( (u-mean(u))^2 )
  lm  <- length(y)*r2
  dof <- ncol(x) - k
  return(list(lm=lm, dof=dof))
}

#----------------------------------------------------------------------------
#  Negative log-likelihood function
#----------------------------------------------------------------------------
neglog <- function(b,y,d,gam) {
  c  <- b[12]
  w  <-  ( 1 + exp(-gam*(trimr(y,6-d,d) - c)))^(-1)
  
  # List of variables in regime 1       
  x1 <- cbind(rep(1, length(y)-6), trimr(y,5,1),   trimr(y,4,2),   trimr(y,3,3),   trimr(y,2,4),   trimr(y,1,5),   trimr(y,0,6))   
  
  # List of variables in regime 2      
  x2 <- cbind(rep(1, length(y)-6),   trimr(y,5,1),   trimr(y,4,2),   trimr(y,3,3))
  
  u  <- trimr(y,6,0) - x1 %*% b[1:7] - x2 %*% b[8:11]*w
  
  sig2 <- as.numeric(t(u) %*% u/length(u))
     lnl  <- - 0.5*log(2*pi) - 0.5*log(sig2) - 0.5*u^2/sig2                                                 
  lf   <- -mean( lnl )
  return(lf)
}

#
#--------------------------- Sunspots  --------------------------------------
#
nlm_sunspots <- function() {
  sunspots <- as.matrix(read.table("sunspots.dat"))
  
  y <- sunspots[,1]
  t <- length(y)
  
  figure()  
  plot(y, type="l",
       main = "",
       ylab = 'Average Monthly Sunspots',
       xlab = "",
       bty = "l")
  
  
  # Perform LST linearity test
  results <- tar_test(y,1)
  lm <- results$lm
  dof <- results$dof
  cat('\nLM statistic (test 1) = ', lm)
  cat('\nDegrees of freedom    = ', dof)
  cat('\np-value               = ', 1-pchisq(lm,dof))
  cat('\n ')
  
  results <- tar_test(y,2)
  lm <- results$lm
  dof <- results$dof        
  cat('\nLM statistic (test 2) = ', lm)
  cat('\nDegrees of freedom    = ', dof)
  cat('\np-value               = ', 1-pchisq(lm,dof))
  cat('\n ')
  
  results <- tar_test(y,3)
  lm <- results$lm
  dof <- results$dof        
  cat('\nLM statistic (test 3) = ', lm)
  cat('\nDegrees of freedom    = ', dof)
  cat('\np-value               = ', 1-pchisq(lm,dof))
  cat('\n ')
  
  results <- tar_test(y,4)
  lm <- results$lm
  dof <- results$dof         
  cat('\nLM statistic (test 4) = ', lm)
  cat('\nDegrees of freedom    = ', dof)
  cat('\np-value               = ', 1-pchisq(lm,dof))
  cat('\n ')
  
  
  # Parameters
  d   <- 2        # Delay parameter     
  gam <- 100        # Adjustment parameter   1, 5, 10, 50, 100       
  
  # Estimate the model
  theta_0 <- c(rnorm(11), 50)
  estResults <- optim(theta_0, neglog, y=y, d=d, gam=gam, method="BFGS", hessian=T)
  theta <- estResults$par
  lnl <- estResults$val
  hess <- estResults$hessian
  
  vc <- (1/t)*inv(hess)
  cat('\n')
  cat('\nLog-likelihood    = ',-lnl, '\n')
  print( cbind(Parameters=theta, 'Std. Errors'=sqrt(diag(vc)) ))  
}





