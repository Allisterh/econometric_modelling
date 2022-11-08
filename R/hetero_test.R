# ============================================================================
#
#    Testing a model of heteroskedasticity 
#
# ============================================================================

rm (list = ls(all=TRUE))
graphics.off()
set.seed(1234567, kind="Mersenne-Twister")

#
# -------------------------- Helper Functions --------------------------------
#
#load required functions - inv, numgrad
source("EMTSUtil.R")

#-----------------------------------------------------------------------------
#   Simulate the data  
#-----------------------------------------------------------------------------
simulatedata <- function(t) {
  beta0 <- 1
  beta1 <- 2
  gam0  <- 0.1
  gam1  <- 0.1

  x <- rnorm(t)                                            
  w <- seq(0.0,0.1*(t-1), 0.1)                                               

  u <- sqrt(exp(gam0 + gam1*w))*rnorm(t)    #  ut is N(0,sig2)                      
  y <- beta0 + beta1*x + u                   #  yt
  
  return(list(y=y, x=x, w=w))
}

#-----------------------------------------------------------------------------
# Negative unconstrained log-likelihood  
#-----------------------------------------------------------------------------
neglog1 <- function(theta,y,x,w) {
  lf <- -mean( lnlt1(theta,y,x,w) )
  return(lf)
}

#-----------------------------------------------------------------------------
# Unconstrained log-likelihood function at each observation
#-----------------------------------------------------------------------------
lnlt1 <- function(b,y,x,w) {
  mu  <- b[1] + b[2]*x
  sig <- sqrt( exp(b[3] + b[4]*w) )
  lnl <- -(1/2)*log(2*pi*sig^2) - (y - mu)^2 /(2*sig^2)
  return(lnl)
}

#-----------------------------------------------------------------------------
# Negative constrained log-likelihood  
#-----------------------------------------------------------------------------
neglog0 <- function(theta,y,x,w) {
  lf <- -mean( lnlt0(theta,y,x,w) )
  return(lf)
}

#-----------------------------------------------------------------------------
# Constrained log-likelihood function at each observation
#-----------------------------------------------------------------------------
lnlt0 <- function(b,y,x,w) {
  mu   <- b[1] + b[2]*x
  sig  <- sqrt( exp(b[3] + 0*w) )
  lnl  <- -(1/2)*log(2*pi*sig^2) - (y - mu)^2 /(2*sig^2)
  return(lnl)
}

#
#------------------------- Model Tests ---------------------------------
#

hetero_test <- function() {
  # Simulate the model   
  t       <- 500
  simResults <- simulatedata(t)
  x <- simResults$x
  y <- simResults$y
  w <- simResults$w

  # Estimate the unconstrained model  
  theta <- c(1, 2, 0.1, 0.1)
  estResults <- optim(theta, neglog1, y=y, x=x, w=w, method="BFGS", hessian=T)  
  theta1 <- estResults$par
  l1     <- estResults$value
  H1     <- estResults$hessian

  l1 <- -l1
  cat('\nLog-likelihood (unconstrained) = ', l1)
  cat('\nUnconstrained parameter estimates\n' )
  sterr <- (1/t)*sqrt( diag(inv(H1) ) )
  print( cbind(theta1, sterr) )
    
  # Estimate the constrained model      
  theta <- c(1, 2, 0.1)
  estResults <- optim(theta, neglog0, y=y, x=x, w=w, method="BFGS", hessian=T)  
  theta0 <- estResults$par
  l0     <- estResults$value
  H0     <- estResults$hessian

  l0 <- -l0
  cat('\nLog-likelihood (constrained) = ', l0)
  cat('\nConstrained parameter estimates\n' )
  sterr <- (1/t)*sqrt( diag(inv(H0) ) )
  print(cbind(theta0, sterr))

  # LR test   
  lr <- -2*t*(l0 - l1)
  cat('\nLR statistic                   = ',lr) 
  cat('\np-value                        = ',1-pchisq(lr,1))

  # Wald test   
  r <- rbind(c(0 , 0 , 0 , 1))
  q <- 0
  wd <- t* t( (r %*% theta1 - q) ) %*% inv(r %*% inv(H1) %*% t(r)) %*% (r %*% theta1 - q)
  cat('\nWald statistic                 = ',wd) 
  cat('\np-value                        = ',1-pchisq(wd,1))

  # LM test      
  theta <- c(theta0,  0)
  gmat  <- numgrad(lnlt1,theta,y,x,w)  
  g     <- cbind(colMeans(gmat))
  j    <- t(gmat) %*% gmat/t
  lm   <- t* t( g ) %*% inv(j) %*% g
  
  cat('\nGradient evaluated at contrained estimates\n')
  print(g)
  cat('\nOuter product of gradients matrix\n')
  print(j)     
  cat('\n\nLM statistic                 = ',lm) 
  cat('\np-value                      = ',1-pchisq(lm,1))

  # LM test (regression form)   
  x <- cbind(rep(1, t) , x)
  # Stage 1 regression
  b <- lm(y ~ x - 1)$coef
  u <- y - x %*% b    
  w <- cbind(rep(1, t) , w)
  v <- u^2
  # Stage 2 regression
  b  <- lm(v ~ w - 1)$coef
  e  <- v - w %*% b
  r2 <- 1 - sum(e^2)/sum( (v-mean(v))^2 )
  lm <- t*r2
  cat('\n\nLM statistic (regression)    = ',lm) 
  cat('\np-value                      = ',1-pchisq(lm,1))  
}   
    

