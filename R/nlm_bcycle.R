#============================================================================
#
#   Estimate Markov Switching model (Hamilton, Econometrica, 1989, 357-384)
#
#============================================================================

rm (list = ls(all=TRUE))
graphics.off()
#
#--------------------------- Helper Functions -------------------------------
#
# load required functions - trimr, inv
source("EMTSUtil.R")

#----------------------------------------------------------------------------
#  Restricted log-likelihood function
#----------------------------------------------------------------------------
neglogr <- function(b,y,flag) {
  # Parameters   
  alpha  <- b[1]
  beta   <- b[2]
  gam    <- b[3]
  delta  <- b[4]
  if (flag) {
    p     <- 1/(1 + exp(-b[5]))  # Ensure p,q [0,1]         
    q     <- 1/(1 + exp(-b[6]))
  }else {
    p <- b[5]
    q <- b[6]
  }   
  # Mean and variance vectors for states 1 and 0    
  m  <- rbind(alpha + beta,alpha)
  s2 <- rbind(gam   + delta, gam)
  
  # Standardized variables for states 1 and 0
  z1 <-  (y - m[1])/sqrt(s2[1])      
  z0 <-  (y - m[2])/sqrt(s2[2])      
  
  
  # Distributions for states 1 and 0
  f1 <- dnorm(z1)/sqrt(s2[1])         
  f0 <- dnorm(z0)/sqrt(s2[2])     
  
  # Define transitional and stationary probabilities    
  p_stat  <- rbind( (1 - q)/(2 - p - q),
                    (1 - p)/(2 - p - q) )
  w <- p_stat                        
  
  
  # Construct likelihood  
  t <- length(y)
  f <- rep(0, t)
  
  for (i in seq(t)) {
    # Conditional distribution of y given lagged y 
    f[i] <- w[1]*f1[i] + w[2]*f0[i]  
    
    if (flag) {
      p <- 1/(1 + exp(-b[5]))           
      q <- 1/(1 + exp(-b[6]))
    }else {
      p <- b[5]
      q <- b[6]
    }
    p_trans <- matrix(c(p, 1 - q,
                        1 - p,     q), nrow=2, byrow=T)
    
    # Update weights of being in each state  
    tmp <- rbind(w[1]*f1[i]/f[i], w[2]*f0[i]/f[i])
    w   <- p_trans %*% tmp  
  }
  lf <- -mean(log(f) )
  return(lf)
}

#-------------------------------------------------------------------------
#  Unrestricted log-likelihood function
#-------------------------------------------------------------------------
neglog <- function(b,y,urate) {
  # Parameters   
  alpha  <- b[1]
  beta   <- b[2]
  gam    <- b[3]
  delta  <- b[4]
  p      <- 1/(1 + exp(-b[5]))       
  q      <- 1/(1 + exp(-b[6]))  
  
  kappa1 <- b[7]   #  cond. prob of being in state 1  
  lam1   <- b[8]   #  cond. prob of being in state 0  
  
  # Mean and variance vectors for states 1 and 0    
  m  <- rbind(alpha + beta, alpha)
  s2 <- rbind(gam   + delta, gam)
  
  # Standardized variables for states 1 and 0
  z1 <-  (y - m[1])/sqrt(s2[1])      
  z0 <-  (y - m[2])/sqrt(s2[2])      
  
  # Distributions for states 1 and 0
  f1 <- dnorm(z1)/sqrt(s2[1])         
  f0 <- dnorm(z0)/sqrt(s2[2])     
  
  # Define transitional and stationary probabilities    
  p_stat  <- rbind((1 - q)/(2 - p - q),     
                   (1 - p)/(2 - p - q))
  w       <- p_stat                        
  
  # Construct likelihood  
  t <- length(y)
  f <- rep(0, t)
  
  for (i in seq(t)) {
    # Conditional distribution of y given lagged y 
    f[i] <- w[1]*f1[i] + w[2]*f0[i]                                               
    p    <- 1/(1 + exp(-b[5]-kappa1*urate[i]))           
    q    <- 1/(1 + exp(-b[6]-lam1*urate[i]))
    
    p_trans <- matrix(c(p,  1 - q,
                        1 - p,    q), nrow=2, byrow=T)
    
    # Update weights of being in each state	
    tmp <- rbind(w[1]*f1[i]/f[i], w[2]*f0[i]/f[i])
    w   <- p_trans %*% tmp  
  }
  lf <- -mean(log(f) )
  return(lf)
}

#
#--------------------------- Business Cycle Model ---------------------------
#
nlm_bcycle <- function( )
{
  # Load data 
  load('nlm_GNP.Rdata')
  gnp <- data[, 'gnp']
  urate <- data[, 'urate']
  
  # Compute percentage growth rate (June 1951 to December 1984)   
  y <- 100*( trimr(log(gnp),1,0) - trimr(log(gnp),0,1) )   
  
  t <- length(y)
  #**********************************************************************
  #***
  #***     Generate graph
  #***
  #**********************************************************************
    
  xData = seqa(1951+2/4,1/4,t)
  figure()
  
  matplot(xData,cbind(y,1.176*rep(1, t), -0.224*rep(1, t)), type="l",
          main = "",
          ylab = expression(y[t]),
          xlab = expression(t),
          bty = "l")  
  
  # Estimate the model without time-varying parameters
  alpha  <- 0.0
  beta   <- 0.0
  gam    <- 0.0
  delta  <- 0.0
  p      <- 0.0
  q      <- 0.0
  # Parameters to control time-variation in the conditional probabilities    
  kappa1 <- 0.0       
  lam1   <- 0.0
  
  start <- c(-0.5, 1, 1.5, 1, 2, 2)
  estResults <- optim(start, neglogr, y=y, flag=T, method="BFGS")
  thetac <- estResults$par
  lf0 <- estResults$val
  
  lf0<- -lf0
  
  # Estimate model again to get standard errros: these are based on the HESSIAN
  start <- thetac   
  start[5] <- 1/(1 + exp(-thetac[5]))          
  start[6] <- 1/(1 + exp(-thetac[6]))
  
  estResults <- optim(start, neglogr, y=y, flag=F, method="BFGS", hessian=T)
  theta <- estResults$par
  hess <- estResults$hess
  
  vc <- (1/t)*inv(hess)
  se <- sqrt(diag(vc))
  
  cat('\nMarkov Switching Model Parameters\n')
  row.names  <- c('alpha','beta','gamma','delta','p','q')
  
  print(matrix(c(theta, se), 6,2, dimnames=list(row.names, c('theta', 'se'))))
  
  cat('\n')
  cat('\nDuration estimate (state=1)     = ',1/(1 - theta[5]))
  cat('\nStandard error                  = ',sqrt(vc[5,5])/(1 - p)^2)
  
  cat('\n ')
  cat('\nDuration estimate (state=0)     = ',1/(1 - theta[6]))
  cat('\nStandard error                  = ',sqrt(vc[6,6])/(1 - p)^2)
  
  # Wald test (delta = 0)
  w <- theta[4]^2/vc[4,4]
  cat('\n')
  cat('\nWald statistic (delta=0)        = ',w)
  cat('\np-value                         = ',1-pchisq(w,1))
  
  # Estimate the model with time-varying probabilities
  start <- c(thetac, 0, 0)
  estResults <- optim(start, neglog, y=y, urate=urate, method="BFGS")
  theta1 <- estResults$par
  lf1 <- estResults$val
  
  lf1 <- -lf1
  theta1[5] <- 1/(1 + exp(-theta1[5]))          
  theta1[6] <- 1/(1 + exp(-theta1[6]))
  
  cat('\nTime-varying Markov Switching Model Parameters\n')
  row.names  <- c('alpha', 'beta','gamma','delta', 'p',  'q', 'kappa', 'lambda')
  print(matrix(theta1, 8, 1, dimnames=list(row.names, c('theta'))))
  
  # LR test  
  lr <- -2*t*(lf0 - lf1)
  cat('\n')
  cat('\nLR statistic (kappa/lambda=0)   = ',lr)
  cat('\np-value                         = ',1-pchisq(lr,2))
}



