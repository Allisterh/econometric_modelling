#============================================================================
#
#   Program to to estimate ex ante real interest rates  
#   from ex post real interest rates
#
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()


#
#--------------------------- Helper Functions -------------------------------
#

# Load required functions - inv
source("EMTSUtil.R")

library(matlab)  # repmat, reshape 

#-------------------------------------------------------------------------
#  Log-likelihood function from Kalman filter
#-------------------------------------------------------------------------
neglog <- function( b,y,flag ) {
  # Define matrices for filter
  alpha <- b[1]      
  Lam   <- 1
  if (flag)
    Phi <- tanh(b[2])
  else
    Phi <- b[2]  
  R     <- b[3]^2
  Q     <- b[4]^2     
  
  lf    <- -mean( lnlt(y,Phi,Lam,R,Q,alpha) )   
  return(lf)
}

#
#--------------------------------------------------------------------------
# Multivariate Kalman filter
#--------------------------------------------------------------------------
lnlt <- function(y,Phi,Lam,R,Q,alpha) {
  # Allocate arrays  
#   t <- nrow(y)
#   n <- ncol(y)
  t <- length(y)
  n <- 1
  k         <- length(Q)
  lnl       <- rep(0, t)  
  
  # Recursions of the Kalman Filter
  # Initialisation following Harvey and Hamilton  
  st <- rep(0, k)

  pt <- reshape(inv(diag(k^2) - Phi^2)*Q, k, k)  
  
  mt <- Lam %*% st + alpha
  vt <- Lam %*% pt %*% t(Lam) + R
  ut <- y[1] - mt
  
  lnl[1] <- - 0.5*n*log(2*pi) - 0.5*log(det(vt)) - 0.5*t(ut) %*% inv(vt) %*% ut
  Kal <- pt %*% t(Lam) %*% inv(vt)
  s0 <- st + Kal %*% ut
  p0 <- pt - Kal %*% Lam %*% pt
  
  # Main loop over observations
  for (i in 2:t) {
    # Prediction 
    st <- Phi %*% s0 
    pt <- Phi %*% p0 %*% cbind(Phi) + Q          
    
    # Observation
    mt <- Lam %*% st + alpha
    vt <- Lam %*% pt %*% t(Lam) + R
    ut <- y[i] - mt
    
    # Construct log-likelihood function
    lnl[i] <- - 0.5*n*log(2*pi) - 0.5*log(det(vt)) - 0.5*t(ut) %*% inv(vt) %*% ut
    
    # Update          
    Kal <- pt %*% t(Lam) %*% inv(vt)
    s0 <- st + Kal %*% ut
    p0 <- pt - Kal %*% Lam %*% pt   
  } 
  return(lnl)
}


#--------------------------------------------------------------------------
#   Extract smoothed factor
#--------------------------------------------------------------------------
kfsmooth <- function(b,y, flag) {
  # Unpack the parameter vector
  alpha <- b[1]      
  Lam   <- 1
  if (flag)
    Phi <- tanh(b[2])
  else
    Phi <- b[2]  
  R     <- b[3]^2
  Q     <- b[4]^2     
  
  # Allocate arrays  
  #t <- nrow(y)
  #n <- ncol(y)
  n <- 1
  t <- length(y)
  k       <- length(Q)
  s10     <- array(0, c(t,k) )    # st|t-1
  s11     <- array(0, c(t,k) )    # st|t  
  p10     <- array(0, c(k,k,t) )
  p11     <- array(0, c(k,k,t) )
  ss      <- rep(0,t)
  
  # Initialisation following Harvey and Hamilton  
  st <- rep(0, k)
  pt <- diag(k)*0.1  
  s10[1] <- st 
  p10[1] <- pt
  
  mt <- Lam %*% st + alpha
  vt <- Lam %*% pt %*% t(Lam) + R
  ut <- y[1] - mt
    
  Kal <- pt %*% t(Lam) %*% inv(vt)
  s0 <- st + Kal %*% ut
  p0 <- pt - Kal %*% Lam %*% pt  

  s11[1] <- s0
  p11[1] <- p0
  
  
  # Main loop over observations  
  for (i in 2:t) {
    # Prediction 
    st <- Phi %*% s0 
    pt <- Phi %*% p0 %*% cbind(Phi) + Q   
    s10[i] <- st 
    p10[,,i] <- pt
    
    # Observation
    mt <- Lam %*% st + alpha
    vt <- Lam %*% pt %*% t(Lam) + R
    ut <- y[i] - mt
    
    # Update          
    Kal <- pt %*% t(Lam) %*% inv(vt)
    s0 <- st + Kal %*% ut
    p0 <- pt - Kal %*% Lam %*% pt
    s11[i] <- s0
    p11[i] <- p0
  }
    
  # Now smooth the factor    
  ss <- s11
  for (j in 1:(t-1)) {
    Jt      <- p11[t-j] %*% t(Phi) %*% inv(p10[t-j])     
    ss[t-j] <- t(s11[t-j]) + Jt %*% ( t(ss[(t-j+1)]) - t(s10[(t-j+1)]) ) 
  }  
  fac <- ss 
  return(list (fac=fac, s11=s11))
}





#
#------------------------- Ex Ante Real Interest Rates ------------------------
#
lfac_exante <- function() {
  # Load data starting Jan 1971 and ending December 2009
  load('exante.Rdata')
  
  interest <- data[,1]
  price <- data[,2]
  
  # Compute ex post real interest rate    
  inflation <- 1200*(trimr(log(price),1,0) - trimr(log(price),0,1))           
  y <- trimr(interest,1,0) - inflation                                
  
  t <- length(y)
  
  # Estimate parameters with restriction imposed
  flag  <- 1
  start <- c(mean(y),0.5,  1, sd(y))
  estResults <- optim(start, neglog, y=y, flag=flag, method="BFGS")
  
  bhat <- estResults$par
  
  # Restimate without restrictions to get standard errors
  flag  <- 0
  start <- c(bhat[1], tanh(bhat[2]),  bhat[3:4])
  estResults <- optim(start, neglog, y=y, flag=flag, method="BFGS", hessian=T)
  bhat <- estResults$par
  lf <- estResults$val
  hess <- estResults$hess
  
  cat('\nLog-likelihood function = ',lf)
  cat('\n')
  cat('\n Mean of ex post real interest rate mean     = ',mean(y))
  cat('\n Variance of ex post real interest rate mean = ',sd(y)^2)
  
  cat('\n' )   
  cat('\n Mean of ex ante real interest rate mean     = ',bhat[1])
  cat('\n Variance of ex ante real interest rate mean = ',bhat[4]^2/(1 - bhat[2]^2)) 
  
  # Real ex ante interest rate
  s_update <- kfsmooth(bhat,y, flag)$s11
  
  # Plot leaves off addition of bhat(1) so you can see two series
  exante <- s_update
  figure()
  
  matplot(seqa(1971+2/12,1/12,t),cbind(exante, y), type="l",
          main = 'Ex Ante Interest Rate',
          xlab = '',
          ylab = '')
}



