#===========================================================================
#
#    Estimate a one factor model of the term structure
#
#===========================================================================

rm (list = ls(all=TRUE))
graphics.off()


#
#--------------------------- Helper Functions -------------------------------
#

# Load required functions - inv
source("EMTSUtil.R")

#----------------------------------------------------------------------------
#   Wrapper function to set up and call the Kalman filter
#----------------------------------------------------------------------------
neglog <- function(b,y,flag) 
{
  # Unpack the parameter vector
    Lam <- b[1:6]
    if (flag)
      Phi <- tanh(b[13])
    else
      Phi <- b[13]
    R <- diag(b[7:12]^2)
    Q <- 1
    lf <- -mean( kalman(y,Phi,Lam,R,Q) )
    return(lf)
}

#--------------------------------------------------------------------------
#   Wrapper function to set up and call the Kalman filter
#--------------------------------------------------------------------------
neglogc <- function(b,y,flag) {
  # Unpack the parameter vector
  Lam <- rep(b[1], 6) 
  if (flag)
    Phi <- tanh(b[8])
  else
    Phi <- b[8]
  R <- diag(b[2:7]^2)
  Q <- 1
    
  lf <- -mean( kalman(y,Phi,Lam,R,Q) )
  return(lf)
}

#----------------------------------------------------------------------------
# Kalman filter
#----------------------------------------------------------------------------
kalman <- function(y,Phi,Lam,R,Q) {
  # Allocate arrays
  t <- nrow(y)
  n <- ncol(y)
  k         <- length(Q)  
  lnl       <- rep(0, t)
  
  # Recursions of the Kalman Filter
  # Initialisation following Harvey and Hamilton  
  st <- rep(0, k)
  pt <- diag(k)*0.1   

  mt <- Lam * st
  vt <- cbind(Lam) %*% pt %*% Lam + R
  ut <- y[1,] - mt
  
  lnl[1] <- - 0.5*n*log(2*pi) - 0.5*log(det(vt)) - 0.5*ut %*% inv(vt) %*% cbind(ut)
  Kal <- pt %*% Lam %*% inv(vt)
  s0 <- st + Kal %*% ut
  p0 <- pt - Kal %*% Lam %*% pt

  # Main loop over observations
  for (i in 2:t) {
    # Prediction 
    st <- Phi %*% s0     
    pt <- Phi %*% p0 %*% cbind(Phi) + Q          
              
    # Observation
    mt <- Lam*st
    vt <- cbind(Lam) %*% pt %*% Lam + R
    ut <- y[i,] - mt

    # Construct log-likelihood function
    lnl[i] <- - 0.5*n*log(2*pi) - 0.5*log(det(vt)) - 0.5*ut %*% inv(vt) %*% cbind(ut)
                
		# Update    			
    Kal <- pt %*% Lam %*% inv(vt)
    s0 <- st + Kal %*% ut
    p0 <- pt - Kal %*% Lam %*% pt    
  }
  return(lnl)
}

#----------------------------------------------------------------------------
#   Extract smoothed factor
#----------------------------------------------------------------------------
kfsmooth <- function(b,y) {
  # Unpack the parameter vector
  Lam <- b[1:6]
  Phi <- b[13]
  R   <- diag(b[7:12]^2)  
  Q   <- 1
  
  # Allocate arrays
  t <- nrow(y)
  n <- ncol(y)
  
  k       <- length(Q)
  s10     <- array(0, c(t,k))    # st|t-1
  s11     <- array(0, c(t,k))    # st|t  
  p10     <- array(0,c(k,k,t))
  p11     <- array(0, c(k,k,t))
  ss      <- rep(0, t)
     
  # Initialisation following Harvey and Hamilton  
  st <- rep(0,k)
  pt <- diag(k)*0.1  
  s10[1,] <- st 
  p10[,,1] <- pt
 
  mt <- Lam*st
  vt <- cbind(Lam) %*% pt %*% Lam + R
  ut <- y[1,] - mt
    
  Kal <- pt %*% Lam %*% inv(vt)
  s0 <- st + Kal %*% ut
  p0 <- pt - Kal %*% Lam %*% pt
  s11[1,] <- s0
  p11[,,1] <- p0
    
  # Main loop over observations
  for (i in 2:t) {
    # Prediction 
    st <- Phi %*% s0                     
    pt <- Phi %*% p0 %*% cbind(Phi) + Q      
    s10[i,] <- st 
    p10[,,i] <- pt
      
    # Observation
    mt <- Lam*st
    vt <- cbind(Lam) %*% pt %*% Lam + R
    ut <- y[i,] - mt
    
    # Update    			
    Kal <- pt %*% Lam %*% inv(vt)
    s0 <- st + Kal %*% ut
    p0 <- pt - Kal %*% Lam %*% pt
    s11[i,] <- s0
    p11[,,i] <- p0
  }
    
  # Now do the backwards recursion to smooth the factor    
  ss[length(ss)] <- s11[length(s11)]
  for (j in rev(seq(t-1))) {
    Jt      <- p11[,,j] %*% t(Phi) %*% inv(p10[,,j])
    ss[j] <- s11[j] + Jt %*% ( ss[j+1] - s10[j+1] )
  }       
  fac <- ss
  return(fac)
}

    
   


#
#--------------------------- Term Structure of Interest Rates ----------------
#
lfac_term <- function () {
    
    # Read the data, rearrange and scale data
    # US daily zero coupon yields starting 10/04/1988 and ending 12/28/2001
    #  	1.	tcm10y
    #		2.	tcm7y
    #		3.	tcm5y
    #		4.	tcm3y
    #		5.	tcm1y
    #		6.	tcm3m
  usdata <- as.matrix(read.table("lfac_usdata.dat"))
  rt <- usdata*100
  yt <- t( apply(rt, 1, '-', colMeans(rt)) ) 
  t <- nrow(yt)
  
  # Estimate the model by MLE  
  start <- c(7.353852219681254,
                7.636371949157256, 
                7.003563287217877, 
                6.293071514715054, 
                5.847811722064529, 
                5.536175149488106, 
                65.74336399439356, 
                43.66727769088171, 
                6.286986479990421, 
                26.59070018994360, 
                36.68045957160562, 
                50.46637407916106, 
                3.979446498229009)
  # Estimate model
  estResults <- optim(start, neglog, y=yt, flag=TRUE, method='BFGS')
  theta1 <- estResults$par
  fval1 <- -estResults$value
  
  theta1[length(theta1)] <- tanh(theta1[length(theta1)])
  
  # Extract and plot the smoothed factor
  fac <- kfsmooth(theta1,yt)
  figure()
  plot(seq(fac),fac, type="l", xlab="", ylab="")
    
  # Estimate restricted model
  start <- c(6.862810483385435,
                    66.30255435166845, 
                    45.85160110029478, 
                    6.589040449849120, 
                    29.53625724264360, 
                    42.59633524199530, 
                    57.62141009765852, 
                    3.988334426025099)
  estResults <- optim(start, neglogc, y=yt, flag=T, method="BFGS")
  theta0 <- estResults$par
  fval0 <- -estResults$val
  
  # Likelihood ratio test of restrictions
  lr  <- -2*t*(fval0 - fval1)
  dof <- length(theta1) - length(theta0)
  
  cat('\nLog-likelihood (unrestricted) = ',fval1)
  cat('\nLog-likelihood (restricted)   = ',fval0)
  cat('\nLR statistic                  = ',lr)
  cat('\np-value                       = ',1-pchisq(lr,dof))
}

