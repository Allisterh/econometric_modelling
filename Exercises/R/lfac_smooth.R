#============================================================================
#
#  	Kalman filter example to compute the smooth factor. 
#       
#
#============================================================================

rm (list = ls(all=TRUE))
graphics.off()
set.seed(123457, kind="Mersenne-Twister")

#
#--------------------------- Helper Functions -------------------------------
#

#load required functions - recserar, figure
source("EMTSUtil.R")


#----------------------------------------------------------------------------
# Univariate Kalman filter
#----------------------------------------------------------------------------
lnlt <- function (b,y) {
  # Unpack parameters
  lam  <- b[1]
  sig  <- b[2]^2
  phi  <- b[3]
  sigv <- b[4]^2    
    
        
  # Allocate arrays
  t   <- length(y)
  zeros <- rep(0, t)
  lnl <- zeros
  spr <- zeros     #     Predicted factor estimates     
  sup <- zeros     #     Updated factor estimates       
  ssm <- zeros     #     Smoothed factor estimates      
  ppr <- zeros     #     Predicted factor cov  
  pup <- zeros     #     Updated factor cov 
  psm <- zeros      #     Smoothed factor cov 
  

  # Recursions of the Kalman Filter
  # Initialisation following Harvey and Hamilton
  st     <- 0.0
  pt     <- 1/(1-phi^2)   
  spr[1] <- st
  ppr[1] <- pt   
  mt     <- lam*st
  vt     <- lam*pt*lam + sig
  ut     <- y[1] - mt
  lnl[1] <- -0.5*log(2*pi) - 0.5*log(det(as.matrix(vt))) - 0.5*ut^2/vt
  kal    <- pt %*% (lam/vt)
  s0     <- st + kal %*% ut
  p0     <- pt - kal %*% lam %*% pt
  sup[1] <- s0
  pup[1] <- p0
       
    
  # Main loop over observations
  for (i in 2:t) {
     # Predict and store predictions
    st      <- phi*s0
    pt     <- phi %*% p0 %*% phi + sigv
    spr[i] <- st
    ppr[i] <- pt
              
    # Observation     
    mt <- lam*st
    vt <- cbind(lam) %*% pt %*% lam + sig
    ut <- y[i] - mt   
    
    
    # Log-likelihood function
    lnl[i] <- -0.5*log(2*pi) - 0.5*log(det(as.matrix(vt))) - 0.5*ut^2/vt
              
    # Update and store updates
    kal <- pt %*% lam/vt
    s0 <- st + kal %*% ut
    p0 <- pt - kal %*% lam %*% pt
    sup[i] <- s0
    pup[i] <- p0
  }
  
  # Backward recursion to smooth and extract the factor
  ssm[t] <- sup[t]
  for (i in rev(seq(t-1))) {
    j      <- (pup[i]*phi)/ppr[i+1]
    ssm[i] <- sup[i] + j*(ssm[i+1] - spr[i+1])*j    
  }
  lfac <- ssm
  return(lfac)
}

#
#---------------------- Kalman Filter (Smooth factor) -----------------------
#

lfac_smooth <- function() 
{
  # Simulate the data       
  t    <- 5
  lam  <- 1.0
  sig  <- 0.5
  phi  <- 0.8
  sigv <- 1

  u <- sig*rnorm(t)
  v <- rnorm(t)
 
  f <- recserar(cbind(v),cbind(0.0), cbind(phi))
  y <- lam*f + u

  # Return the smoothed latent factor     
  theta0 <- c(lam,  sig,  phi,  sigv)
  lfac <- lnlt(theta0,y)
  
  figure()
  par(xaxs="i", yaxs="i")
  matplot(cbind(y, lfac), type="l",
          xlab = "",
          ylab = "",
          bty="l")
}

  

