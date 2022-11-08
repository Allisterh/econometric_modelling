#============================================================================
#    Kalman filter example to show recursive properties of the 
#   algorithm: univariate
#============================================================================

rm (list = ls(all=TRUE))
graphics.off()
set.seed(1234567, kind="Mersenne-Twister")

#
#--------------------------- Helper Functions -------------------------------
#

#load required functions - recserar
source("EMTSUtil.R")


#----------------------------------------------------------------------------
#  Log-likelihood procedure for the kalman filter        
#----------------------------------------------------------------------------
lnlt <- function(b,y) 
{
   # Unpack parameters
  lam  <- b[1]
  sig  <- b[2]^2
  phi  <- b[3]
  sigv <- b[4]^2  
  
  # Allocate arrays
  t   <- length(y)
  lnl <- rep(0, t)
  
  # Recursions of the Kalman Filter
  # Initialisation following Harvey and Hamilton
  st     <- 0.0
  pt     <- 1/(1-phi^2)   
  mt     <- lam*st
  vt     <- lam %*% pt %*% lam + sig
  ut     <- y[1] - mt
  
  lnl[1] <- - 0.5*log(2*pi) - 0.5*log(det(as.matrix(vt))) - 0.5*ut^2/vt
  kal    <- pt %*% (lam/vt)
  s0     <- st + kal %*% ut
  p0     <- pt - kal %*% lam %*% pt 
  
  # Table of results
  cnames <- c("iteration", "yt", "st", "pt", "mt", "vt", "ut", "kal", "lnlt")
  tblResults <- matrix(nrow=t, ncol=length(cnames), dimnames=list(rep("", t), cnames) )
  tblResults[1,] <- c(1, y[1], st, pt, mt, vt, ut, kal, lnl[1])
    
  # Main loop over observations
  for (i in 2:t) {
     # Predict and store predictions
    st      <- phi*s0
    pt     <- phi %*% p0 %*% phi + sigv
              
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
    tblResults[i,] <- c(i, y[i], st, pt, mt, vt, ut, kal, lnl[i])
  }  
  
  # Display the results
  cat('\n')  
  print(t(tblResults))
  cat('\n')  
}

    

#
#---------------------- Univariate Kalman Filter ----------------------------
#

lfac_uni <- function() {
  # Simulate the data
  t <- 5
  
  lam <- 1.0
  sig <- 0.5
  phi <- 0.8
  sige <- 1
  
  u <- sig*rnorm(t)
  v <- rnorm(t)
  
  f <- recserar(cbind(v),cbind(0.0),cbind(phi))  
  y <- lam*f + u
  
  # Call the log-likelihood for the Kalman filter 
  theta0 <- c(lam, sig, phi, sige)
  lnlt(theta0,y)  
}





