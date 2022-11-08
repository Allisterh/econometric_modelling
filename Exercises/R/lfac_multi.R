#=========================================================================
#
#   Recursions of the multivariate Kalman filter
#
#=========================================================================

rm (list = ls(all=TRUE))
graphics.off()

#
#--------------------------- Helper Functions -------------------------------
#

#load required functions - recserar
source("EMTSUtil.R")

# Load required library - reshape, ones, zeros
library("matlab")

#--------------------------------------------------------------------------
# Multivariate Kalman filter
#--------------------------------------------------------------------------
lnlt <- function(y,Phi,Lam,R,Q) {
  # Allocate arrays
  t <- nrow(y)
  n <- ncol(y)
  k         <- nrow(Q)
  lnl       <- rep(0, t)
 
  # Recursions of the Kalman Filter
  # Initialisation following Harvey and Hamilton  
  st <- rep(0, k)
  pt <- t (reshape(inv(eye(k^2) - kronecker(Phi,Phi)) %*% as.vector(Q),k,k) )   

  mt <- Lam %*% st
  vt <- Lam %*% pt %*% t(Lam) + R
  ut <- y[1,] - mt
  
  lnl[1] <- - 0.5*n*log(2*pi) - 0.5*log(det(vt)) - 0.5*t(ut) %*% inv(vt) %*% ut
  
  Kal <- pt %*% t(Lam) %*% inv(vt)
  s0 <- st + Kal %*% ut
  p0 <- pt - Kal %*% Lam %*% pt
  
  # Main loop over observations
  for (i in 2:t) {
    # Prediction 
    st <- Phi %*% s0                     
    pt <- Phi %*% p0 %*% t(Phi) + Q          
              
    # Observation
    mt <- Lam %*% st
    vt <- Lam %*% pt %*% t(Lam) + R
    ut <- y[i,] - mt    
       
    # Construct log-likelihood function
    lnl[i] <- - 0.5*n*log(2*pi) - 0.5*log(det(vt)) - 0.5*t(ut) %*% inv(vt) %*% ut
              
  	# Update    			
    Kal <- pt %*% t(Lam) %*% inv(vt)
    s0 <- st + Kal %*% ut
    p0 <- pt - Kal %*% Lam %*% pt   
    
    cat('\nPredictions, Observations and Updates\n')
    cat('\nst = \n')
    print(st)
    cat('\npt = \n')
    print(pt)
    
    cat('\nvt = \n')
    print(vt)
    cat('\nut = \n')
    print(ut)
    
    cat('\ns0 = \n')
    print(s0)
    cat('\np0 = \n')
    print(p0)    
  }
  
}

#----------------- Recursions of the Multivariate Kalman Filter ------------------
#
lfac_multi <- function(){
  t   <- 5
  # Data from text
  y1 <- c(1.140, 2.315, -0.054, -1.545, -0.576)
  y2 <- c(3.235, 0.552, -0.689,  1.382,  0.718)
  y3 <- c(1.748, 1.472, -1.413, -0.199,  1.481)
  
  y  <- cbind(y1, y2, y3)
  
  # Parameter matrices
  Phi <- matrix(c(0.8,  0.0,
                  0.0,  0.5), nrow=2, byrow=T)
  Lam <- matrix(c(1.0,  0.5,
                  1.0,  0.0,
                  1.0, -0.5), nrow=3, byrow=T)
  R <- matrix(c(0.25,  0.00,  0.00,
                0.00,  0.16,  0.00,
                0.00,  0.00,  0.09), nrow=3, byrow=T)
  Q <- diag(2)

  # Kalman filter     
  lnlt(y,Phi,Lam,R,Q)
}

  

