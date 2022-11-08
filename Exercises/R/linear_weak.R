#==========================================================================
#
#     Program to generate Figure 1a of Stock, Wright and Yogo 
#     weak instrument paper (JBES, 2002)
#
#==========================================================================

rm (list = ls(all=TRUE))
graphics.off()
set.seed(1234, kind="Mersenne-Twister")

#
#--------------------------- Helper Functions -----------------------------------
#
# Load helper functions - inv, figure
source("EMTSUtil.R")
# load ks library - 
library(ks)
#
#--------------------------- Weak Instruments -----------------------------------
#

#---------------------------------------------------------------------------------
#     Generate instruments (scaled) and do Monte Carlo replications
#---------------------------------------------------------------------------------
linear_weak <- function() {
  N     <- 10000          # Number of replications  
  T     <- 5              # Sample size            
  beta  <- 0.0            # Parameter values      
  phi   <- 0.25
  sig11 <- 1.0
  sig22 <- 1.0
  sig12 <- 0.99
  omega <- matrix(c(sig11, sig12,  sig12, sig22), 
                  ncol=2, byrow=T)
  rho   <- sig12/sqrt(sig11*sig22)  
  x <- cbind(runif(T))
  
  biv <- rep(0, N)
  pb <- txtProgressBar(min = 0, max = N, style=3)
  for (i in seq(N)) { 
    u  <- matrix(rnorm(T*nrow(omega)), ncol = 2) %*% chol(omega)
    w <- x*phi  + u[ ,1]                                          # Reduced form equation   
    y <- w*beta + u[,2]                                           # Structural equation
  
    tmp    <- inv(t(x) %*% x)
    tmp1   <- x %*% tmp %*% t(x)
    tmp2   <- t(w) %*% tmp1
    biv[i] <- inv(tmp2 %*% w) %*% (tmp2 %*% y)      # IV estimates
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  #**************************************************************************
  #**
  #**     Generate graph
  #**
  #**************************************************************************
  xi <- seq(-2.0, 2.0, 0.001)
  f  <- kde(biv,h = 0.07, eval.points=xi)
  figure()
  plot(f, mar=c(5, 5, 4, 1), xaxs="i", yaxs="i",
       xlab = expression(paste( hat(beta), ""[IV]*"")),
       ylab =  expression(paste( f,"(",hat(beta), ""[IV]*"", ")")),
       bty="l")

  
}

