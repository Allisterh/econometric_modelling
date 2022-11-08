# ========================================================================
#
#   Estiamte AR(1) coefficient by simulation of a MA(1) model
#
# ========================================================================

rm (list = ls(all=TRUE))
graphics.off()
set.seed(123, kind="Mersenne-Twister")

#
#------------------------- Helper Functions----------------------------------
#

#load required functions -  figure, trimr
source("EMTSUtil.R")

#
#---------------- Autocorrealtion Coefficients by simulation ----------------
#

sim_accuracy <- function() {
  theta <- 0.5        
  t     <- 250        
  h     <- 1         
  n     <- t*h        
  
  # Simulate the data for a MA(1) model    
  u  <- rnorm(n)
  ys <- trimr(u,1,0) - theta*trimr(u,0,1)
  
  # Estimate the first order autocorrelation coefficient    
  ys   <- ys - mean(ys)
  rhos <- lm(trimr(ys,1,0) ~ trimr(ys,0,1) - 1)$coef
  
  cat('\n')
  cat('\nSample size                        = ', t)
  cat('\nh                                  = ', h)
  cat('\nTrue population parameter (theta)  = ', theta)
  cat('\n')
  cat('\nTrue population parameter (rho)    = ', -theta/(1+theta^2))
  cat('\nSimulation estimate (rho)          = ', rhos)
}
