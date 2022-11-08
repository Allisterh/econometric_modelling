#============================================================================
#
#   Estimating GARCH-M Models
#
#============================================================================
rm(list = ls(all=T))
graphics.off()

#
#--------------------------- Helper Functions -------------------------------
# 
# Load required functions - inv
source("EMTSUtil.R")

#----------------------------------------------------------------------------
# Log-likelihood function 
#----------------------------------------------------------------------------
neglog <- function(b,r3m,r10) {
  b <- abs(b) 
  t <- length( r3m )
  m  <- rep(0, t)
  u <- rep(0, t)
  h <- sd( r10 - r3m )^2*rep(1,t)
  
  for (i in 2:t) {
    h[i] <- b[1] + b[2]*u[i-1]^2 + b[3]*h[i-1]
    m[i] <- b[4] + b[5]*r3m[i] + b[6]*(h[i]^(0.5))^b[7]
    u[i] <- r10[i] - m[i]
  } 
  z  <- u/sqrt(h)     
  f  <- -0.5*log(2*pi) - 0.5*log(h) - 0.5*z^2
  lf <- -mean( f )
  return(lf)
}

#
#------------------- Risk Modelling in Term Structure -----------------------
#
garch_m <- function( ) {
  
  # Load data from file
  data <- as.matrix(read.table("yields_us.dat"))
  
  
  r3m <- data[,1]         # Independent variable is the 3 month yield
  r10 <- data[,6]         # Dependent variable is the 10 year yield 
  t   <- length(r10)
  
  # Estimating the GARCH-M(1,1) model
  # Starting values are close to the maximum likelihood estimates         
  start <- c(0.0489928275220468,
             0.9617406702015530, 
             0.2753291619382821,
             2.2590829301737938, 
             0.7763991081157382, 
             0.1167781704457386, 
             1.0354494753290766)
  estResults <- optim(start, neglog, r3m=r3m, r10=r10, method="BFGS", hessian=T)
  theta <- estResults$par
  lf1 <- estResults$val
  hess <- estResults$hess
  
  lf1 <- -lf1
  vc <- (1/t)*inv(hess)
  
  cat('\n')
  cat('\nGARCH-M(1,1) Results')
  cat('\nMean parameters')
  cat('\ngamma0                                  = ',theta[4]) 
  cat('\ngamma1                                  = ',theta[5])
  cat('\nvarphi                                  = ',theta[6])
  cat('\nrho                                     = ',theta[7])
  cat('\nalpha_0                                 = ',theta[1])
  cat('\nalpha_1                                 = ',theta[2]) 
  cat('\nbeta_1                                  = ',theta[3]) 
  cat('\nLog-likelihood function (unrestricted)  = ',lf1)    
  
  # Wald test for risk neutrality    
  w <- (theta[6]-0.0)^2/vc[6,6]
  
  cat('\nWald test of varphi = 0    = ',w) 
  cat('\np-value                    = ',1-pchisq(w,1))
  
  # Wald test for risk-return relationship   
  w <- (theta[7]-1.0)^2/vc[7,7]
  cat('\nWald test of rho = 1       = ',w) 
  cat('\np-value                    = ',1-pchisq(w,1))
}


