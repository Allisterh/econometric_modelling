#============================================================================
#
#      Program to estimate the consumption capm model by GMM
#
#============================================================================

rm (list = ls(all=TRUE))
graphics.off()

#
#------------------------- Helper Functions----------------------------------
#

# Load required functions - trimr, inv
source("EMTSUtil.R")

#----------------------------------------------------------------------------
# Define the moment equations 
#----------------------------------------------------------------------------
meqn <- function(b,cratio,r,zt) {
  cols <- ncol(zt)
  rows <- nrow(zt)
  
  beta <- b[1]       #  Relative risk aversion parameter    
  gam  <- b[2]       #  Discount parameter                  

  ut <- beta*cratio^(-gam) * (1 + r) - 1
  
  k <- 1        
  mt <- array(0, c(rows,cols))  
  for (j in seq(cols)) {
    mt[,k] <- zt[,j]*ut
    k <- k+1    
  }  
  return(mt)  
}
   
    
#----------------------------------------------------------------------------
# Defines the mean of the moment conditions  
#----------------------------------------------------------------------------
meaneqn <- function(b,cratio,r,zt) {
  ret <- colMeans(meqn(b,cratio,r,zt))
  return(ret)
}

#----------------------------------------------------------------------------
# GMM objective function which also computes the optimal w   
#----------------------------------------------------------------------------
q <- function(b,cratio,r,zt,lmax){
  
  d <- meqn(b,cratio,r,zt)
  g <- cbind(colMeans(d))
  
  w  <- t(d) %*% d
  tau <- 1
  while (tau <= lmax) {
    wtau <- t( d[(tau+1):nrow(d),] ) %*% d[1:(nrow(d)-tau),]
    w    <- w + (1.0-tau/(lmax+1))*(wtau + t(wtau))
    tau  <- tau + 1    
  }  
  t <- length(cratio)
  w <- w/t
  ret <- t(g) %*% inv(w) %*% g
  
  return(ret)  
}

#
#----------------------- The Consumption Based on C-CAPM --------------------
#

gmm_ccapm <- function()
{
   xdata <- as.matrix(read.table("ccapm.dat"))
  
  cratio    <- xdata[,1]      # Ratio of real consumption at t+1 and t   
  g_return  <- xdata[,2]      # Real gross return at t+1                        
  r         <- g_return - 1    # Real interest rate at t+1                      
  e         <- xdata[,3]      # Value weighted real returns 
  
  t <- length(cratio)-1
  
  # Instruments <- {const,cratio,rlag,e}
  zt <- cbind(rep(1, t),trimr(cratio,0,1), trimr(r,0,1), trimr(e,0,1))
  
  b0 <- c(0.9, 2.0)
  estResults <- optim(b0, q, cratio=trimr(cratio,1,0), r=trimr(r,1,0), zt=zt, lmax=0, method="BFGS")
  bgmm <- estResults$par  
  qmin <- estResults$value

  # Numerical gradient at optimum
  dg <- numgrad(meaneqn,bgmm,trimr(cratio,1,0), trimr(r,1,0), zt)
    
  # Compute optimal w
  d    <- meqn(bgmm,trimr(cratio,1,0),trimr(r,1,0),zt)
  g    <- cbind(colMeans(d))
  w    <- t(d) %*% d
  tau  <- 1
  lmax <- 0
  while (tau <= lmax) {
    wtau <- t( d[(tau+1):nrow(d),] ) %*% d[1:(nrow(d)-tau),]
    w    <- w + (1.0-tau/(lmax+1))*(wtau + t(wtau))
    tau  <- tau + 1    
  }
  w <- w/t
  v <- t(dg) %*% inv(w) %*% dg

  # Hansen Sargan J Test
  j_stat <- t*q(bgmm,trimr(cratio,1,0),trimr(r,1,0),zt,0)
  
  cat('\nThe value of the objective function (q)  = ', qmin )
  cat('\nJ-test is                                = ', j_stat, '\n' )  
  se <- sqrt(diag(inv(v)/t))
  print(cbind(Estimates=bgmm, se=se, t=bgmm/se))
}


 

