# ===========================================================================
#
#      Program to estimate a linear equation by GMM
#
# ===========================================================================

rm (list = ls(all=TRUE))
graphics.off()
set.seed(1234, kind="Mersenne-Twister")

#
#------------------------- Helper Functions----------------------------------
#

# load required functions - inv
source("EMTSUtil.R")

#----------------------------------------------------------------------------
# Procedure which defines the moment equations  
#----------------------------------------------------------------------------
meqn <- function(b,yt,xt) {
  t  <- length(yt)
  ut <- yt - b[1] - b[2]*xt
  zt <- cbind(rep(1, t), xt)
  cols <- ncol(zt)
  rows <- nrow(zt)
  
  k <- 1        
  mt <- array(0, c(rows,cols))
  for (j in seq(cols) ) {
    mt[,k] <- zt[,j]*ut
    k <- k+1    
  }
  return(mt)
}

    
#----------------------------------------------------------------------------
# Procedure which defines the mean of the moment conditions   
#----------------------------------------------------------------------------
meaneqn <- function(b,yt,xt) {
  return(colMeans(meqn(b,yt,xt)))
}

#----------------------------------------------------------------------------
# GMM objective function which also computes the optimal W
#----------------------------------------------------------------------------
gmmcrit <- function (b,yt,xt,flag,lmax) {
  d <- meqn(b,yt,xt)
  t <- nrow(d)
  g <- cbind(colMeans(d))
  
  if (flag) {    
    w  <- t(d) %*% d
    tau <- 1
    while (tau <= lmax) {
      wtau <- t( d[(tau+1):nrow(d),] ) %*% d[1:(nrow(d)-tau),]    
      w    <- w + (1.0-tau/(lmax+1)) * (wtau + t(wtau))
      tau  <- tau + 1         
    }  
    # Store globally - retrieved later
    w <<- w/t      
  } else {    
     w <- diag(ncol(d))    
  }
  ret <- t(g) %*% inv(w) %*% g
  return(ret)
}

# 
#--------------------- The Relationship Between GMM and OLS------------------
#

gmm_opt <- function() {
    t <- 1000    # Choose the sample size      

  xt <- runif(t)
  ut <- rnorm(t)
  
  yt <- 10 + 5*xt + ut

  lmax <- 2    #  Length of the Newey-West Covariance matrix estimation   

  # OLS estimates
  b0 <- lm(yt ~ cbind(rep(1, length(xt)), xt) - 1)$coef
  cat('\nOLS parameter estimates             = ', b0)

  # First stage estimation 
  
  flag <- FALSE
  estResults <- optim(b0, gmmcrit, yt=yt, xt=xt, flag=flag, lmax=lmax, method="BFGS")
  bgmm <- estResults$par

  # Second stage estimation
  flag <- TRUE
  estResults <- optim(b0, gmmcrit, yt=yt, xt=xt, flag=flag, lmax=lmax, method="BFGS")
  bgmm <- estResults$par  
  
  # computer w
  obj <- gmmcrit(bgmm,yt,xt,flag,lmax)  
  dg      <- numgrad(meaneqn,bgmm,yt,xt)    
  v       <- t(dg) %*% inv(w) %*% dg
  vinv    <- inv(v)/t  
  
  cat('\nThe value of the objective function = ', obj)
  cat('\nJ-test                              = ', t*obj)
  cat('\nGMM estimates                       = ', bgmm)
  cat('\nGMM standard errors                 = ', sqrt(diag(vinv)))
  cat('\nGMM t-statistics                    = ', bgmm/sqrt(diag(vinv)))
}
