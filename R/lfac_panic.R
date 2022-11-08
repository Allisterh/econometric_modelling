#============================================================================
#
#  	Finite sample properties of the Bai and Ng (2004) of the PANIC model
#
#============================================================================

rm (list = ls(all=TRUE))
graphics.off()
set.seed(1234, kind="Mersenne-Twister")
#
#------------------------- Helper Functions----------------------------------
#

#load required functions - trimr, recserar
source("EMTSUtil.R")

#
#------  Sampling Properties of the Principal Components Estimator ----------
#

lfac_panic <- function() 
{
    # Settings
  t <- 200                # Number of observations                   
  n <- 50                 # Number of dependent variables        
  k <- 1                  # Number of factors     
  ndraws <- 5000          # Number of draws   
  lam <- 1 + rnorm(n)    # Choose lambda parameters N(1,1)          
  delta <- 1.0            # Choose delta parameter (same for all u) 

  # Allocate arrays
  lam_sw <- array(0, c(ndraws,n))  # Loadings of the Stock-Watson estimator
  lam_bn <- array(0, c(ndraws,n))  # Loadings of the Bai-Ng estimator

  pb <- txtProgressBar(min=0, max=ndraws, style=3) 
  for (iter in seq(ndraws)) {
    # Simulate the data  
    w  <- matrix(rnorm(t*n), nrow=t)    
    u  <- recserar(w,rbind(rep(0, n)),rbind(delta*rep(1,n)))  # Idiosyncratics                 
    v  <- rnorm(t)
    s  <- recserar(cbind(v),cbind(0),cbind(1))    # Factor                                  
    s  <- (s - mean(s))/apply(s, 2, sd)           # Standardize factor                   
    s0 <- s                                       # Define true factor 
    y  <- t(apply(s, 1, '*', lam)) + u
         
    # First difference the data     
    dy <- trimr(y,1,0) - trimr(y,0,1)
    dy <- dy - mean(dy)
                
    # Principal component analysis using levels  (Stock Watson)                                  
    anResults <- princomp(y)  
    s  <- anResults$scores
    s        <- s[,1:k]                # First k principal components    
    s        <- (s - mean(s))/sd(s)
                
    # Principal component analysis using first differences  (Bai-Ng estimator)
    anResults <- princomp(dy)  
    ds  <- anResults$scores               
    ds       <- ds[,1:k]                # First k principal components
    ds       <- (ds - mean(ds))/sd(ds)
 
    # Estimate factor loadings  
    b <- lm(y ~ cbind(rep(1, nrow(y)), s) - 1)$coef
    #  Stock-Waton: using levels (pick out the slopes)  
    lam_sw[iter,] <- b[2,]    

    # Bai-Ng: using first differences
    tmp            <- lm(dy ~ ds - 1)$coef
    lam_bn[iter,] <- tmp[,]    
    setTxtProgressBar(pb, iter)
  }
  close(pb)
  
  # some plots
  figure()
  par(xaxs="i", yaxs="i")
  matplot(seq(s),cbind(s0, s),type="l",
          main = 'Stock-Watson',
          xlab='',
          ylab='',
          bty = "l")

  figure()
  par(xaxs="i", yaxs="i")
  matplot(seq(s), cbind(s0, cbind(c(s0[1], cumsum(ds)))), type="l",
          main = 'Bai-Ng',
          xlab='',
          ylab='',
          bty = "l")
  # Compute statistics of sampling distributions and print results     
  mse_sw <- colMeans(t(apply(lam_sw, 1, '-', lam))^2)
  mse_bn <- colMeans( t(apply(lam_bn, 1, '-', lam))^2 )

  cat('\n ')
  cat('\nT       = ', t)
  cat('\ndelta   = ', delta)
  cat('\n')

  cat('\nStock-Watson Estimator\n')  
  print( cbind("Mean"=colMeans(lam_sw), "Bias"=colMeans(lam_sw)-lam, "MSE"=mse_sw ))   

  cat('\nBai-Ng Estimator\n')  
  print( cbind("Mean"=colMeans(lam_bn), "Bias"=colMeans(lam_bn)-lam, "MSE"=mse_bn))

  cat('\n ')
  cat('\nStock-Watson Estimator Overall')
  cat('\nMean    = ', mean(lam_sw))
  cat('\nBias    = ', mean(colMeans(lam_sw)-lam))
  cat('\nMSE     = ', mean(mse_sw))

  cat('\n ')
  cat('\nBai-Ng Estimator Overall')
  cat('\nMean    = ', mean(lam_bn))
  cat('\nBias    = ', mean(colMeans(lam_bn)-lam))
  cat('\nMSE     = ', mean(mse_bn))
  
}


