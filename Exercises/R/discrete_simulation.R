#============================================================================
#
#   Program to demonstrate alternative qualitative response models
#
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()
set.seed(123456, kind="Mersenne-Twister")

#
#--------------------------- Helper Functions -------------------------------
#
source("EMTSUtil.R")


#----------------------------------------------------------------------------
#  Probit negative log-likelihood function 
#----------------------------------------------------------------------------
lprobit <- function (b,y,x) {
  f  <- pnorm(x %*% b)
  lf <- -mean( y*log(f) + (1 - y)*log(1 - f) )
  return (lf)
}

#----------------------------------------------------------------------------
#  Tobit negative log-likelihood function 
#----------------------------------------------------------------------------
ltobit <- function(b,y,x) {
  m <- x %*% b
  d <- y != 0  
  lf <- -mean( d*log(dnorm(y - m)) + (1 - d)*log(1 - pnorm(m)) )
  return(lf)
  
}

#----------------------------------------------------------------------------
#  Truncated negative log-likelihood function 
#----------------------------------------------------------------------------

ltrunc <- function(b,y,x) {
  m  <- x %*% b
  lf <- -mean( log(dnorm(y - m)) - log(1 - pnorm(0.0 - m)) )
  return(lf)
}


#
#--------------------------- Simulations  -----------------------------------
#

discrete_simulation <- function() {
  
  # Simulate data set 
  theta0 <- c(1,1)
  t      <- 100
  ndraws <- 3000
  zeros <- array(0, c(ndraws, 2))
  
  theta_full   <- zeros
  theta_probit <- zeros
  theta_tobit  <- zeros
  theta_trunc  <- zeros
  theta_rest   <- zeros
  
  pb <- txtProgressBar(min=0, max=ndraws, style=3)
  for (j in seq(ndraws)) {
    u <- rnorm(t)
    x <- cbind(rep(1, t), rnorm(t))
    y <- x %*% theta0 + u
    
    # Probit data
    y_1 <- rep(1, t)
    ind <- y < 0.0
    y_1[ind] <- 0.0
    
    # Tobit data
    ind <- y > 0.0
    y_2 <- y * ind  
    
    # Truncated data
    y_3 <- y[ind]
    x_3 <- x[ind, ]
    
    # Full data set results  
    theta_full[j,]<- lm(y ~ x - 1)$coef
    
    # Probit data set results  
    theta_probit[j,] <- optim(theta0, lprobit, y=y_1, x=x, method="BFGS")$par
    
    # Tobit data set results  
    theta_tobit[j,] <- optim(theta0, ltobit, y=y_2, x=x, method="BFGS")$par
    
    # Truncated data set results   
    theta_trunc[j,] <- optim(theta0, ltrunc, y=y_3, x=x_3, method="BFGS")$par
    
    # Restricted data set results  
    theta_rest[j,] <- lm(y_3 ~ x_3 - 1)$coef
    setTxtProgressBar(pb, j)
  }
  close(pb)
  
  
  cat('\nTrue parameter values = ',theta0[1],'  ',theta0[2])
  
  cat('\n')
  
  cat('\nMean (full)           = ',colMeans(theta_full))
  cat('\nMean (Probit)         = ',colMeans(theta_probit))
  cat('\nMean (Tobit)          = ',colMeans(theta_tobit))
  cat('\nMean (Trunctated      = ',colMeans(theta_trunc))
  cat('\nMean (Restricted)     = ',colMeans(theta_rest))
  
  cat('\n')
  
  
  cat('\nBias (full)           = ',colMeans(theta_full)-theta0)
  cat('\nBias (Probit)         = ',colMeans(theta_probit)-theta0)
  cat('\nBias (Tobit)          = ',colMeans(theta_tobit)-theta0)
  cat('\nBias (Truncated      = ',colMeans(theta_trunc)-theta0)
  cat('\nBias (Restricted)     = ',colMeans(theta_rest)-theta0)
  
  cat('\n')
  
  tmp <- theta_full - theta0
  cat('\nRMSE (full)           = ',colMeans(tmp)^2)
  
  tmp <- theta_probit - theta0
  cat('\nRMSE (Probit)         = ',colMeans(tmp)^2)
  
  tmp <- theta_tobit - theta0
  cat('\nRMSE (Tobit)          = ',colMeans(tmp)^2)
  
  tmp <- theta_trunc - theta0
  cat('\nRMSE (Truncated)      = ',colMeans(tmp)^2)
  
  tmp <- theta_rest - theta0
  cat('\nRMSE (Restricted)     = ',colMeans(tmp)^2)
}
