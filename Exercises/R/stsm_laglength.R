#============================================================================
#
#   Simulate an AR(3) model and compute the optimal lag length
#
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()
set.seed(1, kind="Mersenne-Twister")

#
#------------------------- Helper Functions----------------------------------
#

#load required functions -  figure
source("EMTSUtil.R")

#
#------------------------ Lag Length Determination --------------------------
#

stsm_laglength <- function() {
  # Parameters
  nreps <- 2000
  t    <- 300                       
  pmax <- 7      
  mu   <- 0.0
  phi1 <- 0.2    
  phi2 <- -0.15    
  phi3 <- 0.05
  
  zeros <- rep(0, pmax)
  acount <- zeros
  scount <- zeros
  hcount <- zeros
  pb <- txtProgressBar(min=0, max=nreps, style=3)
  for (count in seq(nreps)) {
    # Generate data
    yt <- rep(0, t)
    for (i in 4:t) {
      yt[i] <- mu+phi1*yt[i-1]+phi2*yt[i-2]+phi3*yt[i-3]+sqrt(0.5)*rnorm(1)    
    }
    # Set up lags and make sure that T is constant
    # lag the matrix
    ylag <- c(rep(NA, pmax), yt)  
    ylag <- embed(ylag,pmax+1)  
    ylag <- trimr(ylag, 100, 0)
    
    # Loop over the lags (yt is first column of ylag)
    aic <- zeros
    sic <- zeros
    hic <- zeros
    y   <- ylag[,1]
    tt  <- length( y )
    for (j in 2:pmax+1) {
      x <- cbind(rep(1, tt), ylag[,(2:j)])
      k <- ncol(x)
      b <- lm(y ~ x - 1)$coef    
      e <- y-x %*% b
      aic[j-1] <- log(t(e) %*% e/tt) + 2*k/tt
      sic[j-1] <- log(t(e) %*% e/tt) + k*log(tt)/tt
      hic[j-1] <- log(t(e) %*% e/tt) + 2*k*log(log(tt))/tt    
    }
    ind <- which.min(aic)
    acount[ind] <- acount[ind]+1
    
    ind <- which.min(sic)
    scount[ind] <- scount[ind]+1
    
    ind <- which.min(hic)
    hcount[ind] <- hcount[ind]+1
    setTxtProgressBar(pb, count)  
  }
  close(pb)
        
  
  #**************************************************************************
  #**
  #**     Generate graph
  #**
  #**************************************************************************
  
  figure()
  par(mfrow=c(3,1), xaxs="i", yaxs="i")
  barplot(acount, main="AIC", xlab="Lag Order", col="grey")
  barplot(scount, main="SIC", xlab="Lag Order", col="lightgrey")
  barplot(hcount, main="HIC", xlab="Lag Order", col="black")
}




