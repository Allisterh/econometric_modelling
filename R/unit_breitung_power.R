#=========================================================================
#
#   Empirical power of the Breitung test
#
#=========================================================================
rm (list = ls(all=TRUE))
graphics.off()

#
#------------------------- Helper Functions----------------------------------
#

# load required functions -  trimr, recserar
source("EMTSUtil.R")

#
#------------- Unit Root Test Without Lags or Long-Run Variance -------------
#

unit_breitung_power <- function()
{
  # Parameters
  t    <- 1000
  cv   <- c(0,-5,-10,-15,-20)
  reps <- 100000
  tdf <- array(0, c(reps,length(cv)))
  rho <- array(0, c(reps,length(cv)))  
  for (i in seq(length(cv))) {
    set.seed(1234, kind="Mersenne-Twister")  
    pb <- txtProgressBar(min=0, max=reps, style=3)
    for (j in seq(reps)) {    
      y <- recserar(cbind(rnorm(t)),cbind(rnorm(1)),cbind(1+cv[i]/t))
      phihat <- lm(trimr(y,1,0) ~ trimr(y,0,1) - 1)$coef
      vhat <- trimr(y,1,0)-phihat * trimr(y,0,1)
      tdf[j,i] <- (phihat-1)/sqrt((t (vhat) %*% vhat/length(vhat))/(t( trimr(y,0,1)) %*% trimr(y,0,1)))
      rho[j,i] <- sum(cumsum(y)^2)/(sum(y^2)*(t^2))  
      setTxtProgressBar(pb, j)
    } 
    close(pb)
  }  
  print(cbind(cv=cv,tdf=colMeans(tdf< quantile(tdf[,1],0.05)), rho=colMeans(rho < quantile(rho[,1],0.05))  ) )    
}
