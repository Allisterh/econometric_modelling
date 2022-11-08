#=========================================================================
#
#   Program to demonstrate the spurious regression problem 
#   using least squares regression
#
#=========================================================================
rm(list = ls(all=T))
graphics.off()
set.seed(123457, kind="Mersenne-Twister")

#
#------------------------- Helper Functions----------------------------------
#

#load required functions -  recserar, trimr, figure
source("EMTSUtil.R")
#
#-------------------- Spurious Regression Problem ---------------------------
#

nts_spurious2 <- function() {
  t <- 100                            
  r <- 10000
  
  # Generate sampling distributions and plot the histograms
  # I(0) versus I(0)      
  b <- array(0, c(r,2))
  pb <- txtProgressBar(min=0, max=r, style=3)
  for (i in seq(r)) {    
    y1 <- trimr(recserar(cbind(rnorm(t+100)),cbind(0.0),cbind(0.0)),100,0)
    y2 <- trimr(recserar(cbind(rnorm(t+100)),cbind(0.0),cbind(0.0)),100,0)
    
    b[i,] <- lm(y2 ~ cbind(rep(1, t), y1) - 1)$coef
    setTxtProgressBar(pb, i)
  }
  close(pb)
  figure()
  hist(b[,2],breaks=21, col="darkblue",
       main = ' I(0) vs. I(0)',
       xlab="", ylab="")
  
  # I(1) versus I(1)       
  pb <- txtProgressBar(min=0, max=r, style=3)
  for (i in seq(r)) {    
    y1 <- trimr(recserar(cbind(rnorm(t+100)),cbind(0.0),cbind(1.0)),100,0)
    y2 <- trimr(recserar(cbind(rnorm(t+100)),cbind(0.0),cbind(1.0)),100,0)
    
    b[i,] <- lm(y2 ~ cbind(rep(1, t), y1) - 1)$coef
    setTxtProgressBar(pb, i)
  }
  close(pb)
  figure()
  hist(b[,2],breaks=21, col="darkblue",
       main=' I(1) vs. I(1) ',
       xlab="", ylab="")
  # I(1) versus I(2)       
  pb <- txtProgressBar(min=0, max=r, style=3)
  for (i in seq(r)) {    
    y1 <- trimr(recserar(cbind(rnorm(t+100)),cbind(0.0),cbind(1.0)),100,0)
    y2 <- trimr(recserar(cbind(rnorm(t+100)),rbind(0.0, 0.0),rbind(2.0, -1.0)),100,0)
    
    b[i,] <- lm(y2 ~ cbind(rep(1, t), y1) - 1)$coef
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  figure()
  hist(b[,2],21, col="darkblue",
       main=' I(1) vs. I(2) ',
       xlab="", ylab="")

  # I(2) versus I(2)       
  pb <- txtProgressBar(min=0, max=r, style=3)
  for (i in seq(r)) {    
    y1 <- trimr(recserar(cbind(rnorm(t+100)),rbind(0.0, 0.0),rbind(2.0, -1.0)),100,0)
    y2 <- trimr(recserar(cbind(rnorm(t+100)),rbind(0.0, 0.0),rbind(2.0, -1.0)),100,0)
    
    b[i,] <- lm(y2 ~ cbind(rep(1, t), y1) - 1)$coef
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  figure()
  hist(b[,2],21, col="darkblue",
       main = ' I(2) vs. I(2) ',
       xlab="", ylab="")  
}
