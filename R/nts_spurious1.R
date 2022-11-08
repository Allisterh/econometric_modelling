#============================================================================
#
#   Program to demonstrate the spurious regression problem 
#   using correlation coefficients
#
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()
set.seed(15, kind="Mersenne-Twister")

#
#------------------------- Helper Functions----------------------------------
#

#load required functions -  figure, trimr, recserar
source("EMTSUtil.R")


#
#--------------------------- Spurious Regression ----------------------------
#

nts_spurious1 <- function () {
  t <- 100                            
  r <- 10000                       
  
  # I(0) versus I(0)      
  rca <- rep(0, r)
  pb <- txtProgressBar(min=0, max=r, style=3)
  for (i in seq(r)) {
    y1   <- trimr(recserar(cbind(rnorm(t+100)),rbind(0.0),rbind(0.0)),100,0)
    y2   <- trimr(recserar(cbind(rnorm(t+100)),rbind(0.0),rbind(0.0)),100,0)
    rca[i] <- cor(y1, y2) 
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  y <- quantile( rca,probs=c(.01, .025, .05, .10, .50, .90, .95, .975, .99) )
  
  cat('\n ' )
  cat('\nThe empirical distribution (stationary case)' )
  cat('\n********************************************' )
  cat('\n 1.0 per cent      =',y[1])
  cat('\n 2.5 per cent      =',y[2])
  cat('\n 5.0 per cent      =',y[3])
  cat('\n10.0 per cent      =',y[4])
  cat('\n50.0 per cent      =',y[5])
  cat('\n90.0 per cent      =',y[6])
  cat('\n95.0 per cent      =',y[7])
  cat('\n97.5 per cent      =',y[8])
  cat('\n99.0 per cent      =',y[9])
  
   
  # I(1) versus I(1)       
  rcb <- rep(0, r)
  pb <- txtProgressBar(min=0, max=r, style=3)
  for (i in seq(r)) {
    y1   <- trimr(recserar(cbind(rnorm(t+100)),rbind(0.0),rbind(1.0)),100,0)
    y2   <- trimr(recserar(cbind(rnorm(t+100)),rbind(0.0),rbind(1.0)),100,0)
    rcb[i] <- cor(y1,y2)  
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  y <- quantile(rcb, probs=c(.01, .025, .05, .10, .50, .90, .95, .975, .99) )
  
  cat('\n ' )
  cat('\nThe empirical distribution (case 2)' )
  cat('\n********************************************' )
  cat('\n 1.0 per cent      =',y[1]) 
  cat('\n 2.5 per cent      =',y[2]) 
  cat('\n 5.0 per cent      =',y[3])
  cat('\n10.0 per cent      =',y[4])
  cat('\n50.0 per cent      =',y[5])
  cat('\n90.0 per cent      =',y[6])
  cat('\n95.0 per cent      =',y[7])
  cat('\n97.5 per cent      =',y[8])
  cat('\n99.0 per cent      =',y[9])
  
  #
  # I(1) versus I(2)             
  rcc <- rep(0, r)
  pb <- txtProgressBar(min=0, max=r, style=3)
  for (i in seq(r)) {
    y1   <- trimr(recserar(cbind(rnorm(t+100)),rbind(0.0),rbind(1.0)),100,0)
    y2   <- trimr(recserar(cbind(rnorm(t+100)),rbind(0.0,0.0),rbind(2.0,-1.0)),100,0)
    rcc[i] <- cor(y1, y2)
    setTxtProgressBar(pb, i)
  }
  close(pb)
      
  
  y <- quantile(rcc, probs=c(.01, .025, .05, .10, .50, .90, .95, .975, .99) )
  
  cat('\n ' )
  cat('\nThe empirical distribution (case 3)' )
  cat('\n********************************************' )
  cat('\n 1.0 per cent      =',y[1]) 
  cat('\n 2.5 per cent      =',y[2]) 
  cat('\n 5.0 per cent      =',y[3])
  cat('\n10.0 per cent      =',y[4])
  cat('\n50.0 per cent      =',y[5])
  cat('\n90.0 per cent      =',y[6])
  cat('\n95.0 per cent      =',y[7])
  cat('\n97.5 per cent      =',y[8])
  cat('\n99.0 per cent      =',y[9])
   
  # I(2) versus I(2)       
  rcd <- rep(0, r)
  pb <- txtProgressBar(min=0, max=r, style=3)
  for (i in seq(r)) {
    y1   <- trimr(recserar(cbind(rnorm(t+100)),rbind(0.0,0.0),rbind(2.0,-1.0)),100,0)
    y2   <- trimr(recserar(cbind(rnorm(t+100)),rbind(0.0,0.0),rbind(2.0,-1.0)),100,0)
    rcd[i] <- cor(y1, y2)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  y <- quantile(rcd, probs=c(.01, .025, .05, .10, .50, .90, .95, .975, .99))
  
  cat('\n ' )
  cat('\nThe empirical distribution (case 4)' )
  cat('\n********************************************' )
  cat('\n 1.0 per cent      =',y[1]) 
  cat('\n 2.5 per cent      =',y[2]) 
  cat('\n 5.0 per cent      =',y[3])
  cat('\n10.0 per cent      =',y[4])
  cat('\n50.0 per cent      =',y[5])
  cat('\n90.0 per cent      =',y[6])
  cat('\n95.0 per cent      =',y[7])
  cat('\n97.5 per cent      =',y[8])
  cat('\n99.0 per cent      =',y[9])
  
  
  #**********************************************************************
  #***
  #***     Generate graph
  #***
  #**********************************************************************
  
  figure()
  par(xaxs="i", yaxs="i", mfrow=c(2,2), mar=c(5, 5, 4, 1))
      
  #--------------------------------------------------------#
  # Panel (a)
  hist(rca,breaks=21,
       main = '(a) I[0] vs I[0]',
       xlab = expression(hat(rho)),
       ylab = expression(f(hat(rho))),
       bty = "l")
       
  # #--------------------------------------------------------#
  # # Panel (b)
  hist(rcb,breaks=21,
       main = '(a) I[1] vs I[1]',
       xlab = expression(hat(rho)),
       ylab = expression(f(hat(rho))),
       bty = "l")
    
  #--------------------------------------------------------#
  # Panel (c)  
  hist(rcc,breaks=21,
       main = '(a) I[1] vs I[2]',
       xlab = expression(hat(rho)),
       ylab = expression(f(hat(rho))),
       bty = "l")
  #--------------------------------------------------------#
  # Panel (d)
  hist(rcd,breaks=21,
       main = '(a) I[2] vs I[2]',
       xlab = expression(hat(rho)),
       ylab = expression(f(hat(rho))),
       bty = "l")
}
