#============================================================================
#
#     Plot a near unit root process
#
#============================================================================

rm (list = ls(all=TRUE))
graphics.off()
set.seed(123) #mt19937
 
#
#--------------------------- Helper Functions -------------------------------
#

# load utility functions -  recserar, trimr
source("EMTSUtil.R")

#
#--------------------------- Near Unit Root Processes -----------------------
#

t    <- 200                            
sig2 <- 0.1

# Simulate AR1 with phi=1.0 
mue  <- 0.0  
phi1 <- 1.0  
phi2 <- 0.0 
si1  <- 0.0 
si2  <- 0.0

vt  <- sqrt(sig2)*rnorm(t+101)      
yt  <- recserar( cbind(mue + trimr(vt,2,0) + si1*trimr(vt,1,1) + si2*trimr(vt,0,2)) , rbind(0,0) , rbind(phi1, phi2) )   
y1t <- trimr(yt,100,0)                               

# Simulate AR1 with phi=0.99      
mue  <- 0.0  
phi1 <- 0.99  
phi2 <- 0.0 
si1 <- 0.0 
si2 <- 0.0

vt  <- sqrt(sig2)*rnorm(t+101)   
yt  <- recserar( cbind(mue + trimr(vt,2,0) + si1*trimr(vt,1,1) + si2*trimr(vt,0,2)) , rbind(0,0) , rbind(phi1,phi2) )   
y2t <- trimr(yt,100,0)                               


#**********************************************************************
#***
#***     Generate graph
#***
#**********************************************************************

figure()
par(mfrow=c(1,2))
x <- seq(length(y1t))
#--------------------------------------------------------#
# Panel (a)
plot( x,y1t,typ='l', 
      main = expression("(a) " * phi *" = 1.0"),
      xlab = expression(t),
      ylab = expression(y[1]*t),
      ylim = c(-3,3),
      bty = 'l')



x <- seq(length(y2t))
#--------------------------------------------------------#
# Panel (b)
plot( x,y2t,typ='l', 
      main = expression("(b) " * phi *" = 0.99"),
      xlab = expression(t),
      ylab = expression(y[2]*t),
      ylim = c(-3,3),
      bty = 'l')
