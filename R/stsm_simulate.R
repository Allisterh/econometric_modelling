#============================================================================
#
#   Simulate an ARMA(2,2) model and compute the ACF and the PACF
#
#============================================================================

rm (list = ls(all=TRUE))
graphics.off()
set.seed(12, kind="Mersenne-Twister")

#
#------------------------- Helper Functions----------------------------------
#

#load required functions - figure, recserar, trimr
source("EMTSUtil.R")

#
#------------------------  ARMA Model simulation ----------------------------
#
stsm_simulate <- function() {
  t    <- 200                 # Define the sample size      
  lags <- 10                  # Number of lags in the ACF and the PACF      
  
  # Generate error process
  ut <- sqrt(0.1)*rnorm(t+101)      #An additional 100 observations 
  
  # ARMA(0,0) 
  mue  <- 0.0  
  phi1 <- 0.0    
  phi2 <- 0.0     
  si1  <- 0.0      
  si2  <- 0.0   
  
  yt    <- recserar( cbind(mue + trimr(ut,2,0) + si1*trimr(ut,1,1)+ si2*trimr(ut,0,2)) , cbind(c(0.0, 0.0)) , cbind(c(phi1,phi2)))   
  y0t   <- trimr(yt,100,0)                               
  acf0  <- acf( y0t,lags, plot=F)$acf
  pacf0 <- pacf( y0t,lags, plot=F)$acf
  
  
  # ARMA(1,0)
  mue  <- 0.0  
  phi1 <- 0.7    
  phi2 <- 0.0     
  si1  <- 0.0      
  si2  <- 0.0   
  
  yt    <- recserar( cbind(mue + trimr(ut,2,0) + si1*trimr(ut,1,1)+ si2*trimr(ut,0,2)) , cbind(c(0.0, 0.0)) , cbind(c(phi1,phi2)))   
  y1t   <- trimr(yt,100,0)                      # Trim the data  
  acf1  <- acf( y1t,lags, plot=F)$acf           # Compute the ACF                             
  pacf1 <- pacf( y1t,lags, plot=F)$acf           # Compute the PACF                            
  
  
  # ARMA(2,0)
  mue  <- 0.0  
  phi1 <- 0.7    
  phi2 <- -0.5     
  si1  <- 0.0      
  si2  <- 0.0   
  
  yt    <- recserar( cbind(mue + trimr(ut,2,0) + si1*trimr(ut,1,1)+ si2*trimr(ut,0,2)) , cbind(c(0.0, 0.0)) , cbind(c(phi1,phi2)))  
  y2t   <- trimr(yt,100,0)          # Trim the data  
  acf2  <- acf( y2t,lags, plot=F)$acf            # Compute the ACF                             
  pacf2 <- pacf( y2t,lags, plot=F)$acf           # Compute the PACF                            
  
  
  
  # ARMA(0,1)
  mue  <- 0.0  
  phi1 <- 0.0    
  phi2 <- 0.0     
  si1  <- 0.9      
  si2  <- 0.0   
  
  yt    <- recserar( cbind(mue + trimr(ut,2,0) + si1*trimr(ut,1,1)+ si2*trimr(ut,0,2)) , cbind(c(0.0, 0.0)) , cbind(c(phi1,phi2)))   
  y3t   <- trimr(yt,100,0)            
  acf3  <- acf(y3t,lags, plot=F)$acf            # Compute the ACF
  pacf3 <- pacf(y3t,lags, plot=F)$acf           # Compute the PACF 
  
  
  # ARMA(0,2)		
  mue  <- 0.0  
  phi1 <- 0.0    
  phi2 <- 0.0     
  si1  <- -0.2      
  si2  <- 0.7   
  
  yt    <- recserar( cbind(mue + trimr(ut,2,0) + si1*trimr(ut,1,1)+ si2*trimr(ut,0,2)) , cbind(c(0.0, 0.0)) , cbind(c(phi1,phi2))) 
  y4t   <- trimr(yt,100,0)         
  acf4  <- acf(y4t,lags, plot=F)$acf            # Compute the ACF
  pacf4 <- pacf(y4t,lags, plot=F)$acf           # Compute the PACF 
  
  # ARMA(1,1)		
  mue  <- 0.0  
  phi1 <- 0.8    
  phi2 <- 0.0     
  si1  <- 0.7      
  si2  <- 0.0   
  
  yt    <- recserar( cbind(mue + trimr(ut,2,0) + si1*trimr(ut,1,1)+ si2*trimr(ut,0,2)) , cbind(c(0.0, 0.0)) , cbind(c(phi1,phi2)))  
  y5t   <- trimr(yt,100,0)         
  acf5  <- acf(y5t,lags, plot=F)$acf            # Compute the ACF
  pacf5 <- pacf(y5t,lags, plot=F)$acf           # Compute the PACF 
  
  
  #**********************************************************************
  #***
  #***     Plot the series
  #***
  #**********************************************************************
  
  
  figure()
  par(mfrow=c(3,3))
  
  t <- 0:lags  
  
  #--------------------------------------------------------#
  # Panel (a)  
  plot(y2t,type="l", 
       main='(a) ARMA(2,0)',
       xlab = expression(t),
       ylab = expression(y[t]),
       bty="l")    
  #--------------------------------------------------------#
  # Panel (b)
  bp <- barplot(c(acf2),
                main = '(b) ACF ARMA(2,0)',
                xlab = 'lag',
                ylab = 'ACF',
                bty="l")
  axis(1, at=bp, labels=t)
  
  #--------------------------------------------------------#
  # Panel (c)
  bp <- barplot(c(1, pacf2),
                main = '(b) PACF ARMA(2,0)',
                xlab = 'lag',
                ylab = 'PACF',
                bty="l")
  axis(1, at=bp, labels=t)
  
  #--------------------------------------------------------#
  # Panel (d)  
  plot(y4t,type="l",
       main = '(a) ARMA(0,2)',
       xlab = expression(t),
       ylab = expression(y[t]),
       bty="l")
  
  #--------------------------------------------------------#
  # Panel (e)
  bp <- barplot(c(acf4),
                main = '(b) ACF ARMA(0,2)',
                xlab = 'lag',
                ylab = 'ACF',
                bty="l")
  axis(1, at=bp, labels=t)
  
  #--------------------------------------------------------#
  # Panel (f)
  bp <- barplot(c(1, pacf4),
                main = '(b) PACF ARMA(0,2)',
                xlab = 'lag',
                ylab = 'PACF',
                bty="l")
  axis(1, at=bp, labels=t)
  
  #--------------------------------------------------------#
  # Panel (g)  
  plot(y5t,type="l",
       main = '(a) ARMA(1,1)',
       xlab = expression(t),
       ylab = expression(y[t]),
       bty="l")
  
  #--------------------------------------------------------#
  # Panel (h)
  bp <- barplot(c(acf5),
                main = '(b) ACF ARMA(1,1)',
                xlab = 'lag',
                ylab = 'ACF',
                bty="l")
  axis(1, at=bp, labels=t)
  
  #--------------------------------------------------------#
  # Panel (i)
  bp <- barplot(c(1, pacf5),
                main = '(b) PACF ARMA(1,1)',
                xlab = 'lag',
                ylab = 'PACF',
                bty="l")
  axis(1, at=bp, labels=t)
  
  
  figure()
  par(mfrow=c(3,3))
  #--------------------------------------------------------#
  # Panel (a)  
  plot(y0t,type="l", 
       main='(a) ARMA(0,0)',
       xlab = expression(t),
       ylab = expression(y[t]),
       bty="l")    
  #--------------------------------------------------------#
  # Panel (b)
  bp <- barplot(c(acf0),
                main = '(b) ACF ARMA(0,0)',
                xlab = 'lag',
                ylab = 'ACF',
                bty="l")
  axis(1, at=bp, labels=t)
  
  #--------------------------------------------------------#
  # Panel (c)
  bp <- barplot(c(1, pacf0),
                main = '(b) PACF ARMA(0,0)',
                xlab = 'lag',
                ylab = 'PACF',
                bty="l")
  axis(1, at=bp, labels=t)
  
  #--------------------------------------------------------#
  # Panel (d)  
  plot(y1t,type="l",
       main = '(a) ARMA(1,0)',
       xlab = expression(t),
       ylab = expression(y[t]),
       bty="l")
  
  #--------------------------------------------------------#
  # Panel (e)
  bp <- barplot(c(acf1),
                main = '(b) ACF ARMA(1,0)',
                xlab = 'lag',
                ylab = 'ACF',
                bty="l")
  axis(1, at=bp, labels=t)
  
  #--------------------------------------------------------#
  # Panel (f)
  bp <- barplot(c(1, pacf1),
                main = '(b) PACF ARMA(1,0)',
                xlab = 'lag',
                ylab = 'PACF',
                bty="l")
  axis(1, at=bp, labels=t)
  
  #--------------------------------------------------------#
  # Panel (g)  
  plot(y2t,type="l",
       main = '(a) ARMA(2,0)',
       xlab = expression(t),
       ylab = expression(y[t]),
       bty="l")
  
  #--------------------------------------------------------#
  # Panel (h)
  bp <- barplot(c(acf2),
                main = '(b) ACF ARMA(2,0)',
                xlab = 'lag',
                ylab = 'ACF',
                bty="l")
  axis(1, at=bp, labels=t)
  
  #--------------------------------------------------------#
  # Panel (i)
  bp <- barplot(c(1, pacf2),
                main = '(b) PACF ARMA(2,0)',
                xlab = 'lag',
                ylab = 'PACF',
                bty="l")
  axis(1, at=bp, labels=t)
  
  figure()
  par(mfrow=c(3,3))
  
  #--------------------------------------------------------#
  # Panel (a)  
  plot(y3t,type="l", 
       main='(a) ARMA(0,1)',
       xlab = expression(t),
       ylab = expression(y[t]),
       bty="l")    
  #--------------------------------------------------------#
  # Panel (b)
  bp <- barplot(c(acf3),
                main = '(b) ACF ARMA(0,1)',
                xlab = 'lag',
                ylab = 'ACF',
                bty="l")
  axis(1, at=bp, labels=t)
  
  #--------------------------------------------------------#
  # Panel (c)
  bp <- barplot(c(1, pacf3),
                main = '(b) PACF ARMA(0,1)',
                xlab = 'lag',
                ylab = 'PACF',
                bty="l")
  axis(1, at=bp, labels=t)
  
  #--------------------------------------------------------#
  # Panel (d)  
  plot(y4t,type="l",
       main = '(a) ARMA(0,2)',
       xlab = expression(t),
       ylab = expression(y[t]),
       bty="l")
  
  #--------------------------------------------------------#
  # Panel (e)
  bp <- barplot(c(acf4),
                main = '(b) ACF ARMA(0,2)',
                xlab = 'lag',
                ylab = 'ACF',
                bty="l")
  axis(1, at=bp, labels=t)
  
  #--------------------------------------------------------#
  # Panel (f)
  bp <- barplot(c(1, pacf4),
                main = '(b) PACF ARMA(0,2)',
                xlab = 'lag',
                ylab = 'PACF',
                bty="l")
  axis(1, at=bp, labels=t)
  
  #--------------------------------------------------------#
  # Panel (g)  
  plot(y5t,type="l",
       main = '(a) ARMA(1,1)',
       xlab = expression(t),
       ylab = expression(y[t]),
       bty="l")
  
  #--------------------------------------------------------#
  # Panel (h)
  bp <- barplot(c(acf5),
                main = '(b) ACF ARMA(1,1)',
                xlab = 'lag',
                ylab = 'ACF',
                bty="l")
  axis(1, at=bp, labels=t)
  
  #--------------------------------------------------------#
  # Panel (i)
  bp <- barplot(c(1, pacf5),
                main = '(b) PACF ARMA(1,1)',
                xlab = 'lag',
                ylab = 'PACF',
                bty="l")
  axis(1, at=bp, labels=t)
  
  
}
  