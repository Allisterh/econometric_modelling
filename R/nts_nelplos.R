#============================================================================
#
#   Plot and compute ACF of the Nelson-Plosser data (1860 to 1970).
#
#============================================================================

rm (list = ls(all=TRUE))
graphics.off()
set.seed(5, kind="Mersenne-Twister")

#
#------------------------- Helper Functions----------------------------------
#

#load required functions -  figure recserar, trimr, seqa
source("EMTSUtil.R")

#--------------------------------------------------------------------------
# Regression estimates of trend models
#   ntrend = 0 (constant and lag)
#   ntrend = 1 (constant, time trend and lag)
#   ntrend = 2 (constant and time trend)
#--------------------------------------------------------------------------
ar <- function(y,ntrend) {
  t <- length(y)
  if (ntrend == 0.0) {
    x <- cbind(rep(1, t-1), trimr(y,0,1) )                                
    y <- trimr(y,1,0)   
  } else if (ntrend == 1.0) {
    x <- cbind(rep(1, t-1), seqa(0,1,t-1)/100,  trimr(y,0,1))     
    y <- trimr(y,1,0)   
  } else {
    x <- cbind(rep(1, t), seqa(0,1,t)/100 )                          
    y <- trimr(y,0,0)    
  }  
  b <- lm(y ~ x - 1)$coef                                    
  e <- y - x %*%b                                 
  s2 <- t(e) %*% e/length(e)   
  s2 <- as.vector(s2)
  se <- sqrt( diag( s2 * inv(t(x) %*% x) ) )
  
  return(list(b=b, se=se))
}
#--------------------------------------------------------------------------
# Prints the acf and lags to the screen
#--------------------------------------------------------------------------
print.acf <- function(y, lags) {
  y.acf <- acf(y, lags, plot=F)
  cat('\n')
  print(rbind(y.acf$lag, y.acf$acf))
}





#
#------------------------- Nelson-Plosser Data ------------------------------
#

nts_nelplos <- function() {
  
  # Load Nelson-Plosser data set
  load('nelson_plosser.Rdata')
  
  # Variable 'data' constains the following variables
  # Date 1860 - 1970
  # Real GNP, 1909 to 1970  
  # Nominal GNP, 1909 to 1970
  # Real per capita GNP, 1909 to 1970
  # Industrial production, 1860 to 1970
  # Employment, 1890 to 1970
  # Unemployment rate, 1890 to 1970
  # GNP deflator, 1889 to 1970
  # Consumer prices index, 1860 to 1970
  # Wages, 1900 to 1970
  # Real wages, 1900 to 1970 
  # Money stock, 1889 to 1970
  # Velocity, 1869 to 1970
  # Bond yield, 1900 to 1970
  # SP500, 1871 to 1970
  
  # Take logs except for bond yield
  rgnp   <- log(data[50:111,2])
  gnp    <- log(data[50:111,3])
  pcrgnp <- log(data[50:111,4])
  ip     <- log(data[1:111,5])
  emp    <- log(data[31:111,6])
  un     <- log(data[31:111,7])
  prgnp  <- log(data[30:111,8])
  cpi    <- log(data[1:111,9])
  wg     <- log(data[41:111,10])
  rwg    <- log(data[41:111,11])
  m      <- log(data[30:111,12])
  vel    <- log(data[10:111,13])
  bnd    <- data[41:111,14]
  sp500  <- log(data[12:111,15])
  
  
  # Generate deterministic and stochastic variables
  
  t <- 111
  
  y_dtrend <- 0.1 + 0.2*seq(t) + rnorm(t)      
  
  y_strend  <- recserar(cbind(0.3 + rnorm(t)),cbind(0.0), cbind(1.0))          
  
  
  #**********************************************************************
  #***
  #***     Generate graph
  #***
  #**********************************************************************
  
  figure()
  par(xaxs="i", yaxs="i", mfrow=c(4,4))    
  
  xmin <- 1860
  xmax <- 1970
  
  #--------------------------------------------------------#
  # Panel (a)
  plot(seq(1909,1970,1),rgnp,type="l",
       main = 'RGNP',
       xlab = "",
       ylab = "",
       xlim = c(xmin, xmax),
       bty="l")  
  
  
  #--------------------------------------------------------#
  # Panel (b)  
  plot(seq(1909,xmax,1),gnp,type="l",
       main = 'GNP',
       xlab = "",
       ylab = "",
       xlim = c(xmin, xmax),
       bty="l")  
  
  #--------------------------------------------------------#
  # Panel (c)
  plot(seq(1909,xmax,1),pcrgnp,type="l",
       main = 'PCRGNP',
       xlab = "",
       ylab = "",
       xlim = c(xmin, xmax),
       bty="l")
  
  #--------------------------------------------------------#
  # Panel (d)
  plot(seq(xmin,xmax,1),ip,type="l",
       main = 'IP',
       xlab = "",
       ylab = "",
       xlim = c(xmin, xmax),
       bty="l")
  
  #--------------------------------------------------------#
  # Panel (e)
  plot(seq(1890,xmax,1),emp,type="l",
       main = 'Employment',
       xlab = "",
       ylab = "",
       xlim = c(xmin, xmax),
       bty="l")
  
  #--------------------------------------------------------#
  # Panel (f)
  plot(seq(1890,xmax,1),un,type="l",
       main = 'Unemployment',
       xlab = "",
       ylab = "",
       xlim = c(xmin, xmax),
       bty="l")
  
  #--------------------------------------------------------#
  # Panel (g)
  plot(seq(1889,xmax,1),prgnp,type="l",
       main = 'PRGNP',
       xlab = "",
       ylab = "",
       xlim = c(xmin, xmax),
       bty="l")
  
  #--------------------------------------------------------#
  # Panel (h)
  plot(seq(1860,xmax,1),cpi,type="l",
       main = 'CPI',
       xlab = "",
       ylab = "",
       xlim = c(xmin, xmax),
       bty="l")
  
  
  #--------------------------------------------------------#
  # Panel (i)
  plot(seq(1900,xmax,1),wg,type="l",
       main = 'Wages',
       xlab = "",
       ylab = "",
       xlim = c(xmin, xmax),
       bty="l")
  
  #--------------------------------------------------------#
  # Panel (j)
  plot(seq(1900,xmax,1),rwg,type="l",
       main = 'Real Wages',
       xlab = "",
       ylab = "",
       xlim = c(xmin, xmax),
       bty="l") 
  
  #--------------------------------------------------------#
  # Panel (k)
  plot(seq(1889,xmax,1),m,type="l",
       main = 'Money',
       xlab = "",
       ylab = "",
       xlim = c(xmin, xmax),
       bty="l") 
  
  #--------------------------------------------------------#
  # Panel (l)
  plot(seq(1869,xmax,1),vel,type="l",
       main = 'Velocity',
       xlab = "",
       ylab = "",
       xlim = c(xmin, xmax),
       bty="l")
  
  #--------------------------------------------------------#
  # Panel (m)
  plot(seq(1900,xmax,1),bnd,type="l",
       main = 'Bond Yield',
       xlab = "",
       ylab = "",
       xlim = c(xmin, xmax),
       bty="l")
  
  #--------------------------------------------------------#
  # Panel (n)
  plot(seq(1871,xmax,1),sp500,type="l",
       main = 'S & P 500',
       xlab = "",
       ylab = "",
       xlim = c(xmin, xmax),
       bty="l")
  
  #--------------------------------------------------------#
  # Panel (o)
  plot(seq(xmin,xmax,1),y_dtrend,type="l",
       main = 'Deterministic Trend',
       xlab = "",
       ylab = "",
       xlim = c(xmin, xmax),
       bty="l")
  
  #--------------------------------------------------------#
  # Panel (p)
  plot(seq(xmin,xmax,1),y_strend,type="l",
       main = 'Stochastic Trend',
       xlab = "",
       ylab = "",
       xlim = c(xmin, xmax),
       bty="l")  
  
  # Compute autocorrelations for the first 6 lags
  lags <- 6
  
  # Nelson and Plosser (1982) Table 2, p.147  
  cat('\nAutocorrelation functions')
  cat('\nNominal GNP')
  print.acf(gnp, lags)
  
  cat('\n')
  cat('\nPer capita GNP')
  print.acf( pcrgnp,lags)
  cat('\n')
  cat('\nIndustrial production')
  print.acf(ip,lags)
  cat('\n')
  cat('\nUnemployment')
  print.acf(emp,lags)
  cat('\n')
  cat('\nUnemployment rate')
  print.acf(un,lags)
  cat('\n')
  cat('\nGNP deflator')
  print.acf( prgnp,lags)
  cat('\n')
  cat('\nConsumer price index')
  print.acf(cpi,lags)
  cat('\n')
  cat('\nWages')
  print.acf(wg,lags)
  cat('\n')
  cat('\nReal wages')
  print.acf( rwg,lags)
  cat('\n')
  cat('\nMoney stock')
  print.acf( m,lags) 
  cat('\n')
  cat('\nVelocity')
  print.acf(vel,lags)
  cat('\n')
  cat('\nBond yield')
  print.acf( bnd,lags)
  cat('\n')
  cat('\nStock price')
  print.acf( sp500,lags)
  cat('\n')
  cat('\nDeterministic trend')
  print.acf(y_dtrend,lags)
  cat('\n')
  cat('\nStochastic trend')
  print.acf( y_strend,lags) 
  cat('\n')
  
  # Estimate the stochastic trend model 
  # an AR(1) model with a constant and no time tend
  
  ar.obj <- ar(rgnp,0)
  cat('\nReal GNP')
  cat('\nconstant         ',ar.obj$b[1],'  (',ar.obj$se[1],')')
  cat('\nAR1 parameter    ',ar.obj$b[2],'  (',ar.obj$se[2],')')
  
  ar.obj <- ar(gnp,0)
  cat('\n ')
  cat('\nNominal GNP')
  cat('\nconstant         ',ar.obj$b[1],'  (',ar.obj$se[1],')')
  cat('\nAR1 parameter    ',ar.obj$b[2],'  (',ar.obj$se[2],')')
  ar.obj <- ar(pcrgnp,0)
  cat('\n ')
  cat('\nPer capita GNP')
  cat('\nconstant         ',ar.obj$b[1],'  (',ar.obj$se[1],')')
  cat('\nAR1 parameter    ',ar.obj$b[2],'  (',ar.obj$se[2],')')
  ar.obj <- ar(ip,0)
  cat('\n ')
  cat('\nIndustrial production')
  cat('\nconstant         ',ar.obj$b[1],'  (',ar.obj$se[1],')')
  cat('\nAR1 parameter    ',ar.obj$b[2],'  (',ar.obj$se[2],')')
  ar.obj <- ar(emp,0)
  cat('\n ')
  cat('\nEmployment')
  cat('\nconstant         ',ar.obj$b[1],'  (',ar.obj$se[1],')')
  cat('\nAR1 parameter    ',ar.obj$b[2],'  (',ar.obj$se[2],')')
  ar.obj <- ar(un,0)
  cat('\n ')
  cat('\nUnemployment rate')
  cat('\nconstant         ',ar.obj$b[1],'  (',ar.obj$se[1],')')
  cat('\nAR1 parameter    ',ar.obj$b[2],'  (',ar.obj$se[2],')')
  ar.obj <- ar(prgnp,0)
  cat('\n ')
  cat('\nGNP deflator')
  cat('\nconstant         ',ar.obj$b[1],'  (',ar.obj$se[1],')')
  cat('\nAR1 parameter    ',ar.obj$b[2],'  (',ar.obj$se[2],')')
  ar.obj <- ar(cpi,0)
  cat('\n ')
  cat('\nConsumer price index')
  cat('\nconstant         ',ar.obj$b[1],'  (',ar.obj$se[1],')')
  cat('\nAR1 parameter    ',ar.obj$b[2],'  (',ar.obj$se[2],')')
  ar.obj <- ar(wg,0)
  cat('\n ')
  cat('\nWages')
  cat('\nconstant         ',ar.obj$b[1],'  (',ar.obj$se[1],')')
  cat('\nAR1 parameter    ',ar.obj$b[2],'  (',ar.obj$se[2],')')
  ar.obj <- ar(rwg,0)
  cat('\n ')
  cat('\nReal wages')
  cat('\nconstant         ',ar.obj$b[1],'  (',ar.obj$se[1],')')
  cat('\nAR1 parameter    ',ar.obj$b[2],'  (',ar.obj$se[2],')')
  ar.obj <- ar(m,0)
  cat('\n ')
  cat('\nMoney stock')
  cat('\nconstant         ',ar.obj$b[1],'  (',ar.obj$se[1],')')
  cat('\nAR1 parameter    ',ar.obj$b[2],'  (',ar.obj$se[2],')')
  ar.obj <- ar(vel,0)
  cat('\n ')
  cat('\nVelocity')
  cat('\nconstant         ',ar.obj$b[1],'  (',ar.obj$se[1],')')
  cat('\nAR1 parameter    ',ar.obj$b[2],'  (',ar.obj$se[2],')')
  ar.obj <- ar(bnd,0)
  cat('\n ')
  cat('\nBond yield')
  cat('\nconstant         ',ar.obj$b[1],'  (',ar.obj$se[1],')')
  cat('\nAR1 parameter    ',ar.obj$b[2],'  (',ar.obj$se[2],')')
  ar.obj <- ar(sp500,0)
  cat('\n ')
  cat('\nStock price')
  cat('\nconstant         ',ar.obj$b[1],'  (',ar.obj$se[1],')')
  cat('\nAR1 parameter    ',ar.obj$b[2],'  (',ar.obj$se[2],')')
  ar.obj <- ar(y_dtrend,0)
  cat('\n ')
  cat('\nDeterministic trend')
  cat('\nconstant         ',ar.obj$b[1],'  (',ar.obj$se[1],')')
  cat('\nAR1 parameter    ',ar.obj$b[2],'  (',ar.obj$se[2],')')
  ar.obj <- ar(y_strend,0)
  cat('\n ')
  cat('\nStochastic trend')
  cat('\nconstant         ',ar.obj$b[1],'  (',ar.obj$se[1],')')
  cat('\nAR1 parameter    ',ar.obj$b[2],'  (',ar.obj$se[2],')')
  
  # Estimate the combined model: 
  # an AR(1) model with a constant and a time tend 
  ar.obj <- ar(rgnp,1)
  cat('\nReal GNP')
  cat('\nconstant         ',ar.obj$b[1],'  (',ar.obj$se[1],')')
  cat('\ntrend            ',ar.obj$b[2],'  (',ar.obj$se[2],')')    
  cat('\nAR1 parameter    ',ar.obj$b[3],'  (',ar.obj$se[3],')') 
  ar.obj <- ar(gnp,1)
  cat('\n ')
  cat('\nNominal GNP')
  cat('\nconstant         ',ar.obj$b[1],'  (',ar.obj$se[1],')')
  cat('\ntrend            ',ar.obj$b[2],'  (',ar.obj$se[2],')')    
  cat('\nAR1 parameter    ',ar.obj$b[3],'  (',ar.obj$se[3],')')
  ar.obj <- ar(pcrgnp,1)
  cat('\n ')
  cat('\nPer capita GNP')
  cat('\nconstant         ',ar.obj$b[1],'  (',ar.obj$se[1],')')
  cat('\ntrend            ',ar.obj$b[2],'  (',ar.obj$se[2],')')    
  cat('\nAR1 parameter    ',ar.obj$b[3],'  (',ar.obj$se[3],')')
  ar.obj <- ar(ip,1)
  cat('\n ')
  cat('\nIndustrial production')
  cat('\nconstant         ',ar.obj$b[1],'  (',ar.obj$se[1],')')
  cat('\ntrend            ',ar.obj$b[2],'  (',ar.obj$se[2],')')    
  cat('\nAR1 parameter    ',ar.obj$b[3],'  (',ar.obj$se[3],')')
  ar.obj <- ar(emp,1)
  cat('\n ')
  cat('\nEmployment')
  cat('\nconstant         ',ar.obj$b[1],'  (',ar.obj$se[1],')')
  cat('\ntrend            ',ar.obj$b[2],'  (',ar.obj$se[2],')')    
  cat('\nAR1 parameter    ',ar.obj$b[3],'  (',ar.obj$se[3],')')
  ar.obj <- ar(un,1)
  cat('\n ')
  cat('\nUnemployment rate')
  cat('\nconstant         ',ar.obj$b[1],'  (',ar.obj$se[1],')')
  cat('\ntrend            ',ar.obj$b[2],'  (',ar.obj$se[2],')')    
  cat('\nAR1 parameter    ',ar.obj$b[3],'  (',ar.obj$se[3],')')
  ar.obj <- ar(prgnp,1)
  cat('\n ')
  cat('\nGNP deflator')
  cat('\nconstant         ',ar.obj$b[1],'  (',ar.obj$se[1],')')
  cat('\ntrend            ',ar.obj$b[2],'  (',ar.obj$se[2],')')    
  cat('\nAR1 parameter    ',ar.obj$b[3],'  (',ar.obj$se[3],')')
  ar.obj <- ar(cpi,1)
  cat('\n ')
  cat('\nConsumer price index')
  cat('\nconstant         ',ar.obj$b[1],'  (',ar.obj$se[1],')')
  cat('\ntrend            ',ar.obj$b[2],'  (',ar.obj$se[2],')')    
  cat('\nAR1 parameter    ',ar.obj$b[3],'  (',ar.obj$se[3],')')
  ar.obj <- ar(wg,1)
  cat('\n ')
  cat('\nWages')
  cat('\nconstant         ',ar.obj$b[1],'  (',ar.obj$se[1],')')
  cat('\ntrend            ',ar.obj$b[2],'  (',ar.obj$se[2],')')    
  cat('\nAR1 parameter    ',ar.obj$b[3],'  (',ar.obj$se[3],')')
  ar.obj <- ar(rwg,1)
  cat('\n ')
  cat('\nReal wages')
  cat('\nconstant         ',ar.obj$b[1],'  (',ar.obj$se[1],')')
  cat('\ntrend            ',ar.obj$b[2],'  (',ar.obj$se[2],')')    
  cat('\nAR1 parameter    ',ar.obj$b[3],'  (',ar.obj$se[3],')')
  ar.obj <- ar(m,1)
  cat('\n ')
  cat('\nMoney stock')
  cat('\nconstant         ',ar.obj$b[1],'  (',ar.obj$se[1],')')
  cat('\ntrend            ',ar.obj$b[2],'  (',ar.obj$se[2],')')    
  cat('\nAR1 parameter    ',ar.obj$b[3],'  (',ar.obj$se[3],')')
  ar.obj <- ar(vel,1)
  cat('\n ')
  cat('\nVelocity')
  cat('\nconstant         ',ar.obj$b[1],'  (',ar.obj$se[1],')')
  cat('\ntrend            ',ar.obj$b[2],'  (',ar.obj$se[2],')')    
  cat('\nAR1 parameter    ',ar.obj$b[3],'  (',ar.obj$se[3],')')
  ar.obj <- ar(bnd,1)
  cat('\n ')
  cat('\nBond yield')
  cat('\nconstant         ',ar.obj$b[1],'  (',ar.obj$se[1],')')
  cat('\ntrend            ',ar.obj$b[2],'  (',ar.obj$se[2],')')    
  cat('\nAR1 parameter    ',ar.obj$b[3],'  (',ar.obj$se[3],')')
  ar.obj <- ar(sp500,1)
  cat('\n ')
  cat('\nStock price')
  cat('\nconstant         ',ar.obj$b[1],'  (',ar.obj$se[1],')')
  cat('\ntrend            ',ar.obj$b[2],'  (',ar.obj$se[2],')')    
  cat('\nAR1 parameter    ',ar.obj$b[3],'  (',ar.obj$se[3],')')
  ar.obj <- ar(y_dtrend,1)
  cat('\n ')
  cat('\nDeterministic trend')
  cat('\nconstant         ',ar.obj$b[1],'  (',ar.obj$se[1],')')
  cat('\ntrend            ',ar.obj$b[2],'  (',ar.obj$se[2],')')    
  cat('\nAR1 parameter    ',ar.obj$b[3],'  (',ar.obj$se[3],')')
  ar.obj <- ar(y_strend,1)
  cat('\n ')
  cat('\nStochastic trend')
  cat('\nconstant         ',ar.obj$b[1],'  (',ar.obj$se[1],')')
  cat('\ntrend            ',ar.obj$b[2],'  (',ar.obj$se[2],')')    
  cat('\nAR1 parameter    ',ar.obj$b[3],'  (',ar.obj$se[3],')')
  
  # Estimate the trend model: 
  # a model with a constant and a time tend 
  ar.obj <- ar(rgnp,2)
  cat('\nReal GNP')
  cat('\nconstant         ',ar.obj$b[1],'  (',ar.obj$se[1],')')
  cat('\ntrend            ',ar.obj$b[2],'  (',ar.obj$se[2],')')       
  ar.obj <- ar(gnp,2)
  cat('\n ')
  cat('\nNominal GNP')
  cat('\nconstant         ',ar.obj$b[1],'  (',ar.obj$se[1],')')
  cat('\ntrend            ',ar.obj$b[2],'  (',ar.obj$se[2],')')    
  ar.obj <- ar(pcrgnp,2)
  cat('\n ')
  cat('\nPer capita GNP')
  cat('\nconstant         ',ar.obj$b[1],'  (',ar.obj$se[1],')')
  cat('\ntrend            ',ar.obj$b[2],'  (',ar.obj$se[2],')')    
  ar.obj <- ar(ip,2)
  cat('\n ')
  cat('\nIndustrial production')
  cat('\nconstant         ',ar.obj$b[1],'  (',ar.obj$se[1],')')
  cat('\ntrend            ',ar.obj$b[2],'  (',ar.obj$se[2],')')    
  ar.obj <- ar(emp,2)
  cat('\n ')
  cat('\nEmployment')
  cat('\nconstant         ',ar.obj$b[1],'  (',ar.obj$se[1],')')
  cat('\ntrend            ',ar.obj$b[2],'  (',ar.obj$se[2],')')    
  ar.obj <- ar(un,2)
  cat('\n ')
  cat('\nUnemployment rate')
  cat('\nconstant         ',ar.obj$b[1],'  (',ar.obj$se[1],')')
  cat('\ntrend            ',ar.obj$b[2],'  (',ar.obj$se[2],')')    
  ar.obj <- ar(prgnp,2)
  cat('\n ')
  cat('\nGNP deflator')
  cat('\nconstant         ',ar.obj$b[1],'  (',ar.obj$se[1],')')
  cat('\ntrend            ',ar.obj$b[2],'  (',ar.obj$se[2],')')    
  ar.obj <- ar(cpi,2)
  cat('\n ')
  cat('\nConsumer price index')
  cat('\nconstant         ',ar.obj$b[1],'  (',ar.obj$se[1],')')
  cat('\ntrend            ',ar.obj$b[2],'  (',ar.obj$se[2],')')    
  ar.obj <- ar(wg,2)
  cat('\n ')
  cat('\nWages')
  cat('\nconstant         ',ar.obj$b[1],'  (',ar.obj$se[1],')')
  cat('\ntrend            ',ar.obj$b[2],'  (',ar.obj$se[2],')')    
  ar.obj <- ar(rwg,2)
  cat('\n ')
  cat('\nReal wages')
  cat('\nconstant         ',ar.obj$b[1],'  (',ar.obj$se[1],')')
  cat('\ntrend            ',ar.obj$b[2],'  (',ar.obj$se[2],')')    
  ar.obj <- ar(m,2)
  cat('\n ')
  cat('\nMoney stock')
  cat('\nconstant         ',ar.obj$b[1],'  (',ar.obj$se[1],')')
  cat('\ntrend            ',ar.obj$b[2],'  (',ar.obj$se[2],')')    
  ar.obj <- ar(vel,2)
  cat('\n ')
  cat('\nVelocity')
  cat('\nconstant         ',ar.obj$b[1],'  (',ar.obj$se[1],')')
  cat('\ntrend            ',ar.obj$b[2],'  (',ar.obj$se[2],')')    
  ar.obj <- ar(bnd,2)
  cat('\n ')
  cat('\nBond yield')
  cat('\nconstant         ',ar.obj$b[1],'  (',ar.obj$se[1],')')
  cat('\ntrend            ',ar.obj$b[2],'  (',ar.obj$se[2],')')    
  ar.obj <- ar(sp500,2)
  cat('\n ')
  cat('\nStock price')
  cat('\nconstant         ',ar.obj$b[1],'  (',ar.obj$se[1],')')
  cat('\ntrend            ',ar.obj$b[2],'  (',ar.obj$se[2],')')    
  ar.obj <- ar(y_dtrend,2)
  cat('\n ')
  cat('\nDeterministic trend')
  cat('\nconstant         ',ar.obj$b[1],'  (',ar.obj$se[1],')')
  cat('\ntrend            ',ar.obj$b[2],'  (',ar.obj$se[2],')')    
  ar.obj <- ar(y_strend,2)
  cat('\n ')
  cat('\nStochastic trend')
  cat('\nconstant         ',ar.obj$b[1],'  (',ar.obj$se[1],')')
  cat('\ntrend            ',ar.obj$b[2],'  (',ar.obj$se[2],')')    
  
  
}
