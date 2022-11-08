#============================================================================
#
#   Unit root tests: Nelson-Plosser data (1860 to 1970).
#
#============================================================================          
rm (list = ls(all=TRUE))
graphics.off()

#
#------------------------- Helper Functions----------------------------------
#

# load required functions -  trimr
source("EMTSUtil.R")

#-----------------------------------------------------------------------------
#   Compute acf and Mtests for Nelson Plosser data  
#-----------------------------------------------------------------------------
np <- function(y,lags,flag) {
  t <- length(y)
  
  if (!flag) {
    # acf of y
    acy <- acf(y,lags, plot=F)$acf[-1]
    
    # acf acf of differences
    acdy <- acf(trimr(y,1,0)-trimr(y,0,1),lags, plot=F)$acf 
    
    # acf of residuals from a linear trend
    x    <- cbind(rep(1, length(y)), seqa(1,1,t)) 
    uols <- glsdetrend(y,x,-t)
    acu  <- acf(uols,lags, plot=F)$acf[-1]
    # MZt tests
    cbar   <- -13.5 
    k      <- kmaic(uols)
    ugls   <- glsdetrend(y,x,cbar)
    
    Mztols <- trimr(cbind(Mtests(uols,k)),2,0)
    Mztgls <- trimr(cbind(Mtests(ugls,k)),2,0)
  } else { # Unemployment flagged
    # acf of y
    acy <- acf(y,lags, plot=F)$acf[-1]
    
    # acf acf of differences
    acdy <- acf(trimr(y,1,0)-trimr(y,0,1),lags, plot=F)$acf[-1] 
    
    # acf of residuals from a constant only
    x    <- matrix(1, nrow=length(y))
    uols <- glsdetrend(y,x,-t)
    acu  <- acf(uols,lags, plot=F)$acf[-1]
    
    # MZt tests
    cbar   <- -7 
    k      <- kmaic(uols)
    ugls   <- glsdetrend(y,x,cbar)
    
    Mztols <- trimr(cbind(Mtests(uols,k)),2,0)
    Mztgls <- trimr(cbind(Mtests(ugls,k)),2,0)
  }   
  return (list(acy=acy, acdy=acdy,acu=acu,Mztols=Mztols,Mztgls=Mztgls))
}

#----------------------------------------------------------------------------
#  Detrending function: 
#       cbar = -7 constant 
#       cbar = -13.5 linear trend  
#       cbar = -T for OLS detrending
#----------------------------------------------------------------------------

glsdetrend <- function(y,x,cbar) {
  
  t <- length(y)
  yc <- cbind(c(y[1], (trimr(y,1,0)-(1+cbar/t)*trimr(y,0,1))) )
  xc <- rbind(x[1,], cbind(trimr(x,1,0)-(1+cbar/t)*trimr(x,0,1)))    
  b <- lm(yc ~ xc - 1)$coef
  u <- y - x %*% b  
  return(u)
}


##---------------------------------------------------------------------------- 
# Autoregressive long run variance estimator 
#---------------------------------------------------------------------------- 
ar1rvar <- function(u,k)  {
  
  du <- cbind(trimr(u,1,0)-trimr(u,0,1) )
  x <- trimr(u,0,1) 
  
  if (k > 0)  { 
    x <- cbind(x,lag.matrix(du,1:k)) 
    x <- cbind(x[-(1:k), ])
  }
  x <- cbind(x)  
  b <- lm(trimr(du,k,0) ~ x - 1)$coef
  
  e <- trimr(du,k,0)- x %*% b 
  s2 <- t(e) %*% e/ length(e)  
  if (k > 0)  { 
    s2 <- s2/ (1-sum(trimr(b,1,0)))^2
  }   
  return(s2)
} 


#----------------------------------------------------------------------------
#  M tests
#----------------------------------------------------------------------------
Mtests <- function( u,k ){
  
  t  = length(u)
  s2 = ar1rvar( u,k )
  u2 = sum( u[1:(t-1)]^2/t^2 )
  
  cat('\nSum of u(1) to u(t-1) = ',sum(u[1:t-1]^2))
  cat('\nLast value: u(t)      = ',u[t])
  cat('\n ')   
  
  mza = (u[t]^2/t-s2)/(2*u2)
  msb = sqrt(u2/s2)
  mzt = msb*mza
  
  mz = c(mza,  msb, mzt)
  return(mz)
}

#-------------------------------------------------------------------------
#  Select ADF lag length by MAIC: u should be OLS de-trended residuals
#-------------------------------------------------------------------------
kmaic <- function(u) {
  kmax <- floor(12*(length(u)/100)^0.25) 
  maic <- rep(0, kmax+1)
  
  # Set up lagged regressors
  du <- cbind(trimr(u,1,0)-trimr(u,0,1))
  x  <- cbind(trimr(u,0,1), lag.matrix(du,1:kmax))
  x <- cbind(x[-(1:kmax), ])
  
  
  for (j in 0:kmax) {
    b <- lm(trimr(du,kmax,0) ~ x[,1:(1+j)] - 1)$coef
    e <- trimr(du,kmax,0)- as.matrix(x[,1:(1+j)]) %*% b 
    s2 <- t(e) %*% e/length(e)
    
    maic[j+1] <- log(s2) + 2*(j+b[1]^2*sum(x[,1]^2)/s2)/length(e)    
  }
  k <- which.min(maic)
  k     <- k-1
  return(k)
}

#
#------------------------- Nelson-Plosser ----------------------------------
#
unit_nelplos <- function() {
  # Load Nelson-Plosser data set
  data <- as.matrix(read.table("nelson_plosser.dat"))
  
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
  
  # Choose data and lag length and flag unemployment)
  y    <- un
  flag <- 1
  lags <- 6
  
  np.data <- np(y,lags,flag)
  acy <- np.data$acy
  acdy <- np.data$acdy
  acu <- np.data$acu
  Mztols <- np.data$Mztols
  Mztgls <- np.data$Mztgls
  
  
  cat('\n         Autocorrelation Functions\n')
  print( cbind(Lag=seqa(1,1,lags), Levels=acy, "1st Diff"=acdy, Detrend=acu))
  
  cat('\n     Unit Root Tests\n')
  print( cbind(Mztols, Mztgls) )
  
  # Union of rejections
  if (!flag) {
    cvgls <- -2.867 
    cvols <- -3.130 
    tau <- 1.038
  } else {
    cvgls <- -2.867 
    cvols <- -3.130 
    tau <- 1.038
  }
  cat('\nUOR Mztols classification is I(',Mztols > tau*cvols,')')
  cat('\nUOR Mztgls classification is I(',Mztgls > tau*cvols,')')  
  
}

