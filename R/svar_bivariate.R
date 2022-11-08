#==========================================================================
#
#     Program to estimate a bivariate svar model and demonstrate
#   	alternative restrictions
#
#==========================================================================

rm (list = ls(all=TRUE))
graphics.off()
set.seed(123, kind="Mersenne-Twister")

#
#------------------------- Helper Functions----------------------------------
#

#load required functions -  figure, trimr, inv, reshapeg
source("EMTSUtil.R")

#
#-------------------- SVAR Bivariate Model ----------------------------------
#

svar_bivariate <- function() {
    # Load data
  # Data from 1970Q1 to 2002Q2    
  #1. Oil price 
	#2. Output EMU
	#3. CPI EMU
	#4. Interest rate EMU
	#5. Output US
	#6. CPI US
  #7. Interest rate US
  
  load('peersman_data.Rdata')
      
  # Choose identification method
  # 0 = non-orthogonal shocks
  # 1 = orthogonal shocks with short-run restrictions
  # 2 = orthogonal shocks with long-run restrictions
  # 3 = orthogonal shocks with short- and long-run restrictions
  # 4 = orthogonal shocks with sign restrictions
  itype <- 4      

  # Data transformations
  loil <- log( ytdata[,1] )     ##ok<NASGU,NODEF> # log price of oil    
  lo   <- log( ytdata[,5] )     # log output          
  lp   <- log( ytdata[,6] )     # log price level     
  r    <- ytdata[,7]            ##ok<NASGU> # interest rate       

  # Construct variables for use in VAR: 
  # # growth rates in oil, output and price, and level of interest rate
 
  yact <- cbind(lo, lp)
  yact <- trimr(yact,36,0)       # Sample period 1970Q1 to 2002Q2

  y <- 400*( trimr(yact,1,0) - trimr(yact,0,1) ) 

  # Set SVAR parameters
  p <- 3      # Order of VAR lags
  q <- 40     # Order of VMA lags      

  # Estimate the VAR with p lags and a constant     
  ylag   <- rep(1, nrow(y))
  nconst <- 1

  for (i in seq(p)) {    
    ylag <- cbind(trimr(ylag,1,0), trimr(y,0,i))
  }
  
  # OLS Estimates
  bar    <- lm(trimr(y,p,0) ~ ylag - 1)$coef
  v      <- trimr(y,p,0) - ylag %*% bar
  omegav <- t(v) %*% v/nrow(v)

  cat('\nCovariance matrix of VAR residuals\n')
  print(unname(omegav))

  # Constuct A(L) (VAR) matrices and hence the C(1) long-run matrix
  bar <- trimr(bar,nconst,0)
  k   <- ncol(v)
  a   <- array(0, c(k^2,p))
  a1  <- diag(k)  # I
  
  for (i in seq(p)) {
    tmp    <- bar[(1+k*(i-1)):(i*k),]
    a[,i] <- c(tmp)
	  a1   <- a1 - matrix(a[,i],k,k)    
  }

  # Invert A(1) matrix needed for the MLE
  a1inv <- inv(a1)            
  
  # Identification
  if (itype == 0) {       # Non-orthogonal one unit impulse responses 
    s <- matrix(c(1,0, 0,1), nrow=2, byrow=T) 
  } else if(itype == 1) { # Orthogonal one-standard deviation impulse response based on short-run restrictions
    s <- matrix(c(3, 0, -1, 2), nrow=2, byrow=T)    
  }else if (itype == 2) { #   # Orthogonal one-standard deviation impulse response based on long-run restrictions
    f <- matrix(c(3, 0, -1, 2), nrow=2, byrow=T)
    s <- a1 %*% f
  }else if (itype == 3) {  # Orthogonal impulse response based on short-run and long-run restrictions
    s11 <- 2 
    s12 <- 1 
    s21 <- 0 
    s22 <- (-a1inv[1,1]/a1inv[1,2])*s12 
    s   <- matrix(c(s11,  s12, 
                    s21,  s22), nrow=2, byrow=T)
  }else if (itype == 4) { # Orthogonal impulse response based on sign short-run restrictions 
    s <- matrix(c(3, 0,  
                  -1, 2), nrow=2, byrow=T)
    th <- pi*0.4
    qmat <- matrix(c(cos(th),  -sin(th),
                   sin(th),   cos(th)), nrow=2, byrow=T)

    sq <- s %*% t(qmat)
    cat('\n------------------------')
    cat('\ns\n' )
    print( s )
    cat('\n------------------------')
    cat('\nqmat\n')
    print(qmat)
    cat('\n------------------------')
    cat('\nsq\n')
    print( sq )
    s <- sq
    
  }    
  cat('\n------------------------')
  cat('\na1\n' )
  print( a1 )
  cat('\n------------------------')
  cat('\na1inv\n' )
  print( a1inv )
  cat('\n------------------------')
  cat('\ns\n' )
  print( s )
  cat('\n------------------------')

  # Construct C(L) matrices (vector moving average)   
  c <- cbind(c(diag(k)))
  
  for (i in seq(q)) {    
    ss <- array(0, c(k,k))
    j   <- 1.0
    
    while (j <= min(c(p, i))) {      
      ss <- ss + reshapeg(a[,j],k,k) %*% reshapeg(c[,i-j+1],k,k)
      j   <- j + 1
    }
    tmp <- t(ss)
    c  <- cbind(c,  c(tmp))
  }      
  
   # Construct orthogonal impulse response functions
  impulse <- c(t(s))  

  for (i in 2:(q+1)) {     
    tmp      <- t( (reshapeg( c[,i],k,k ) %*% s ) )
    impulse  <- cbind(impulse, c(tmp))    
  }
  
	impulse <- t(impulse)
 
  # Compute cumulative sum of impulse responses for levels   
  impulse <- apply(impulse, 2, cumsum)

  #**********************************************************************
  #***
  #***     Generate graph
  #***
  #**********************************************************************
  
  figure()
  par(mfrow=c(2,2))
  t <- seq(0, q, 1)
  
  #--------------------------------------------------------#
  # Panel (a)
  if (itype == 0) {
    title <- expression(paste("v",""[1]*",t", " shock"))    
  }else {
    title <- expression(paste("z",""[1]*",t", " shock"))    
  }
  plot(t,impulse[1:(q+1),1],type="l",
          main = title,
          ylab = expression(paste(y[1],",t")),
          xlab = "Quarter",
          bty="l")
  
  #--------------------------------------------------------#
  # Panel (b)
  if (itype == 0) {
    title <- expression(paste("v",""[2]*",t", " shock"))    
  }else {
    title <- expression(paste("z",""[2]*",t", " shock"))    
  }
  plot(t,impulse[1:(q+1),2],type="l",
          main = title,
          ylab = expression(paste(y[1],",t")),
          xlab = "Quarter",
          bty="l")


  #--------------------------------------------------------#
  # Panel (c)
  plot(t,impulse[1:(q+1),3],type="l",
          main = title,
          ylab = expression(paste(y[2],",t")),
          xlab = "Quarter",
          bty="l")


  #--------------------------------------------------------#
  # Panel (d)    
  plot(t,impulse[1:(q+1),4],type="l",
          main = title,
          ylab = expression(paste(y[2],",t")),
          xlab = "Quarter",
          bty="l")
}

