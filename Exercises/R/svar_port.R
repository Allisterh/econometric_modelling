#============================================================================
#
#     Program to estimate a Structural VAR for the share example 
#       using Blanchard-Quah long-run restrictions
#		with present value long-run restriction.
#
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()

#
#------------------------- Helper Functions----------------------------------
#

# Load required functions -  trimr, reshapeg, figure
source("EMTSUtil.R")

# Load required library - repmat
library("matlab")

#--------------------------------------------------------------------------
# Log-likelihood function for restricted SVAR
#--------------------------------------------------------------------------
neglogl <- function (b,v,a1 ) {
  t <- nrow(v)
  n <- ncol(v)
  lf      <- rep(0, t)
  f   <-     matrix(c(b[1],   0,       0,      0,     0,
                      b[2],   b[3],    b[4],   0,  0,
                      b[5],  -b[3],   b[6],   0,  b[12],
                      b[7],   b[8],    b[9],  b[10],  0,
                      0,     0,       0,       0,      b[11]), nrow=n, byrow=T)
  s      <- a1 %*% f		# Long-run restrictions
	omegav <- s %*% t(s)  
  for (i in seq(t)) {
    lf[i] <- -0.5*n*log(2*pi) - 0.5*log(det(omegav)) - 0.5*v[i,] %*% inv(omegav) %*% cbind(v[i,])    
  }
  lnl <- -mean(lf)  
  return(lnl)
}

#--------------------------------------------------------------------------
# Log-likelihood function for unrestricted SVAR
#--------------------------------------------------------------------------
neglogu <- function( b,v,a1 ) {
  t <- nrow(v)
  n <- ncol(v)
  lf      <- rep(0, t)
  f   <-     matrix(c(b[1],   0,       0,      0,    0,
                      b[2],   b[3],   b[4],   0,     0,
                      b[5],   b[13],  b[6],   0,     b[12],
                      b[7],   b[8],   b[9],   b[10], 0,
                      0,     0,       0,       0,      b[11]), nrow=n, byrow=T)
  s      <- a1 %*% f  	# Long-run restrictions
	omegav <- s %*% t(s)
  for (i in seq(t)) {
    lf[i] <- -0.5*n*log(2*pi) - 0.5*log(det(omegav)) - 0.5*v[i,] %*% inv(omegav) %*% cbind(v[i,])    
  }
  lnl <- -mean(lf)  
  return(lnl)
}

#
#------------------------------- A Portfolio SVAR Model----------------------
#
svar_port <- function () {
    # Load data
  # Data from 1980:1 to 2002:2, taken from the RBA Bulletin database    
  #1. US Share price index, average over the quarter, S&P500
	#2. US/AU exchange rate, average over the quarter
	#3. CPI, all groups, 1989-1990 =100, seasonally adjusted
	#4. Nominal GDP, expenditure, total, s.a.
	#5. Real GDP, expenditure, total, 2000/20001 $m, s.a.
	#6. Australian share price index, average over the quarter, S&P/ASX 200
  #7. 90 day bank accepted bill interest rate, average over the quarter

  load('portfolio_data.Rdata')
  ytdata <- as.matrix(ytdata)
  # Data transformations
  ipd   <- ytdata[,4]/ytdata[,5]		# Implicit price deflator
  r3    <- ytdata[,7]/100
  lipd  <- log(ipd)
  lcpi  <- log(ytdata[,3])
  lrao  <- log(ytdata[,6]/ipd)
  lr3   <- log(r3)
  lrgdp <- log(ytdata[,5])
  lrdj  <- log(ytdata[,1]/(ytdata[,2]*ipd))	# US share price converted into Australian dollars and then expressed in reals 
  inf   <- 1*(trimr(lipd,4,0) - trimr(lipd,0,4))

  yact  <- 100*( cbind(lrgdp, lr3, lrao, lcpi, lrdj))
  y     <- trimr(yact,1,0) - trimr(yact,0,1)
  # Create dummy variables
  d87        <- rep(0, 102)
  d87[32]    <- 1.0			# Stock market crash dummy variable
  d2000      <- rep(0,102)
  d2000[83]  <- 1.0           # GST introduction in 2000:3c

  # SVAR Parameters
  p <- 2		# VAR order
  q <- 40		# VMA order

  # Estimate the VAR with p lags
  ylag     <- cbind(rep(1, nrow(y)), trimr( cbind(d87, d2000),1,0))
  nconst   <- ncol(ylag)

  for (i in seq(p)) {
    ylag <- cbind(trimr(ylag,1,0), trimr(y,0,i))
  }

  reg <- lm(trimr(y,p,0) ~ ylag - 1)
  bar   <- reg$coef
  mu    <- bar[1:nconst,]
  v     <- reg$residuals
  omega <- t(v) %*% v/nrow(v)

  # Constuct A(L) (VAR) matrices and hence the C(1) long-run matrix
  bar <- trimr(bar,nconst,0)
  k   <- ncol(v)
  a   <- array(0, c(k^2,p))
  a1  <- diag(k)
  for (i in seq(p)) {
    tmp    <- bar[(1+k*(i-1)):(i*k),]
    a[,i] <- as.vector(tmp)
    a1     <- a1 - reshapeg(a[,i],k,k)
  }

  bstart <- c(1.0389699473962597,
                4.3007505638494345, 
                5.1966616592138655, 
                10.0996788841097480, 
                2.7076795082063354, 
                1.4379097951392796, 
                -0.7123778562836982, 
                0.7750813995549357, 
                1.6190786936559989, 
                1.3494397321299636,
                10.7524502796372890,
                3.5142132747486432, 
                -5.1782573061530375)
  estResults <- optim(bstart, neglogu, v=v, a1=a1, method="BFGS")
  b <- estResults$par
  fvalu <- -estResults$val

  cat('\nUnrestricted Log-likelihood function = ', fvalu)


  # Estimate restritced SVAR by maximum likelihood 
  bstart <- c(1.038918123051897,
                4.301175420613745, 
                5.180760075404540, 
                10.10913871549879, 
                2.707544494664958, 
                1.428884030946225, 
               -0.7120443107363003, 
                0.7725526070345229, 
                1.620606390960305, 
                1.349442082771462, 
                10.75208272049359, 
                3.513290773831895)  

  estResults <- optim(bstart, neglogl, v=v, a1=a1, method="BFGS")
  b <- estResults$par
  fvalr <- -estResults$val
  cat('\nRestricted Log-likelihood function = ',fvalr )
    
  # Likelihood ratio test
  lr <- -2*nrow(v)*(fvalr-fvalu)
  cat('\nLR test        = ',lr)
  cat('\np-value        = ',1-pchisq(lr,1))
  
  # Construct C(L) matrices (vector moving average) 
  cl <- diag(k)
  cl <- cbind(as.vector(cl))

  for (i in seq(q)) {
    ss <- array(0, c(k,k))    
    j   <- 1.0    
    while (j <= min(c(p,i))) {
      ss <- ss + reshapeg(a[,j],k,k) %*% reshapeg(cl[,(i-j+1)],k,k)      
      j   <- j + 1    
    }
    tmp <- t(ss)
    cl  <- cbind(cl,  as.vector(tmp))    
  }


  # Construct orthogonal impulse response functions
  # Long run restrictions at optimal parameters
  f   <- matrix(c(b[1],   0,       0,      0,     0,
                    b[2],   b[3],    b[4],   0,  0,
                    b[5],  -b[3],   b[6],   0,  b[12],
                    b[7],   b[8],    b[9],  b[10],  0,
                    0,     0,       0,       0,      b[11]), nrow=ncol(v), byrow=T)

	s      <- a1 %*% f		# Long-run restrictions

  impulse <- as.vector(t(s))
  for(i in 2:(q+1)) {    
    tmp      <- t( (reshapeg( cl[,i],k,k) %*% s) )
    impulse  <- cbind(impulse,  as.vector(tmp))
  }
  impulse <- t(impulse)
  
  # Compute cumulative sum of impulse responses for levels
  impulse <- apply(impulse, 2, cumsum)
  
  # Construct variance decompositions 
  csum <- apply(trimr(impulse, 0, 1)^2, 2, cumsum)
  tmp0   <- reshapeg(csum, q*k,k )
  tmp1   <- repmat( rowSums(tmp0),1,k )
  decomp <- reshapeg( 100*tmp0/tmp1 , q , k^2 )

  #**********************************************************************
  #***
  #***     Generate graph
  #***
  #**********************************************************************

  figure()
  par(xaxs="i", yaxs="i", mfrow=c(5,5))
    
  # Plot impulse responses
  #--------------------------------------------------------#
  matplot(seqa(0,1,q+1), cbind(impulse[1:(q+1),1], rep(0, (q+1))), type="l",
          main = 'Agg. supply shock',
          ylab = 'Real Output',
          xlab = 'Quarter',
          bty = "l")


  #--------------------------------------------------------#
  matplot(seqa(0,1,q+1), cbind(impulse[1:(q+1),2], rep(0, q+1)), type="l",
          main = 'Agg. demand shock',
          ylab = '',
          xlab = 'Quarter',
          bty = "l")

  #--------------------------------------------------------#
  matplot(seqa(0,1,q+1), cbind(impulse[1:(q+1),3], rep(0, q+1)), type="l",
          main = 'Aust. portfolio shock',
          ylab = '',
          xlab = 'Quarter',
          bty = "l")

 
  #--------------------------------------------------------#
  matplot(seqa(0,1,q+1), cbind(impulse[1:(q+1),4], rep(0, q+1)), type="l",
          main = 'Nominal shock',
          ylab = '',
          xlab = 'Quarter',
          bty = "l")   
	
    
  #--------------------------------------------------------#
  matplot(seqa(0,1,q+1), cbind(impulse[1:(q+1),5], rep(0, q+1)), type="l",
          main = 'US portfolio shock',
          ylab = '',
          xlab = 'Quarter',
          bty = "l")   
    
  #--------------------------------------------------------#
  matplot(seqa(0,1,q+1), cbind(impulse[1:(q+1),6], rep(0, q+1)), type="l",
          main = '',
          ylab = 'Interest rate',
          xlab = 'Quarter',
          bty = "l")


  #--------------------------------------------------------#
  matplot(seqa(0,1,q+1), cbind(impulse[1:(q+1),7], rep(0, q+1)), type="l",
          main = '',
          ylab = '',
          xlab = 'Quarter',
          bty = "l")
  
    
  #--------------------------------------------------------#
  matplot(seqa(0,1,q+1), cbind(impulse[1:(q+1),8], rep(0, q+1)), type="l",
          main = '',
          ylab = '',
          xlab = 'Quarter',
          bty = "l")
 
  #--------------------------------------------------------#
  matplot(seqa(0,1,q+1), cbind(impulse[1:(q+1),9], rep(0, q+1)), type="l",
          main = '',
          ylab = '',
          xlab = 'Quarter',
          bty = "l")

 
  #--------------------------------------------------------#
  matplot(seqa(0,1,q+1), cbind(impulse[1:(q+1),10], rep(0, q+1)), type="l",
          main = '',
          ylab = '',
          xlab = 'Quarter',
          bty = "l")  

  #--------------------------------------------------------#
  matplot(seqa(0,1,q+1), cbind(impulse[1:(q+1),11], rep(0, q+1)), type="l",
          main = '',
          ylab = 'Real Aust. equity',
          xlab = 'Quarter',
          bty = "l")

  #--------------------------------------------------------#
  matplot(seqa(0,1,q+1), cbind(impulse[1:(q+1),12], rep(0, q+1)), type="l",
          main = '',
          ylab = '',
          xlab = 'Quarter',
          bty = "l")
  #--------------------------------------------------------#
  matplot(seqa(0,1,q+1), cbind(impulse[1:(q+1),13], rep(0, q+1)), type="l",
          main = '',
          ylab = '',
          xlab = 'Quarter',
          bty = "l")
  
  #--------------------------------------------------------#
  matplot(seqa(0,1,q+1), cbind(impulse[1:(q+1),14], rep(0, q+1)), type="l",
          main = '',
          ylab = '',
          xlab = 'Quarter',
          bty = "l")

 #--------------------------------------------------------#
  matplot(seqa(0,1,q+1), cbind(impulse[1:(q+1),15], rep(0, q+1)), type="l",
          main = '',
          ylab = '',
          xlab = 'Quarter',
          bty = "l")


  #--------------------------------------------------------#
  matplot(seqa(0,1,q+1), cbind(impulse[1:(q+1),16], rep(0, q+1)), type="l",
          main = '',
          ylab = 'Price',
          xlab = 'Quarter',
          bty = "l")

  #--------------------------------------------------------#
  matplot(seqa(0,1,q+1), cbind(impulse[1:(q+1),17], rep(0, q+1)), type="l",
          main = '',
          ylab = '',
          xlab = 'Quarter',
          bty = "l")

    
  #--------------------------------------------------------#
  matplot(seqa(0,1,q+1), cbind(impulse[1:(q+1),18], rep(0, q+1)), type="l",
          main = '',
          ylab = '',
          xlab = 'Quarter',
          bty = "l")
    
     
  #--------------------------------------------------------#
  matplot(seqa(0,1,q+1), cbind(impulse[1:(q+1),19], rep(0, q+1)), type="l",
          main = '',
          ylab = '',
          xlab = 'Quarter',
          bty = "l")

  #--------------------------------------------------------#
  matplot(seqa(0,1,q+1), cbind(impulse[1:(q+1),20], rep(0, q+1)), type="l",
          main = '',
          ylab = '',
          xlab = 'Quarter',
          bty = "l") 
	
  #--------------------------------------------------------#
  matplot(seqa(0,1,q+1), cbind(impulse[1:(q+1),21], rep(0, q+1)), type="l",
          main = '',
          ylab = 'Real US equity',
          xlab = 'Quarter',
          bty = "l")
  #--------------------------------------------------------#
  matplot(seqa(0,1,q+1), cbind(impulse[1:(q+1),22], rep(0, q+1)), type="l",
          main = '',
          ylab = '',
          xlab = 'Quarter',
          bty = "l")

  #--------------------------------------------------------#
  matplot(seqa(0,1,q+1), cbind(impulse[1:(q+1),23], rep(0, q+1)), type="l",
          main = '',
          ylab = '',
          xlab = 'Quarter',
          bty = "l")
	
  #--------------------------------------------------------#
  matplot(seqa(0,1,q+1), cbind(impulse[1:(q+1),24], rep(0, q+1)), type="l",
          main = '',
          ylab = '',
          xlab = 'Quarter',
          bty = "l")

 
  #--------------------------------------------------------#
  matplot(seqa(0,1,q+1), cbind(impulse[1:(q+1),25], rep(0, q+1)), type="l",
          main = '',
          ylab = '',
          xlab = 'Quarter',
          bty = "l")
}


