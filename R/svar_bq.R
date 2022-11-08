#============================================================================
#
#      Estimate the Blanchard-Quah SVAR model
#
#============================================================================

rm (list = ls(all=TRUE))
graphics.off()

#
#------------------------- Helper Functions----------------------------------
#

# Load required functions -  trimr, reshapeg, figure, lag.matrix
source("EMTSUtil.R")

# Load required library - repmat
library("matlab")


#----------------------------------------------------------------------------
# Log-likelihood function for SVAR with long-run restrictions
#----------------------------------------------------------------------------
neglog <- function(b,v,a1) {
  t <- nrow(v)
  n <- ncol(v)
  lnl  <- rep(0, t)
  f   <-   matrix(c(b[1],  0,
                    b[2],  -abs(b[3])), nrow=2, byrow=T)
  s   <- a1 %*% f  	# Long-run restrictions		 
	vc  <- s %*% t(s)
  ivc <- inv(vc)    
	for (i in seq(t)) {
    lnl[i] <- -0.5*n*log(2*pi) - 0.5*log(det(vc)) - 0.5*v[i,] %*% ivc %*% cbind(v[i,])
	}
  lf <- -mean(lnl)
  return(lf)
}

#----------------------------------------------------------------------------
# Return SVAR matrices
#----------------------------------------------------------------------------
svarmat <- function(b,a1) {
  f <- matrix(c(b[1],   0,  
            b[2],   -abs(b[3])), nrow=2, byrow=T)   

  s   <- a1 %*% f		# Long-run restrictions
  return(s)
}

#----------------------------------------------------------------------------
# Log-likelihood function for SVAR with Okuns Law and long-run restrictions
#----------------------------------------------------------------------------
neglogc <- function(b,v,a1) {
  t <- nrow(v)
  n <- ncol(v)
  lf  <- rep(0, t)
  f   <-  matrix(c(b[1],      0,
                      -3*b[1], -abs(b[2])), nrow=2, byrow=T)
                    
  s      <- a1 %*% f    # Long-run restrictions
	vc <- s %*% t(s)  
  ivc <- inv(vc)
  for (i in seq(t)) {
    lf[i] <- -0.5*n*log(2*pi) - 0.5*log(det(vc)) - 0.5*v[i,] %*% ivc %*% cbind(v[i,])    
  }
  lnl <- -mean(lf)    
  return(lnl)
}

#
#-------------------- The Blanchard-Quah Model and Okun's Law ---------------
#
svar_bq <- function() {
    # Load data (1950Q1 - 2009Q3)
  #    (1) Real GDP, s.a.
  #    (2) Unemployment rate, s.a.
  load('bq_us.Rdata')
  #load('bq_jp.Rdata')
  #load('bq_uk.Rdata')

  ytdata <- ytdata_us
  
  # Construct variables for use in VAR
  lrgdp  <- log(ytdata[,1])
  urate <- ytdata[,2]

  y <- cbind(400*(trimr(lrgdp,1,0) - trimr(lrgdp,0,1)),  trimr(urate,1,0))
  t <- dim(y)[1]
  n <- dim(y)[2]
    
  cat('\nCovariance matrix of the data\n')
  print(cov(y))

  # Estimate the VAR with p lags, a constant and a time trend    
  p <- 8		# Order of VAR lags	 
  q <- 40		# Order of VMA lags		 

  ylag <- lag.matrix(y,1:p)
  ylag <- cbind(rep(1,t), ylag)
  ylag <- ylag[-(1:p),]

  yt  <- trimr(y,p,0)
  reg <- lm(yt ~ ylag - 1)
  bar <- reg$coef
  v   <- reg$residuals
 
  vc  <- t(v) %*% v/nrow(v)

  cat('\nCovariance matrix of VAR residuals\n')
  print(vc)

  # Constuct A(L) (VAR) matrices and hence the C(1) long-run matrix
  bar <- trimr(bar,1,0)
  k   <- ncol(v)
  a   <- array(0, c(k^2,p))
  a1  <- diag(k)

  for (i in seq(p)) {
    tmp    <- bar[(1+k*(i-1)):(i*k),]
    a[,i]  <- as.vector(tmp)
    a1    <- a1 - reshapeg(a[,i],k,k)    
  }
  # Invert A(1) matrix needed for the MLE
  a1inv <- inv(a1) 
  cat('\nInverse of a1 matrix\n')    
  print(a1inv)

  # Esimate Blanchard-Quah SVAR model with long-run restrictions
  bstart <- vech(t(chol(vc)))
  estResults <- optim(bstart, neglog, v=v, a1=a1, method="BFGS", hessian=T)
  b <- estResults$par
  fval <- estResults$val
  H <- estResults$hessian
  
  cat('\nLog-likelihood function  = ',-fval)
  vcov <- inv(H)/nrow(v)

  # Test of Okun's Law: long-run version 
  cc <- c(b[2]/b[1] + 3)              # Aggregate supply shock parameters 
  d <- c(-b[2]/b[1]^2,   1/b[1],  0)
  w <- cbind(cc) %*% inv(d %*% vcov %*% cbind(d)) %*% cc
  
  

  cat('\nWald test of Okuns Law (aggregate supply version)')
  cat('\nMean            = ', b[2]/b[1])
  cat('\nWald statistic  = ', w)
  cat('\np-value         = ',1-pchisq(w,1))
  
  cat('\n ' )
  cat('\nWald test of Okuns Law (aggregate demand version)')
  cat('\nAggregate demand shocks have no affect on output in the long-run')
  cat('\n')
    
  # Impuse responses
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
  s <- svarmat(b,a1)

  impulse <- as.vector(t(s))
  for(i in 2:(q+1)) {    
    tmp      <- t( (reshapeg( cl[,i],k,k) %*% s) )
    impulse  <- cbind(impulse,  as.vector(tmp))
  }

  impulse <- t(impulse)
    
  # Compute cumulative sum of impulse responses for levels 
  csum <- apply(impulse[, c(1,2)], 2, cumsum)
	impulse <- cbind(csum, impulse[,c(3,4)])

  # Construct variance decompositions
  csum <- apply(trimr(impulse, 0, 1)^2, 2, cumsum)
  tmp0   <- reshapeg(csum, q*k,k )
  tmp1   <- repmat( rowSums(tmp0),1,k )
  decomp <- reshapeg( 100*tmp0/tmp1 , q , k^2 )

  # SVAR with long-run restriction and Okuns Law restriction  
  bstart <- b[c(1,3)]
  estResults <- optim(bstart, neglogc, v=v, a1=a1, method="BFGS")
  fval <- estResults$val
  cat('\nLog-likelihood function  = ',-fval)
  cat('\n ')

  #**********************************************************************
  #***
  #***     Generate graph
  #***
  #**********************************************************************

  figure()
  par(xaxs="i", yaxs="i", mfrow=c(2,2))
    
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
          ylab = 'Real Output',
          xlab = 'Quarter',
          bty = "l")

  #--------------------------------------------------------#
  matplot(seqa(0,1,q+1), cbind(impulse[1:(q+1),3], rep(0, q+1)), type="l",
          main = 'Agg. supply shock',
          ylab = 'Unemployment',
          xlab = 'Quarter',
          bty = "l")
  #--------------------------------------------------------#
  matplot(seqa(0,1,q+1), cbind(impulse[1:(q+1),4], rep(0, q+1)), type="l",
          main = 'Agg. demand shock',
          ylab = 'Unemployment',
          xlab = 'Quarter',
          bty = "l")
}


