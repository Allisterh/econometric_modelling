#============================================================================
#
#     Program to estimate the Peerman SVAR model
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
# Wrapper function for log-liklihood
#--------------------------------------------------------------------------
neglogl <- function( b,v,a1 ) {
   lf <- - mean( loglt( b,v,a1 ) )  
   return(lf)
}

#--------------------------------------------------------------------------
# Log-likelihood function for SVAR
#--------------------------------------------------------------------------
loglt <- function(b,v,a1inv) {
  t <- nrow(v)
  n <- ncol(v)
  lf      <- rep(0, t)
  
  # Long-run restrictions
  s43 <- -(a1inv[2,2]*b[8]+a1inv[2,3]*b[9])/a1inv[2,4]  
  s34 <- -(a1inv[2,4]*b[10])/a1inv[2,3]     

  s   <-   matrix(c(b[1],    0,      0,     0,
                    b[2],   b[5],   b[8],   0,
                    b[3],   b[6],   b[9],   s34,
                    b[4] ,  b[7],    s43,  b[10]), nrow=n, byrow=T)
		  
	omegav <- s %*% t(s)

  for (i in seq(t)) {
    lf[i] <- -0.5*n*log(2*pi) - 0.5*log(det(omegav)) - 0.5*v[i,] %*% inv(omegav) %*% cbind(v[i,])    
  }
  return(lf)  
}

#
#---------------The Peersman SVAR Model of Oil Price Shocks------------------
#
svar_peersman <- function() {
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
  ytdata <- as.matrix(ytdata)
 

  # Data transformations
  loil <- log( ytdata[,1] )     ##ok<NODEF> # log price of oil    
  lo   <- log( ytdata[,5] )     # log output          
  lp   <- log( ytdata[,6] )     # log price level     
  r    <- ytdata[,7]            # interest rate       

  # Construct variables for use in VAR: 
  # # growth rates in oil, output and price, and level of interest rate
  yact <- cbind(loil, lo, lp, r)
  yact <- trimr(yact,36,0)       # Sample period 1970Q1 to 2002Q2
  y <- cbind(100*(trimr(yact[,1:3],1,0) - trimr(yact[,1:3],0,1)),  trimr(yact[,4],1,0))
  # Set SVAR parameters
  p <- 3      # Order of VAR lags
  q <- 40     # Order of VMA lags      

  # Estimate the VAR with p lags, a constant and a time trend     
  ylag   <- cbind(rep(1, nrow(y)),  seq(nrow(y)))
  nconst <- ncol(ylag)
  for (i in seq(p)) {
     ylag <- cbind(trimr(ylag,1,0), trimr(y,0,i))
  }

  # OLS Estimates
  reg <- lm(trimr(y,p,0) ~ ylag - 1)
  bar    <- reg$coef
  mue    <- bar[1:nconst,]
  v      <- reg$residuals
  omegav <- t(v) %*% v/nrow(v)

  # Constuct A(L) (VAR) matrices and hence the C(1) long-run matrix
  bar <- trimr(bar,nconst,0)
  k   <- ncol(v)
 	a   <- array(0, c(k^2,p))
  a1  <- diag(k)
  for (i in seq(p)) {
    tmp    <- bar[(1+k*(i-1)):(i*k),]
    a[,i] <- c(tmp)
    a1     <- a1 - reshapeg(a[,i],k,k)   
  }
    
  # Invert A(1) matrix needed for the MLE
  a1inv <- inv(a1)            

  bstart <- c(12.9662120731029020,
            0.0133377541945983, 
            0.0917703235318815, 
            0.1492472650864925, 
            0.3093300004256599, 
            -0.1277889658083315, 
            -0.1312058376209482, 
            0.4307222392355360,
            0.1183146060546325, 
            0.6944925202098823)
  estResults <- optim(bstart, neglogl, v=v, a1=a1inv, method="BFGS")
  b <- estResults$par
  fval <- -estResults$val

  cat('\nLog-likelihood function', fval)

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
  s43 <- -(a1inv[2,2]*b[8]+a1inv[2,3]*b[9])/a1inv[2,4]  
  s34 <- -(a1inv[2,4]*b[10])/a1inv[2,3]     
  s   <-  matrix(c(b[1],   0,      0,       0,      
                  b[2],   b[5],    b[8],    0,       
	                b[3],   b[6],    b[9],    s34,     
   		            b[4],   b[7],    s43,     b[10]), nrow=4, byrow=T) 
  cat('\nS matrix at optimal parameters\n')            
  print(s)

  impulse <- as.vector(t(s))
  for(i in 2:(q+1)) {    
    tmp      <- t( (reshapeg( cl[,i],k,k) %*% s) )
    impulse  <- cbind(impulse,  as.vector(tmp))
  }

  impulse <- t(impulse)
  # Compute cumulative sum of impulse responses for levels ie first three series 
  csum <- apply(impulse[,1:12], 2, cumsum)
  impulse <- cbind(csum, impulse[,13:16])

  # Construct variance decompositions 
  csum <- apply(trimr(impulse, 0, 1)^2, 2, cumsum)
  tmp0   <- reshapeg(csum,q*k,k)
  tmp1   <- repmat( rowSums(tmp0),1,k )
  decomp <- reshapeg( 100*tmp0/tmp1 , q , k^2 )
 

  #**********************************************************************
  #***
  #***     Generate graph
  #***
  #**********************************************************************

  figure()
  par(xaxs="i", yaxs="i", mfrow=c(4,4))
  t <- 0:q

  #--------------------------------------------------------#
  # Panel (a)
  matplot(t,cbind(impulse[1:(q+1),1], rep(0, (q+1))), type="l",
          main = 'Oil Price Shock',
          ylab = 'Oil price',
          xlab = 'Quarter',
          bty = "l")
    
  #--------------------------------------------------------#
  # Panel (b)
  matplot(t,cbind(impulse[1:(q+1),2], rep(0, (q+1))), type="l",
          main = 'Supply Shock',
          ylab = '',
          xlab = 'Quarter',
          bty = "l")

  #--------------------------------------------------------#
  # Panel (c)
  matplot(t,cbind(impulse[1:(q+1),3], rep(0, (q+1))), type="l",
          main = 'Demand Shock',
          ylab = '',
          xlab = 'Quarter',
          bty = "l")


  #--------------------------------------------------------#
  # Panel (d)
  matplot(t,cbind(impulse[1:(q+1),4], rep(0, (q+1))), type="l",
          main = 'Money Shock',
          ylab = '',
          xlab = 'Quarter',
          bty = "l")
  #--------------------------------------------------------#
  # Panel (e)
  matplot(t,cbind(impulse[1:(q+1),5], rep(0, (q+1))), type="l",
          main = '',
          ylab = 'Output',
          xlab = 'Quarter',
          bty = "l")
   

  #--------------------------------------------------------#
  # Panel (f)
  matplot(t,cbind(impulse[1:(q+1),6], rep(0, (q+1))), type="l",
          main = '',
          ylab = '',
          xlab = 'Quarter',
          bty = "l")
   
  #--------------------------------------------------------#
  # Panel (g)
  matplot(t,cbind(impulse[1:(q+1),7], rep(0, (q+1))), type="l",
          main = '',
          ylab = '',
          xlab = 'Quarter',
          bty = "l")

  #--------------------------------------------------------#
  # Panel (h)
  matplot(t,cbind(impulse[1:(q+1),8], rep(0, (q+1))), type="l",
          main = '',
          ylab = '',
          xlab = 'Quarter',
          bty = "l")

  #--------------------------------------------------------#
  # Panel (i)
  matplot(t,cbind(impulse[1:(q+1),9], rep(0, (q+1))), type="l",
          main = '',
          ylab = 'Price',
          xlab = 'Quarter',
          bty = "l")

  #--------------------------------------------------------#
  # Panel (j)
  matplot(t,cbind(impulse[1:(q+1),10], rep(0, (q+1))), type="l",
          main = '',
          ylab = '',
          xlab = 'Quarter',
          bty = "l")
    
  #--------------------------------------------------------#
  # Panel (k)
  matplot(t,cbind(impulse[1:(q+1),11], rep(0, (q+1))), type="l",
          main = '',
          ylab = '',
          xlab = 'Quarter',
          bty = "l")
    
  #--------------------------------------------------------#
  # Panel (l)
  matplot(t,cbind(impulse[1:(q+1),12], rep(0, (q+1))), type="l",
          main = '',
          ylab = '',
          xlab = 'Quarter',
          bty = "l")
    
    
  #--------------------------------------------------------#
  # Panel (m)
  matplot(t,cbind(impulse[1:(q+1),13], rep(0, (q+1))), type="l",
          main = '',
          ylab = 'Interest rate',
          xlab = 'Quarter',
          bty = "l")
     
  #--------------------------------------------------------#
  # Panel (n)
  matplot(t,cbind(impulse[1:(q+1),14], rep(0, (q+1))), type="l",
          main = '',
          ylab = '',
          xlab = 'Quarter',
          bty = "l")

  #--------------------------------------------------------#
  # Panel (o)
  matplot(t,cbind(impulse[1:(q+1),15], rep(0, (q+1))), type="l",
          main = '',
          ylab = '',
          xlab = 'Quarter',
          bty = "l")
 
  #--------------------------------------------------------#
  # Panel (p)
  matplot(t,cbind(impulse[1:(q+1),16], rep(0, (q+1))), type="l",
          main = '',
          ylab = '',
          xlab = 'Quarter',
          bty = "l")
}
