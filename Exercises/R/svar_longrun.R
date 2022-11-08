#============================================================================
#
#   Estimate a Structural VAR based on long-run restrictions
#
#============================================================================

rm (list = ls(all=TRUE))
graphics.off()

#
#------------------------- Helper Functions----------------------------------
#

#load required functions -  trimr, inv
source("EMTSUtil.R")

#----------------------------------------------------------------------------
#      Negative log-likelihood function   
#----------------------------------------------------------------------------
neglog <- function(b,v,a1) {
  t <- nrow(v)
  n <- ncol(v)
  f <-  matrix(c(b[1], 0,  0,
                 b[4],  b[2], 0,
                 b[5],  0,  b[3]), nrow=n, byrow=T)
  s  <- (diag(3) - a1) %*% f
	vc <- s %*% t(s)
  
  lnl <- rep(0, t)
  for (i in seq(t)) {
    lnl[i] <- -0.5*n*log(2*pi) - 0.5*log(det(vc)) - 0.5*v[i,] %*% inv(vc) %*% cbind(v[i,])    
  }
  lf = -mean(lnl)
  return(lf)
}

#----------------------------------------------------------------------------
#   Return SVAR matrices   
#----------------------------------------------------------------------------
func <- function(b,a1)
{    
  f <-  matrix(c(b[1], 0,  0,
                 b[4],  b[2], 0,
                 b[5],  0,  b[3]), nrow=3, byrow=T)
  s  <- (diag(3) - a1) %*% f
  vc0 <- s %*% t(s) 
  return(list(f=f, s=s, vc0=vc0))  
}


#
#--------- A Model of Interest Rates with Long-Run Restrictions ------------
#
svar_longrun <- function() 
{
  load('mcnew.Rdata')
  rt <- as.matrix(rt)
  yt   <- trimr(rt[,1:3],1,0)-trimr(rt[,1:3],0,1)

  # Estimate a VAR with one lag
  y  <- trimr(yt,1,0)
  x  <- cbind(rep(1, nrow(y)),  trimr(yt,0,1) )
  a  <- lm(y ~ x - 1)$coef
  v  <- y - x %*% a                            
  vc <- t(v) %*% v/nrow(v)

  # Lag 1 autoregressive parameters   
  a1 <- trimr(a,1,0)   

  # Esimate SVAR
  theta0   <- rep(1,5)
  estResults <- optim(theta0, neglog, v=v, a1=a1, method="BFGS")
  b <- estResults$par
  fval <- estResults$val

  cat('\nLog-likelihood value = ', -fval)
  cat('\nVAR parameter estimates\n')
  print(b)
  cat('\n')
  cat('\nVAR variance-covariance matrix (unconstrained)\n')
  print(vc)
    
  # Compute SVAR matrices
  smat <- func(b, a1)
  f <- smat$f
  s <- smat$s
  vc0 <- smat$vc0

  cat('\nVAR variance-covariance matrix (constrained)\n')
  print( vc0 )
  cat('\nDeterminant of omegav   = ', det(vc0))
  cat('\nI-A1\n')
  print( diag(3) - a1)
  cat('\nf\n')     
  print(f)
  cat('\ns\n') 		
  print(s)  
}

  
