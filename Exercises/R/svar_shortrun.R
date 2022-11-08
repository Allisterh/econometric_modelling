#============================================================================
#
#   Estimate a Structural VAR based on short-run restrictions
#
#============================================================================

rm (list = ls(all=TRUE))
graphics.off()

#
#------------------------- Helper Functions----------------------------------
#

#load required functions -  trimr, inv, diagrv
source("EMTSUtil.R")

#----------------------------------------------------------------------------
#      Negative log-likelihood function   
#----------------------------------------------------------------------------
neglog <- function(b,v)
{
  t <- nrow(v)
  n <- ncol(v)
  b0 <-   matrix(c(1,    0,  0,
            -b[1], 1, 0, 
            -b[2],  0, 1), nrow=n, byrow=T)

  # Structural residual variances 
	d  <- diagrv(diag(n),abs(b[3:5]))		       
  s  <- inv(b0) %*% sqrt(d)
	vc <- s %*% t(s)
  lnl <- rep(0,t)
  for (i in seq(t)) {
    lnl[i] <- -0.5*n*log(2*pi) - 0.5*log(det(vc))- 0.5*v[i,] %*% inv(vc) %*% cbind(v[i,])    
  }  
  lf <- -mean(lnl)
  
  return(lf)  
}

#-----------------------------------------------------------------------
#   Return SVAR matrices   
#-----------------------------------------------------------------------
func <- function(b,v)
{
  n <- ncol(v)
  b0 <-   matrix(c(1,   0,   0, 
            -b[1],    1,  0,
	          -b[2],  0,  1), nrow=n, byrow=T)
  # Structural residual variances 
	d   <- diagrv(diag(n),abs(b[3:5]))		           
  s   <- inv(b0) %*% sqrt(d)
	vc0 <- s %*% t(s)
  return(list(b0=b0, d=d, s=s, vc0=vc0))
}

#
#--------- A Model of Interest Rates with Short-Run Restrictions ------------
#
svar_shortrun <- function() 
{
  # U.S.zero coupon yields Dec 1946 - Feb 2002 (0,1,3,6,9 mth and 10 yrs)
  load('mcnew.Rdata')
  yt   <- as.matrix(rt[,1:3])
 
  # Estimate a VAR with one lag
  y  <- trimr(yt,1,0)
  x  <- cbind(rep(1, nrow(y)),  trimr(yt,0,1))
  a  <- lm(y ~ x - 1)$coef
  v  <- y - x %*% a                            
  vc <- t(v) %*% v/nrow(v)

  # Esimate SVAR  
  theta0   <- 0.1*rep(1,5)
  estResults <- optim(theta0, neglog, v=v, method="BFGS")
  b <- estResults$par
  fval <- estResults$value
  
  cat('\nLog-likelihood value = ', -fval)
  cat('\nVAR parameter estimates = ', b)
  cat('\nVAR variance-covariance matrix (unconstrained)\n')
  print(vc)

  # Compute SVAR matrices
  smat <- func(b,v)
  b0 <- smat$b0
  d <- smat$d
  vc0 <- smat$vc0
  s <- smat$s
    
  cat('\nVAR variance-covariance matrix (constrained)\n')
  print( vc0 )
  cat('\nDeterminant of omegav   = ', det(vc0))
  cat('\nb0\n')
  print(b0)
  cat('\nd\n')
  print(d)
  cat('\ns\n')
  print(s)
}
