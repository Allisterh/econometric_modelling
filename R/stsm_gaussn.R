#===========================================================================
#
#   Program to estimate an ARMA(1,1) using the GAUSS-NEWTON algorithm
#
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()
set.seed(12345, kind="Mersenne-Twister")
#
#------------------------- Helper Functions----------------------------------
#

#load required functions - trimr, recserar
source("EMTSUtil.R")

#
#------------------------- The Gauss-Newton Algorithm -----------------------
#
stsm_gaussn <- function()
{
  # Alternative starting values
theta <- c(0.2, 0.2, 0.2)       
#theta <- c(0.0, 0.0, 0.0)       
#theta <- c(0.0,0.3, 1.1 ) 
#theta <- c(0.0, 1.1, 0.3 )

# Alternative sample sizes
t <- 200
#t <- 500

# Simulate the data (discard first 100 simulated observations)
mu   <- 0.0            
phi1 <- 0.8           
si1  <- 0.3
vt   <- rnorm(t+101)
yt   <- recserar( cbind(mu + trimr(vt,1,0) + si1*trimr(vt,0,1)) , cbind(0.0) , cbind(phi1))
yt   <- trimr(yt,100,0)   

# Gauss Newton algorithm
crit  <- 0.00001                # Convergence criterion          
maxit <- 20                      # Maximum number of iterations    


i <- 1
while (i <= maxit)
{
  # vt is being redefined in the gauss-newton algorithm   
  vt  <- recserar( cbind(trimr(yt,1,0) - theta[1] - theta[2]*trimr(yt,0,1)) , cbind(0.0) , cbind(-theta[3]))
  z1t <- recserar( cbind(rep(1, length(vt))), cbind(0.0), cbind(-theta[3]))
  z2t <- recserar( cbind(trimr(yt,0,1)) , cbind(0.0) , cbind(-theta[3]) )
  z3t <- recserar( cbind(c(0.0, trimr(vt,0,1))) , cbind(0.0) , cbind(-theta[3]))
  zt  <- cbind(z1t,  z2t,  z3t)  

  # Remove all starting zeros   
  vt  <- trimr(vt,2,0)       
  zt  <- trimr(zt,2,0) 
  cat('\nIteration        = ', i) 
  cat('\nParameter vector = ', theta)
  cat('\n ')

  dtheta <- cbind(lm(vt ~ zt - 1)$coef)

  # Check convergence and update parameters
  if (t(dtheta) %*% dtheta < crit)
      break
  else {
    theta   <- theta  + dtheta 
    
    if (i == maxit)
      cat('\nFailed to converge after iteration   ', maxit)
  }  
  i <- i + 1                 
}
 
sig2 <- (t(vt) %*% vt)/t
vc   <- as.vector(sig2) * inv(t(zt) %*% zt)
se   <- sqrt(diag(vc))

cat('\n ')
print( cbind("Parameters"=theta,  "Std. errors"=se,  "t-stats"=theta/se))

cat('\n ')
cat('\nEstimated asymptotic covariance matrix\n')
print( vc )

# Wald test    
w <- theta[3]^2/vc[3,3]

cat('\n ')
cat('\nWald statistic (si1 = 0.0)     = ', w)
cat('\np-value                        = ', 1-pchisq(w,1))
cat('\n ')

# LM test   
# First regression
y <- trimr(yt,1,0)   
x <- cbind(rep(1, length(y)), trimr(yt,0,1))
b <- lm(y ~ x - 1)$coef
e <- y - x %*% b

# Second regression
y   <- trimr(e,1,0)    
x   <- cbind(rep(1, length(y)), trimr(yt,1,1), trimr(e,0,1))
b <- lm(y ~ x - 1)$coef
e   <- y - x %*% b
y <- y - mean(y)

rsq <- 1 - (t(e) %*% e)/(t(y) %*% y)     
lm  <- length(yt)*rsq

cat('\n ')
cat('\nLM statistic (si1 = 0.0)       = ', lm) 
cat('\np-value                        = ', 1-pchisq(lm,1))
cat('\n ')

}
