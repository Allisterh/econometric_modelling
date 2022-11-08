#=========================================================================
#
#   Compute GMM estimates of the parameters of a gamma distribution
#
#=========================================================================
rm (list = ls(all=TRUE))
graphics.off()
set.seed(123, kind="Mersenne-Twister")

#
#------------------------- Helper Functions----------------------------------
#

# Load required functions - inv, 
source("EMTSUtil.R")
# Load required library - grad, hessian
library("numDeriv")

#-------------------------------------------------------------------------
# GMM objective function - 2 moment conditions
#-------------------------------------------------------------------------
gmmcrit2 <- function(alpha,y) {
   m1 <- y - alpha
   m2 <- y^2 - alpha*(alpha+1)
   m  <- cbind(m1,  m2)
   w  <- t(m) %*% m/length(y)
   q  <- 0.5*colMeans(m) %*% inv(w) %*% cbind(colMeans(m))
   return(q)  
}

#-------------------------------------------------------------------------
# GMM objective function - 3 moment conditions
#-------------------------------------------------------------------------
gmmcrit3 <- function(alpha,y) {
  m1 <- y - alpha
  m2 <- y^2 - alpha*(alpha+1)
  m3 <- 1/y - 1/(alpha-1)
  m  <- cbind(m1,  m2, m3)
  w  <- t(m) %*% m/length(y)
  q  <- 0.5*colMeans(m) %*% inv(w) %*% cbind(colMeans(m))
  return(q)  
}

#-------------------------------------------------------------------------
# GMM objective function - 2 parameter gamma function
#-------------------------------------------------------------------------
gmmcrit <- function(theta,y) {
  alpha <- theta[1]
  beta  <- theta[2]
  m1    <- y - alpha/beta
  m2    <- y^2 - alpha*(alpha+1)/beta^2
  m3    <- 1/y - beta/(alpha-1)
  m     <- cbind(m1,  m2,  m3)
  w     <- t(m) %*% m/length(y)
  q     <- 0.5*colMeans(m) %*% inv(w) %*% cbind(colMeans(m))  
  return (q)
}


#
#----------------- Estimating a Gamma Distribution --------------------------
#
gmm_gamma <- function ()
{
    t <- 10 
  # Simulate data
   y <- round(5 + 2*rnorm(t))

  # Load GAUSS data
  #y <- as.matrix(read.table("table.dat"))

  # Using 2 moment conditions 
  # Zero iteration
  alpha0 <- mean(y)  
  g <- grad(gmmcrit2, alpha0, y=y)  
  h <- hessian(gmmcrit2,alpha0,y=y)

  cat('\nZero iteration')
  cat('\n Parameter estimate = ', alpha0)
  cat('\n Objective function = ',  gmmcrit2(alpha0,y))
  cat('\n Gradient           = ', g)
  cat('\n Hessian            = ', h)

  # Newton Raphson update
  alpha1 <- as.vector(alpha0 - inv(h)*g)

  # First iteration
  g <- grad(gmmcrit2,alpha1,y=y)           
  h <- hessian(gmmcrit2,alpha1,y=y)

  cat('\n ')
  cat('\nFirst iteration')
  cat('\n Parameter estimate = ', alpha1)
  cat('\n Objective function = ', gmmcrit2(alpha1,y))
  cat('\n Gradient           = ', g)
  cat('\n Hessian            = ', h)
  
   # Newton Raphson update
   alpha2 <- as.vector(alpha1 - inv(h)*g)
  
  # Second iteration
  g <- grad(gmmcrit2,alpha2,y=y)           
  h <- hessian(gmmcrit2,alpha2,y=y)

  cat('\n ')
  cat('\nSecond iteration')
  cat('\n Parameter estimate = ', alpha2)
  cat('\n Objective function = ',  gmmcrit2(alpha2,y))
  cat('\n Gradient           = ', g)
  cat('\n Hessian            = ', h)

  v <- inv(h)/t
  cat('\nVariance of alpha   = ', v)
  cat('\nStd error of alpha  = ', sqrt(v))

  # Iterative solution
  estResults <- optim(alpha0, gmmcrit2, y=y, method="BFGS", hessian=T)
  alphahat <- estResults$par
  fc <- estResults$value
  H <- estResults$hessian

  cat('\n ')
  cat('\nTwo moment conditions')
  cat('\n Parameter estimate = ', alphahat)
  cat('\n Objective function = ',  fc )
  v <- inv(H)/t
  cat('\nVariance of alpha   = ', v)
  cat('\nStd error of alpha  = ', sqrt(v))

  # Using 3 moment conditions
  estResults <- optim(alpha0, gmmcrit3, y=y, method="BFGS", hessian=T)
  alphahat <- estResults$par
  fc <- estResults$value
  H <- estResults$hessian

  cat('\n ')
  cat('\nThree moment conditions')
  cat('\n Parameter estimate = ', alphahat)
  cat('\n Objective function = ', fc )
  v <- inv(H)/t
  cat('\nVariance of alpha   = ', v)
  cat('\nStd error of alpha  = ', sqrt(v))

  # Estimating two-parameter gamma distribution
  theta0                 <- c(8,2)
  estResults <- optim(theta0, gmmcrit, y=y, method="BFGS", hessian=T)
  thetahat <- estResults$par
  fc <- estResults$value
  H <- estResults$hessian
  
  cat('\n ')
  cat('\nTwo parameter gamma distribution')
  cat('\n Estimate of alpha  = ', thetahat[1])
  cat('\n Estimate of beta   = ', thetahat[2])
  cat('\n Objective function = ',  fc )
  v <- inv(H)/t
  cat('\nCovariance matrix\n')
  print(v)
}
