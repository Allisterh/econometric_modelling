#============================================================================
#
#   Compute GMM estimates of the parameters of a student t distribution
#
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()
set.seed(123, kind="Mersenne-Twister")

#
#------------------------- Helper Functions----------------------------------
#

# load required functions - inv, numhess, numgrad
source("EMTSUtil.R")
# Load required library - grad, hessian
library("numDeriv")

#----------------------------------------------------------------------------
# GMM objective function - 2 moment conditions
#----------------------------------------------------------------------------
gmmcrit2 <- function (theta,y) {
  mu <- theta[1]
  nu <- theta[2]
  m1 <- y - mu
  m2 <- (y - mu)^2 - nu/(nu - 2)
  m  <- cbind(m1,  m2)
  
  w  <- t(m) %*% m/length(y)
  q  <- 0.5*colMeans(m) %*% inv(w) %*% cbind(colMeans(m))  
  return(q)
}
#----------------------------------------------------------------------------
# GMM objective function - 3 moment conditions
#----------------------------------------------------------------------------
gmmcrit3 <- function(theta,y) {
   mu <- theta[1]
   nu <- theta[2]
   m1 <- y - mu
   m2 <- (y - mu)^2 - nu/(nu - 2)
   m3 <- (y - mu)^4 - 3*nu^2/((nu - 2)*(nu - 4))
   m  <- cbind(m1,  m2, m3)
   w  <- t(m) %*% m/length(y)
   q  <- 0.5*colMeans(m) %*% inv(w) %*% cbind(colMeans(m))  
   return(q)
}


gmm_student <- function() {
  t <- 10                             
  y <- as.matrix(read.table("table.dat"))

  # Using 2 moment conditions 
  # Zero iteration
  theta0 <- c(mean(y),  5)
  g <- cbind((grad(gmmcrit2,theta0,y=y)))
  h <- hessian(gmmcrit2,theta0,y=y)

  cat('\nZero iteration')
  cat('\n Parameter estimate = ', theta0 )
  cat('\n Objective function = ',  gmmcrit2(theta0,y) )
  cat('\n Gradient           = ', g )
  cat('\n Hessian            = \n')
  print(h)
  
  # Newton Raphson update
  theta1 <- theta0 - inv(h) %*% g   
  theta1 <- as.vector(theta1)

  # First iteration
  g <- cbind( grad(gmmcrit2,theta1,y=y) )
  h <- hessian(gmmcrit2,theta1,y=y)

  cat('\n ')
  cat('\nFirst iteration')
  cat('\n Parameter estimate = ', theta1 )
  cat('\n Objective function = ',  gmmcrit2(theta1,y) )
  cat('\n Gradient           = ', g )
  cat('\n Hessian            = \n')
  print(h)

  # Newton Raphson update
  theta2 <- as.vector(theta1 - inv(h) %*% g)
  
  # Second iteration
  g <- t( grad(gmmcrit2,theta2,y=y) )
  h <- hessian(gmmcrit2,theta2,y=y)
  
  cat('\n ')
  cat('\nSecond iteration')
  cat('\n Parameter estimate = ', theta2 )
  cat('\n Objective function = ',  gmmcrit2(theta2,y) )
  cat('\n Gradient           = ', g )
  cat('\n Hessian            = \n')
  print( h )

  v <- inv(h)/t
  cat('\nCovariance matrix   = \n')
  print( v )

  # Iterative solution
  estResults <- optim(theta0, gmmcrit2, y=y, method="BFGS", hessian=T)
  thetahat <- estResults$par
  fc <- estResults$value
  H <- estResults$hessian
  
  cat('\n ')
  cat('\nTwo moment conditions')
  cat('\n Parameter estimate = ', thetahat )
  cat('\n Objective function = ',  fc  )
  v <- inv(H)/t
  cat('\nCovariance matrix   = \n')
  print( v )

  # Using 3 moment conditions
  estResults <- optim(theta0, gmmcrit3, y=y, method="BFGS", hessian=T)
  thetahat <- estResults$par
  fc <- estResults$value
  H <- estResults$hessian

  cat('\n ')
  cat('\nThree moment conditions')
  cat('\n Parameter estimate = ', thetahat )
  cat('\n Objective function = ',  fc  )
  v = inv(H)/t
  cat('\nCovariance matrix   = \n')
  print( v )
}
  
