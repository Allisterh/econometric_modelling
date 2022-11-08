#=========================================================================
#
#     Program to estimate the parameters of a Weibull distribution using
#     the Newton-Raphson and BHHH algorithms
#
#=========================================================================

rm (list = ls(all=TRUE))
graphics.off()

set.seed(123457, kind="Mersenne-Twister")

max_weibull <- function () { 
  # Sample size
  t <- 20

  # Generate Weibull random numbers
  alpha <- 1.0
  beta  <- 2.0
  
   # Use built in function
   y <- rweibull(t, scale = (1/alpha)^(1/beta), shape = beta)  # Generate Weibull random numbers  

    # or load data
#   y <- c(0.293, 0.589, 1.374, 0.954, 0.608, 1.199, 1.464, 
#           0.383, 1.743, 0.022, 0.719, 0.949, 1.888, 0.754,
#           0.873, 0.515, 1.049, 1.506, 1.090, 1.644)
  
  cat('\n')
  print(cbind(1:t, y))
  
  theta0 <- cbind(c(0.5,1.5))             #     Starting estimates 
  alpha <- theta0[1]
  beta  <- theta0[2]
  
  g <- gradvec(theta0,y)
  h <- hessmat(theta0,y)
  j <- opgmat(theta0,y)
  
  cat('\nEstimates at iteration 0\n')
  print(theta0)
  
  cat('\nValue of log-likelihood at iteration 0    = ', -lnl(theta0, y) )
  
  cat('\n\nGradient vector at iteration 0\n')
  print(g)
  
  cat('\n\nHessian at iteration 0\n')
  print(h)
  
  cat('\n\nOPG matrix at iteration 0\n')  
  print(j)
  
  # Newton Raphson update
  theta1 <- theta0 - solve(h, g)
  
  cat('\n\nEstimates at iteration 1 (Newton Raphson) \n')
  print(theta1)
  
  cat('\nValue of log-likelihood at iteration 1    = ', -lnl(theta1, y) )
  
  # BHHH
  theta1 <- theta0 + solve(j, g)
  
  cat('\n\nEstimates at iteration 1 (BHHH) \n' ) 
  print(theta1 )
  
  cat('\n\nValue of log-likelihood at iteration 1    = ', -lnl(theta1, y) )
  
  # Call optim to optimize function
  theta <- optim(theta1, lnl, y = y, method = "BFGS")$par;
  
  cat('\nMLE of theta\n' )
  print(theta)

  g <- gradvec(theta,y)
  h <- hessmat(theta,y)
  j <- opgmat(theta,y)

  cat('\nCovariance matrix (Hessian)\n' ) 
  print ( (1/t)*solve(h) )

  cat('\nCovariance matrix (OPG)\n' )
  print( (1/t)*solve(j) )

  # Call optim to optimize the concentrated likelihood  
  beta <- optim(theta1[2], lnlc, y = y, method = "BFGS")$par  
    
  alpha <- 1/mean(y^beta)
  cat('\nResults for concentrated likelihood\n' )
  cat('\nMLE of alpha   = ', alpha )
  cat('\nMLE of beta    = ', beta )
    
  # Call optim to optimize the transformed likelihood
  results <- optim(theta1, lnlt, y = y, method = "BFGS", hessian=TRUE)
  theta <- results$par
  h <- results$hessian

  cat('\nResults for transformed likelihood\n' )
  cat('\nMLE of lambda   = ', theta[1])
  cat('\nMLE of beta     = ', theta[2]) 
  
  # Std error by delta method
  d = c( -(1/theta[2])*theta[1]^(-1/theta[2]-1),
        (log(theta[1])/theta[2]^2)*theta[1]^(-1/theta[2]) )
  
  cat('\nStandard error of lambda by delta method' ) 
  cat('\n', -d %*% solve(h) %*% cbind(d) )  
} 

#------------------------- Functions ----------------------------- #

# Log-likelihood 
lnl <-  function(theta,y) {
  a  <- theta[1]
  b  <- theta[2]
  f  <- log(a) + log(b) + (b-1)*log(y) - a*y^b
  lf <- -mean( f )
  return (lf)          
} 

# Concentrated log-likelihood 
lnlc <- function(b,y) {
  a  <- 1/mean(y^b)    
  f  <- log(a) + log(b) + (b-1)*log(y) - a*y^b
  lf <- -mean( f )
  return(lf)
}

# Transformed log-likelihood function
lnlt <- function(theta,y) {
  l  <- theta[1]
  b  <- theta[2]
  f  <- log(b) - log(l ) + (b-1)*log(y/l) - (y/l)^b
  lf <- -mean(f)
  return(lf)
}

# Return gradient vector 
gradvec <- function(theta,y) {
  alpha <- theta[1]
  beta  <- theta[2]
  
  g1 <- 1/alpha - mean(y^beta)
  g2 <- 1/beta + mean(log(y)) - alpha*mean(log(y)*y^beta )
  g  <- cbind(c(g1, g2))  
  return(g)
}

# Return Hessian matrix 
hessmat <- function(theta,y) {
  alpha <- theta[1]
  beta  <- theta[2]

  h11 <- - 1/alpha^2 
  h12 <- - mean(log(y)*y^beta)
  h21 <- h12
  h22 <- -1/beta^2 - alpha*mean( (log(y)^2)*(y^beta) )
  h <- matrix(c(h11, h12, h21, h22), nrow=2, byrow =TRUE)
  return(h)
}

# Return OPG matrix 
opgmat <- function(theta,y) {
  alpha <- theta[1]
  beta  <- theta[2]

  g22 <- 1/beta + log(y) - alpha*log(y)*y^beta
  j11 <- mean( ( 1/alpha - y^beta)^2 )
  j12 <- mean( ( 1/alpha - y^beta) * g22 )
  j21 <- j12
  j22 <- mean( g22^2 )
  j <- matrix(c(j11, j12, j21, j22), nrow=2, byrow =TRUE)  
  return(j)
}
