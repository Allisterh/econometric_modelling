#==========================================================================
#
#    Estimate the exponential model by maximum likelihood.
#     
#==========================================================================

rm (list = ls(all=TRUE))
graphics.off()
set.seed(123457, kind="Mersenne-Twister")

#
#--------------------------- Helper Functions -----------------------------------
#
# load required functions - numhess and inv
source("EMTSUtil.R")


#--------------------------------------------------------------------------------
# Simulate data an return the results as a list of x and y
#--------------------------------------------------------------------------------
simulatedata <- function(t) {
  b0  <- 1.0
  b1  <- 0.05
  sig <- 0.5
  
  u <- sig*rnorm(t)
  x <- 1:t  
  
  y <- b0*exp( b1*x ) + u  
  return(list(y=y, x=x))
}

#--------------------------------------------------------------------------
#  Log-likelihood function
#--------------------------------------------------------------------------    
neglog <- function(b,y,x) {
  lf <- -mean( lnlt(b,y,x) )
  return(lf)
}
 
#--------------------------------------------------------------------------
#   Log-likelihood (concentrated) at each observation
#--------------------------------------------------------------------------    
lnlt <- function(b,y,x) {
  e  <- y - b[1]*exp( b[2]*x )  
  s2 <- e %*% e/length(e)                        
  lf <- - 0.5*log(2*pi) - 0.5*log(s2) - 0.5*e^2/s2
  return(lf)
}

#
#----------------------- Exponential Model by MLE  ------------------------------
#
nls_exponential <- function() {
  t <- 50
  make.new.random <- FALSE
  if (make.new.random) {
    simResults <- simulatedata(t)
    x <- simResults$x
    y <- simResults$y    
  } else {
    simResults <- read.table("nls_expdata.dat")
    y <- simResults[,1]
    x <- simResults[,2]
  }
  # Estimate the model using Newton Raphson algorithm
  start <- c(0.1, 0.1)
  estResults <- optim(start, neglog, y=y, x=x, method = "BFGS")
  theta <- estResults$par
  
  ht <- numhess(neglog, theta, y, x)
  vcov <- 1/t * inv(ht)
  cat('\nParameter estimates = ', theta)
  cat('\nNegative Hessian matrix\n')
  print(ht)
  cat('\nCovariance matrix\n')
  print(vcov)
}

