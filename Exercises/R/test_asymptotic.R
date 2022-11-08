#==========================================================================
#
#   Program to generate asymptotic distribution of the Wald test
#   applied to the regression model.
#
#==========================================================================

rm (list = ls(all=TRUE))
graphics.off()
set.seed(1234, kind="Mersenne-Twister")

#
#--------------------------- Helper Functions -----------------------------------
# 
#-------------------------------------------------------------------------
# Wrapper function for matrix inverse
#-------------------------------------------------------------------------
inv <- function(M) {
  return(solve(M))
}

#-------------------------------------------------------------------------
# Log-likelihood function of unconstrained model
#-------------------------------------------------------------------------

neglog <- function(theta,y,x1,x2,x3) {
  m  <- theta[1] + theta[2]*x1 + theta[3]*x2 + theta[4]*x3  
  s2 <- abs(theta[5])                                       
  lf <- -mean( -0.5*log(2*pi) - 0.5*log(s2) - 0.5*(y - m)^2/s2  )
}

# load required functions
source("EMTSUtil.R")
source("EMTSUtil.R")

#
#--------------------------- Exercises  -----------------------------------
#
test_asymptotic <- function() {
  # Set parameter values
  beta0  <- 1.0                                           
  beta1  <- 0.0                                                     
  beta2  <- 0.0                                                       
  beta3  <- 0.0                                                         
  sig2   <- 0.1                 
  sig    <- sqrt(sig2)

  t      <- 1000                                                                     
  ndraws <- 10000

  x1 <- runif(t)                                                                           
  x2 <- rnorm(t)                                                                           
  x3 <- rnorm(t)^2                                                                        

  # Arrays to hold results
  zeros <- array(0, c(ndraws,1))
  wd1 <- zeros                                                 
  wd2 <- zeros                                                
  wd3 <- zeros

  theta0  <- cbind(c(beta0, beta1, beta2, beta3, sig2 ))

  pb <- txtProgressBar(min=0, max=ndraws, style=3)
 

  # Loop over number of replications 
  for (i in seq(ndraws)) {    
    u <- sig*rnorm(t)                                  
    y <- beta0 + beta1*x1 + beta2*x2  + beta3*x3 + u   
      
    results <- optim(theta0, neglog, y = y, x1 = x1, x2 = x2, x3 = x3)
    p <- results$par
       
    H <- numhess(neglog,p,y,x1,x2,x3)
    theta <- p
    cov   <- inv(H)
        
    # One restriction
    R <- rbind(c(0,  1,  0,  0,  0))
    Q <- 0
    wd1[i] <- t* t( (R %*% theta - Q) ) %*% inv(R %*% cov %*% t(R)) %*% (R %*% theta - Q)

    # Two restrictions
    R <- rbind(c(0,  1,  0,  0,  0), 
               c(0,  0,  1,  0,  0))
    Q <- cbind(c(0,0))
    wd2[i] <- t* t( (R %*% theta - Q) ) %*% inv(R %*% cov %*% t(R)) %*% (R %*% theta - Q)
        
    # Three restrictions
    R <- rbind( c(0,  1,  0,  0,  0), 
         c(0,  0,  1,  0,  0),
         c(0,  0,  0,  1,  0) )      
    Q <- cbind(c(0,0,0))
    
    wd3[i] <- t* t( (R %*% theta - Q) ) %*% inv(R %*% cov %*% t(R)) %*%(R %*%theta - Q)
    
    setTxtProgressBar(pb, i)      
  }
  close(pb)
  
  #********************************************************************
  #***
  #***     Generate graph
  #***
  #********************************************************************
  
  figure()
  
  par(xaxs="i", yaxs="i")  
  hist(wd3, freq=FALSE, breaks = 50, 
       main = "Distribution of Wald Statistics",
       ylab = expression(f(W)),
       xlab = expression(W), bty="l",
      xlim = c(0.0001, 15))
  ygrid <- seq(0.0001, 15, 0.1)
  
  cpdf <- dchisq(ygrid,3)
  lines(ygrid, cpdf, type="l", lwd=2)
}
    

    


