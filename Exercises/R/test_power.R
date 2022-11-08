#=========================================================================
#
#    Simulating the power of the Wald test using exponential regression
#    model
#
#=========================================================================

rm (list = ls(all=TRUE))
graphics.off()
set.seed(123457, kind="Mersenne-Twister")

#
#--------------------------- Helper Functions -----------------------------------
# 

#-------------------------------------------------------------------------
# Define the unconstrained log of the likelihood at each observations
#-------------------------------------------------------------------------

lnlt <- function(theta,y,x) {
  beta0 <- theta[1]
  beta1 <- theta[2]
  mue   <- exp(beta0 + beta1*x)
  ln1   <- -log(mue) - y/mue
  return(ln1)
}

#-------------------------------------------------------------------------
# Define the unconstrained log of the likelihood 
#-------------------------------------------------------------------------
lnl <- function(theta,y,x) {
  lnlt1 <- lnlt(theta,y,x)
  ln <- -mean(lnlt1)
  return(ln)
}

#-------------------------------------------------------------------------
# Define the analytic Hessian
#-------------------------------------------------------------------------

analytic_hessian <- function(theta,y,x) {
  beta0 <- theta[1]
  beta1 <- theta[2]
  mue  <- exp(beta0 + beta1*x)

  H      <- array(0, dim = c(2,2))
  H[1,1] <- -sum(y/mue)
  H[1,2] <- -sum((x*y)/mue )
  H[2,1] <- H[1,2]
  H[2,2] <- -sum( ((x^2)*y)/mue )
  return (H)
}


#-------------------------------------------------------------------------
# Wrapper function for matrix inverse
#-------------------------------------------------------------------------
inv <- function(M) {
  return(solve(M))
}

#
#--------------------------- Wald Statistics Tests -----------------------------------
#
test_power <- function() {
  beta0       <- 1            # Intercept parameter                                             
  t           <- 5            # Sample size                                 
  ndraws      <- 10000
  wd          <- array(0, c(ndraws,1))
  beta1_range <- c(-4, -3, -2, -1, 0, 1, 2, 3, 4)
  power       <- rep(0,length(beta1_range))

  for (k in seq(beta1_range)) {
    beta1 <- beta1_range[k]
    
    j  <- 0 
    x  <- rnorm(t)      # Explanatory variable (fixed in repeated samplea) 
    
     # Main do loop to generate Monte Carlo results
     pb <- txtProgressBar(min=0, max=ndraws, style=3)
     for (i in seq(ndraws)) {
        mue     <- exp(beta0 + beta1*x)  #mean
        u       <- runif(t)
        y       <- -mue*log(1 - u)
          
        theta0  <- c(beta0, beta1)
        results <- optim(theta0, function(theta) {lnl(theta, y, x)}, method="BFGS")
        theta <- results$par
          
        H       <- analytic_hessian(theta,y,x)
          
        # One restriction
        R       <- rbind(c(0,1))
        Q       <- 0
        wd[i,1] <- t* t( (R %*% cbind(theta) - Q)) %*% inv( R %*% inv(-H) %*% t(R) ) %*% (R %*% cbind(theta) - Q)
        if (wd[i] > 4.288) {
          j <- j+1
        }
        # Update progress
        setTxtProgressBar(pb, i)
    }
    power[k] <- j/ndraws      
    close(pb)
      
  }
  print(matrix(c(beta1_range, power), 
               nrow=2, byrow=T,
               dimnames = list(c("Beta1", "Power"), 
                               c(rep("", length(beta1_range))) ) ))
}


