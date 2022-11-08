#=========================================================================
#
#    Simulating the size of the wald test using exponential regression
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
# Unconstrained log-likelihood function at each observation
#-------------------------------------------------------------------------
lnlt <- function(theta,y,x) {
  beta0 <- theta[1]
  beta1 <- theta[2]
  mue   <- exp(beta0 + beta1*x)
  ln1   <- -log(mue) - y/mue  
  return(ln1)
}


#-------------------------------------------------------------------------
# Unconstrained log-likelihood function
#-------------------------------------------------------------------------
lnl <-function(theta,y,x) {
  lnlt1 <- lnlt(theta,y,x)
  ln <- -mean(lnlt1)
  return (ln)
}

#-------------------------------------------------------------------------
# The Hessian
#-------------------------------------------------------------------------

analytic_hessian <- function (theta,y,x) {
  beta0  <- theta[1]
  beta1  <- theta[2]
  mue    <- exp(beta0 + beta1*x)
  H      <- array(0, c(2,2))
  H[1,1] <- -sum(y/mue)
  H[1,2] <- -sum((x*y)/mue )
  H[2,1] <- H[1,2]
  H[2,2] <- -sum( ((x^2)*y)/mue ) 
  return(H)
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

test_size <- function( ) {
  beta0  <- 1                                   
  beta1  <- 0                                         
  t      <- c(5,10,25,100)
  ndraws <- 10000
  wd     <- array(0, c(ndraws,1))
  chi2   <- qchisq(0.95,1)

  # table array to store the results
  rnames <- c("Sample T", "Size", "Critical Value (5%)")
  cnames <- c(rep("", length(t)))
  res <- array(dim = c(length(rnames), length(cnames)), 
                dimnames = list(rnames, cnames))
  for (l in seq(t)) {
    j  <- 0
    x  <- rnorm(t[l])     # Explanatory variable (fixed in repeated samples)
      
    # Main do loop to generate Monte Carlo results
    pb <- txtProgressBar(min=0, max=ndraws, style=3)
    for (i in seq(ndraws)) {
      mue     <- exp(beta0 + beta1*x)  #mean
      u       <- runif(t[l])        
      y       <- -mue*log(1 - u)
        
      theta0  <- c(beta0, beta1)
      results <- optim(theta0, function(theta) {lnl(theta, y, x)}, method="BFGS")
      theta <- results$par
        
      H       <- analytic_hessian(theta,y,x)      
        
      #   One restriction
      R       <- rbind(c(0, 1))
      Q       <- 0
      wd[i,1] <- t[l]* t( (R %*% cbind(theta) - Q)) %*% inv( R %*% inv(-H) %*% t(R) ) %*% (R %*% cbind(theta) - Q)
      if (wd[i] > chi2) {
        j <- j+1         
      } 
      # Update progress
      setTxtProgressBar(pb, i)
    }
    
    wd_sort <- sort(wd)          
    
    # Store the results
    res[1,l] <- t[l]
    res[2,l] <- j/ndraws
    res[3,l] <- wd_sort[ndraws*0.95]
    
    close(pb)
  }
  print(res)  
}





