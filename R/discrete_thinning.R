#============================================================================
  #
#   Program to simulate the Poisson autoregressive model 
#   using the binomial thinning operator.
#
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()
set.seed(123, kind="Mersenne-Twister")

#
#------------------------ Binomial thinning model  --------------------------
#
discrete_thinning <- function() {
  # Parameters
  t   <- 10                                                    
  rho <- 0.3                                 
  lam <- 3.5
  
  # Generate poisson random numbers
  u <- rpois(lambda=lam, n=t)
  
  # Initialize y by choosing the median of draws from the unconditional distribution 
  # to ensure that the starting value is positive)  
  y <- rep(1,t)*median(rpois(lambda=lam/(1-rho),n=t))
  
  for (i in 2:t) {
    # Generate y_t-1 uniform random numbers   
    e <- runif(y[i-1])             
    
    # Sum bernoulli random draws  
    b_thin <- sum( e < rho )      
    
    # Generate realisation of y at time t   
    y[i] <- b_thin + u[i]           
    
    cat('\nIteration     = ',i)
    cat('\ny(i-1)        = ',y[i-1])
    cat('\nb_thin        = ',b_thin)
    cat('\nu(i)          = ',u[i])
    cat('\ny(i)          = ',y[i])
    cat('\n')
  }  
}
