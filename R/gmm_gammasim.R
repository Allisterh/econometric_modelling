# ===========================================================================
#
#      Monte Carlo properties of GMM estimator of the gamma distribution
#
# ===========================================================================

rm (list = ls(all=TRUE))
graphics.off()
set.seed(123, kind="Mersenne-Twister")

#
#------------------------- Helper Functions----------------------------------
#

# Load required functions - inv, trimr, reshapeg
source("EMTSUtil.R")

#-------------------------------------------------------------------------
# Negative log-likelihood function for gamma distribution with beta = 1
#-------------------------------------------------------------------------

neglog <- function( a,y ) {
  logl <-  (a-1)*mean(log(y)) - mean(y) - lgamma(a)
  return(-logl)
}

#-------------------------------------------------------------------------
# GMM criterion function for gamma distribution with beta = 1
#-------------------------------------------------------------------------
gmmcrit <- function(a,y) {
  m <- cbind((y-a),  (y^2 - a*(a+1)) )
  w <- t(m) %*% m/length(y)
  q <- 0.5*colMeans(m) %*% inv(w) %*% cbind(colMeans(m))
  return(q) 
}

#-------------------------------------------------------------------------
#  Monte Carlo experiment for different parameter inputs
#-------------------------------------------------------------------------

simexp <- function(t,reps,dist,a0v) {
  cols <- length(a0v)
    
    # Allocate arrays
  zeros <- array(0, c(reps, cols))
  aML  <- zeros
  a1   <- zeros 
  a2   <- zeros
  bML  <- zeros
  b1   <- zeros
  b2   <- zeros
  seML <- zeros
  se1  <- zeros
  se2  <- zeros
  seQ2 <- zeros
  tML  <- zeros
  t1   <- zeros  
  t2   <- zeros
  tQ2  <- zeros
  j2   <- zeros
  qmin <- zeros
  
  
  for (ac in seq(cols)) {
    for (rep in seq(reps)) {
      if (dist == "gam") {
         y <- rgamma(t, shape=a0v[ac],scale=1)
      } else if (dist == "exp") {
        y <- -log(runif(t))*a0v[ac]        
      } else {
         cat('\nWrong distribution')
      }
      # Maximum likelihood
      estResults <- optim(mean(y), neglog, y=y, method="BFGS", hessian=T)
      aML[rep, ac] <- estResults$par
      H <- estResults$hessian
      
      bML[rep,ac]  <- aML[rep,ac]-a0v[ac]
      seML[rep,ac] <- sqrt(1/(t*H))
      tML[rep,ac]  <- bML[rep,ac]/seML[rep,ac]
      
      # GMM using the first moment
      a1[rep,ac]  <- mean(y)           
      b1[rep,ac]  <- a1[rep,ac]-a0v[ac]
      se1[rep,ac] <- sd(y)/sqrt(t)
      t1[rep,ac]  <- b1[rep,ac]/se1[rep,ac]
      
      # GMM using two moments
      estResults <- optim(mean(y), gmmcrit, y=y, method="BFGS", hessian=T)
      a2[rep, ac] <- estResults$par
      qmin[rep, ac] <- estResults$value
      H <- estResults$hessian
      
      b2[rep,ac]  <- a2[rep,ac]-a0v[ac]      
      se2[rep,ac] <- sqrt(1/(t*H))
      t2[rep,ac]  <- b2[rep,ac]/se2[rep,ac]
      
      m            <- cbind((y-a2[rep,ac]),  (y^2-a2[rep,ac]*(a2[rep,ac]+1)))
      WTi          <- inv(t(m) %*% m/nrow(m))
      D            <- cbind(c(1, (2*a2[rep,ac]+1)))
      seQ2[rep,ac] <- 1/sqrt(t* t(D) %*% WTi %*% D)
      tQ2[rep,ac]  <- b2[rep,ac]/seQ2[rep,ac]          
      j2[rep,ac]   <- 2*t*qmin[rep,ac]
    }
  
  }
  cat('\n Distribution    ', dist) 
  cat('\nResults for Maximum Likelihood\n')  
  print(cbind(
    "Bias"=colMeans(bML), 
    "Std err"=colMeans(seML), 
    " t > 1.96"=colMeans(abs(tML)>1.96)))
  
  cat('\n')
  cat('\nResults for GMM (1 moment condition)')
  cat('\nResults for Maximum Likelihood\n')  
  print(cbind(
    "Bias"=colMeans(b1), 
    "Std err"=colMeans(se1), 
    " t > 1.96"=colMeans(abs(t1)>1.96)))
  
  cat('\n ')
  cat('\nResults for GMM (2 moment conditions)')
  cat('\nResults for Maximum Likelihood\n')  
  print(cbind("Bias"=colMeans(b2), 
              "Std err"=colMeans(se2), 
              " t > 1.96"=colMeans(abs(t2)>1.96), 
              "Std err(Q)"=colMeans(seQ2), 
              " t(Q) > 1.96"=colMeans(abs(tQ2)>1.96), 
              "J > 3.84"=colMeans(j2>3.84)))  
  }



#
#--------------- Monte Carlo Evidence for the Gamma model -------------------
#
gmm_gammasim <- function() {
  # Set the parameters of the Monte Carlo experiment  
  reps <- 200
    
  # Now uncomment only one of the following options
  # Estimation    
  dist <- 'gam' 
  t <- 100 
  a0v <- c(1,2,3,4,5)
  simexp(t,reps,dist,a0v)
}
