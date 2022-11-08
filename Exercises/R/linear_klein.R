#==========================================================================
#
#        Program to estimate Klein's macroeconomic model by full 
#	  	information maximum likelihood.
#
#==========================================================================

rm (list = ls(all=TRUE))
graphics.off()
set.seed(1234567, kind="Mersenne-Twister")

#
#--------------------------- Helper Functions -----------------------------------
#

# load utility functions - inv 
source("EMTSUtil.R")

#-----------------------------------------------------------------
# Log-likelihood function
#-----------------------------------------------------------------
        
neglog <- function(theta,y,x) {
   lf <- -mean(lnlt(theta,y,x))
   return(lf)
}

#-----------------------------------------------------------------
# Log-likelihood function at each observation 
#-----------------------------------------------------------------
lnlt <- function (theta,y,x) {
  t <- nrow(y)
  n <- ncol(y)
  b <- matrix(c(1-theta[2],  	-theta[6],		-theta[10],
    	         -theta[2], 		1-theta[6], 	-theta[10],
              theta[2]-theta[4], 	-theta[6], 		1), byrow = T, ncol = 3)

  a <- matrix(c(-theta[1], 	-theta[5], -theta[9],
                  -theta[2], -theta[6], -theta[10],
                   theta[2], theta[6],   0,
                   -theta[4], 0,        0,
                     0,       0,     -theta[12],
                  -theta[3], -theta[7],	0,
                     0, 		  0, 		  -theta[11],
                     0, 	 -theta[8],  0), byrow = T, ncol = 3)
  
  u <- array(0, c(t,n))
  for (i in seq(t)) {
     u[i,] <- y[i,] %*% b + x[i,] %*% a    
  }
  omega <- t(u) %*% u/t                         # Concentrate out resid var-covar matrix  
  lnl <- array (0, c(t,1))

  for (i in seq(t)) {
    lnl[i] <- -n*0.5*log(2*pi) + log(abs(det(b))) - 0.5*log(det(omega)) - 0.5*u[i,] %*% inv(omega) %*% cbind(u[i,])
  }
  return(lnl)  
}

#-----------------------------------------------------------------
# Instrumental variable estimation
#-----------------------------------------------------------------
iv <- function(y,w,x) {
  
  tmp0   <- inv(t(w) %*% x %*% inv(t(x) %*% x) %*% t(x) %*% w)
  tmp1   <- (t(w) %*% x %*% inv( t(x) %*% x) %*% t(x) %*% y)
  b      <- tmp0 %*% tmp1                                             # IV estimates 
  # Standard error of regression
  e       <- y - w %*% b  
  t       <- length(y)
  sigma   <- sqrt(t(e) %*% e/t)
  
  # Variance-covariance matrix
  vcov    <- c(sigma^2) *tmp0
  sterr   <- sqrt(diag(vcov))
  tstats  <- b/sterr
  return(b)
}

#-----------------------------------------------------------------
# Create a results table
#-----------------------------------------------------------------
create.empty.table <- function() {
  rnames <-  c("alpah0", "alpah1","alpah2","alpah3",
           "beta0", "beta1","beta2","beta3", 
            "gam0", "gam1","gam2","gam3")
  cnames <- c("OLS", "IV", "FIML")
  results <- matrix(nrow = length(rnames), ncol = length(cnames), 
                    dimnames = list(rnames, cnames))
  return(results)
}
        

#
#---------------------------  Klein's Model -----------------------------------
#

linear_klein <- function () {
  t <- 22

  # Read in the data
  klein_data <- read.table("klein.dat") # [22,10]   


  ct    <- klein_data[,1]
  p      <- klein_data[,2]
  pw     <- klein_data[,3]
  i      <- klein_data[,4]
  klag   <- klein_data[,5]
  d      <- klein_data[,6]
  income <- klein_data[,7]
  gw     <- klein_data[,8]
  g      <- klein_data[,9]
  tax    <- klein_data[,10]
  
  trend <- seq(-11, -11+t-1, 1)

  # Estimate consumption function by OLS 
  y <- cbind(ct[-1])
  pwgw <- pw + gw

  x <- matrix(c(rep(1, t-1), p[-1], p[-length(p)], pwgw[-1]), ncol=4)
  alpha_ols <- lm(y ~ x - 1)$coef
  
  # Create a table of estimates to compare
  results <- create.empty.table()
  results[1:4,"OLS"] <- alpha_ols
  
  # Estimate investment function by OLS
  y <- i[-1]
  x <- matrix(c(rep(1, t-1), p[-1], p[-length(p)], klag[-1]), ncol = 4)
  beta_ols <- lm(y ~ x - 1)$coef
  results[5:8,"OLS"] <- beta_ols


  # Estimate private wage function by OLS 
  y <- pw[-1]
  x <- matrix(c(rep(1, t-1), d[-1], d[-length(d)], trend[-1]), ncol = 4)
  
  gam_ols <- lm(y ~ x - 1)$coef
  results[9:12,"OLS"] <- gam_ols
  
  # Define the endogenous variables for the system
  cipw <- matrix(c(ct, i, pw), ncol = 3)
  y <- cipw[-1,]

  # Define exogenous/predetermined variables (instruments) for the system    
  x <- matrix(c(rep(1, t-1), g[-1], tax[-1], gw[-1], trend[-1],
              p[-length(p)], d[-length(d)], klag[-1]), ncol = 8)

  # Estimate consumption function by IV
  pwgw <- pw + gw

  tmpy <- ct[-1]
  tmpw <- matrix(c(rep(1, t-1), p[-1], p[-length(p)], pwgw[-1]), ncol = 4)

  alpha_iv <- iv(tmpy, tmpw , x)
  results[1:4,"IV"] <- alpha_iv

  # Estimate investment function by IV
  tmpy <- i[-1]
  tmpw <- matrix(c(rep(1, t-1), p[-1], p[-length(p)], klag[-1]), ncol = 4)
  beta_iv <- iv(tmpy, tmpw , x)
  results[5:8,"IV"] <- beta_iv

  # Estimate private wage function by IV 
  tmpy <- pw[-1]
  tmpw <- matrix(c(rep(1, t-1), d[-1], d[-length(p)], trend[-1]), ncol = 4)
  gam_iv <- iv(tmpy, tmpw, x)
  results[9:12,"IV"] <- gam_iv

  t <- nrow(y)   	# Redefine the sample to allow for lags in the system
  theta_0 <- c(alpha_iv, beta_iv, gam_iv)  # Use IV estimates as starting values     

  # Estimate the model using BFGS algorithm and compute QMLE se 
  estResults <- optim(theta_0, neglog, y = y, x = x, method = "BFGS")
  theta <- estResults$par
  a0 <- estResults$value 
  
  results[1:12,"FIML"] <- theta
  
  print(results)  
  
  cat('\n\nLog-likelihood = ', -a0)

  t <- nrow(y)
  n <- ncol(y) 
  
   b <- matrix(c(1-theta[2],    -theta[6],		-theta[10],
    	         -theta[2], 		1-theta[6], 	-theta[10],
              theta[2]-theta[4], 	-theta[6], 		1), byrow = T, ncol = 3)

  a <- matrix(c(-theta[1], 	-theta[5], -theta[9],
                  -theta[2], -theta[6], -theta[10],
                   theta[2], theta[6],   0,
                   -theta[4], 0,        0,
                     0,       0,     -theta[12],
                  -theta[3], -theta[7],	0,
                     0, 		  0, 		  -theta[11],
                     0, 	 -theta[8],  0), byrow = T, ncol = 3) 
  
  u <- array(0, c(t,n))
  for (i in seq(t)) {
     u[i,] <- y[i,] %*% b + x[i,] %*% a    
  }
  rvc <- t(u) %*% u/t                         # Concentrate out resid var-covar matrix
        
  cat('\nResidual variance-covariance matrix\n')
  print(rvc)
  cat('\n\n')
}





