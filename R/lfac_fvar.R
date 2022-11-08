# ===========================================================================
#
#   Algorithm to estimate the FVAR model of Stock and Watson 
#   using 30 yields for the US.
#   US daily yields starting Jan 10th 1990 and ending May 31 2005
# ===========================================================================
rm (list = ls(all=TRUE))
graphics.off()

#
#--------------------------- Helper Functions -------------------------------
#

# Load required functions - trimr
source("EMTSUtil.R")

#
#--------------------------- F-VAR Model -----------------------------------
#

lfac_fvar <- function() {
  # Load the data 
  # Yields 1yr, 2yr ... to 30yr (30 series) 
  load('daily_finance.Rdata')
  
  # Compute spreads
  yields <- yields[,1:30]
    
  #Standardize the data 
  ry <- nrow(yields)
  cy <- ncol(yields)
    
  y <- array(0, c(ry,cy))
  for (i in seq(cy)) {
     y[,i] <- (yields[,i] - mean( yields[,i] ))/sd(yields[,i] )  
  }
  
  # Number of dependent variables
  rr <- nrow(y)
  n <- ncol(y)  # Number of dependent variables  
  k <- 3        # Number of factors (chosen to be 3 here)                 
  
  lam <- array(0, c( n,k ))
  gam <- array(0, c( n,n ))
    
  lamnew <- lam
  gamnew <- gam
  
  maxiter <- 100       # Maximum number of iterations   
  tol     <- 1e-3      # Convergence criterion
  
  # Estimate F-VAR parameters using the Stock-Watson algorithm  
  iter <- 1
  pb <- txtProgressBar(min=0, max=maxiter, style=3)
  while (iter <= maxiter) {
    # Principal component analysis       
    va <- eigen(cor(y))$vectors
    
    s  <- y %*% va[, 1:k]    # First k principal components
    
    # Regressions on measurement equations
    for (i in seq(n)) {
      ylag <- trimr(y[,i],0,1)
      b <-  lm(trimr(y[,i],1,0) ~ cbind(trimr(s,1,0), ylag) - 1)$coef
      lam[i,] <- b[1:k]
      gam[i,i] <- b[k+1]
      y[,i] <-  y[,i] - gam[i,i]*( c(0.0,ylag) ) 
    }   
    
    # Check for convergence
    chk1 <- c(lam)-c(lamnew)
    chk2 <- c(gam)-c(gamnew)
    
    if (sqrt(sum(chk1)^2) < tol 
        && sqrt(sum(chk2)^2) < tol ) {      
      setTxtProgressBar(pb, 100)    
      break
    }
    else {
      lamnew <- lam
      gamnew <- gam    
    }
    iter <- iter + 1
    setTxtProgressBar(pb, iter)
  }  
  close(pb)
  # Regressions on factor equations 
  phi <- lm(trimr(s,1,0) ~ trimr(s,0,1) - 1)$coef
  
  # Scale lam so shock is a positive shock to interest rates
  lam <- -lam
  cat('\nConvergence achieved after iteration = ', iter)
  
  cat('\nLambda\n')
  print(lam)
  
  cat('\nGamma\n' )
  print(cbind(diag(gam)))
  
  cat('\nPhi\n')
  print(unname(phi))
  cat('\n')

  
  #**********************************************************************
  #***
  #***     Generate graphs
  #***
  #**********************************************************************
  
  figure()
  par(xaxs="i", yaxs="i")
  
  matplot(seq(n),lam,type="l",
       main = "Three Factors",
       xlab = "Maturity",
       ylab = "",
       ylim = c(-0.5, 0.5),
       col = c(4:6),
       lty = c(1:3),
       bty = "l") 
  legend("topright",
         legend=c("Level", "Slope", "Curvature"),
         lty=c(1:3),
         lwd=c(1,1,1),
         col=c(4:6))
}

