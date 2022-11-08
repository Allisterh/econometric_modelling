#============================================================================
#
#   Program to estimate the GENTS model using ftse equity returns
#
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()

#
#--------------------------- Helper Functions -------------------------------
#
# load required functions - figure. trimr
source("EMTSUtil.R")

#----------------------------------------------------------------------------
#  Log-likelihood function
#----------------------------------------------------------------------------
neglog <- function(b,y,x) {
  t   <- length(y)
  gam <- b[1]                          
  t1  <- b[2]*x                      
  t2  <- -0.5*(1+gam^2)*rep(1, t)
  t3  <- b[3]*x
  t4  <- -0.5*rep(1, t)
  
  eta <- rep(0, t)
  lnl <- rep(0, t)  
  # Range of integration includes the max and min values of yt  
  for (i in seq(t)){    
    eta[i]  <- log( integrate(gst, lower=-15, upper=15, gam,t1[i],t2[i],t3[i],t4[i],rel.tol=1e-8,abs.tol=1e-12)$value )  
    lnl[i] <- t1[i]*atan(y[i]/gam) + t2[i]*log(gam^2+y[i]^2) + t3[i]*y[i] + t4[i]*y[i]^2 - eta[i]    
  }
  
  lf <- -mean( lnl )
  
  return (lf)
}

#-------------------------------------------------------------------------
#  Compute eta for the generalised Student t distribution
#-------------------------------------------------------------------------
gst <- function(s,gam,t1,t2,t3,t4) {
  f <- exp( t1*atan(s/gam) + t2*log(gam^2 + s^2) + t3*s + t4*s^2 )
  return(f)
}

#-------------------------------------------------------------------------
#  Compute the conditional mean using the generalised Student t distribution
#-------------------------------------------------------------------------
gst1 <- function(s,gam,t1,t2,t3,t4,eta){
  f <- s*exp( t1*atan(s/gam) + t2*log(gam^2 + s^2) + t3*s + t4*s^2  - eta)
  return(f)
}

#-------------------------------------------------------------------------
#  Compute the conditional variance using the generalised Student t distribution
#-------------------------------------------------------------------------
gst2 <- function(s,gam,mm,t1,t2,t3,t4,eta) {
  tmp <- (s - mm)^2
  f   <- tmp*exp( t1*atan(s/gam) + t2*log(gam^2 + s^2) + t3*s + t4*s^2  - eta)
  return(f)
}

#-------------------------------------------------------------------------
#  Compute the conditional skewness using the generalised Student t distribution
#-------------------------------------------------------------------------
gst3 <- function(s,gam,mm,t1,t2,t3,t4,eta) {
  tmp <- (s - mm)^3
  f   <- tmp*exp( t1.*atan(s/gam) + t2*log(gam^2 + s^2) + t3*s + t4*s^2 - eta)
  return(f)
}

#-------------------------------------------------------------------------
#  Compute the conditional kurtosis using the generalised Student t distribution
#-------------------------------------------------------------------------
gst4 <- function(s,gam,mm,t1,t2,t3,t4,eta) {
  tmp <- (s - mm)^4
  f   <- tmp*exp( t1*atan(s/gam) + t2*log(gam^2 + s^2) + t3*s + t4*s^2 - eta)
  return(f)
}


nlm_gents <- function( ) {
  
  # Load data: daily returns on ftse 20 Nov. 1973 to 23 July 2001
  # Scale by 100 to convert to percentages   
  data <- as.matrix(read.table("ftse.dat"))
  y <- 100*data
  
  #Treat outliers at observations 222, 233, 3527, 3528 and 3529
  dum1       <- rep(0, length(y))
  dum1[222]  <- 1
  dum2       <- rep(0, length(y))
  dum2[233]  <- 1
  dum3       <- rep(0, length(y))
  dum3[3527] <- 1
  dum4       <- rep(0, length(y))
  dum4[3528] <- 1
  dum5       <- rep(0, length(y))
  dum5[3529] <- 1
  
  # Rename y as the residuals from a dummy variable regression
  d <- cbind(rep(1, length(y)), dum1,   dum2,   dum3,   dum4,   dum5)
  y <- residuals(lm(y ~ d - 1))
  
  # Current and lagged returns
  x <- trimr(y,0,1)                                                 
  y <- trimr(y,1,0)                              
  t <- length(y)
  
  # Estimate model
  theta_0 <- c(0.928834745978081,
               -0.700513899637738,
               0.695374988151895)
  estResults <- optim(theta_0, neglog, y=y, x=x, method="BFGS", hessian=T)
  theta <- estResults$par
  hess <- estResults$hess
  
  gam <- theta[1]
  vc  <- (1/t)*inv(hess)
  
  # Perform Wald test of skewness   
  r <- matrix(c(0,   1,   0,
                0,   0,   1), nrow=2, ncol=3, byrow=T)
  q <- rbind(0,0)
  w <- t( (r %*% theta - q) ) %*% inv(r %*% vc %*% t(r)) %*% (r %*% theta - q)
  
  cat('\n')
  cat('\nWald test of skewness = ',w)
  cat('\np-value               = ',1-pchisq(w,2))
  
  # Compute conditional moments    
  y_mean <- rep(0, t)
  y_var  <- rep(0, t)
  y_skew <- rep(0, t)
  y_kurt <- rep(0, t)
  
  gam <- theta[1]
  t1  <- theta[2]*x                      
  t2  <- -0.5*(1+gam^2)*rep(1, t)
  t3  <- theta[3]*x
  t4  <- -0.5*rep(1, t)
  eta <- rep(0, t)
  for (i in seq(t)){
    eta[i]  <- log( integrate(gst, lower=-15, upper=15, gam,t1[i],t2[i],t3[i],t4[i],rel.tol=1e-8,abs.tol=1e-12)$value )
    y_mean[i] <- integrate(gst1,-15, 15, gam,t1[i],t2[i],t3[i],t4[i],eta[i])$value
    y_var[i]  <- integrate(gst2, -15, 15, gam,y_mean[i],t1[i],t2[i],t3[i],t4[i],eta[i])$value
  }
  
  figure()
  plot(seq(t),y_mean, type="l",
       main = "", xlab="", ylab="",
       bty = "l") 
  
  figure()
  
  plot(seq(t),y_var, type="l", 
       main="",
       xlab = "",
       ylab = "",
       bty = "l")
}


