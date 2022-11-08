#=============================================================================
#
#     Program to estimate a symmetric BEKK MGARCH model of US yields
#       US daily zero coupon yields are expressed in percentages, 
#       starting 10/04/1988 and ending 12/28/2001
#
#=============================================================================

rm(list = ls(all=T))
graphics.off()

#
#--------------------------- Helper Functions  ------------------------------
# 

# Load require functions - inv, trimr
source("EMTSUtil.R")

#-----------------------------------------------------------------------------
# Log-liklihood function for asymmetric BEKK uncontrained
#-----------------------------------------------------------------------------
negloga <- function( b,y ) {
  t <- nrow(y)
  n <- ncol(y)
  
  f <- rep(0, t)
  
  h <- cov(y)     # Initialise  conditional covariance matrix
  u <- apply(y, 2, sd)     # Initialise the disturbance vector
  #u <- t(u)
  
  c  <- matrix(c(b[1], 0.0 ,
                 b[2],   b[3]), nrow=2, byrow=T)
  
  a <-  matrix(c(b[4],b[6],
                 b[5],  b[7]), nrow=2, byrow=T)
  
  d <-  matrix(c(b[8],  b[10],
                 b[9], b[11]), nrow=2, byrow=T)
  for (i in seq(t)) {
    m    <- b[c(12, 13)]                   # Update conditional mean
    u    <- y[i,]- m                   # Update residuals                                
    f[i] <- -0.5*n*log(2*pi) - 0.5*log(det(h)) - 0.5*t(u) %*% inv(h) %*% u        
    h    <- c %*% t(c) + a %*% (u %*% t(u)) %*% t(a) + d %*% h %*% t(d)  # Update conditional covariance matrix
  }
  lf <- -mean( f )
  return(lf)
}


#--------------------------------------------------------------------------
# Log-liklihood function for symmetric BEKK model
#--------------------------------------------------------------------------
neglogs <- function( b,y ) {
  t <- nrow(y)
  n <- ncol(y)
  f      <- rep(0, t)
    
  h <- cov(y)     # Initialise  conditional covariance matrix
  u <- apply(y, 2, sd)     # Initialise the disturbance vector
  #u <- t(u)
    
  c  <- matrix(c(b[1],0.0,
                 b[2], b[3]), nrow=2, byrow=T)
  
  a <-  matrix(c(b[4],b[5],
               b[5], b[6]), nrow=2, byrow=T)
  
  d <-  matrix(c(b[7],  b[8],
               b[8], b[9]), nrow=2, byrow=T)
    
  for (i in seq(t)) {
    m    <- b[c(10, 11)]                   # Update conditional mean
    u    <- y[i,]- m                   # Update residuals   
    f[i] <- -0.5*n*log(2*pi) - 0.5*log(det(h)) - 0.5*t(u) %*% inv(h) %*% u       
    h    <- c %*% t(c) + a %*% (u %*% t(u)) %*% t(a) + d %*% h %*% t(d)  # Update conditional covariance matrix     
  }
  lf <- -mean( f )    
  return(lf)
}

#--------------------------------------------------------------------------
# Log-liklihood function for constant covariance BEKK model
#--------------------------------------------------------------------------
neglogc <- function( b,y ) {
  t <- nrow(y)
  n <- ncol(y) 
  f      <- rep(0, t)
  
  h <- cov(y)     # Initialise  conditional covariance matrix
  u <- apply(y, 2, sd)     # Initialise the disturbance vector
  #u <- t(u)
  
  c  <- matrix(c(b[1],0.0,
                 b[2], b[3]), nrow=2, byrow=T)
  
  a <-  matrix(c(b[4],0.00,
                 0.00, b[5]), nrow=2, byrow=T)
  
  d <-  matrix(c(b[6],  0.00,
                 0.00, b[7]), nrow=2, byrow=T)
  
  for (i in seq(t)) {
    m    <- b[c(8,9)]                   # Update conditional mean
    u    <- y[i,]- m                   # Update residuals   
    f[i] <- -0.5*n*log(2*pi) - 0.5*log(det(h)) - 0.5*t(u) %*% inv(h) %*% u       
    h    <- c %*% t(c) + a %*% (u %*% t(u)) %*% t(a) + d %*% h %*% t(d)  # Update conditional covariance matrix
  }
  lf <- -mean( f )
  return(lf)
}



                          
#
#--------------------------- BEKK Model of US Yields  -----------------------
#                          
mgarch_bekk <- function( ) {
  # Load data
  load('yields_us.Rdata')
  
  # Choose variables
  r  <- rdata[,c(1, 2)]
  y  <- 100*(trimr( r,1,0 ) - trimr( r,0,1 )) 
  t  <- nrow( y )
  
  # Estimate the BEKK Aymmetric MGARCH(1,1) model             
  start <- c(1.0711572152073034,
             0.2847761383915719,
             0.5666734919888730, 
             0.3849272367534938, 
             0.0079982965205830, 
             0.0492176759790856, 
             0.2235088706616157, 
             0.8810887055460768, 
             -0.0064483299253031, 
             0.0271203983383643, 
             0.9730134169147020, 
             0.0603486110799331, 
             0.0587667747700864)
  estResults <- optim(start, negloga, y=y, method="BFGS")
  thetaa <- estResults$par
  lfa <- estResults$val
  
  lfa <- -lfa
  
  cat('\nLikelihood function (asymmetric) = ',lfa)
  
  # Estimate the BEKK Symmetric MGARCH(1,1) model  
  start <- c(1.0827864786474668,
             0.5341723472962869,
             0.6326817767441521,
             0.4018538161228361,
             0.0241620281075959,
             0.2103118107333996,
             0.8982242699069298, 
             0.0048863543542688, 
             0.9609684149405032, 
             0.0713182051952838, 
             0.0474958290682067)
  estResults <- optim(start, neglogs, y=y, method="BFGS")
  thetas <- estResults$par
  lfs <- estResults$val
  
  lfs <- -lfs
  cat('\nLikelihood function (symmetric)   = ',lfs)
  
  # Compute conditional variance at optimal parameters
  h <- cov(y)     
  u <- apply(y, 2, sd)
  
  
  c  <- matrix(c(thetas[1],  0.0,
                 thetas[2],  thetas[3]), nrow=2, byrow=T)
  
  a <-  matrix(c(thetas[4], thetas[5],
                 thetas[5], thetas[6]), nrow=2, byrow=T)
  
  d <-  matrix(c(thetas[7], thetas[8],
                 thetas[8], thetas[9]), nrow=2, byrow=T)
  m <- thetas[c(10, 11)]
  
  cov_mat <- array(0, c( t,4 ))
  for (i in seq(t)) {
    u <- y[i,]- m                  # Update residuals    
    h <- c %*% t(c) + a %*% (u %*% t(u)) %*% t(a) + d %*% h %*% t(d) # Update covariance matrix
    cov_mat[i,] <- as.vector(h)
  }
  
  #*********************************************************************
  #**     Generate graph of conditional variances
  #*********************************************************************
  
  vec <- seq(t)
  
  figure()
  par(mfrow=c(2,2))
  
  #--------------------------------------------------------#
  # Panel (a)
  plot(vec,cov_mat[,1],type="l",
       main =' (a) Variance of 3-month rate',
       ylab = expression(hat(h[t])),
       xlab = 't',
       bty = "l")
  
  #--------------------------------------------------------#
  # Panel (b)
  plot(vec,cov_mat[,4],type="l",
       main =' (b) Variance of 1-year rate',
       ylab = expression(hat(h[t])),
       xlab = 't',
       bty = "l")
  
  #--------------------------------------------------------#
  # Panel (c)
  plot(vec,cov_mat[,2],type="l",
       main =' (c) Conditional covariance',
       ylab = expression(hat(h[t])),
       xlab = 't',
       bty = "l")
  
  #--------------------------------------------------------#
  ccor <- cov_mat[,2]/( sqrt( cov_mat[,1]*cov_mat[,4] ) )
  # Panel (d)
  plot(vec,ccor,type="l",
       main =' (d) Conditional correlation',
       ylab = expression(hat(h[t])),
       xlab = 't',
       bty = "l")
  # Estimate the constant cov BEKK  MGARCH(1,1) model  
  start <- c(1.215169555478360,
                  0.781175687461714,
                  0.400847056654352,
                  0.401280406108247,
                  0.213852259080924,
                  0.899426993215549,
                  0.963984125877018,
                  0.073850907078739,
                  0.040752440808561)
  
  estResults <- optim(start, neglogc, y=y, method="BFGS")
  thetac <- estResults$par
  lfc <- estResults$val 
  
  lfc <- -lfc
  
  
  # LR test of symmetry    
  lr  <- -2*t*(lfs - lfa)
  
  cat('\n')
  cat('\nLR test of symmetry   = ',lr)
  cat('\np-value               = ',1-pchisq(lr,4))
  
  # LR test of constant covariance
  lr  <- -2*t*(lfc - lfa)
  
  cat('\n')
  cat('\nLR test of const. cov.   = ',lr)
  cat('\np-value                  = ',1-pchisq(lr,4))
}


