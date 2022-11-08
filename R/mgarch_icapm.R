#============================================================================
#
#   Program to estimate a international capital asset pricing model 
#   with time-varying beta risk using the symmetric BEKK specification.
#
#============================================================================
rm(list = ls(all=T))
graphics.off()

#
#--------------------------- Helper Functions  ------------------------------
# 

# Load require functions - inv, trimr, figure
source("EMTSUtil.R")

#----------------------------------------------------------------------------
# Log-liklihood function for symmetric BEKK model
#----------------------------------------------------------------------------
neglog <- function( b,y ) {
  t <- nrow(y)
  n <- ncol(y)
  f      <- rep(0, t)
  
  h <- cov(y)     # Initialise  conditional covariance matrix
  u <- apply(y, 2, sd)     # Initialise the disturbance vector
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

#----------------------------------------------------------------------------
# Calculates the beta risk from symmetric BEKK model
#----------------------------------------------------------------------------
bekk.beta <- function( b,y ) {
  t <- nrow(y)
  n <- ncol(y) 
  br <- rep(0, t)
  
  h <- cov(y)     # Initialise  conditional covariance matrix
  u <- apply(y, 2, sd)     # Initialise the disturbance vector
  
  c  <- matrix(c(b[1],0.0,
                 b[2], b[3]), nrow=2, byrow=T)
  
  a <-  matrix(c(b[4],b[5],
                 b[5], b[6]), nrow=2, byrow=T)
  
  d <-  matrix(c(b[7],  b[8],
                 b[8], b[9]), nrow=2, byrow=T)
  
  for (i in seq(t)) {
    m    <- b[c(10, 11)]                   # Update conditional mean
    u    <- y[i,]- m                   # Update residuals   
    h    <- c %*% t(c) + a %*% (u %*% t(u)) %*% t(a) + d %*% h %*% t(d)  # Update conditional covariance matrix     
    br[i] <- h[1,2]/h[2,2] 
  }  
  return(br)
}

#
#--------------------------- MGARCH ICAPM model  ----------------------------
#

mgarch_icapm <- function ()
{
  # Daily data starting 3-Feb-1988 to 29-Dec-1995
  #  1.	date
  #	2.	interest rate (risk free)
  #	3.	nyse returns (NYSE stock index)
  #	4.	world returns (from MSCI index)
  
  rdata <- read.table("icapm.dat")
  r <- 100*(rdata[,3] - rdata[,2])      # NYSE excess return          
  m <- 100*(rdata[,4] - rdata[,2])      # World excess return         
  t <- length(r)                                    
  y <- cbind(r,  m )
  
  
  # Constant beta risk estimate (based on OLS)      
  x     <- cbind(rep(1, t), m)
  
  b_ols <- lm(r ~ x - 1)$coef
  
  cat('\nConstant estimate of beta risk <- ', b_ols[2])
  
  # Estimate the BEKK Symmetric MGARCH(1,1) model             
  start <- c(0.065125621828620,
             -0.082834944191625,
             0.054312771634006,
             0.141289906762252,
             -0.037479056981814,
             0.257666483252922,
             0.974472477018916,
             0.024314746700112,
             0.946834227119610,
             0.063159015564876,
             0.040511737867333 )
  estResults <- optim(start, neglog, y=y, method="BFGS")
  theta <- estResults$par
  
  # Compute beta risk at optimal parameters
  br <- bekk.beta( theta,y )
  
  # Plot conditional beta risk
  figure()
  
  plot( seqa(1988+2/12,1/257,t) , br, type="l", 
        main = 'Conditional Beta Risk',              
        xlab = "t",
        ylab = expression(beta[t]),
        bty = "l")
  
}

