#================================================================================
#
#     Vuoung's Nonnested Test
#
#================================================================================

rm (list = ls(all=TRUE))
graphics.off()

#
#--------------------------- Helper Functions -----------------------------------
#
# load required functions - inv
source("EMTSUtil.R")

#--------------------------------------------------------------------------------
#   Likelihood function model 1
#--------------------------------------------------------------------------------
model1 <- function( b,y,x ) {
  u  <- y-b[1]*x[,1]-b[2]*x[,2]-b[3]*x[,3]
  lf <- -0.5*log(2*pi) - 0.5*log( b[4] ) - u^2/(2*b[4]) + log(abs(1.0))
  return(lf)  
}
#--------------------------------------------------------------------------------
#   Likelihood function model 1
#--------------------------------------------------------------------------------
model2 <- function( b,y,x ) {
  u  <- y - b[1]*x[,1] - b[2]*x[,2] - b[3]*x[,3]
  lf <- -0.5*log(2*pi) - 0.5*log( b[4] ) - u^2/(2*b[4])  + log( abs( 1.0/exp( y ) ) )
  return(lf)
}

#
#--------------------------- Vuong's model -----------------------------------
#

nls_money <- function() {
  # Load data
  # [cpi fedfunds gdp m2 tbill]
  load('moneydemand.Rdata')
  
  mt <- money[,4]/money[,1]
  yt <- money[,3]/money[,1]
  rt <- money[,5]/100

  t  <- length( mt )         # Define the sample size      

  # Estimate model 1
  y    <- mt                    # Estimate Model 1    
  x    <- cbind(rep(1, t), rt, yt)
  b    <- lm(y ~ x - 1)$coef
  sig2 <- (t( (y - x %*% b)) %*% (y - x %*% b))/t
  lf1  <- model1( c(b, sig2),y,x)

  # Estimate Model 2 
  y    <- log(mt)              	   
  x    <- cbind(rep(1,t), log(rt), log(yt))
  b    <- lm(y ~ x - 1)$coef
  sig2 <- (t( (y - x %*% b)) %*% (y - x %*% b))/t
  lf2  <- model2( c(b, sig2),y,x)

  # Perform Vuong's test
  dt   <- lf1 - lf2
  dbar <- mean( dt )
  sbar <- sd(dt)

  V <- sqrt(t)*(dbar/sbar)
    
  cat('\nVuongs test and p-value')
  cat('\n------------------------\n')
  cat(V, pnorm(V,0,1), '\n')      
}


