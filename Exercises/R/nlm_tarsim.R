#============================================================================
#
#   Program to demonstrate the sampling properties of TAR models
#
#============================================================================

rm (list = ls(all=TRUE))
graphics.off()
set.seed(12, kind="Mersenne-Twister")

#
#--------------------------- Helper Functions -------------------------------
#
# load required functions - trimr
source("EMTSUtil.R")

#----------------------------------------------------------------------------
#  Simulate a Threshold Autoregressive models
#----------------------------------------------------------------------------
tar <- function (gtype,nobs,phi1,beta1,gam,c,sig)
{
  y <- rep(0, nobs+100)
  u <- rnorm(nobs+100)
  
  for  (t in 2:(nobs+100)) {
    # SETAR
    if (gtype == 1)
      g <- y[t-1] > c 
                                                                     
    # STAR
    if (gtype == 2)
      g <- pnorm( gam*(y[t-1] - c) ) 
                                                                       
    # LSTAR
    if (gtype == 3)
      g <- 1/( 1 + exp(-gam*(y[t-1] - c)) ) 
                                                                      
    # ESTAR
    if (gtype == 4)
      g <- 1 - exp(-gam*(y[t-1] - c)^2) 
                                                                      
    y[t] <- phi1*y[t-1] + beta1*y[t-1]*g + sig*u[t] 
  }
  y <- trimr(y,100,0)  
  return(y)  
}

#----------------------------------------------------------------------------
#  Compute the LM statistic to test a threshold autoregressive model 
#  assuming one lag in the auxiliary model
#----------------------------------------------------------------------------
tar_test <- function(yvar,p) {
  # First stage regression      
  y <- trimr(yvar,1,0)
  x <- cbind(rep(1, length(y)) , trimr(yvar,0,1))
  k <- ncol(x)
  u <- lm(y ~ x - 1)$residuals                                       
  
  # Second stage regression      
  if (p == 1)
    x <- cbind(x,  x[,2]*(x[,2]))
  else if (p == 2)
    x <- cbind(x,  x[,2]*(x[,2])^1, x[,2]*(x[,2]^2))
  else if (p == 3)
    x <- cbind(x , x[,2]*(x[,2]^1), x[,2]*x[,2]^2 , x[,2]*x[,2]^3)
  else if (p == 4)
    x <- cbind(x , x[,2]*(x[,2]^1), x[,2]*(x[,2]^3))
  e <- lm(u ~ x - 1)$residuals                                     
    
  # Compute LM statistic        
  r2  <- 1 - sum(e^2)/sum( (u-mean(u))^2 )
  lm  <- length(y)*r2
  dof <- ncol(x) - k  
  return(list (lm=lm, dof=dof))
}


nlm_tarsim <- function( )
{
  # Choose model to simulate
  # 1 = SETAR, 2 = STAR, 3 = LSTAR, 4 = ESTAR
  type <- 3 
  
  # Choose test type
  # p = 1: test based on regressing yt on {constant ylag, ylag*ylag}       
  # p = 2: test based on regressing yt on {constant ylag, ylag*ylag, ylag*ylag^2}   
  # p = 3: test based on regressing yt on {constant ylag, ylag*ylag, ylag*ylag^2, ylag*ylag^3}  
  # p = 4: test based on regressing yt on {constant ylag, ylag*ylag,              ylag*ylag^3}  
  p <- 3
  
  # Parameters
  c      <- 0.0                                                     
  phi1   <- -0.5                                           
  beta1  <- c(0.0, 0.5, 1.0)      
  gam    <- c(0.5, 1.0, 1.5)
  sig    <- 5                        
  nobs   <- 100                                                     
  ndraws <- 10000                                       
  
  # Initialise arrays
  lm       <- rep(0, ndraws)
  lm_power <- rep(0, length(beta1))
  
  for (j in seq(beta1)) 
  {
    for (k in seq(ndraws)) {
      #  Simulate a lstar model
      y <- tar(type,nobs,phi1,beta1[j],gam[1],c,sig)    
      res <- tar_test(y,p)
      lm[k] <- res$lm
      dof <- res$dof
    } 
    
    # Size
    if (j==1) {
      
      cve <- quantile(lm,0.95)
      cva <- qchisq(0.95,dof)
      pow <- 100*mean(lm>cve)   
      nom <- 100*mean(lm>cva)
      
      cat('\n')
      cat('\nSize of test: beta1             = ', beta1[j])
      cat('\n5%  critical value (empirical)  = ', cve) 
      cat('\n5%  critical value (asymptotic) = ', cva)
      cat('\nNominal size                    = ', nom)
      
      
      lm_power[j] <- pow
    }  
    else  {
      
      cat('\n')
      cat('\nPower of test: beta1        = ', beta1[j])
      cat('\nUnadjusted                  = ', 100*mean(lm>qchisq(0.95,dof)))
      cat('\nSize adjusted               = ', 100*mean(lm>cve))
      
      lm_power[j] <- 100*mean(lm>cve)
    } 
  }
}


