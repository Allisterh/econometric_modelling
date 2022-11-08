#============================================================================
#
#   Estimate the ACH model of Hamilton-Jorda for AMR on Aug01, 2006 
#
#   This program can't reproduce the GAUSS results reported in the book
#   because the likelihood function is flat relative to 
#   the delta1 parameter in particular. Consequently the results
#   are particularly sensitive to choice of algorithm and starting values.
#
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()
#
#--------------------------- Helper Functions -------------------------------
# 
# Load required functions -  trimr, recserar
source("EMTSUtil.R")

#----------------------------------------------------------------------------
#  The log-likelihood of the restricted ACH model
#----------------------------------------------------------------------------
neglog0 <- function(b,y,u,ubar) {
  # ACD(1,1) model
  si <- recserar( cbind(b[1] + b[2]*trimr(c(0.0, u),0,1)),cbind(ubar),cbind(b[3]))     
  
  # Hazard rate   
  h <- 1/(1 + exp(si))      
  
  lf <- -mean( y*log(h) + (1 - y)*log(1 - h) )
  return(lf)  
}


#----------------------------------------------------------------------------
#  Return h and psi for the restricted ACH model
#----------------------------------------------------------------------------
 hazard0 <- function(b,u,ubar) {
   # ACD(1,1) model
   si <- recserar( cbind(b[1] + b[2]*trimr(c(0.0, u),0,1)),cbind(ubar),cbind(b[3]))     
   
   # Hazard rate   
   h <- 1/(1 + exp(si))
   return(h)
 }

#----------------------------------------------------------------------------
#  The log-likelihood of the unrestricted ACH model
#----------------------------------------------------------------------------
neglog <- function (b,y,u,ubar,firm_lag,macro_lag) {
  # ACD(1,1) model
  si <- recserar( cbind(b[1] + b[2]*trimr(c(0.0, u),0,1)),cbind(ubar),cbind(b[3]))     
  
  # Hazard rate   
  h <- 1/(1 + exp(si + b[4]*firm_lag + b[5]*macro_lag))      
  
  lf <- -mean(  y*log(h) + (1 - y)*log(1 - h) )
  return(lf)
}

#----------------------------------------------------------------------------
#  Return h and si for the unrestricted ACH model
#----------------------------------------------------------------------------
hazard <- function (b,u,ubar,firm_lag,macro_lag) {
  # ACD(1,1) model
  si <<- recserar( cbind(b[1] + b[2]*trimr(c(0.0, u),0,1)),cbind(ubar),cbind(b[3]))     
  
  # Hazard rate   
  h <- 1/(1 + exp(si + b[4]*firm_lag + b[5]*macro_lag))
  return(h)
}
#
#--------------------------- ACH model of Trade -----------------------------
#

discrete_trade <- function() {
  # Load data for AMR on Aug 1, 2006.
  # Order of the variables:
  #   1.  Hour
  #   2.  Minute
  #   3.  Second
  #   4.	Time (time effectively starts at 9.30am and ends at 4.00pm with the time interval being one second)
  #   5.	y (a binary variable giving a one if a trade has occured, and 0 otherwise
  #   6.	N (a counter variable which increases by one when a trade occurs)
  #   7.	u (gives the duration between last trade)
  #   8.	Firm News dummy (AMR)
  #   9.	Macro News dummy
  
  data <- as.matrix(read.table("amr_aug1.dat"))
  
  # Trim he first 32 observations and exclude the last observation
  # Ensure postive durations
  data <- trimr(data,32,1)
  
  # Define variables
  y     <- data[,5]
  nt    <- data[,6]
  u     <- data[,7]
  firm  <- data[,8]
  macro <- data[,9]
  
  ubar      <- mean(u)           # Average duration time                      
  firm_lag  <- trimr(firm,0,1)   # Lagged firm news variable                              
  macro_lag <- trimr(macro,0,1)  # Lagged macroeconomic news variable                         
  y         <- trimr(y,1,0)      # Lineup other data  
  u         <- trimr(u,1,0)
  t         <- length(y)
  
  # Descriptive statistics
  cat('\nTotal number of observations             = ',t) 
  cat('\nTotal number of trades                   = ',sum(y))
  cat('\nAverage number of trades (per second)    = ',mean(y))
  cat('\nAverage duration time (in seconds)       = ',mean(u))
  cat('\nMinimum duration time (in seconds)       = ',min(u))
  cat('\nMaximum duration time (in seconds)       = ',max(u))
  cat('\nUnconditional Hazard rate                = ',1/mean(u))
  cat('\n ')
  
  
  # Estimate the restricted model (without news variables)
  theta_0 <-   c(0.0196493517339459,
                 0.0010464918437577, 
                 0.9840278397571428 )
  estResults <- optim(theta_0, neglog0, y=y, u=u, ubar=ubar, method="BFGS")
  theta0 <- estResults$par
  l0 <- estResults$val
  
  l0 <- -l0   
  cat('\nRestricted parameter estimates = ', theta0)
  cat('\nLog-likelihood (restricted) = ',l0) 
  cat('\nTxLog-likelihood        = ',t*l0)
  cat('\n ')
  
  
  # Compute h at optimal parameters
  h <- hazard0(theta0,u,ubar) 
  cat('\nACH Model without annoucement variables')
  cat('\nHazard rate                         = ',mean(h)) 
  cat('\nTime to the next trade (in seconds) = ',1/mean(h))
  cat('\n')
  
  # Estimate the unrestricted model (with news variables)  
  theta_0 <-   c(theta_0[1],
                 theta_0[2], 
                 theta_0[3], 
                 0.1,
                 0.1)
  estResults <- optim(theta_0, neglog, y=y,u=u,ubar=ubar,firm_lag=firm_lag, macro_lag=macro_lag, method="BFGS")
  theta1 <- estResults$par
  l1 <- estResults$val
    
  l1 <- -l1  
  cat('\nUnrestricted parameter estimates = ', theta1)
  cat('\nLog-likelihood (unrestricted) = ',l1)
  cat('\nTxLog-likelihood        = ',t*l1)
  cat('\n')
  
  # # Compute h at optimal parameters
  h <- hazard(theta1,u,ubar,firm_lag,macro_lag) 
  cat('\nACH Model without annoucement variables')
  cat('\nHazard rate                         = ',mean(h))
  cat('\nTime to the next trade (in seconds) = ',1/mean(h))
  cat('\n ')
  
  
  
  cat('\nACH Model without announcement variables: general model')
  h_adj <- 1/(1 + exp(mean(si)))         
  cat('\nHazard rate                         = ',h_adj)
  cat('\nTime to the next trade (in seconds) = ',1/h_adj)
  cat('\n')
  
  cat('\nACH Model without firm news variables')
  h_adj <- 1/(1 + exp(mean(si) + theta1[5]))    
  cat('\nHazard rate                         = ',h_adj)
  cat('\nTime to the next trade (in seconds) = ',1/h_adj)
  cat('\n ')
  
  
  cat('\nACH Model without macro news variables')
  h_adj <- 1/(1 + exp(mean(si) + theta1[4]))    
  cat('\nHazard rate                         = ',h_adj)
  cat('\nTime to the next trade (in seconds) = ',1/h_adj)
  cat('\n ')
  
  cat('\nACH Model with both firm and macro news variables')
  h_adj <- 1/(1 + exp(mean(si) + theta1[4] + theta1[5]))    
  cat('\nHazard rate                         = ',h_adj)
  cat('\nTime to the next trade (in seconds) = ',1/h_adj)
  cat('\n ')
}

