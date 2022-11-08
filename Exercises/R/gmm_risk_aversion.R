# ===========================================================================
#
# GMM estimation of risk aversion using the Mehra-Prescott data
#    
# ===========================================================================

rm (list = ls(all=TRUE))
graphics.off()

#
#------------------------- Helper Functions----------------------------------
#

# Load required functions - inv, numgrad, figure
source("EMTSUtil.R")
# Load required library - repmat
library("matlab")


#----------------------------------------------------------------------------#
# GMM objective function   
#----------------------------------------------------------------------------#   
q <- function(theta,rc,rb,rs,inst,flag) {
  
  beta <- theta[1]
  gam  <- theta[2]
  m1   <- (beta*(1 + rc)^(-gam)*(1 + rb) - 1) * inst
  m2   <- (beta*(1 + rc)^(-gam)*(1 + rs) - 1) * inst
  m    <- cbind(m1, m2)
  g    <- colMeans(m)
  if (flag)
    w <- t(m) %*% m/nrow(m)
  else
    w <- eye(ncol(m))
  
  ret  <- 0.5*t(g) %*% inv(w) %*% g
  
  return(ret)
}

#
#--------------- Risk Aversion and the Equity Premium Puzzle ----------------
#
gmm_risk_aversion <- function()
{
  # Load Mehra-Prescott data (annual for the period 1889 to 1978)
  load('equity_mp.RData')
  
  yr <- data[,1]      #     Year                                            
  r  <- data[,2]      #     Nominal risk free rate, yield percentage p.a.   
  p  <- data[,3]      #     Price of consumption goods                      
  c  <- data[,4]      #     Real per capita consumption                     
  d  <- data[,5]      #     Real dividend                                   
  s  <- data[,8]      #     Real stock price                                
  
  
  rs <- (s[2:91] + d[1:90] - s[1:90])/s[1:90]
  
  rb <- (1 + r[1:90]/100)*((p[1:90])/p[2:91]) -1 
  
  rc <- (c[2:91] - c[1:90])/c[1:90] 
  
  # Instruments
  const <- ones(length(rs)-1,1)
  inst  <- cbind(const,   trimr(rc,0,1))
  #inst <- [ const   trimr(rc,0,1)   trimr(rb,0,1) ]
  #inst <- [ const   trimr(rc,0,1)   trimr(rb,0,1)   trimr(rs,0,1) ] 
  
  # Adjust sample size
  rs <- trimr(rs,1,0)
  rb <- trimr(rb,1,0)
  rc <- trimr(rc,1,0)
  t  <- length(rs)
  
  # Identity weighting matrix
  flag        <- 0
  
  theta0      <- c(0.9, 1.0)
  estResults <- optim(theta0, q, rc=rc, rb=rb, rs=rs, inst=inst,flag=flag, method="BFGS")
  theta <- estResults$par
  qmin <- estResults$val
  
  dof   <- ncol(inst)-length(theta)
  Jstat <- 2*t*qmin
  
  cat('\nResults using identity weighting matrix ')
  cat('\nRisk aversion parameter estimate  = ', theta[2])
  cat('\nJ test statistic                  = ', Jstat)
  
  if (dof > 0.0) 
    cat('\np-value                       = ', 1-pchisq(Jstat,dof))
  
  
  
  # Optimal weighting matrix
  flag        <- 1
  estResults <- optim(theta, q, rc=rc, rb=rb, rs=rs, inst=inst,flag=flag, method="BFGS")
  theta <- estResults$par
  qmin <- estResults$val
  
  dof   <- ncol(inst)-length(theta)
  Jstat <- 2*t*qmin
  
  cat('\n ')
  cat('\nResults using optimal weighting matrix ')
  cat('\nRisk aversion parameter estimate  = ', theta[2])
  cat('\nJ test statistic                  = ', Jstat)
  
  if (dof > 0.0) 
    cat('\np-value                       = ', 1-pchisq(Jstat,dof))
  
}



