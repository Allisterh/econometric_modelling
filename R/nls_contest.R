#=========================================================================
#
#     Program to perform tests on a nonlinear consumption function
#     
#     U.S. data on real consumption and real disposable income(2005 $)
#     1960:Q1 to 2009:Q4 200 observations
#
#=========================================================================
rm (list = ls(all=TRUE))
graphics.off()

#
#--------------------------- Helper Functions -----------------------------------
#

# load required functions - inv, numhess, numgrad
source("EMTSUtil.R")

#-------------------------------------------------------------------------
# Log-likelihood function of constrained model
#-------------------------------------------------------------------------
neglog0 <- function(b,ct,yt) {
  t  <- length( ct )
  e  <- ct - b[1] - b[2]*yt                                 
  s2 <- t(e) %*% e/t                               # Concentrate the variance      
  f  <- - 0.5*log(2*pi) - 0.5*log(s2) - 0.5*e^2/s2 
  lf <- -mean(f)
  return(lf)
}

#-------------------------------------------------------------------------
# Log-likelihood function of unconstrained model
#-------------------------------------------------------------------------
neglog1 <- function(b,ct,yt) {
  t  <- length(ct)
  e  <- ct - b[1] - b[2]*yt^b[3]                                
  s2 <- t(e) %*% e/t                              # Concentrate the variance   
  f  <- - 0.5*log(2*pi) - 0.5*log(s2) - 0.5*e^2/s2
  lf <- -mean(f)
  return(lf)
}

     
#-------------------------------------------------------------------------
# Log-likelihood function of constrained model at each observation
#-------------------------------------------------------------------------
lnlt1 <- function(b,ct,yt) {
  t  <- length(ct)
  e  <- ct - b[1] - b[2]*yt^b[3]                               
  s2 <- t(e) %*% e/t                              # Concentrate the variance   
  lf <- - 0.5*log(2*pi) - 0.5*log(s2) - 0.5*e^2/s2
  return(lf)
}


#
#--------------------------- Nonlinear Consumption Test------------------------------
#

nls_contest <- function() {
  # Load data
  # [cons, inc]
  load('USdata.Rdata')
  cons <- USdata[,1]
  inc <- USdata[,2]  
  
  yt <- inc
  ct <- cons
  t <- length(yt)

  # Estimate the constrained model 
  b0 <- c(-228,  0.9)
  estResults <- optim(b0, neglog0, ct=ct, yt=yt, method = "Nelder-Mead")
  theta0 <- estResults$par
  f0     <- estResults$value

  cat('\nRestricted Parameter Estimates = ')
  theta0 <- c(theta0, 1.000)
  cat(theta0, '\n')

  # Estimate the unconstrained model
  b0 <- theta0
  estResults <- optim(b0, neglog1, ct=ct, yt=yt, method = "Nelder-Mead")
  theta1 <- estResults$par
  f1     <- estResults$value

  f0 <- -f0
  f1 <- -f1
    
  # Compute relevant matrices
  g    <- numgrad(lnlt1,theta0,ct,yt) 
  G    <- rbind(colMeans(g))
  J    <- t(g) %*% g/t
  H    <- numhess(neglog1,theta1,ct,yt)
  invH <- inv(H)

  cat('\nUnestricted Parameter Estimates = ')
  cat(theta1, '\n')
  cat('\nCovariance Matrix of the Parameters\n')
  print(invH/t)
  

  # Perform likelihood ratio test     
  lr <- -2*t*(f0 - f1)               
  cat('\nRestricted Likelihood Function\n')
  cat(t*f0 )
  cat('\nUnrestricted Likelihood Function\n')
  cat(t*f1)
  cat('\n')
  cat('\nLR test and p-value\n')
  cat(lr, 1-pchisq(lr,1))

  # Perform Wald test           
  R <- rbind(c(0, 0, 1))
  theta1 <- rbind(theta1)

  Q <- 1
  W <- t* t( (R %*% t(theta1) - Q) ) %*% inv(R %*% invH %*% t(R)) %*% (R %*% t(theta1) - Q)
  cat('\n')
  cat('\nWald statistic and p-value\n')
  cat(W,  1-pchisq(W,1))

  # Perform Lagrange multiplier test (OPG)    
  cat('\n')
  cat('\nOuter product of gradients matrix\n')
  print(J)   
  LM <- t*(G %*% inv(J) %*% t(G))                          
  
  cat('\n ')  
  cat('\nLM statistic and p-value\n')
  cat(LM,  1-pchisq(LM,1))

  # Perform 2-step LM test    
  x <- cbind(rep(1,t),inc)
  y <- cons 
  
  # Estimate constrained model
  b <- lm(y ~ x - 1)$coef
  e <- y - x %*% b             
  b <- c(b, 1.0)

  # Evaluate derivatives at constrained estimates
  z1 <- rep(1,t)                           
  z2 <- inc^b[3]
  z3 <- b[2]*log(inc)*inc^b[3]
  z  <- cbind(z1 , z2 , z3)

  # Second step regression
  zereg <- lm(e ~ z - 1)$coef
  v  <- e - z %*% (zereg)  
  r2 <- 1 - lm((t(v) %*% v) ~ (t(e) %*% e) -1)$coef

  # LM statistic
  LM <- t*r2
  cat('\n\n2-step LM statistic and p-value\n')
  cat(LM, 1-pchisq(LM,1))
  
}
 


