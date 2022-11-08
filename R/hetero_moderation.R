#=========================================================================
#
#    Program to test the "Great Moderation" hypothesis, using 
#    real annual US GDP data per capita from 1946 to 2006.
#                                                         
#=========================================================================
rm (list = ls(all=TRUE))
graphics.off()

#
#------------------------- Helper Functions-------------------------------
#
#load required functions - inv, figure
source("EMTSUtil.R")

#-----------------------------------------------------------------------
#      Negative log-likelihood function (unconstrained)   
#-----------------------------------------------------------------------
neglog <- function(b,y,x,flag ) {
  if (flag) {
    mu   <- b[1] + b[2]*x
    sig2 <- exp(b[3] + b[4]*x)    
  } else {
    mu <- b[1]
    sig2 <- exp(b[2] + b[3]*x)  
  }
  lnl  <- -(1/2)*log(2*pi*sig2) - (y - mu)^2 /(2*sig2) 
  lf <- -mean( lnl )
  return(lf)  
}

#-----------------------------------------------------------------------
#      Negative log-likelihood function (constrained)   
#-----------------------------------------------------------------------
neglog0 <- function(b,y,x,flag ) {
  if (flag) {
    mu   <- b[1] + b[2]*x
    sig2 <- exp(b[3] + 0*x)    
  } else {
    mu   <- b[1]
    sig2 <- exp(b[2] + 0*x)    
  }
  lnl  <- -(1/2)*log(2*pi*sig2) - (y - mu)^2 /(2*sig2)
  lf   <- -mean( lnl )
  return(lf)
}

#
#----------------------- The Great Moderation ---------------------------
#
hetero_moderation <- function() {
  # read gdp data
  gdp <-  c(11241,
            10924,
            11206,
            10957,
            11717,
            12412,
            12668,
            13032,
            12719,
            13389,
            13410,
            13435,
            13088,
            13783,
            13840,
            13933,
            14552,
            14971,
            15624,
            16420,
            17290,
            17532,
            18196,
            18573,
            18392,
            18771,
            19555,
            20485,
            20195,
            19961,
            20822,
            21565,
            22526,
            22982,
            22666,
            23007,
            22347,
            23146,
            24593,
            25382,
            26024,
            26664,
            27514,
            28221,
            28429,
            28007,
            28556,
            28941,
            29741,
            30128,
            30880,
            31886,
            32833,
            33904,
            34755,
            34645,
            34837,
            35361,
            36300,
            37052,
            37752)
  
  # Compute growth rate
  y <- 100*( log(gdp[-1]) - log(gdp[-length(gdp)] ) )
  
  #**************************************************************************
  #**
  #**     Generate graph
  #**
  #**************************************************************************
  figure()
  par(xaxs="i", yaxs="i")
  gr <- seq(1947, 2006, 1)
  plot(gr,y, type="l",
       xlab = "Years",
       ylab = "Growth Rate of US Per Capital GDP",
       bty = "l")
  # Generate descriptive statistics   
  cat('\nSample mean (1947-1983)     = ', mean(y[1:37]))
  cat('\nSample mean (1984-2006)     = ', mean(y[38:60]))
  cat('\nSample variance (1984-2006) = ', mean((y[1:37]-mean(y[1:37]))^2) )
  cat('\nSample variance (1984-2006) = ', mean((y[38:60]-mean(y[38:60]))^2))
  
  d <- cbind(c(rep(0, 37),rep(1, 23)))       #     Construct dummy variable   
  t <- length(y)
  
  # Estimate the unconstrained model 
  theta <- 0.1*rep(1, 4)
  flag  <- TRUE
  estResults <- optim(theta, neglog, y=y, x=d, flag=flag, method="BFGS", hessian=T)
  theta1 <- estResults$par
  a1 <- estResults$value
  h1 <- estResults$hessian      
  
  omega1 <- inv(h1)
  lnl1 <- -a1                         
  
  cat('\n')
  cat('\nML estimate of the variance (1946-1983) = ',exp(theta1[3]))
  cat('\nML estimate of the variance (1984-2006) = ',exp(theta1[3] + theta1[4]) )
  
  # Estimate the constrained model    
  theta <- 0.1*rep(1, 3)
  estResults <- optim(theta, neglog0, y=y, x=d, flag=flag, method="BFGS")  
  a0 <- estResults$value
  lnl0 <- -a0
  
  # LR test for heteroskedasticity
  lr <- -2*t*(lnl0 - lnl1)
  
  cat('\n')
  cat('\nLR statistic            = ',lr) 
  cat('\np-value                 = ',1-pchisq(lr,1))
  
  # Wald test for heteroskedasticity
  r <- rbind(c(0 , 0 , 0 , 1))
  q <- 0
  wd <- t* t( (r %*% theta1 - q) ) %*% inv(r %*% omega1 %*% t(r)) %*% (r %*% theta1 - q)
  cat('\nWald statistic          = ',wd) 
  cat('\np-value                 = ',1-pchisq(wd,1))
  
  # LM test (regression form)    
  x <- cbind(rep(1, t), d)   
  
  # Stage 1 regression
  b <- lm(y ~ x - 1)$coef
  u <- y - x %*% b    
  w <- cbind(rep(1,t) , d)
  v <- u^2
  
  # Stage 2 regression
  b  <- lm(v ~ w - 1)$coef
  e  <- v - w %*% b
  r2 <- 1 - sum(e^2)/sum( (v-mean(v))^2 )
  lm <- t*r2
  
  cat('\nLM statistic            = ',lm) 
  cat('\np-value                 = ',1 - pchisq(lm,1))
  
  # Wald test that beta1 = 0  
  r  <- rbind(c(0 , 1 , 0 , 0))
  q  <- 0
  wd <- t( (r %*% theta1 - q) ) %*% inv(r %*% omega1 %*% t(r)) %*% (r %*% theta1 - q)
  
  cat('\n')
  cat('\nWald statistic of beta1 = 0   = ',wd)
  cat('\np-value                       = ',1 - pchisq(wd,1))
  
  #-----------------------------------------------------------------------
  # Tests based on the assumption that beta1 = 0.
  #-----------------------------------------------------------------------
  
  # Estimate the unconstrained model but with beta1 = 0     
  theta <- c(0.1, 0.1, 0.1)
  flag <- 0
  estResults <- optim(theta, neglog, y=y, x=d, flag=flag, method="BFGS", hessian=T)
  theta1 <- estResults$par
  a1 <- estResults$value
  h1 <- estResults$hessian
  
  omega1 <- inv(h1)
  lnl1 <- -a1      
  
  # Estimate the constrained model but with beta1 = 0   
  theta <- c(0.1, 0.1)
  estResults <- optim(theta, neglog0, y=y, x=d, flag=flag, method="BFGS")  
  a0 <- estResults$value
  
  omega1 <- inv(h1)
  lnl0 <- -a0
  
  # LR test for heteroskedasticity
  cat('\n')
  lr <- -2*t*(lnl0 - lnl1)
  cat('\nLR statistic            = ',lr) 
  cat('\np-value                 = ',1-pchisq(lr,1))
  
  # Wald test   
  r <- rbind(c(0 , 0 , 1))
  q <- 0
  wd <- t* t( (r %*% theta1 - q)) %*% inv(r %*% omega1 %*% t(r)) %*% (r %*% theta1 - q)
  cat('\nWald statistic          = ',wd) 
  cat('\np-value                 = ',1-pchisq(wd,1))
  
  # LM test (regression form)   
  x <- cbind(rep(1, t))
  
  # Stage 1 regression
  b <- lm(y ~ x - 1)$coef
  u <- y - x %*% b    
  w <- cbind(rep(1, t), d)
  v <- u^2
  # Stage 2 regression
  b  <- lm(v ~ w - 1)$coef
  e  <- v - w %*% b
  r2 <- 1 - sum(e^2)/sum( (v-mean(v))^2 )
  lm <- t*r2
  cat('\nLM statistic            = ', lm) 
  cat('\np-value                 = ', 1-pchisq(lm, 1))  
}
