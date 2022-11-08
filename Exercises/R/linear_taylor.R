#=========================================================================
#
#   Program to estimate parameters of a Taylor Rule using various methods.
#
#=========================================================================

rm (list = ls(all=TRUE))
graphics.off()

#
#--------------------------- Helper Functions -----------------------------------
#

# load utility functions - inv 
source("EMTSUtil.R")

#
#--------------------------- Linear Taylor Rule -----------------------------------
#

linear_taylor <- function () {  
  load("taylor.RData")
  
  # Choose data 1987:Q1 to 1999:Q4
  infl <- taylor[104:nrow(taylor), 2]  
  ygap <- taylor[104:nrow(taylor), 3]
  ffr  <- taylor[104:nrow(taylor), 4]
  
  t <-length(ffr)
  tmp <- matrix(c(t, sum(infl), sum(ygap), 
                  sum(infl), sum(infl*infl), sum(infl*ygap),
                  sum(ygap), sum(infl*ygap), sum(ygap*ygap)), ncol = 3, byrow = T)
  print(tmp)
  tmp1 <- cbind(c(sum(ffr), sum(ffr*infl), sum(ffr*ygap)))
    
  print(tmp1)

  # OLS estimates long-hand
  betahat1 <- inv(tmp) %*% tmp1

  # OLS estimates
  ones <- rep(1, t)
  x <- matrix( c(ones, infl,  ygap), ncol = 3)
  
  y <- ffr
  betahat  <- lm(y ~ x - 1)$coef
    
  cat('\nML/OLS estimates - both methods\n' )
  print(cbind(betahat1, betahat) )
    
  e       <- y - x %*% betahat              
  sig2hat <- t(e) %*% e/t
  cat('\nML/OLS estimate of variance\n' )
  cat('\n', sig2hat)

  vcov <- c(sig2hat) * inv( t(x) %*% x)
  cat('\nCovariance matrix of parameters\n' )
  print(vcov)
     
  # Wald test of restrictions b(2)=1.5 and b(3)=0.5
  R <- matrix(c(0, 1, 0,
         0, 0, 1), nrow = 2, byrow = T)
  Q <- cbind(c(1.5,
              0.5))

  W <- t( (R %*% betahat - Q) ) %*% inv(R %*% vcov %*% t(R)) %*% (R %*% betahat - Q)
  pv <- 1 - pchisq(W,2)

  cat('\nWald test results\n')
  cat('\nWald statistic    = ', W)
  cat('\np value           = ', pv)  
  
  #***************************************************************************
  # Plot the graph
  #**************************************************************************
  figure()
  # plot of the entire data
  xs <- as.Date(taylor[,1], origin="0000-01-01")
  plot (xs, taylor[,4], type="l", col="green", xaxs="i", yaxs="i",
        main = "US data of various rates",
        xlab = "Years",
        ylab = "Percent",
        ylim = c(-10,20),
        bty = "l", lty = 1)
  lines(xs, taylor[,2], col="blue", lty = 3)
  lines(xs, taylor[,3], col="red", lty = 5)
  legend("topright", 
         c("Federal funds rate","Inflation rate", "Output gap rate"), 
         lty = c(1,3,5),
         col = c("green", "blue", "red"), lwd = 1)
}






