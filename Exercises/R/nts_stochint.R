#============================================================================
#
#   Program to simulate the distribution of a stochastic integral
#
#============================================================================

rm (list = ls(all=TRUE))
graphics.off()
set.seed(15, kind="Mersenne-Twister")

#
#------------------------- Helper Functions----------------------------------
#

#load required functions -  figure, trimr, recserar
source("EMTSUtil.R")


#
#---------- Simulating the Distribution of Stochastic Integrals -------------
#

nts_stochint <- function() {
  imoment <- 1
  t     <- 500
  nreps <- 10000      
  delta <- 0.0
  phi   <- 1.0
  sig2  <- 1.0
  sig   <- sqrt(sig2)
  y0    <- 0.0  
  m <- rep(0, nreps)
  
  pb <- txtProgressBar(min=0, max=nreps, style=3)
  
  for (i in seq(nreps)) {
    # Generate random walk and discard first element
    v  <- sqrt(sig2)*rnorm(t+1)
    y  <- trimr( recserar(cbind(delta + v) , cbind(y0) , cbind(phi)),1,0)    
    v  <- trimr(v,1,0)                                      
    if (imoment == 1) {                                                                              
        # Sample mean with standardization based on t^(-0.5)
        m[i] <- (1/sig^2)*sum( trimr(y,0,1)*trimr(v,1,0))*t^(-1)                                  
    }
    else if (imoment == 2) {
      # Standardized sum of y(t-1)^2 * v(t) with standardization based on sig^(-3)xt^(-3/2)                         
      m[i] <- (1/sig^3)*sum( trimr(y^2,0,1)*trimr(v,1,0))*t^(-3/2)       
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  cat('\nSample size = ',t )
  cat('\n')
  
  
  if (imoment == 1) {
    xi <- seq(-0.99, length.out=201, by=0.1)
    fhat  <- density(m, from=xi[1], to=xi[length(xi)], n=length(xi))$y
    
    tmp <- 2*xi + 1                                                    
    jacobian <- 2
    ftrue <- jacobian*(tmp^(-0.5))*exp(-tmp/2)/(gamma(0.5)*sqrt(2))
    
    cat('\nSample mean of m           = ',mean(m) )
    cat('\nTheoretical mean of m      = ',0.0 )
    cat('\nSample variance of m       = ',sd(m)^2)
    cat('\nTheoretical variance of m  = ',0.5)  
  } else if (imoment == 2 ) {
    xi <- seq(-5, length.out=101, by=0.1)
    fhat <- density(m, from=xi[1], to=xi[length(xi)], n=length(xi))$y  
    ftrue <- dnorm(xi)
  
    cat('\nSample mean of m           = ',mean(m))
    cat('\nSample variance of m       = ',sd(m)^2)  
  }
  
  
  #**********************************************************************
  #***
  #***     Generate graph
  #***
  #**********************************************************************
  
  figure()
  par(xaxs="i", yaxs="i")
  #--------------------------------------------------------#
  matplot(xi,cbind(ftrue,fhat), type="l",
          ylim = c(0,2),
          xlim = c(-1, 4),
          xlab = "m",
          ylab = expression(f(m)),
          bty = "l") 
}

