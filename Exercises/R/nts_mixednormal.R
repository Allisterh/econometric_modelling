#============================================================================
#
#   Program to simulate a mixed normal distribution and 
#   compare it to the standard normal distribution
#
#============================================================================

rm (list = ls(all=TRUE))
graphics.off()
set.seed(1234, kind="Mersenne-Twister")

#
#------------------------- Helper Functions----------------------------------
#

#load required functions -  figure, trimr, recserar
source("EMTSUtil.R")


#
#--------------- Simulating a Mixed Normal Distribution ---------------------
#

nts_mixednormal <- function() {
  t    <- 4000
  nreps <- 20000
  b    <- rep(0, nreps)
  
  pb <- txtProgressBar(min=0, max=nreps, style=3)
  for (i in seq(nreps)) {  
    # Independent random numbers used for y2 and u
    y2 <- cumsum(rnorm(t))
    u  <- rnorm(t)
  
    # regress y2 on u and scale the slope estimate by t  
    b[i] <- t*lm(u ~ y2 - 1)$coef  
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  cat('\nMean     = ',mean(b)) 
  cat('\nStd.dev. = ',sd(b))
  cat('\nSkewness = ',mean((b-mean(b))^3)/sd(b)^3)
  cat('\nKurtosis = ',mean((b-mean(b))^4)/sd(b)^4)
  
  minx <- -10 
  maxx <- 10
  x    <- seq(from=minx, by=(maxx-minx)/200,length.out=201)
  
  ftrue <- dnorm((x-mean(b))/sd(b) )/sd(b)
  fhat  <- density(b,from=x[1], to=x[length(x)], n=length(x))$y
  
  #**********************************************************************
  #***
  #***     Generate graph
  #***
  #**********************************************************************
  figure()
  par(xaxs="i", yaxs="i")
  
  #--------------------------------------------------------#
  matplot(x,cbind(ftrue,fhat), type="l",
          main = "Comparison of mixed normal and normal distribuion",
          ylab = expression(f(m)),
          xlab = "m",
          xlim = c(-10, 10),
          ylim = c(0, 0.3),
          lty = c(1,3),
          col = 2:3,        
          bty = "l")
    legend("topright", 
           legend=c("Normal", "Mixed Normal"),
           lty=c(1,3), col=2:3)  
}

