#============================================================================
#
#   Monte Carlo experiment on the Breitung test
#
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()

#
#------------------------- Helper Functions----------------------------------
#

# load required functions -  trimr, recserar
source("EMTSUtil.R")


#
#------------- Unit Root Test Without Lags or Long-Run Variance -------------
#

unit_breitung_size <- function()
{
  reps <- 10000
  Tv   <- c(100,200)
  arv  <- c(0,0.3,0.6,0.9)
  mav  <- c(-0.8,-0.4,0,0.4,0.8)
  rej  <- array(0, c(length(arv),length(mav)))
  rho  <- rep(0, reps)
 
  cat('\n---------------------------------------------')
  totRows = length(Tv)*length(arv)*length(mav)
  rnames <- seq(totRows)
  cnames <- c('T', 'arv','mav', 'Rej Freq')
  tabResults <- matrix(nrow=totRows, ncol=length(arv), dimnames=list(rnames, cnames))
  row <- 1
  for (i in seq(Tv)) {
    for (j in seq(arv)) {
      for (m in seq(mav)) {
        set.seed(1234, kind="Mersenne-Twister")
        pb <- txtProgressBar(min=0, max=reps, style=3)
        for (k in seq(reps)) {
          v      <- rnorm(Tv[i]+1)
          y      <- cumsum(recserar(cbind( trimr(v,1,0)+mav[m]*trimr(v,0,1)),cbind(v[1]),cbind( arv[j])))
          rho[k] <- sum(cumsum(y)^2)/(sum(y^2)*(length(y)^2)) 
          setTxtProgressBar(pb, k)
        }
        close(pb)
        rej[j,m] <- mean(rho < 0.02)
        tabResults[row, ] <- c(Tv[i], arv[j], mav[m], rej[j,m])          
        row <- row+1
      }
    }   
  }
  print(tabResults)
}
