#=========================================================================
#
#   Program to demonstrate Slutsky's theorem by simulation 
#
#=========================================================================

rm (list = ls(all=TRUE))
graphics.off()

set.seed(12, kind="Mersenne-Twister")


# Choose parameters
t  <- 10
n  <- 50000
mu <- 2

zeros <- rep(0, n)
m1 <- zeros        
m2 <- zeros        

pb <- txtProgressBar( min = 0, max = n, style = 3)
for (i in seq(n)) {
  # Generate exponential random numbers from a Gamma distribution
  y <- mu * rgamma(t, 1)
  m1[i] <- ( mean(y) / sd(y) )^2                                                     
  m2[i] <- ( sqrt(t) * mean(y) )^2
  setTxtProgressBar(pb, i)
}
close(pb)

cat('\nSample size    = ', t)
cat('\n')
cat( '\nMoment results for m1 ')
cat('\nMean of m1     = ' , mean(m1))
cat('\nVariance of m1 = ' , sd(m1)^2)
cat('\n')
cat('\nMoment results for m2')
cat('\nMean of m2     = ' , mean(m2))
cat('\nVariance of m2 = ' , sd(m2)^2, '\n')


