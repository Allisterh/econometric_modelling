#=========================================================================
#
#   Program to demonstrate by simulation the necessary and sufficient 
#   conditions of the weak law of large numbers
#
#=========================================================================

rm (list = ls(all=TRUE))
graphics.off()

set.seed(12, kind="Mersenne-Twister")

# Choose parameters
t <- 100
n <- 50000

zeros <- rep(0, n)

m1 <- zeros        
m2 <- zeros       
m3 <- zeros       
m4 <- zeros       

pb <- txtProgressBar(min = 0, max = n, style = 3)
for(i in seq(n))
{
  y <- runif(t) - 0.5          # Simulate y from a uniform distribution (-0.5,0.5).
  #y <- 2.0 + rt(t, 3)     # Simulate Student t (mue,dof=3) random numbers            

  m1[i] <- sum(y^1)/t
  m2[i] <- sum(y^2)/t                                                       
  m3[i] <- sum(y^3)/t                                                        
  m4[i] <- sum(y^4)/t
  setTxtProgressBar(pb, i)
}
close(pb)


cat('\nSample size    = ', t,  '\n' )
cat( '\n' )

cat('\nMean of m1     = ' , mean(m1),  '\n' )
cat('\nVariance of m1 = ' , sd(m1)^2,  '\n' )

cat( '\n' )
cat('\nMean of m2     = ' , mean(m2),  '\n' )
cat('\nVariance of m2 = ' , sd(m2)^2,  '\n' )

cat( '\n' )
cat('\nMean of m3     = ' , mean(m3),  '\n' )
cat('\nVariance of m3 = ' , sd(m3)^2,  '\n' )

cat( '\n' )
cat('\nMean of m4     = ' , mean(m4),  '\n' )
cat('\nVariance of m4 = ' , sd(m4)^2,  '\n' )

